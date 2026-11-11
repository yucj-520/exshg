"""Second-harmonic generation with excitons from an ELK Bethe-Salpeter run.

Length-gauge density matrix to second order (Aversa and Sipe, PRB 52, 14636 (1995)),
hbar = 1, charge q = +1, z = w + i eta (retarded):

    rho1^b_nm(z)   = -f_nm r^b_nm g1_nm,                      g1 = 1/(w_nm - z)
    S^bc_nm(z)     =  [r^c, rho1^b]_nm + i (rho1^b_nm);c       (second-order source)
    rho2^bc_nm(2z) =  S^bc_nm g2_nm,                           g2 = 1/(w_nm - 2z)
    chi^abc(2z)    = (i / 2z) (1/Omega) sum_k w_k Tr( v^a rho2^bc )

With v = i w_nm r, the pole part of (i/2z) v rho2 is r_mn S_nm g2_nm (polarization form).
On the BSE transition space a = (k, v, c) this independent-particle 2w propagator is
replaced by the Tamm-Dancoff exciton propagator,

    rho2_cv(a) = sum_a' [ sum_S X_aS X*_a'S / ( Omega_S - 2z) ] S_cv(a'),
    rho2_vc(a) = sum_a' [ sum_S X*_aS X_a'S / (-Omega_S - 2z) ] S_vc(a').

The first-order step (rho1) and all pairs outside the BSE window stay independent-particle.
The same substitution at first order must rebuild ELK's own EPSILON_BSE; this is checked
in every run and fixes the X / X* convention.
"""
import time
from pathlib import Path

import numpy as np

from .elkio import (input_value, read_avec_omega, read_complex_spectrum, read_eigval, read_evbse,
                    read_kpoints, read_pmat, read_wplot)

HA_EV = 27.211386
AU_PMV = 24.44          # chi(2) in Hartree atomic units -> pm/V
BOHR_ANGSTROM = 0.529177
W_MIN = 0.004           # Ha; frequencies below this are skipped (the response has a 1/2z prefactor)

# files read from the ELK run directory, and the ELK task that writes them
REQUIRED_FILES = {
    "elk.in": "the ELK input file",
    "KPOINTS.OUT": "ground state (task 0 or 1)",
    "EIGVAL.OUT": "ground state (task 0 or 1)",
    "INFO.OUT": "ground state (task 0 or 1)",
    "LATTICE.OUT": "ground state (task 0 or 1)",
    "PMAT.OUT": "momentum matrix elements (task 120)",
    "EVBSE.OUT": "BSE eigenvectors (task 186)",
}


class ExshgError(Exception):
    """A problem with the input that the user can fix."""


def rv_k(e, es, f, p, tol):
    """Position, velocity and derivative matrices of one k-point.

    e, es : bare and scissor-shifted eigenvalues; f : occupations; p[i, j, a] : momentum.
    Returns r[a], v[a], Delta[n, m, a], the generalised derivative rd[b, a] = r^b_nm;a,
    f_nm and the scissored w_nm. r and r;a use bare energies (they are wavefunction
    properties); pairs with |w_nm| <= tol count as degenerate (no r, v = p).
    """
    n = len(e)
    w, ws = e[:, None] - e[None, :], es[:, None] - es[None, :]
    nd = np.abs(ws) > tol
    winv = np.where(nd, 1.0 / np.where(nd, w, 1.0), 0.0)
    r = np.stack([-1j * p[:, :, a] * winv for a in range(3)])
    vd = p[np.arange(n), np.arange(n), :].real
    dlt = vd[:, None, :] - vd[None, :, :]
    rd = np.empty((3, 3, n, n), complex)
    for a in range(3):
        for b in range(3):
            t = -(r[b] * dlt[:, :, a] + r[a] * dlt[:, :, b])
            s = r[a] @ (w * r[b]) - (w * r[b]) @ r[a]
            rd[b, a] = (t + 1j * s) * winv
    v = np.stack([np.where(nd, 1j * ws * r[a], p[:, :, a]) for a in range(3)])
    fnm = f[:, None] - f[None, :]
    return r, v, dlt, rd, fnm, ws


def film_factor(bse_dir, thickness):
    """Factor that converts chi(2) per unit-cell volume to a film of `thickness` (Angstrom).

    For a layer in a supercell with in-plane vectors a1, a2: Omega / (|a1 x a2| thickness).
    """
    avec, omega = read_avec_omega(Path(bse_dir))
    area = np.linalg.norm(np.cross(avec[:, 0], avec[:, 1]))
    return omega / area * BOHR_ANGSTROM / thickness


def _input(d, key, default=None):
    try:
        return input_value(d, key)
    except KeyError:
        if default is not None:
            return default
        raise ExshgError(f"the block '{key}' must be given explicitly in {d / 'elk.in'}") from None


def excitonic_shg(bse_dir, component="222", scissor=0.0, wmax_ev=1.8, every=4, verbose=True):
    """chi^(2)_abc(2w; w, w) with and without excitons from a finished ELK BSE run.

    bse_dir   : directory of the ELK run (tasks 120, 121, 180, 185, 186, 187; reducek 0; scissor 0)
    component : Cartesian indices 'abc' (1 = x, 2 = y, 3 = z), e.g. '222' for yyy
    scissor   : rigid shift of the conduction bands in Hartree, applied consistently: all
                energy denominators and the exciton energies are shifted, the matrix
                elements r and r;k are not
    wmax_ev   : largest fundamental photon energy in eV
    every     : take every n-th point of the ELK wplot grid

    Returns a dict with the fundamental frequency w (Ha), chi_ip and chi_x (retarded,
    atomic units, per unit-cell volume), the per-exciton resonant contributions, the
    exciton energies and the first-order check (eps_x, eps_ipb, eps_elk, linear_corr,
    linear_ratio) on the same grid.

    Raises ExshgError if the run directory is incomplete or not usable.
    """
    d = Path(bse_dir)
    say = print if verbose else (lambda *a, **k: None)
    component = str(component)
    if len(component) != 3 or any(ch not in "123" for ch in component):
        raise ExshgError(f"component must be three digits out of 1, 2, 3 (e.g. 222 for yyy), got '{component}'")
    if every < 1:
        raise ExshgError("every must be 1 or larger")
    ia, ib, ic = (int(ch) - 1 for ch in component)
    if not d.is_dir():
        raise ExshgError(f"{d} is not a directory")
    eps_name = f"EPSILON_BSE_{ia + 1}{ib + 1}.OUT"
    needed = dict(REQUIRED_FILES)
    needed[eps_name] = f"BSE dielectric function (task 187, with the pair '{ia + 1} {ib + 1}' in optcomp)"
    missing = [f"  {name:22s} {what}" for name, what in needed.items() if not (d / name).is_file()]
    if missing:
        raise ExshgError(f"files missing in {d}:\n" + "\n".join(missing))

    if float(_input(d, "scissor", ["0.0"])[0]) != 0.0:
        raise ExshgError("the ELK run must have scissor 0; apply the scissor shift in exshg instead (--scissor)")
    eta = float(_input(d, "swidth")[0])
    nv, nc = int(_input(d, "nvbse")[0]), int(_input(d, "ncbse")[0])
    _input(d, "wplot")
    _, wk = read_kpoints(d)
    e, f = read_eigval(d)
    nk, nst = e.shape
    try:
        _, p = read_pmat(d, nk, nst)
    except ValueError as err:
        raise ExshgError(f"{err}; rerun task 120 with the k-points and bands of the ground state") from None
    _, omega = read_avec_omega(d)
    if not np.allclose(wk, wk[0]):
        raise ExshgError("the run uses symmetry-reduced k-points; rerun ELK with 'reducek 0'")
    fmax = f.max()
    nocc = int(round(f.sum() / (nk * fmax)))
    # occupations carry the smearing of the ELK run; every k must have the same nocc filled bands
    if f[:, :nocc].min() < 0.5 * fmax or f[:, nocc:].max() > 0.5 * fmax:
        raise ExshgError("the number of occupied bands differs between k-points (metal); "
                         "exshg needs an insulator or semiconductor")
    if nv > nocc or nocc + nc > nst:
        raise ExshgError(f"nvbse = {nv}, ncbse = {nc} do not fit the {nocc} occupied and {nst - nocc} empty bands")
    es = e.copy()
    es[:, nocc:] += scissor
    vst, cst = np.arange(nocc - nv, nocc), np.arange(nocc, nocc + nc)
    om, x = read_evbse(d)
    nm = len(om)
    if nm != nk * nv * nc:
        raise ExshgError(f"EVBSE.OUT has {nm} states, but nk x nvbse x ncbse = {nk * nv * nc}. "
                         "Check that nvbse and ncbse in elk.in are those of the BSE run and that "
                         "the BSE used the Tamm-Dancoff approximation (the ELK default)")

    wg = read_wplot(d)
    w = wg[(wg > W_MIN) & (wg < wmax_ev / HA_EV)][::every]
    if len(w) < 2:
        raise ExshgError("fewer than two frequencies selected; check wplot in elk.in, --wmax and --every")
    z = w + 1j * eta
    nz = len(z)
    say(f"  Run directory : {d}")
    say(f"  Component     : chi_{''.join('xyz'[i] for i in (ia, ib, ic))} ({component})")
    say(f"  k-points      : {nk}")
    say(f"  Bands         : {nst} ({nocc} occupied)")
    say(f"  BSE space     : {nm} transitions = {nk} k x {nv} valence x {nc} conduction")
    say(f"  Broadening    : {eta:g} Ha ({eta * HA_EV:.3f} eV)")
    say("  Scissor       : " + ("none" if scissor == 0 else f"{scissor:g} Ha ({scissor * HA_EV:.3f} eV)"))
    say(f"  Frequencies   : {nz} points, 2hw = {2 * w[0] * HA_EV:.2f} to {2 * w[-1] * HA_EV:.2f} eV")
    say("")

    t0 = time.time()
    say("  [1/3] Independent-particle response", flush=True)
    # BSE basis: a = k*nv*nc + v*nc + c, v upward to the valence-band maximum
    vi = np.repeat(vst[None, :], nc, 0).T.ravel()
    ci = np.tile(cst, nv)
    chi_ip = np.zeros(nz, complex)
    s_cv = np.zeros((nz, nm), complex)
    s_vc = np.zeros((nz, nm), complex)
    v_vc = np.zeros(nm, complex)
    v_cv = np.zeros(nm, complex)
    r_cv1 = np.zeros(nm, complex)
    r_vc1 = np.zeros(nm, complex)
    va_vc1 = np.zeros(nm, complex)
    va_cv1 = np.zeros(nm, complex)
    ip_lin = np.zeros(nz, complex)
    pairs = {(ib, ic), (ic, ib)}                      # E^b E^c symmetric
    step = max(1, nk // 10)
    for k in range(nk):
        r, v, dlt, rd, fnm, ws = rv_k(e[k], es[k], f[k], p[k], eta)
        g1 = 1.0 / (ws[None] - z[:, None, None])
        g2 = 1.0 / (ws[None] - 2 * z[:, None, None])
        src = 0
        for b, c in pairs:
            rho1 = -fnm * r[b] * g1
            src = src + (r[c] @ rho1 - rho1 @ r[c]
                         + 1j * (-fnm * (rd[b, c] * g1 - r[b] * dlt[:, :, c] * g1 ** 2)))
        src = src / len(pairs)
        chi_ip += wk[k] * np.einsum("mn,znm->z", v[ia], src * g2)
        sl = slice(k * nv * nc, (k + 1) * nv * nc)
        s_cv[:, sl] = src[:, ci, vi]
        s_vc[:, sl] = src[:, vi, ci]
        v_vc[sl] = r[ia][vi, ci]                      # r^a: polarization form on the BSE space
        v_cv[sl] = r[ia][ci, vi]
        # first-order ingredients for the EPSILON_BSE check
        r_cv1[sl] = -fnm[ci, vi] * r[ib][ci, vi]
        r_vc1[sl] = -fnm[vi, ci] * r[ib][vi, ci]
        va_vc1[sl] = r[ia][vi, ci]
        va_cv1[sl] = r[ia][ci, vi]
        ip_lin += wk[k] * (np.einsum("a,za->z", va_vc1[sl], r_cv1[sl] * g1[:, ci, vi])
                           + np.einsum("a,za->z", va_cv1[sl], r_vc1[sl] * g1[:, vi, ci]))
        if (k + 1) % step == 0 or k + 1 == nk:
            say(f"        k-point {k + 1:5d} of {nk}   {time.time() - t0:6.0f} s", flush=True)

    de = (es[:, ci] - es[:, vi]).ravel()              # w_cv per basis state, k-major
    xc = x.conj()
    om = om + scissor                                 # Tamm-Dancoff: rigid shift

    def x_resp(src_cv, src_vc, vcv, vvc, zz):
        """v_vc rho_cv + v_cv rho_vc on the BSE space with exciton propagators at zz."""
        a = xc.T @ src_cv.T
        rcv = x @ (a / (om[:, None] - zz[None, :]))
        b_ = x.T @ src_vc.T
        rvc = xc @ (b_ / (-om[:, None] - zz[None, :]))
        return vvc @ rcv + vcv @ rvc

    # first order: must rebuild ELK's EPSILON_BSE
    say(f"  [2/3] Check against {eps_name} of ELK", flush=True)
    lin_x = x_resp(np.repeat(r_cv1[None], nz, 0), np.repeat(r_vc1[None], nz, 0), va_cv1, va_vc1, z) * wk[0]
    eps_x = 1 + 4 * np.pi * lin_x / omega
    eps_ipb = 1 + 4 * np.pi * ip_lin / omega
    we, epse = read_complex_spectrum(d / eps_name)
    eps_elk = np.interp(w - scissor, we, epse.real) + 1j * np.interp(w - scissor, we, epse.imag)
    m = np.abs(eps_elk.imag) > 0.1 * np.abs(eps_elk.imag).max()
    lin_corr = float(np.corrcoef(eps_x.imag, eps_elk.imag)[0, 1])
    lin_ratio = float(np.median(eps_x.imag[m] / eps_elk.imag[m]))
    if scissor != 0:
        verdict = "approximate comparison (scissor shift applied)"
    elif lin_corr > 0.999 and abs(lin_ratio - 1) < 0.02:
        verdict = "passed"
    else:
        verdict = "FAILED"
    say(f"        correlation {lin_corr:.4f}, amplitude ratio {lin_ratio:.4f}: {verdict}")
    if verdict == "FAILED":
        say("        WARNING: the BSE dielectric function of ELK is not reproduced.")
        say("        The files were not read correctly or belong to different runs.")
        say("        Do not use the result below.")

    # second order: replace the 2w propagator on the BSE space
    say("  [3/3] Excitons at the second-harmonic frequency", flush=True)
    g2ip_cv = 1.0 / (de[None] - 2 * z[:, None])
    g2ip_vc = 1.0 / (-de[None] - 2 * z[:, None])
    ip_space = (np.einsum("a,za->z", v_vc, s_cv * g2ip_cv) + np.einsum("a,za->z", v_cv, s_vc * g2ip_vc)) * wk[0]
    xs = np.array([x_resp(s_cv[i:i + 1], s_vc[i:i + 1], v_cv, v_vc, 2 * z[i:i + 1])[0]
                   for i in range(nz)]) * wk[0]
    # resonant contribution of each exciton S to chi_x
    contrib = ((v_vc @ x)[:, None] * (xc.T @ s_cv.T) / (om[:, None] - 2 * z[None, :])) * wk[0] / omega
    chi_ip = (1j / (2 * z)) / omega * chi_ip
    chi_x = chi_ip + (xs - ip_space) / omega
    say(f"        done   {time.time() - t0:6.0f} s")
    return dict(w=w, chi_ip=chi_ip, chi_x=chi_x, contrib=contrib, Om=om, eps_x=eps_x, eps_ipb=eps_ipb,
                eps_elk=eps_elk, linear_corr=lin_corr, linear_ratio=lin_ratio)
