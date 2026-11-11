"""End-to-end tests on a synthetic ELK run directory written here in the documented file formats.

The reference is analytic: if the electron-hole kernel is zero, the BSE eigenvectors are the
identity and the exciton energies are the band-energy differences, so the response with
"excitons" must equal the independent-particle response, with and without a scissor shift.
"""
import struct

import numpy as np
import pytest

from exshg import ExshgError, excitonic_shg, film_factor
from exshg.cli import main

NK, NST, NOCC, NV, NC = 4, 6, 2, 2, 2
OMEGA, ETA = 250.0, 0.01


def _record(payload):
    n = struct.pack("<i", len(payload))
    return n + payload + n


def _spectrum(path, w, y):
    with open(path, "w") as fh:
        for part in (y.real, y.imag):
            for wi, yi in zip(w, part):
                fh.write(f"{wi:18.10G}{yi:18.10G}\n")
            fh.write("\n")


@pytest.fixture()
def run_dir(tmp_path):
    rng = np.random.default_rng(7)
    e = np.sort(np.concatenate([rng.uniform(-0.40, -0.15, (NK, NOCC)),
                                rng.uniform(0.10, 0.60, (NK, NST - NOCC))], axis=1), axis=1)
    f = np.zeros((NK, NST))
    f[:, :NOCC] = 2.0
    p = rng.normal(size=(NK, NST, NST, 3)) + 1j * rng.normal(size=(NK, NST, NST, 3))
    p = 0.5 * (p + np.conj(np.swapaxes(p, 1, 2)))                      # Hermitian momentum
    kpt = rng.uniform(0, 1, (NK, 3))

    (tmp_path / "elk.in").write_text(
        "tasks\n  0\n\nswidth\n  %g\n\nwplot\n  200  100  0\n  0.0  0.5\n\nnvbse\n  %d\n\nncbse\n  %d\n" % (ETA, NV, NC))
    (tmp_path / "KPOINTS.OUT").write_text(
        f"{NK:6d} : nkpt; k-point, vkl, wkpt, nmat below\n"
        + "".join(f"{i + 1:6d}{k[0]:18.10G}{k[1]:18.10G}{k[2]:18.10G}{1 / NK:18.10G}{100:8d}\n"
                  for i, k in enumerate(kpt)))
    with open(tmp_path / "EIGVAL.OUT", "w") as fh:
        fh.write(f"{NK:6d} : nkpt\n{NST:6d} : nstsv\n")
        for i in range(NK):
            fh.write(f"\n{i + 1:6d}{kpt[i, 0]:18.10G}{kpt[i, 1]:18.10G}{kpt[i, 2]:18.10G} : k-point, vkl\n")
            fh.write(" (state, eigenvalue and occupancy below)\n")
            for j in range(NST):
                fh.write(f"{j + 1:6d}{e[i, j]:22.14G}{f[i, j]:18.10G}\n")
    (tmp_path / "INFO.OUT").write_text(f"Unit cell volume      :   {OMEGA:.6f}\n")
    (tmp_path / "LATTICE.OUT").write_text(
        "vector a1 :   5.0   0.0   0.0\nvector a2 :   0.0   5.0   0.0\nvector a3 :   0.0   0.0  10.0\n")
    rec = np.zeros(NK, dtype=[("vkl", "<f8", 3), ("nst", "<i4"), ("p", "<c16", NST * NST * 3)])
    rec["vkl"], rec["nst"] = kpt, NST
    rec["p"] = p.transpose(0, 3, 2, 1).reshape(NK, -1)                 # Fortran pmat(i, j, a)
    rec.tofile(tmp_path / "PMAT.OUT")
    # zero kernel: X = identity, Omega = e_c - e_v in the ELK order a = k*nv*nc + v*nc + c
    om = np.array([e[k, c] - e[k, v] for k in range(NK)
                   for v in range(NOCC - NV, NOCC) for c in range(NOCC, NOCC + NC)])
    nm = len(om)
    (tmp_path / "EVBSE.OUT").write_bytes(
        _record(struct.pack("<i", nm)) + _record(om.tobytes())
        + _record(np.eye(nm, dtype=complex).tobytes()))
    w = 0.5 * np.arange(200) / 200
    _spectrum(tmp_path / "EPSILON_BSE_22.OUT", w, 1 + 1 / (0.3 - w - 0.02j))   # placeholder, replaced below
    return tmp_path


def test_zero_kernel_gives_independent_particles(run_dir):
    res = excitonic_shg(run_dir, "222", verbose=False)
    assert np.abs(res["chi_ip"]).max() > 0
    assert np.allclose(res["chi_x"], res["chi_ip"], rtol=1e-9, atol=1e-12 * np.abs(res["chi_ip"]).max())
    assert np.allclose(res["eps_x"], res["eps_ipb"], rtol=1e-9)


def test_zero_kernel_with_scissor(run_dir):
    res = excitonic_shg(run_dir, "222", scissor=0.05, verbose=False)
    assert np.allclose(res["chi_x"], res["chi_ip"], rtol=1e-9, atol=1e-12 * np.abs(res["chi_ip"]).max())
    bare = excitonic_shg(run_dir, "222", verbose=False)
    assert not np.allclose(res["chi_ip"], bare["chi_ip"], rtol=1e-3)   # the shift does change the spectrum


def test_linear_check_reads_the_elk_spectrum(run_dir):
    first = excitonic_shg(run_dir, "222", every=1, verbose=False)
    _spectrum(run_dir / "EPSILON_BSE_22.OUT", first["w"], first["eps_x"])
    res = excitonic_shg(run_dir, "222", every=1, verbose=False)
    assert res["linear_corr"] == pytest.approx(1.0, abs=1e-8)
    assert res["linear_ratio"] == pytest.approx(1.0, abs=1e-8)


def test_film_factor(run_dir):
    # cell 5 x 5 x 10 bohr: height 10 bohr = 5.29177 Angstrom
    assert film_factor(run_dir, 5.29177) == pytest.approx(OMEGA / 25.0 / 10.0, rel=1e-6)


def test_command_line_writes_files(run_dir, capsys):
    assert main([str(run_dir), "222", "--thickness", "3.0"]) == 0
    out = capsys.readouterr().out
    assert "Files written" in out
    table = np.loadtxt(run_dir / "XSHG_BSE_222.txt")
    data = np.load(run_dir / "XSHG_BSE_222.npz")
    assert table.shape == (len(data["w"]), 6)
    assert "film of thickness 3 Angstrom" in (run_dir / "XSHG_BSE_222.txt").read_text()


@pytest.mark.parametrize("component", ["22", "224", "yyy"])
def test_bad_component(run_dir, component):
    with pytest.raises(ExshgError):
        excitonic_shg(run_dir, component, verbose=False)


def test_missing_file_is_reported(run_dir, capsys):
    (run_dir / "PMAT.OUT").unlink()
    assert main([str(run_dir)]) == 1
    err = capsys.readouterr().err
    assert "PMAT.OUT" in err and "task 120" in err


def test_reduced_k_grid_is_rejected(run_dir):
    lines = (run_dir / "KPOINTS.OUT").read_text().splitlines()
    parts = lines[1].split()
    parts[4] = "0.5"
    lines[1] = "  ".join(parts)
    (run_dir / "KPOINTS.OUT").write_text("\n".join(lines) + "\n")
    with pytest.raises(ExshgError, match="reducek"):
        excitonic_shg(run_dir, "222", verbose=False)


def test_metal_is_rejected(run_dir):
    txt = (run_dir / "EIGVAL.OUT").read_text().splitlines()
    i = [j for j, line in enumerate(txt) if "k-point, vkl" in line][0]
    parts = txt[i + 2 + NOCC].split()                      # first empty band of the first k-point
    txt[i + 2 + NOCC] = f"{parts[0]:>6}{parts[1]:>22}{2.0:18.10G}"
    parts = txt[i + 2 + NOCC - 1].split()                  # and empty its last filled band at the next k-point
    j = [j for j, line in enumerate(txt) if "k-point, vkl" in line][1]
    parts = txt[j + 2 + NOCC - 1].split()
    txt[j + 2 + NOCC - 1] = f"{parts[0]:>6}{parts[1]:>22}{0.0:18.10G}"
    (run_dir / "EIGVAL.OUT").write_text("\n".join(txt) + "\n")
    with pytest.raises(ExshgError, match="insulator"):
        excitonic_shg(run_dir, "222", verbose=False)

