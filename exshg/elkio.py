"""Readers for the ELK output files used by exshg.

All quantities are returned in atomic units, as ELK writes them.
"""
import re
from pathlib import Path

import numpy as np


def input_value(d, key):
    """Tokens of the line that follows the block name `key` in <d>/elk.in."""
    lines = (Path(d) / "elk.in").read_text().split("\n")
    hits = [j for j, line in enumerate(lines) if line.strip() == key]
    if not hits:
        raise KeyError(f"block '{key}' not found in {Path(d) / 'elk.in'}")
    return lines[hits[0] + 1].split()


def read_wplot(d):
    """Frequency grid of the ELK run (block wplot): w_i = w0 + (w1 - w0) i / nwplot."""
    lines = (Path(d) / "elk.in").read_text().split("\n")
    j = [i for i, line in enumerate(lines) if line.strip() == "wplot"][0]
    nw = int(lines[j + 1].split()[0])
    w0, w1 = map(float, lines[j + 2].split()[:2])
    return w0 + (w1 - w0) * np.arange(nw) / nw


def read_kpoints(d):
    """KPOINTS.OUT: k-points in lattice coordinates (nk, 3) and their weights (nk,)."""
    lines = (Path(d) / "KPOINTS.OUT").read_text().splitlines()
    nk = int(lines[0].split()[0])
    rows = np.array([[float(x) for x in line.split()[1:5]] for line in lines[1:1 + nk]])
    return rows[:, :3], rows[:, 3]


def read_eigval(d):
    """EIGVAL.OUT: second-variational eigenvalues and occupations, both (nk, nstsv)."""
    txt = (Path(d) / "EIGVAL.OUT").read_text().splitlines()
    nk, nst = int(txt[0].split()[0]), int(txt[1].split()[0])
    e, f = np.zeros((nk, nst)), np.zeros((nk, nst))
    i, ik = 2, 0
    while ik < nk:
        if "k-point, vkl" in txt[i]:
            vals = np.array([[float(x) for x in txt[i + 2 + j].split()[1:3]] for j in range(nst)])
            e[ik], f[ik] = vals[:, 0], vals[:, 1]
            ik += 1
            i += 2 + nst
        else:
            i += 1
    return e, f


def read_pmat(d, nk, nst):
    """PMAT.OUT (task 120): momentum matrix elements p[k, i, j, a] = <i| p_a |j>.

    Direct-access records: vkl(3) real*8, nstsv int*4, pmat(nstsv, nstsv, 3) complex*16
    in Fortran order.
    """
    rec = np.dtype([("vkl", "<f8", 3), ("nst", "<i4"), ("p", "<c16", nst * nst * 3)])
    raw = np.fromfile(Path(d) / "PMAT.OUT", dtype=rec)
    if len(raw) != nk or not np.all(raw["nst"] == nst):
        raise ValueError("PMAT.OUT does not match EIGVAL.OUT (k-points or number of states)")
    # C-order (k, a, j, i) is Fortran pmat(i, j, a)
    p = raw["p"].reshape(nk, 3, nst, nst).transpose(0, 3, 2, 1)
    return raw["vkl"], p


def read_avec_omega(d):
    """Lattice vectors (columns) from LATTICE.OUT and the unit-cell volume from INFO.OUT."""
    info = (Path(d) / "INFO.OUT").read_text()
    omega = float(re.search(r"Unit cell volume\s*:\s*(\S+)", info).group(1))
    lat = (Path(d) / "LATTICE.OUT").read_text()
    a = [list(map(float, re.search(rf"vector a{i}\s*:\s*(\S+)\s+(\S+)\s+(\S+)", lat).groups()))
         for i in (1, 2, 3)]
    return np.array(a).T, omega


def fortran_records(path):
    """Records of a gfortran sequential unformatted file, as uint8 views.

    Records longer than 2 GiB are stored as subrecords: the leading marker is negative
    while more subrecords follow. They are joined here.
    """
    m = np.memmap(path, np.uint8, mode="r")
    recs, parts, pos = [], [], 0
    while pos < m.size:
        n = int(np.frombuffer(m[pos:pos + 4], np.int32)[0])
        parts.append(m[pos + 4:pos + 4 + abs(n)])
        pos += abs(n) + 8
        if n >= 0:
            recs.append(parts[0] if len(parts) == 1 else np.concatenate(parts))
            parts = []
    return recs


def read_evbse(d):
    """EVBSE.OUT (task 186): exciton energies (nm,) and eigenvectors X[a, S] (nm, nm)."""
    recs = fortran_records(Path(d) / "EVBSE.OUT")
    nm = int(np.frombuffer(recs[0], np.int32)[0])
    om = np.frombuffer(recs[1], np.float64).copy()
    x = np.frombuffer(recs[2], np.complex128).reshape(nm, nm, order="F")
    return om, x


def read_complex_spectrum(path):
    """ELK two-block spectrum file (real part, blank line, imaginary part)."""
    blocks = [s for s in Path(path).read_text().strip().split("\n\n") if s.strip()]
    re_, im_ = np.loadtxt(blocks[0].split("\n")), np.loadtxt(blocks[1].split("\n"))
    return re_[:, 0], re_[:, 1] + 1j * im_[:, 1]
