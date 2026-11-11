#!/usr/bin/env python3
"""Plot |chi_yyy| of monolayer MoS2 with and without excitons.

usage: python plot.py [table.txt] [out.png]
The table is the --txt output of exshg (default: the reference spectrum shipped here).
Values are rescaled from the 30 A supercell to a film of thickness 6.15 A.
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

here = Path(__file__).resolve().parent
table = sys.argv[1] if len(sys.argv) > 1 else here / "reference" / "XSHG_BSE_222.txt"
out = sys.argv[2] if len(sys.argv) > 2 else here / "mos2_shg.png"
CELL, FILM = 30.0, 6.15                                # Angstrom

d = np.loadtxt(table)
e2 = d[:, 1]
ip = np.abs(d[:, 2] + 1j * d[:, 3]) * CELL / FILM
bse = np.abs(d[:, 4] + 1j * d[:, 5]) * CELL / FILM

plt.rcParams.update({"font.size": 13})
fig, ax = plt.subplots(figsize=(5.2, 5.2))
ax.plot(e2, ip, "--", color="0.45", lw=1.8, label="independent particles")
ax.plot(e2, bse, "-", color="#CC3311", lw=2.4, label="with excitons (BSE)")
ax.set_xlabel(r"$2\hbar\omega$ (eV)")
ax.set_ylabel(r"$|\chi^{(2)}_{yyy}|$ (pm/V)")
ax.set_xlim(0.2, 3.6)
ax.set_ylim(bottom=0)
ax.set_box_aspect(1)
ax.legend(frameon=False)
fig.tight_layout()
fig.savefig(out, dpi=200)
print("wrote", out)
