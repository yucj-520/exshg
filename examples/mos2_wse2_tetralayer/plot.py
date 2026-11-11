#!/usr/bin/env python3
"""Plot |chi_yyy| of the MoS2/WSe2 tetralayer with and without excitons, both stackings.

usage: python plot.py [H_table.txt R_table.txt] [out.png]
The tables are --txt outputs of exshg (default: the reference spectra shipped here).
Values are rescaled from the 40 A supercell to a film of thickness 27 A.
"""
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

here = Path(__file__).resolve().parent
tables = sys.argv[1:3] if len(sys.argv) > 2 else [here / "reference" / "XSHG_BSE_222_H.txt",
                                                  here / "reference" / "XSHG_BSE_222_R.txt"]
out = sys.argv[3] if len(sys.argv) > 3 else here / "tetralayer_shg.png"
CELL, FILM = 40.0, 27.0                                # Angstrom

plt.rcParams.update({"font.size": 13})
fig, axs = plt.subplots(1, 2, figsize=(9.6, 5.0), sharey=True)
for ax, table, title in zip(axs, tables, ("H-like stacking", "R-like stacking")):
    d = np.loadtxt(table)
    e2 = d[:, 1]
    ip = np.abs(d[:, 2] + 1j * d[:, 3]) * CELL / FILM
    bse = np.abs(d[:, 4] + 1j * d[:, 5]) * CELL / FILM
    ax.plot(e2, ip, "--", color="0.45", lw=1.8, label="independent particles")
    ax.plot(e2, bse, "-", color="#CC3311", lw=2.4, label="with excitons (BSE)")
    ax.set_title(title, fontsize=13)
    ax.set_xlabel(r"$2\hbar\omega$ (eV)")
    ax.set_xlim(1.0, 3.4)
    ax.set_box_aspect(1)
    w = (e2 >= 1.8) & (e2 <= 2.25)
    i = np.argmax(np.where(w, bse, 0))
    print(f"{title}: mean |chi| for 2hw in 1.8-2.25 eV: IP {ip[w].mean():.1f}, BSE {bse[w].mean():.1f} pm/V; "
          f"BSE maximum {bse[i]:.1f} pm/V at {e2[i]:.2f} eV")
axs[0].set_ylabel(r"$|\chi^{(2)}_{yyy}|$ (pm/V)")
axs[0].set_ylim(bottom=0)
axs[0].legend(frameon=False, loc="upper left")
fig.tight_layout()
fig.savefig(out, dpi=200)
print("wrote", out)
