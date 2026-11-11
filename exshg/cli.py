"""Command line: exshg <bse_dir> [abc] [--scissor Ha] [--wmax eV] [--every n] [--out prefix]."""
import argparse
from pathlib import Path

import numpy as np

from .core import AU_PMV, HA_EV, excitonic_shg


def main(argv=None):
    ap = argparse.ArgumentParser(
        prog="exshg",
        description="Second-harmonic generation with excitons from a finished ELK BSE run.")
    ap.add_argument("bse_dir", help="directory of the ELK run (PMAT.OUT, EVBSE.OUT, EPSILON_BSE_*.OUT, ...)")
    ap.add_argument("component", nargs="?", default="222", help="Cartesian component abc, default 222 (yyy)")
    ap.add_argument("--scissor", type=float, default=0.0, help="consistent scissor shift in Hartree (default 0)")
    ap.add_argument("--wmax", type=float, default=1.8, help="largest fundamental photon energy in eV (default 1.8)")
    ap.add_argument("--every", type=int, default=4, help="use every n-th point of the ELK wplot grid (default 4)")
    ap.add_argument("--out", default=None, help="output prefix (default <bse_dir>/XSHG_BSE_<abc>[_sc<D>])")
    ap.add_argument("--txt", action="store_true", help="also write a text table (eV, pm/V)")
    a = ap.parse_args(argv)

    res = excitonic_shg(a.bse_dir, a.component, scissor=a.scissor, wmax_ev=a.wmax, every=a.every)
    tag = f"XSHG_BSE_{a.component}" + ("" if a.scissor == 0 else f"_sc{a.scissor:g}")
    prefix = Path(a.out) if a.out else Path(a.bse_dir) / tag
    npz = prefix.with_name(prefix.name + ".npz")
    np.savez(npz, **res)
    print(f"  wrote {npz}")
    if a.txt:
        txt = prefix.with_name(prefix.name + ".txt")
        ip, x = res["chi_ip"] * AU_PMV, res["chi_x"] * AU_PMV
        np.savetxt(txt, np.c_[res["w"] * HA_EV, 2 * res["w"] * HA_EV, ip.real, ip.imag, x.real, x.imag],
                   fmt="%12.6f",
                   header="hw(eV) 2hw(eV) Re_chi_IP Im_chi_IP Re_chi_BSE Im_chi_BSE   "
                          "(pm/V, per unit-cell volume, retarded)")
        print(f"  wrote {txt}")


if __name__ == "__main__":
    main()
