"""Command line interface: exshg <bse_dir> [abc] [options]."""
import argparse
import sys
from pathlib import Path

import numpy as np

from . import __version__
from .core import AU_PMV, HA_EV, ExshgError, excitonic_shg, film_factor

EPILOG = """examples:
  exshg .                          chi_yyy of the ELK run in the current directory
  exshg run_dir 222 --plot         also draw the spectrum (needs matplotlib)
  exshg run_dir 222 --thickness 6.15 --plot
                                   values for a film of thickness 6.15 Angstrom (2D material)
  exshg run_dir 111 --scissor 0.02646
                                   chi_xxx with a scissor shift of 0.72 eV

documentation: https://yucj-520.github.io/exshg/
"""


def build_parser():
    ap = argparse.ArgumentParser(
        prog="exshg",
        description="Second-harmonic generation with excitons from a finished ELK Bethe-Salpeter run.",
        epilog=EPILOG, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("bse_dir", help="directory of the ELK run")
    ap.add_argument("component", nargs="?", default="222",
                    help="Cartesian component abc with 1 = x, 2 = y, 3 = z (default: 222 = yyy)")
    ap.add_argument("--scissor", type=float, default=0.0, metavar="HA",
                    help="scissor shift of the conduction bands in Hartree (default: 0)")
    ap.add_argument("--wmax", type=float, default=1.8, metavar="EV",
                    help="largest fundamental photon energy in eV (default: 1.8, i.e. 2hw up to 3.6 eV)")
    ap.add_argument("--every", type=int, default=4, metavar="N",
                    help="use every N-th point of the ELK wplot grid (default: 4)")
    ap.add_argument("--thickness", type=float, default=None, metavar="ANGSTROM",
                    help="layer thickness of a 2D material: the table and the plot are given for a film "
                         "of this thickness instead of per unit-cell volume")
    ap.add_argument("--plot", action="store_true", help="also write a figure <prefix>.png (needs matplotlib)")
    ap.add_argument("--out", default=None, metavar="PREFIX",
                    help="output prefix (default: <bse_dir>/XSHG_BSE_<abc>)")
    ap.add_argument("-q", "--quiet", action="store_true", help="print only errors")
    ap.add_argument("--version", action="version", version=f"exshg {__version__}")
    ap.add_argument("--txt", action="store_true", help=argparse.SUPPRESS)   # accepted for compatibility
    return ap


def write_table(path, res, scale, norm):
    ip, x = res["chi_ip"] * AU_PMV * scale, res["chi_x"] * AU_PMV * scale
    header = (f"exshg {__version__}: chi(2) in pm/V, {norm}, retarded response\n"
              "hw(eV)       2hw(eV)      Re_chi_IP    Im_chi_IP    Re_chi_BSE   Im_chi_BSE")
    np.savetxt(path, np.c_[res["w"] * HA_EV, 2 * res["w"] * HA_EV, ip.real, ip.imag, x.real, x.imag],
               fmt="%12.6f", header=header)


def write_plot(path, res, scale, component, norm):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    e2 = 2 * res["w"] * HA_EV
    plt.rcParams.update({"font.size": 13})
    fig, ax = plt.subplots(figsize=(5.2, 5.2))
    ax.plot(e2, np.abs(res["chi_ip"]) * AU_PMV * scale, "--", color="0.45", lw=1.8, label="independent particles")
    ax.plot(e2, np.abs(res["chi_x"]) * AU_PMV * scale, "-", color="#CC3311", lw=2.4, label="with excitons (BSE)")
    ax.set_xlabel(r"$2\hbar\omega$ (eV)")
    ax.set_ylabel(r"$|\chi^{(2)}_{%s}|$ (pm/V)" % "".join("xyz"[int(ch) - 1] for ch in component))
    ax.set_xlim(e2[0], e2[-1])
    ax.set_ylim(bottom=0)
    ax.set_box_aspect(1)
    ax.set_title(norm, fontsize=11)
    ax.legend(frameon=False)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def main(argv=None):
    a = build_parser().parse_args(argv)
    say = (lambda *args, **kw: None) if a.quiet else print
    try:
        if a.thickness is not None and a.thickness <= 0:
            raise ExshgError("--thickness must be positive")
        say(f"exshg {__version__}: second-harmonic generation with excitons from ELK BSE\n")
        res = excitonic_shg(a.bse_dir, a.component, scissor=a.scissor, wmax_ev=a.wmax, every=a.every,
                            verbose=not a.quiet)
        if a.thickness is None:
            scale, norm = 1.0, "per unit-cell volume"
        else:
            scale, norm = film_factor(a.bse_dir, a.thickness), f"film of thickness {a.thickness:g} Angstrom"
    except ExshgError as err:
        print(f"exshg: error: {err}", file=sys.stderr)
        return 1

    tag = f"XSHG_BSE_{a.component}" + ("" if a.scissor == 0 else f"_sc{a.scissor:g}")
    prefix = Path(a.out) if a.out else Path(a.bse_dir) / tag
    files = [prefix.with_name(prefix.name + ".npz"), prefix.with_name(prefix.name + ".txt")]
    np.savez(files[0], **res)
    write_table(files[1], res, scale, norm)
    if a.plot:
        try:
            png = prefix.with_name(prefix.name + ".png")
            write_plot(png, res, scale, a.component, norm)
            files.append(png)
        except ImportError:
            print("exshg: warning: matplotlib is not installed, no figure written "
                  "(pip install matplotlib)", file=sys.stderr)

    e2 = 2 * res["w"] * HA_EV
    say(f"\n  Result ({norm})")
    for label, key in (("independent particles", "chi_ip"), ("with excitons", "chi_x")):
        y = np.abs(res[key]) * AU_PMV * scale
        i = int(y.argmax())
        say(f"    {label:22s}: largest |chi| = {y[i]:8.1f} pm/V at 2hw = {e2[i]:.2f} eV")
    say("\n  Files written")
    for path in files:
        say(f"    {path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
