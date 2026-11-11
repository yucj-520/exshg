---
layout: default
title: "Reference"
nav_order: 6
---

# Reference

## Requirements on the ELK run

| Block in `elk.in` | Requirement |
|---|---|
| `tasks` | `0` (or `1`), `120`, `121`, `180`, `185`, `186`, `187` |
| `scissor` | `0.0` |
| `reducek` | `0` unless the cell has no symmetry |
| `swidth` | broadening η in Hartree; must be given explicitly |
| `wplot` | frequency grid; must be given explicitly |
| `nvbse`, `ncbse` | number of valence and conduction bands in the BSE; must be given explicitly |
| `optcomp` | must contain the pair `a b` of the requested component `abc`, so that `EPSILON_BSE_ab.OUT` exists |

The material must be an insulator or semiconductor, and the BSE must use the Tamm–Dancoff
approximation (the ELK default).

Files read from the run directory: `elk.in`, `KPOINTS.OUT`, `EIGVAL.OUT`, `INFO.OUT`,
`LATTICE.OUT`, `PMAT.OUT`, `EVBSE.OUT`, `EPSILON_BSE_ab.OUT`.

## Command line

```
exshg <bse_dir> [abc] [--scissor HA] [--wmax EV] [--every N]
                      [--thickness ANGSTROM] [--plot] [--out PREFIX] [-q]
```

| Argument | Meaning | Default |
|---|---|---|
| `bse_dir` | directory of the ELK run | |
| `abc` | Cartesian component, 1 = x, 2 = y, 3 = z | `222` |
| `--scissor HA` | scissor shift of the conduction bands in Hartree | 0 |
| `--wmax EV` | largest fundamental photon energy ħω in eV | 1.8 |
| `--every N` | use every N-th point of the ELK `wplot` grid | 4 |
| `--thickness ANGSTROM` | layer thickness of a two-dimensional material; the table and the figure are given for a film of this thickness | per unit-cell volume |
| `--plot` | also write a figure (needs Matplotlib) | off |
| `--out PREFIX` | output prefix | `<bse_dir>/XSHG_BSE_<abc>` |
| `-q`, `--quiet` | print only errors | |
| `--version` | print the version | |
| `-h`, `--help` | print the help text | |

`python -m exshg ...` is equivalent. The exit status is 0 on success and 1 if the run
directory cannot be used.

## Python

```python
from exshg import excitonic_shg, film_factor, HA_EV, AU_PMV

res = excitonic_shg("path/to/elk_run", "222", scissor=0.0, wmax_ev=1.8, every=4)
two_hw = 2 * res["w"] * HA_EV                       # eV
chi = abs(res["chi_x"]) * AU_PMV                    # pm/V, per unit-cell volume
chi_film = chi * film_factor("path/to/elk_run", 6.15)   # film of 6.15 Angstrom
```

`excitonic_shg` raises `exshg.ExshgError` with an explanation if the run directory cannot
be used. Pass `verbose=False` to suppress the printed report.

## Output

`<prefix>.txt`, six columns: ħω (eV), 2ħω (eV), Re χ_IP, Im χ_IP, Re χ_BSE, Im χ_BSE
(pm/V). The first line states the normalisation (per unit-cell volume, or the film
thickness given with `--thickness`).

`<prefix>.png` (with `--plot`): the absolute value of χ⁽²⁾ with and without excitons
against 2ħω.

`<prefix>.npz`, always in atomic units and per unit-cell volume:

| Key | Content |
|---|---|
| `w` | fundamental frequency ω (Hartree) |
| `chi_ip` | independent-particle χ⁽²⁾ |
| `chi_x` | χ⁽²⁾ with excitons |
| `contrib` | resonant contribution of each exciton to `chi_x`, shape (excitons, ω) |
| `Om` | exciton energies (Hartree), including the scissor shift |
| `eps_x` | dielectric function with excitons rebuilt by `exshg` |
| `eps_ipb` | independent-particle dielectric function restricted to the BSE window |
| `eps_elk` | `EPSILON_BSE` of ELK on the same grid |
| `linear_corr`, `linear_ratio` | result of the built-in check |

## Units and conventions

- The `.npz` file is in atomic units. Multiply χ⁽²⁾ by `AU_PMV` = 24.44 for pm/V and
  energies by `HA_EV` = 27.211386 for eV.
- χ⁽²⁾ is normalised to the unit-cell volume Ω. For a layer in a supercell with in-plane
  lattice vectors a₁ and a₂, the value for a film of thickness d is larger by the factor
  Ω / (|a₁ × a₂| d). This is what `--thickness` and `film_factor` apply.
- The response is retarded. ELK task 125 writes the advanced response, which is the
  complex conjugate.
- Frequencies below 0.004 Ha (0.11 eV) are skipped.

## Resources

The BSE eigenvectors are held in memory: 16 N² bytes for N = n_k × `nvbse` × `ncbse`
transitions.

| Run | N | Memory | Time on 20 cores |
|---|---|---|---|
| monolayer MoS₂, 12 × 12, 4 + 4 bands | 2304 | 0.6 GB | 20 seconds |
| tetralayer, 6 × 6, 16 + 16 bands | 9216 | 8 GB | 3 minutes |

NumPy uses the number of threads given by `OMP_NUM_THREADS`.

## Troubleshooting

| Message | Cause and remedy |
|---|---|
| `files missing in ...` | the listed ELK tasks have not been run in this directory; the message names the task for each file |
| `the block '...' must be given explicitly in elk.in` | add `swidth`, `wplot`, `nvbse` and `ncbse` to `elk.in` with the values used in the run |
| `the ELK run must have scissor 0` | rerun ELK with `scissor 0.0` and apply the shift with `--scissor` |
| `the run uses symmetry-reduced k-points` | rerun ELK with `reducek 0` |
| `PMAT.OUT does not match EIGVAL.OUT` | rerun task 120 with the k-points and bands of the ground state |
| `EVBSE.OUT has ... states, but nk x nvbse x ncbse = ...` | `nvbse`/`ncbse` in `elk.in` differ from those of the BSE run, or the BSE was not run in the Tamm–Dancoff approximation |
| `the number of occupied bands differs between k-points (metal)` | `exshg` treats insulators and semiconductors only |
| check reports `FAILED` | the files were not read correctly (compiler, ELK version) or belong to different runs; do not use the result |
