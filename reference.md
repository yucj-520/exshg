# Reference

[Back to contents](index.md)

## Requirements on the ELK run

| Block in `elk.in` | Requirement |
|---|---|
| `tasks` | `0` (or `1`), `120`, `121`, `180`, `185`, `186`, `187` |
| `scissor` | `0.0` |
| `reducek` | `0` unless the cell has no symmetry |
| `swidth` | broadening η in Hartree; read by `exshg` |
| `wplot` | frequency grid; read by `exshg` |
| `nvbse`, `ncbse` | number of valence and conduction bands in the BSE; read by `exshg` |
| `optcomp` | must contain the pair `a b` of the requested component `abc`, so that `EPSILON_BSE_ab.OUT` exists |

Files read from the run directory: `elk.in`, `KPOINTS.OUT`, `EIGVAL.OUT`, `INFO.OUT`,
`LATTICE.OUT`, `PMAT.OUT`, `EVBSE.OUT`, `EPSILON_BSE_ab.OUT`.

## Command line

```
exshg <bse_dir> [abc] [--scissor D] [--wmax E] [--every n] [--out prefix] [--txt]
```

| Argument | Meaning | Default |
|---|---|---|
| `bse_dir` | directory of the ELK run | |
| `abc` | Cartesian component, 1 = x, 2 = y, 3 = z | `222` |
| `--scissor D` | consistent scissor shift in Hartree | 0 |
| `--wmax E` | largest fundamental photon energy ħω in eV | 1.8 |
| `--every n` | use every n-th point of the ELK `wplot` grid | 4 |
| `--out prefix` | output prefix | `<bse_dir>/XSHG_BSE_<abc>` |
| `--txt` | also write a text table in eV and pm/V | off |

`python -m exshg ...` is equivalent.

## Python

```python
from exshg import excitonic_shg, HA_EV, AU_PMV

res = excitonic_shg("path/to/elk_run", "222", scissor=0.0, wmax_ev=1.8, every=4)
two_hw = 2 * res["w"] * HA_EV            # eV
chi = abs(res["chi_x"]) * AU_PMV         # pm/V
```

## Output

`<prefix>.npz`:

| Key | Content |
|---|---|
| `w` | fundamental frequency ω (Hartree) |
| `chi_ip` | independent-particle χ⁽²⁾ (atomic units) |
| `chi_x` | χ⁽²⁾ with excitons (atomic units) |
| `contrib` | resonant contribution of each exciton to `chi_x`, shape (excitons, ω) |
| `Om` | exciton energies (Hartree), including the scissor shift |
| `eps_x` | dielectric function with excitons rebuilt by `exshg` |
| `eps_ipb` | independent-particle dielectric function restricted to the BSE window |
| `eps_elk` | `EPSILON_BSE` of ELK on the same grid |
| `linear_corr`, `linear_ratio` | result of the built-in check |

`<prefix>.txt` (with `--txt`), six columns:
ħω (eV), 2ħω (eV), Re χ_IP, Im χ_IP, Re χ_BSE, Im χ_BSE (pm/V).

## Units and conventions

- Atomic units in the `.npz` file. Multiply χ⁽²⁾ by `AU_PMV` = 24.44 for pm/V and energies
  by `HA_EV` = 27.211386 for eV.
- χ⁽²⁾ is normalised to the unit-cell volume. For a layer in a supercell of height c,
  multiply by c / d to obtain the value for a film of thickness d.
- The response is retarded. ELK task 125 writes the advanced response, which is the
  complex conjugate.

## Resources

The BSE eigenvectors are held in memory: 16 N² bytes for N = n_k × `nvbse` × `ncbse`
transitions.

| Run | N | Memory | Time |
|---|---|---|---|
| monolayer MoS₂, 12 × 12, 4 + 4 bands | 2304 | below 2 GB | a few minutes on 4 cores |
| tetralayer, 6 × 6, 16 + 16 bands | 9216 | 7.5 GB | 2 minutes on 20 cores |

NumPy uses the threads given by `OMP_NUM_THREADS`.

## Troubleshooting

| Message or symptom | Cause |
|---|---|
| `block 'swidth' not found in elk.in` | `swidth`, `wplot`, `nvbse` and `ncbse` must be given explicitly in `elk.in` |
| `the BSE run must use the full uniform k-grid (reducek 0)` | the run used symmetry-reduced k-points; rerun ELK with `reducek 0` |
| `PMAT.OUT does not match EIGVAL.OUT` | task 120 was run with different k-points or bands than the ground state |
| `EVBSE.OUT has ... states, expected ...` | `nvbse`/`ncbse` in `elk.in` differ from those of the BSE run |
| `EPSILON_BSE_ab.OUT` not found | add the pair `a b` to `optcomp` and rerun task 187 |
| linear check correlation not 1 | files not read correctly (compiler, ELK version) or inconsistent runs; do not use the result |
