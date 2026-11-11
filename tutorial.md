---
layout: default
title: "Tutorial: monolayer MoS₂"
nav_order: 3
---

# Tutorial: monolayer MoS₂

This tutorial runs the complete workflow for χ⁽²⁾_yyy of free-standing monolayer MoS₂ and
compares the result with the reference shipped in the repository. All files are in
[`examples/mos2_monolayer`](https://github.com/yucj-520/exshg/tree/main/examples/mos2_monolayer).

You need a working ELK 11.2.3 (see [Installation](install.md)). The ELK part takes about
half an hour on 40 cores; `exshg` itself takes less than a minute.

## Step 1. Get the example

```bash
git clone https://github.com/yucj-520/exshg.git
pip install "./exshg[plot]"
cd exshg/examples/mos2_monolayer
```

## Step 2. Set the species path

Open `elk.in` and replace the path in the `sppath` block by the `species` directory of
your ELK copy (keep the trailing slash):

```
sppath
  '/path/to/elk/species/'
```

Nothing else needs to be changed. The input asks ELK for

| Task | What it does | File used by `exshg` |
|---|---|---|
| `0` | ground state (PBE with spin–orbit coupling) | `EIGVAL.OUT`, `KPOINTS.OUT`, `INFO.OUT`, `LATTICE.OUT` |
| `120` | momentum matrix elements | `PMAT.OUT` |
| `121` | independent-particle dielectric function | |
| `180` | screened interaction for the BSE | |
| `185` | BSE Hamiltonian | |
| `186` | BSE eigenvalues and eigenvectors | `EVBSE.OUT` |
| `187` | BSE dielectric function | `EPSILON_BSE_22.OUT` |

## Step 3. Run ELK

```bash
export OMP_NUM_THREADS=2
mpirun -np 20 elk > elk.out
```

(or simply `elk > elk.out` for an OpenMP-only build). On a cluster, put these lines in a
job script. The run is finished when `EPSILON_BSE_22.OUT` exists.

## Step 4. Run exshg

```bash
exshg . 222 --thickness 6.15 --plot
```

- `222` selects the component χ_yyy.
- `--thickness 6.15` gives the values for a film of 6.15 Å, the thickness of one MoS₂
  layer, instead of the value per volume of the 30 Å supercell.
- `--plot` also writes a figure.

The output should be

```
exshg 0.3.0: second-harmonic generation with excitons from ELK BSE

  Run directory : .
  Component     : chi_yyy (222)
  k-points      : 144
  Bands         : 148 (26 occupied)
  BSE space     : 2304 transitions = 144 k x 4 valence x 4 conduction
  Broadening    : 0.003 Ha (0.082 eV)
  Scissor       : none
  Frequencies   : 62 points, 2hw = 0.23 to 3.55 eV

  [1/3] Independent-particle response
        k-point    14 of 144        1 s
        ...
        k-point   144 of 144       13 s
  [2/3] Check against EPSILON_BSE_22.OUT of ELK
        correlation 1.0000, amplitude ratio 1.0000: passed
  [3/3] Excitons at the second-harmonic frequency
        done       14 s

  Result (film of thickness 6.15 Angstrom)
    independent particles : largest |chi| =    784.9 pm/V at 2hw = 3.55 eV
    with excitons         : largest |chi| =    631.0 pm/V at 2hw = 2.03 eV

  Files written
    XSHG_BSE_222.npz
    XSHG_BSE_222.txt
    XSHG_BSE_222.png
```

Step `[2/3]` is the built-in check: `exshg` rebuilds the BSE dielectric function of ELK
from the same files and compares. `passed` means that the files were read correctly.

## Step 5. Compare with the reference

Compare your `XSHG_BSE_222.png` with the reference figure
`reference/XSHG_BSE_222.png`:

![chi2 of monolayer MoS2](https://raw.githubusercontent.com/yucj-520/exshg/main/examples/mos2_monolayer/reference/XSHG_BSE_222.png)

and your table `XSHG_BSE_222.txt` with `reference/XSHG_BSE_222.txt`. Small differences in
the last digits between machines are normal.

The peaks of the curve with excitons are at 2ħω = 1.05, 2.03, 2.73 and 3.06 eV. The first
two are the two-photon resonances with the A exciton (1.07 eV) and the C exciton
(2.01 eV) of the same BSE run (`EIGVAL_BSE.OUT`).

## Notes

- The calculation uses PBE energies without a quasiparticle correction and a 12 × 12
  k-mesh. It tests the workflow; it is not a converged prediction for MoS₂.
- By symmetry, χ_xxx vanishes for this orientation. `exshg . 111` gives a curve about 2 %
  of χ_yyy, which is the size of the symmetry breaking of the ELK BSE itself on this mesh.
- For your own material, start from this input, replace `avec` and `atoms`, and choose
  `ngridk`, `nvbse`, `ncbse` and `nempty`. Keep `reducek 0` and `scissor 0.0`.
