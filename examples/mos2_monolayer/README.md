# Monolayer MoS₂

ELK input and reference result for χ⁽²⁾_yyy of free-standing monolayer MoS₂
(PBE with spin–orbit coupling, 12 × 12 k-mesh, 4 valence and 4 conduction bands in the BSE).

```bash
# 1. ELK: ground state, PMAT.OUT and BSE (set sppath in elk.in first)
elk > elk.out

# 2. second-harmonic response with excitons
exshg . 222 --txt

# 3. plot
python plot.py XSHG_BSE_222.txt
```

Expected output of step 2:

```
linear check Im eps_22 vs ELK EPSILON_BSE: corr 1.0000, median ratio 1.0000
```

`reference/XSHG_BSE_222.txt` is the spectrum obtained with ELK 11.2.3, and `mos2_shg.png`
the corresponding plot. The peaks of the curve with excitons lie at 2ħω = 1.05 eV and
2.03 eV, the energies of the A (1.07 eV) and C (2.01 eV) excitons of the same BSE run.

The calculation uses PBE energies without a quasiparticle correction and a coarse k-mesh,
so the exciton energies are too low and the binding energy is too large. The example is
meant as a test of the workflow, not as a converged result. The xxx component is forbidden
by symmetry and can be used as a further check (`exshg . 111`).

The ELK steps take about half an hour on 20 cores; `exshg` takes a few minutes.
