---
layout: default
title: "Method"
nav_order: 5
---

# Method

## Independent particles

The independent-particle response is obtained from the length-gauge density matrix to
second order (Aversa and Sipe, 1995). With z = ω + iη, occupations f_nm = f_n − f_m and
transition energies ω_nm:

```
ρ1^b_nm(z)    = − f_nm r^b_nm / (ω_nm − z)
S^bc_nm(z)    =   [r^c, ρ1^b]_nm + i (ρ1^b_nm);c
ρ2^bc_nm(2z)  =   S^bc_nm / (ω_nm − 2z)
χ^abc(2z)     =   (i / 2z) (1/Ω) Σ_k w_k Tr( v^a ρ2^bc )
```

The position matrix elements are r_nm = −i p_nm / ω_nm for non-degenerate states, with the
momentum matrix elements p_nm from `PMAT.OUT`. The generalised derivative r_nm;k is
evaluated with the sum rule over bands. Ω is the unit-cell volume and η is the ELK
broadening `swidth`.

## Excitons at 2ω

The BSE of ELK is solved in the Tamm–Dancoff approximation on the transitions
a = (k, v, c) between `nvbse` valence and `ncbse` conduction bands. It gives the exciton
energies Ω_S and eigenvectors X_S (`EVBSE.OUT`).

With v = i ω_nm r, the resonant part of the last line above is r_mn S_nm / (ω_nm − 2z).
On the BSE transition space, the independent-particle propagator 1/(ω_cv − 2z) is replaced
by the exciton propagator:

```
ρ2_cv = Σ_S X_S  (X_S† S_cv) / ( Ω_S − 2z)
ρ2_vc = Σ_S X_S* (X_Sᵀ S_vc) / (−Ω_S − 2z)
```

All transitions outside the BSE window keep their independent-particle form. If the
electron–hole kernel is zero, X is the identity and Ω_S = ω_cv, and the result is the
independent-particle χ⁽²⁾.

## Built-in check

The same substitution at first order gives the dielectric function with excitons. ELK
writes this quantity itself (`EPSILON_BSE_ab.OUT`), so the two can be compared. Every run
prints

```
  [2/3] Check against EPSILON_BSE_22.OUT of ELK
        correlation 1.0000, amplitude ratio 1.0000: passed
```

The correlation tests the reading of the files and the phase convention of the
eigenvectors, and the ratio tests the normalisation. The check passes if the correlation
is above 0.999 and the ratio is within 2 % of 1. In our runs the correlation is 1.0000 and
the ratio is between 0.994 and 1.000. If the check reports `FAILED`, the files were not
read correctly or belong to different runs, and the second-order result must not be used.

With a scissor shift (`--scissor`) the comparison is made with the ELK spectrum shifted by
the same amount. It is then only approximate, because the shift also changes the
amplitudes, and it is reported as such. Run once without the shift to make the exact check.

## Scissor shift

The ELK run itself must use `scissor 0`. A scissor shift D can be applied in `exshg`
(`--scissor`, in Hartree). It is applied consistently:

- conduction-band energies are shifted by D in all energy denominators,
- exciton energies are shifted by D (rigid shift; eigenvectors and kernel unchanged),
- the matrix elements r and r;k are computed from the unshifted energies, because they are
  properties of the wavefunctions.

## Approximations

- Excitons enter only in the second-harmonic step (2ω). The first-order step at ω is
  independent-particle, and transitions between excitons are not included.
- Tamm–Dancoff approximation.
- The result depends on the k-mesh and the band window of the BSE run, as any BSE result.

## References

- C. Aversa and J. E. Sipe, Phys. Rev. B **52**, 14636 (1995).
- M. Grüning and C. Attaccalite, Phys. Rev. B **89**, 081102(R) (2014).
