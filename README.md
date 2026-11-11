# exshg

**Second-harmonic generation with excitons from ELK Bethe–Salpeter calculations**

[![tests](https://github.com/yucj-520/exshg/actions/workflows/tests.yml/badge.svg)](https://github.com/yucj-520/exshg/actions/workflows/tests.yml)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE)
![Python 3.9+](https://img.shields.io/badge/python-3.9%2B-blue.svg)

`exshg` computes the second-order optical susceptibility χ⁽²⁾(2ω; ω, ω) of a crystal
**including excitons**. It is a post-processing tool for the all-electron code
[ELK](https://elk.sourceforge.io):

- ELK solves the Bethe–Salpeter equation (BSE) for the *linear* optical response.
- ELK computes second-harmonic generation (SHG) only for *independent particles*.
- `exshg` reads the output of a finished ELK BSE run and returns the SHG spectrum with
  the electron–hole interaction included at the second-harmonic frequency.

The ELK source is not modified, and no additional ELK run is needed.

![chi2 of monolayer MoS2 with and without excitons](examples/mos2_monolayer/mos2_shg.png)

*Monolayer MoS₂: the excitons move the SHG resonances to the exciton energies and change
their strength ([tutorial](docs/tutorial.md)).*

## Installation

```bash
pip install git+https://github.com/yucj-520/exshg.git
```

Requirements: Python 3.9 or newer and NumPy. ELK 11.2.3 (built with gfortran) is needed
to produce the input data. See [Installation](docs/install.md).

## Quick start

1. Run ELK with the ground state, the momentum matrix elements and the BSE
   (`elk.in`):

   ```
   tasks
     0
     120
     121
     180
     185
     186
     187

   scissor
     0.0

   reducek
     0
   ```

2. Run `exshg` in the directory of that run:

   ```bash
   exshg . 222 --txt
   ```

   `222` selects the component χ_yyy.

3. Read the result: `XSHG_BSE_222.txt` contains the photon energy and χ⁽²⁾ with and
   without excitons in pm/V.

A complete example with all files and the expected output is in the
[tutorial](docs/tutorial.md).

## Documentation

| Page | Content |
|---|---|
| [Installation](docs/install.md) | `pip install`, requirements, ELK |
| [Tutorial: monolayer MoS₂](docs/tutorial.md) | the full workflow, step by step, with the expected output |
| [Example: MoS₂/WSe₂ tetralayer](docs/tetralayer.md) | a 12-atom heterostructure with spin–orbit coupling |
| [Method](docs/method.md) | equations, approximations and the built-in check |
| [Reference](docs/reference.md) | command line, Python function, output files, units, troubleshooting |

## Examples

| Example | System | Files |
|---|---|---|
| [`examples/mos2_monolayer`](examples/mos2_monolayer) | monolayer MoS₂, 12 × 12 k-mesh | input, reference spectrum, plot |
| [`examples/mos2_wse2_tetralayer`](examples/mos2_wse2_tetralayer) | MoS₂/WSe₂ tetralayer, H-like and R-like stacking, 6 × 6 k-mesh | inputs, reference spectra, plot |

## Method in brief

The independent-particle χ⁽²⁾ is obtained from the length-gauge density matrix to second
order (Aversa and Sipe). On the transitions included in the BSE, the propagator of the
second-harmonic step, 1/(ω_cv − 2ω), is replaced by the exciton propagator built from the
BSE eigenvalues and eigenvectors of ELK. All other transitions keep their
independent-particle form. Details are in [Method](docs/method.md).

**Built-in check.** The same replacement at first order must reproduce the BSE dielectric
function that ELK writes itself. Every run prints this comparison, so a wrongly read file
or an inconsistent run is detected immediately.

**Approximations.** Excitons enter only at 2ω; the step at ω is independent-particle and
transitions between excitons are not included. The BSE is solved in the Tamm–Dancoff
approximation.

## Validation

- The independent-particle part agrees with ELK task 125 at scissor 0 to about 1 %.
- The first-order replacement reproduces ELK's `EPSILON_BSE` (correlation 1.0000).
- Monolayer MoS₂ without spin–orbit coupling, 21 × 21 k-mesh, scissor 0.72 eV: the main
  peak is at ħω = 1.45 eV with 530 pm/V (thickness 6.15 Å). The real-time calculation of
  Grüning and Attaccalite, Phys. Rev. B **89**, 081102(R) (2014), gives 1.56 eV and about
  650 pm/V. Their method also includes the excitons at ω.

## Citation

If you use `exshg`, please cite this repository (see [`CITATION.cff`](CITATION.cff)).
A publication describing the method and its first application is in preparation.

## Licence and author

`exshg` is distributed under the GNU General Public License, version 3 or later.

Author: Cuiju Yu, Aalto University. Questions and bug reports are welcome in the
[issue tracker](https://github.com/yucj-520/exshg/issues).

Parts of the code were developed with the assistance of an AI coding tool (Claude). All
results were validated against ELK as described above.
