---
layout: default
title: "Installation"
nav_order: 2
---

# Installation

## exshg

```bash
pip install git+https://github.com/yucj-520/exshg.git
```

or, from a copy of the repository,

```bash
git clone https://github.com/yucj-520/exshg.git
cd exshg
pip install .
```

Requirements: Python 3.10 or newer and NumPy. Matplotlib is needed only for the example
plots (`pip install "exshg[plot] @ git+https://github.com/yucj-520/exshg.git"`).

Check the installation:

```bash
exshg --help
```

## ELK

`exshg` reads the output of ELK. It was developed and tested with **ELK 11.2.3** compiled
with **gfortran**. The ELK source is not modified.

1. Download ELK from <https://elk.sourceforge.io> and build it as described in its manual
   (`make` after choosing the compiler options in `make.inc`). MPI and OpenMP builds both
   work.
2. Note the path of the `species` directory of your ELK copy. It goes into the `sppath`
   block of every input file.

Two points matter for `exshg`:

- **Compiler.** `EVBSE.OUT` is a Fortran unformatted file. `exshg` reads the record layout
  written by gfortran. A file written by another compiler with a different record-marker
  convention will not be read correctly.
- **Version.** The file formats of `PMAT.OUT` and `EVBSE.OUT` are those of ELK 11.2.3. Other
  versions have not been tested; the built-in check (see [Method](method.md)) tells you
  whether the files were read correctly.
