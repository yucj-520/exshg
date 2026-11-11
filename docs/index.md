# exshg documentation

`exshg` computes the second-order susceptibility χ⁽²⁾(2ω; ω, ω) with excitons from a
finished [ELK](https://elk.sourceforge.io) Bethe–Salpeter (BSE) run. It is a
post-processing tool: ELK is used unchanged.

| Page | Content |
|---|---|
| [Installation](install.md) | `pip install`, requirements, ELK |
| [Tutorial: monolayer MoS₂](tutorial.md) | the full workflow, step by step, with the expected output |
| [Example: MoS₂/WSe₂ tetralayer](tetralayer.md) | a 12-atom heterostructure with spin–orbit coupling |
| [Method](method.md) | equations, approximations and the built-in check |
| [Reference](reference.md) | command line, Python function, output files, units, troubleshooting |

Quick start:

```bash
pip install git+https://github.com/yucj-520/exshg.git
exshg path/to/elk_bse_run 222 --txt
```

Source code and examples: <https://github.com/yucj-520/exshg>
