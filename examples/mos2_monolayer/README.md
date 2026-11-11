# Monolayer MoS₂

ELK input and reference result for χ⁽²⁾_yyy of free-standing monolayer MoS₂
(PBE with spin–orbit coupling, 12 × 12 k-mesh, 4 valence and 4 conduction bands in the BSE).

```bash
# set sppath in elk.in to your ELK species directory
mpirun -np 20 elk > elk.out
exshg . 222 --txt
python plot.py XSHG_BSE_222.txt my_result.png
```

Expected output of `exshg`:

```
linear check Im eps_22 vs ELK EPSILON_BSE: corr 1.0000, median ratio 1.0000
```

Reference result: `reference/XSHG_BSE_222.txt` and `mos2_shg.png`.

![chi2 of monolayer MoS2](mos2_shg.png)

The step-by-step description is in the [tutorial](../../docs/tutorial.md).
