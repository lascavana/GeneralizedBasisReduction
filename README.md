## Generalized Basis Reduction
Implementation of the [generalized lattice basis reduction](https://www.jstor.org/stable/3689761)  of Lovász and Scarf, following [Cook et al.](https://pubsonline.informs.org/doi/abs/10.1287/ijoc.5.2.206?casa_token=oKwfy2a-aZAAAAAA%3AQLvaYZRbumvZAv6_87WR9poBZ7n0W4yoYRu7Xuz5focL9tX48gb-BcHtLtOXujBmkQwneywLbQ&journalCode=ojoc).

### Requirements
[Gurobi](https://www.gurobi.com/)

### Input
The input file must be an ```.lp``` or ```.mps``` file.

### Output
A basis of $\mathbb{Z}^n$ that is reduced for the given polytope.

### Instructions
This repository provides a shell script to compile and run the program. Modify the filename in this script (```run.sh```) and then run
```
sh run.sh
```

