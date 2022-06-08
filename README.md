## Lovász-Scarf reduction
Implementation of the [lattice basis reduction](https://www.jstor.org/stable/3689761)  of Lovász and Scarf, following [Cook et al.](https://pubsonline.informs.org/doi/abs/10.1287/ijoc.5.2.206?casa_token=oKwfy2a-aZAAAAAA%3AQLvaYZRbumvZAv6_87WR9poBZ7n0W4yoYRu7Xuz5focL9tX48gb-BcHtLtOXujBmkQwneywLbQ&journalCode=ojoc).

### Requirements
[Gurobi](https://www.gurobi.com/)

### Input
The input file must be a ```.txt``` file in the following form:
```shell
n m
A
b
lb
ub
```
where n is the number of variables, m is the number of constraints, A is the constraint matrix (m x n), b is the right-hand side (m), and lb and ub are the vectors containing the variable lower and upper bounds, respectively, (n).

### Instructions
This repository provides a shell script to compile and run the program. Modify the filename in this script (```run.sh```).
