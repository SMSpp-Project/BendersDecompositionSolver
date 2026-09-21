# test

A tester for `BendersDecompositionSolver`.

It writes the instance it runs on, so that the module knows of no Block but
the ones it builds: a small two-stage linear program, in two equivalent
models. One is monolithic, the first-stage `Variable` and the second-stage
ones in a single `AbstractBlock`, which a `:MILPSolver` solves to give the
reference optimum; the other is structured, a master `Block` carrying the
first-stage `Variable` and one nested sub-`Block` per scenario carrying the
second-stage ones, whose `Constraint` couple the master, which
`BendersDecompositionSolver` solves. What is checked is that the two optima
agree within a relative tolerance, whichever way the Solver is asked to work.

Every way it has of writing a cut is run on that instance, one
`BlockSolverConfig` each: the multi-cut and the single-cut master, the
cutting-plane and the bundle one, the Pareto-optimal cut, the two
normalizations of the feasibility cut and its phase one, the combinatorial cut
of a binary master, the cut that serves for both feasibility and optimality,
what keeping a subproblem in the master costs, and what the Block is given
back as at the end of a `compute()`. The instances that make each of them
matter are built on purpose, i.e., with the subproblems feasible everywhere
for the optimality cuts and infeasible at some first-stage point for the
feasibility ones, and in the nested shape a subproblem has when it is a model
of its own rather than a bare set of rows.

The comparisons that are posed on the Block of another module are not here:
they live in the test suite of that Block, where its instances and its
batteries are.


## Authors

- **Antonio Frangioni**  
  Dipartimento di Informatica  
  Università di Pisa

- **Donato Meoli**  
  Dipartimento di Informatica  
  Università di Pisa


## License

This code is provided free of charge under the [GNU Lesser General Public
License version 3.0](https://opensource.org/licenses/lgpl-3.0.html) -
see the [LICENSE](LICENSE) file for details.
