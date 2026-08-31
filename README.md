# BendersDecompositionSolver

Definition and implementation of the `BendersDecompositionSolver` class, which
implements the `CDASolver` interface within the SMS++ framework to automate the
Benders decomposition of a block-structured problem.

Attached to a suitably structured Block (B), `BendersDecompositionSolver` keeps
the "complicating" (first-stage) `Variable` of (B) in a master problem and
projects out the `Variable` of each sub-`Block`, representing the resulting
value function as a `BendersBFunction`. The Benders optimality and feasibility
cuts are produced by the `BendersBFunction` through the standard `C05Function`
interface, so the solver concentrates on the automation layer: detecting the
complicating `Variable`, building one `BendersBFunction` per sub-`Block`,
assembling the master problem and driving the cut loop.

The Block (B) is assumed to have the following structure:

- the root (B) plays the role of the master: its `ColVariable` are the
  complicating `Variable`, and its `Objective` and `Constraint` are the
  first-stage objective and feasible region;

- each sub-`Block` is a subproblem, with its own `Variable`, `Objective` and
  `Constraint`, coupled to the master only through the appearance of (some of)
  the complicating `Variable` in (some of) its `Constraint`.

Two regimes are supported for solving the master problem, selected by a
parameter:

- a convex regime, in which the value functions enter the master objective as
  `C05Function` and the master is solved by a bundle-type `CDASolver`, which
  performs the (stabilized) cutting-plane loop internally;

- a MILP regime, in which the master is a MILP with epigraph `Variable` and the
  Benders cuts are added as dynamic `Constraint`, the solver driving the outer
  cut loop until no violated cut remains.

The master and the subproblem `Solver` are instantiated through the `Solver`
factory and configured via `Configuration`, so that how the master and the
subproblems are solved is not hard-wired. See the class documentation for the
full description of the assumptions, parameters and solution recovery.

Both kinds of cut are duals of the subproblem, so the `Solver` given to it has
to be asked for a *vertex* solution: an interior-point one has optimal but
non-basic duals, giving a valid yet weaker optimality cut, and it proves
infeasibility without producing the Farkas certificate that the feasibility cut
is. The certificate only exists if the infeasibility is proved by the simplex,
hence the presolve of the subproblem `Solver` has to be switched off as well;
if it is missing, an exception is thrown rather than silently converging to a
wrong optimum. A subproblem that is always feasible, e.g. because the
constraints the master can make unsatisfiable carry a slack with a large cost,
needs no feasibility cut in the first place and is therefore the robust choice
whenever the model allows it.

The method is described in

W. van Ackooij, A. Frangioni, W. de Oliveira "Inexact Stabilized Benders'
Decomposition Approaches, with Application to Chance-Constrained Problems with
Finite Support" *Computational Optimization and Applications* 65(3), 637-669,
2016

and in

D. Baena, J. Castro, A. Frangioni "Stabilized Benders Methods for Large-scale
Combinatorial Optimization, with Application to Data Privacy" *Management
Science* 66(7), 3051-3068, 2020


## Getting started

These instructions will let you build `BendersDecompositionSolver`.


### Requirements

- [SMS++ core library](https://gitlab.com/smspp/smspp)

It's not a build requirement, but you will need a SMS++ `Solver` capable of
solving the master problem, such as
[BundleSolver](https://gitlab.com/smspp/bundlesolver) for the convex regime or
[MILPSolver](https://gitlab.com/smspp/milpsolver) for the MILP regime.


### Build and install with CMake

Configure and build the library with:

```sh
mkdir build
cd build
cmake ..
cmake --build .
```

The library has the same configuration options of
[SMS++](https://gitlab.com/smspp/smspp-project/-/wikis/Customize-the-configuration).

Optionally, install the library in the system with:

```sh
cmake --install .
```

### Usage with CMake

After the library is built, you can use it in your CMake project with:

```cmake
find_package(BendersDecompositionSolver)
target_link_libraries(<my_target> SMS++::BendersDecompositionSolver)
```

### Build and install with makefiles

Carefully hand-crafted makefiles have also been developed for those unwilling
to use CMake. Makefiles build the executable in-source (in the same directory
tree where the code is) as opposed to out-of-source (in the copy of the
directory tree constructed in the build/ folder) and therefore it is more
convenient when having to recompile often, such as when developing/debugging
a new module, as opposed to the compile-and-forget usage envisioned by CMake.

Each executable using `BendersDecompositionSolver` has to include a "main
makefile" of the module, which typically is either [makefile-c](makefile-c)
including all necessary libraries comprised the "core SMS++" one, or
[makefile-s](makefile-s) including all necessary libraries but not the "core
SMS++" one (for the common case in which this is used together with other
modules that already include them). The makefiles in turn recursively include
all the required other makefiles, hence one should only need to edit the "main
makefile" for compilation type (C++ compiler and its options) and it all should
be good to go. In case some of the external libraries are not at their default
location, it should only be necessary to create the `../extlib/makefile-paths`
out of the `extlib/makefile-default-paths-*` for your OS `*` and edit the
relevant bits (commenting out all the rest).

Check the [SMS++ installation wiki](https://gitlab.com/smspp/smspp-project/-/wikis/Customize-the-configuration#location-of-required-libraries)
for further details.


## Getting help

If you need support, you want to submit bugs or propose a new feature, you can
[open a new issue](https://gitlab.com/smspp/bendersdecompositionsolver/-/issues/new).


## Contributing

Please read [CONTRIBUTING.md](CONTRIBUTING.md) for details on our code of
conduct, and the process for submitting merge requests to us.


## Authors

### Current Lead Authors

- **Antonio Frangioni**  
  Dipartimento di Informatica  
  Università di Pisa

- **Donato Meoli**  
  Dipartimento di Informatica  
  Università di Pisa

### Contributors


## License

This code is provided free of charge under the [GNU Lesser General Public
License version 3.0](https://opensource.org/licenses/lgpl-3.0.html) -
see the [LICENSE](LICENSE) file for details.

## Disclaimer

The code is currently provided free of charge under an open-source license.
As such, it is provided "*as is*", without any explicit or implicit warranty
that it will properly behave or it will suit your needs. The Authors of
the code cannot be considered liable, either directly or indirectly, for
any damage or loss that anybody could suffer for having used it. More
details about the non-warranty attached to this code are available in the
license description file.
