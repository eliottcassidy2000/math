# Session 324: prescribed targets, Busch inverse cuts, and the seven–Fermat rail

**Chan–Pak remains open.** This packet proves an infinite strong target family, identifies a restricted but exact inverse recursion, and separates arithmetic propagation of a gap from the geometry needed to create or eliminate one.

## Main construction

[The seven-times-dyadic rail](FERMAT_SEVEN_RAIL.md) strongly realizes every 7(2^s+1) for s≥2. At s=1 it would be 21, which is impossible. The infinite tail is a single repeatable physical cut with output 3584·2^m+7, m≥1; the residual 7 is a signed boundary correction 29·13−5·74, not a reservoir count. Eight literal initial witnesses bridge to that tail. It follows in particular that both F_r and 7F_r occur for every Fermat number index r≥1, whether or not F_r is prime. [The boundary analysis](RAIL_BOUNDARY.md) proves why this particular target curve cannot begin at any of the first three phase shifts under a single growing transitive module. [Small companion rigidity](SMALL_FERMAT_RIGIDITY.md) rules out all nontrivial automorphisms at 35, 119 and 1799; larger Fermat companion counts can only have automorphism group order 1 or 7.

## Inverting rather than forecasting

[Ordered Busch overlap](BUSCH_INVERSE.md) gives the exact surplus over H(A)H(B), with a complete classification of the zero and unit defects. Attaching a nearly transitive block yields N=(2^j+1)H+2^jD. The prescribed N therefore determines a finite list of smaller source counts and exact requested defects. A zero-defect root has a unique exit and a reducible deletion. A unit-defect root has a unique exit to a vertex that ends exactly one Hamilton path of its strong deletion. The latter points to [a complete Fibonacci grammar](UNIT_END_GRAMMAR.md) for rooted unit-end sources, with exact prefix transfer and a smaller two-sided-unit subgrammar. The naive three-digit version already fails at H=11: no zero-defect source supplies 33 by the small ternary digit, even though 33 occurs.

The right root on the actual count-63 seed has defect zero. Moving to the new terminal root after each two-vertex step gives strong counts 63·3^k indefinitely. The one-vertex closure tree is different: [the gap packet](GAP_STRUCTURE.md) gives exactly three consecutive one-vertex tripling moves, and no fourth from any strong 63 seed. [The commutative inverse](COMMUTATIVE_INVERSE.md) removes factorization from the earlier two-module solver, gives literal halving from either side, and bounds the reachable language rather than assuming it is dense.

## How extra holes could behave

[Gap structure](GAP_STRUCTURE.md) distinguishes the inherited prime/7p frontier reduction from new propagation constraints. A sole extra gap p would force infinitely many *strong witnesses* at other targets pq, not infinitely many gaps. A sole extra gap 7p would force strong 7p² and 21p. When there are no additional missing primes and only finitely many candidate 7p axes, an infinite exception set is equivalent to a wholly missing ray 7p^j on one of those axes. More generally a finite residual prime-exponent cage gives necessary and sufficient tail and axis tests for cofiniteness. Multiplicative closure alone permits more than one finite extra hole. [An S-unit obstruction](ACYCLIC_CLOCKS.md) shows why a fixed quotient with independent transitive lengths at acyclic marks cannot hit infinitely many odd targets of bounded prime support, even asynchronously.

## Busch, Neményi, Fermat

[The mechanics comparison](BUSCH_NEMENYI_FERMAT.md) makes the inverse-fiber analogy precise: for a Cartesian stress net the five free parameters are real and compatible with the equilibrium equations; four ordered Fermat-count SCCs produce five discrete neutral gap coordinates. Neither dimension statement supplies native endpoint readiness. Busch's size floor and the multiplicative order of 2 show that **all** Fermat-count realizers with r≥1 are rigid, including composite Fermat numbers. [The fixed-count class series](FIXED_COUNT_FIBER.md) gives the exact growth polynomial and converse quasipolynomial for each attained h from its finitely many strong skeletons.

## Reproduce

Run each `replay_*.py` here with Python 3. Replays test literal masks, the ordered seam and both root choices, the unit-end decoder through order five, the phase-shift boundary, the circulant orbit profiles, and the commutative inverse and deletion laws. These are controls for the proved statements, not extrapolations to the full spectrum.
