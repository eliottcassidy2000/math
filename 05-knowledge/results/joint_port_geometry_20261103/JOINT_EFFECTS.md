# Joint attachments enlarge the blocks, but preserve an effect cube

**PROVED for the specified linear observation language, with exact rational
matrix certificates and actual-source rank witnesses.** Local labels JE1--5.
Same-agent analytic review; independent graph replay and rank constructions.
This is not a smallest arbitrary encoding, a full tournament invariant, or
closure under every native complement and every re-marking.

## JE1. The language, cap, and physical source space

Keep ordered boundary ports (r,b). A path cover partitions all vertices into
unordered directed paths. Its skeleton retains which live port starts or
ends each path; an unmarked endpoint is written -1. There are eighteen
two-port skeletons and an unbounded number m of wholly unmarked paths.
Retain states with at most two unmarked terminal ends:

`m + #{(a,b) in skeleton : b=-1} <= 2`.

This gives46 states. A birth of polarity0 beats every unmarked old vertex;
its incoming arcs from live ports are specified by a mask. Such a birth can
consume an old terminal endpoint only if that endpoint is marked. Forgetting
a mark can increase the number of unmarked terminal ends, but cannot reduce
it. Thus discarded states cannot later contribute to the four queries
H, s_r, E_r(1), S_r(2). The cap is exact for every word in this language.
Opposite-polarity births do not satisfy this argument.

Generators are:

| Name | Actual operation | Permanent births | Native moves |
|---|---|---:|---:|
| N | NR: birth x beaten only by r, then y beaten only by x; output (x,b) |2|0|
| C3,C5 | inherited native oldest-child cycle, with replacement working port |6,8|4,8|
| J_m, m=0,1,2,3 | one polarity0 birth with incoming mask m on (r,b), forget that newborn as a port, then N |3|0|

All original vertices survive these generators. Only N and J_m preserve b
as a live boundary mark. C3/C5 replace that mark at their first birth and
finish with new working marks. [literal_interfaces.py](literal_interfaces.py)
checks all four J_m and both native cycles against independent full graph
observers. J0 may have a weak intermediate graph; its sealed endpoint is
strong. Temporary cycle helpers are separate from permanent birth counts.

Initial physical sources are NR-sealed nonempty strong tournaments, with
the second output mark either another old vertex, the old root, or the
latest newborn. Take the Q-linear span of their capped cover vectors. The
three exact sealing matrices, with the impossible empty-input row removed,
give a31-dimensional invariant upper bound B. The31 actual sources saved
in [rank_transfer.json](rank_transfer.json), of orders5--7, are independent
modulo1009 and hence over Q. Therefore they span B exactly.

**Ghost boundary.** If the empty-input coordinate is allowed, the joint seal
image has dimension32. Its extra output is the formal two-singleton-port
cover of a graph with no other vertices. No nonempty input realizes it.
Its removal is a source-domain condition, not an optional simplification.

## JE2. Minimal future-query dimensions

Start with H,s,E(1),S(2) and close their columns under left multiplication
by every generator. Exact elimination terminates; the resulting spaces are
invariant. Restriction to B gives these ranks:

| Extra joint births, besides C3,C5,N | Minimum Q-linear observer dimension |
|---|---:|
| J0 only |6|
| J1 only |6|
| J2 only |28|
| J3 only |30|
| all J0,J1,J2,J3 |31|

Each lower bound is the rank of actual source responses to listed future
words and queries, modulo1009. Each matching upper bound is an exact
rational invariant basis. Consequently these are all-word minima among
linear observers of actual sources, not merely ranks of formal variables.
The matrix files include the words furnishing each query column.

The asymmetric jump J1 versus J2 is rooted: the mask bit0 refers to the
root, bit1 to the spare. They are not interchangeable while the queried
root is fixed. Joint attachments detect information absent from the old
six-dimensional rooted count/golden observer.

## JE3. Exact algebra and complete invariant-quotient classification

On B, filter by support on at least two, at least one, or at least zero
unmarked terminal ends. The dimensions are14,26,31, giving graded sizes
14,12,5. The generators preserve this filtration.

[filtration.py](filtration.py) constructs an adapted rational basis. Two
common Sylvester systems, solved in [split_algebra.py](split_algebra.py),
split the off-diagonal extensions simultaneously for every generator.
This is an explicit exact conjugation to three diagonal blocks, not a
deduction from the dimensions alone.

The diagonal matrix words contain365 independent elements modulo1009.
All displayed rational denominators are invertible at1009. Their Q-rank is
therefore at least365, while the block upper bound is
14^2+12^2+5^2=365. Hence the full unital word algebra is exactly

`M14(Q) direct-sum M12(Q) direct-sum M5(Q)`.                 (JE)

Its action is the direct sum of the three natural modules, each once.
Central block idempotents split every invariant subspace; the only
submodules of a natural full-matrix module are0 and the whole block.
There are therefore exactly eight invariant subspaces and eight invariant
linear quotients, with quotient dimensions

`0,5,12,14,17,19,26,31`.

Only the full31-dimensional quotient preserves every required future query
in JE1. Each of the seven smaller invariant quotients forgets some of those
observations; invariance alone does not make it a sufficient carrier.

This extends the previous M2(Q) direct-sum Q direct-sum M3(Q) effect cube:
the sizes change, but the Boolean lattice of three independently retained
blocks survives. The identification is a lattice isomorphism of block
choices, not an identification of the physical sources or their algebras.

## JE4. A genuinely persistent second port

To keep the original spare through C3/C5, three live working ports suffice:
keep that spare while the other two perform the oldest-child program.
Its virtual insertion history has the same final graph as the native
program, with the proved native oldest-child permutation. The persistent
spare remains an extra supplied mark throughout. This is a different
operation language from replacing it and freely recovering its name later.

Compile its46-state boundary operator with the three-port kernels in
[persistent_port.py](persistent_port.py). Start only with the N seal keeping
an old spare distinct from the old root. Its invariant physical span and
minimal future observer both have dimension27;27 actual independent source
witnesses establish the lower bound. The J_m generators remain allowed.
The corresponding unital algebra has dimension95 modulo1009. Its exact
rational algebra classification is OPEN;95 is not silently promoted to a
rational dimension or a decomposition theorem.

The earlier phase report's sentence claiming C3/C5 already preserved the
second live port has been corrected in place and in the mistakes ledger.
The old matrices and clock orders are unchanged. A retained vertex and a
retained port are different resources.

## JE5. Recurrences and exact weighted program series

The complete31-dimensional state, for any supplied source, obeys a uniform
order-seven recurrence under repeated C3. Its minimum common annihilator is

`x(x-4)(x^2-61x+18)(x^3-472x^2+3936x-3072)`.

For C5 it is

`x(x-16)(x^2-421x+66)(x^3-4984x^2+147744x-431616)`.

These are exact matrix identities, established blockwise and checked on the
full matrices in [recurrences.py](recurrences.py). That receipt also gives
minimum polynomials for N and all J_m. J2/J3 have repeated nonzero factors:
the whole algebra is semisimple, but individual generators can have Jordan
blocks. These are different assertions.

Let A_L be the31-state matrix of letter L, v an actual source row, and q a
query column. With z counting permanent births and t counting native moves,
the exact generating function summing query outcomes over all supplied
program words is

`v [I-z^2 A_N-z^6 t^4 A_C3-z^8 t^8 A_C5-z^3 sum_m A_Jm]^-1 q`. (GF)

Its constant term is the empty word. Expanding the geometric matrix series
proves (GF), and all entries are rational functions. Letter-specific weights
can retain full word multiplicities. This enumerates programs with weighted
outcomes, not isomorphism classes: multiple histories can produce the same
tournament. Deduplication, unrestricted complement/re-marking, and a full
A000568 reconstruction remain additional problems.
