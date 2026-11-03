# Finite joint responses from an acyclic variable interface

**PROVED specialization of the inherited substitution kernel.**

## R1. A finite polynomial for all positive module sizes

Fix a tournament Q, f individually named fixed vertices, and d other vertices
whose induced tournament is transitive. Replace the latter by nonempty
transitive modules of sizes k1,...,kd. Fix q ordered, possibly empty path bins,
and any simultaneous start/end requirements at the fixed named vertices.
Then the number of spanning covers obeying these requirements is

    sum_r c_r product_i r_i! S(ki,r_i),
    0 <= r_i <= f+q,   c_r nonnegative integers.             (R1)

The coefficients depend on the complete quotient relation, bin order and
joint endpoint requirements. They are not products of endpoint marginals.

Proof. Split each path at every change of module and at every fixed vertex.
Each resulting variable run is increasing inside its transitive module.
Between successive fixed vertices or bin boundaries, the variable-only word
is a path in an acyclic tournament; it visits each variable module at most
once. There are f+q such intervals in total. Hence r_i<=f+q. A legal compressed
word with r_i occurrences of module i lifts in exactly r_i! S(ki,r_i) ways:
assign its labelled vertices surjectively to the ordered occurrences, then
sort each occurrence. This assignment retains every endpoint and component
incidence. Compression and lifting are inverse. Empty-bin conventions are
handled by the separators. QED.

`run_polynomial.py` implements the finite word recursion with all joint guards
checked at bin boundaries. It rejects cyclic variable quotients. In that case
the general August rational walk kernel still applies, but the bounded-degree
conclusion above need not: C3 already permits arbitrarily many variable-only
cycles. Splitting a distinguished member out of a module gives another fixed
port when its incidences are supplied; silently re-marking an anonymous member
is not an operation of this compressed observer.

For one variable block and scalar H, with e fixed vertices, c_j<=e!*binom(e+1,j).
The functions j!S(k,j), j>=1, are linearly independent on positive k by their
triangular expansion in 1^k,...,j^k. Thus equality of finite coefficient vectors
is an all-size identity. No finite extrapolation hypothesis is needed.

## R2. Complete two-ear boundary classification

Use the [prior word atlas](../ready_return_pentagon_20261103/ATLAS.md):
ports a->b, root r with in-neighbours a,b, residual chain word 11^k w, k>=1.
All fixed vertices (a,b,r and every letter of w) may be distinguished separately
by each new vertex. Both births must be uniform on the initial 11 block.
The edge between the two newborns is arbitrary. This includes every literal
two-ear attachment in this language, including those using more than the
old three ports. It permits different choices at different k.

All sixteen neutral C4 tails and all thirty-one C6 chart tails were tested.
The C4 list is the prior all-size classification of neutral defect-four roots;
the C6 list is only the transitive-residual chart, not all defect-six roots.
`boundary.json` gives every accepted ordered mask and its exact coefficients.

| Tail | C | Uniform strong 3H+2 programs | Uniform 3H+4 |
|---|---:|---:|---:|
| 00,01 | 4 | 6 | 0 |
| 01,11 | 4 | 2 | 0 |
| 01,11,11 | 6 | 2 | 0 |
| Every other tail | 4 or6 | 0 | 0 |

There are also four weak +2 programs on tail 10,10,01. They are retained in
the data, not counted as strong returns. The six programs in the first row
all return to the actual word 11^k,00,00,00,01, up to their displayed newborn
ordering; its H is 27*2^k-1. This extends 9*2^k-1 by a uniform +2 step.
The other strong rows append 00,01 to the old word. None of these conclusions
says the resulting source is ready for both subsequent branches.

Exhaustiveness: at k=1 enumerate every two-birth incidence mask. Any all-size
identity must occur there. For each candidate compute its complete R1 run
polynomial and compare coefficients with 3H+b. This is a finite certificate
for the entire operation alphabet, not a search only over successful samples.
Strongness for all k follows from strongness of its one-vertex quotient and
nonempty transitive substitution. Converse gives the transported identical
classification with all incidences reversed.

## R3. State dependence cannot conceal an eventual +4 return

For each two-ear program the discrepancy Hnew(k)-3Hold(k)-4 is an integer
linear combination of j^k, j<=e+1, where e<=9 in this atlas. The sum of the
absolute coefficients is at most

    B=9!*3^10 + 3*7!*3^8 +4 = 21526903444.                 (R2)

This follows from the single-block coefficient bound above and the binomial
expansion of j!S(k,j). If J>=2 is its highest nonzero base, its term dominates
all others whenever (J/(J-1))^k>B. For J=1 the discrepancy is a nonzero constant.
The exact integer inequality 10^226>B*9^226 proves domination for k>=226.
R2's finite classification excludes an identically zero +4 discrepancy for
every program. Therefore **no selection from this entire alphabet**, even one
depending on k, realizes 3H+4 at k>=226. This is not a restriction on ears that
distinguish different vertices inside the long block, on larger gadgets, or
on complement moves that change old edges.

## R4. What the birth evaluator retains

The separate three-old-port search has 512 ordered two-birth commands and a
165-coordinate finite-horizon input. Only covers with at most three components
can contribute: each birth joins at most two old components and reduces their
number by at most one. The named-birth lift retains the new port as well as
all old ports. Its anonymous-component action is a zeroth-plus-first derivative,
because a newborn uniform on anonymous vertices cannot both receive from and
send to an anonymous endpoint. `dynamics.py` checks its nonnegative derivative
and compares it with the independent directed-cut evaluator.

The finite horizon is sufficient for these two births. It is not a minimal
carrier for arbitrarily many births, arbitrary joint tensors, or complements.
