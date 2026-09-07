# Independent referee: third native wedge and thirteen connected-word deletions

**Status: PASS, independently audited analytic proof and FINITE-EXACT certificate.**
No repair to the producer was required. This audit accepts the complete
`(24,18,30)` strict-atlas wedge and its two additional word deletions, relative
to the pinned inherited full-profile and minimum-edge suppliers. It does not
close clock 7200 or identify either remaining word with an unsafe row.

The reviewed primary is `continuing10_20260907_lrc_third_wedge.md`, with
pre-promotion SHA-256
`c96b35df14eb583c4431512f170281ce5dcadc739f43fc1bdf6be349d8f57146`.
The first two wedge theorems and their initial eleven-word topology were
already independently audited; their existing proofs are dependencies here.
Our new numerical paths import no producer implementation.

## 1. Types and physical consumer

The row consists of thirteen distinct positive integer speeds with total gcd
one. Six chosen labels have gcd 7200; the other seven labels have margins
`d_i=gcd(7200,u_i)`. The new wedge theorem assumes an actual two-edge atlas
path with the ordered margins `(24,18,30)`. It does not assume the complement
is connected, the decoder span is exact, or a finite speed box.

For an actual pair `(u,v)=D(p,q)` with coprime `p,q`, put
`e=gcd(7200,D)` and `n=7200/e`. Its danger intersection on every translated
7200-grid is `e` times a corresponding translated n-grid intersection for
the primitive pair. The quotient phase is multiplied by `D/e`; it is not
discarded pointwise. An all-phase minimum is therefore a valid uniform
credit. The pair minima in this proof use strict danger `<1/14`, so their
complement supplies weak safety `>=1/14` even at grid walls.

If a forest has k dangerous vertices at a grid point, its induced dangerous
edge count is at most k-1 for k>0 and zero for k=0. Summing the resulting
pointwise union bound gives a safe-grid lower bound
`7200-sum(single danger counts)+sum(forest pair intersections)`.
Each single danger count is at most `d_i ceil(7200/(7d_i))`.
Thus a forest credit strictly above the full-word excess E gives a safe
grid point for all phases. The inherited proper-six phase supplier supplies
the starting phase, whose 7200 lifts remain safe for all selected six labels.
This accepts the claimed physical weak-safety consequence, not just a local
pair statistic. No new literature or stronger proper-subset supplier is used.

## 2. Complete third-wedge computation and valuation obstruction

The independent engine enumerates all 5855 strict coprime atlas ratios and
both orientations, then tests the two exact margin identities. It recovers
231 arms over the middle margin 18 from margin 24 and 141 from margin 30.
Every one of their 32,571 products is considered. Exactly 32,541 have arm
credit sum greater than 107; the remaining 30 have credits 0 and 66.

For a low product, normalize `(a/b,1,c/d)` to a primitive integer triple R.
Every integral realization has the form hR: the gcd-one coordinates force
its common rational scale to be integral. Its prescribed margin gcd forces
`gcd(7200,h)=6`. For every integer R_i,
`gcd(7200,hR_i)=gcd(7200,6R_i)`, by the prime valuation formula for a gcd.
This proves necessity of the joint test. Choosing h=6 proves sufficiency
for this three-label margin realization, without asserting a thirteen-label
completion. All 30 low triples have distinct coordinates.

Exactly 26 fail this joint test, only at prime 3. We verify every rejected
triple and its actual margins. For example `(9588,8109,7685)` has actual
margins `(72,18,30)` at scale gcd 6 and cannot realize the proposed triple.
The four accepted primitive triples are

```
(3476,21093,445), (5372,19671,415),
(14852,12561,265), (16748,11139,235).
```

Each endpoint pair has quotient clock 1200 and uniform count minimum 16,
hence actual credit 96. Combining that endpoint with the actual 66-credit
arm gives a two-edge tree of credit 162, exceeding 107. It is not a sum of
all three edges of a cycle. The endpoint alone would be insufficient.

The raw-wall referee constructs spatial danger intersections by testing
literal rational midpoints between all walls, then counts strict grid
departures before arrivals. This differs from the producer's intersection
construction. All 376 saved profiles agree in extrema, wall cardinality,
component normalization and literal extremal phases. Each referee minimum
and maximum also has its own independently checked literal grid owner.
The hostile n=7, p=q=1, phase=1/2 has strict count zero but closed count two,
explicitly testing the weak-endpoint convention used by the proof.

We also recompute the conditioned word ceiling directly from the inherited
profile table. The complete 26-value divisor alphabet gives 23,751 four-slot
multiset completions. Literal subset gcd evaluation checks all 126 proper
positional profiles of a retained word. Exactly 162 full words survive;
their maximum excess is 107, attained by `(9,16,18,24,25,30,32)`. The entire
list agrees with the first producer's pinned conditioned table.

## 3. Exact topological quantifiers

Connectedness of the actual complementary seven-label graph is required
only from this point onward. Under hypothetical failure, at each margin-24
vertex either all its 18-neighbor edges are absent or all its 4/16-neighbor
edges are absent. The third theorem additionally forces, at each relevant
margin-18 vertex, absence of all 24-neighbor edges or all 30-neighbor edges.
These disjunctions follow from the presence of the forbidden actual paths;
they are not constraints on a chosen convenient tree alone.

For some simultaneous branch choice the whole actual graph is contained in
the generously allowed graph. Every actual edge credit is bounded below by
the inherited minimum for its two margin classes. An actual spanning tree
therefore has credit at least the **minimum** weighted spanning-tree cost of
the generous graph. If the generous graph is disconnected, no connected
actual graph lies in it. This proves the universal branch deletion rule.
It would be invalid to use a maximum tree assembled from possible edges.

Starting from the pinned complete fifteen-word universe and minimum-edge
table, our verifier reconstructs every positional branch, retaining distinct
vertices even when margins coincide. It enumerates all `7^5=16807` labeled
Pruefer trees and tests containment directly. This uses neither a greedy
minimum-tree procedure nor the producer's cut-condition check.

All 37 branch graphs and minima agree. The new costs at indices 11 and 14
are `[162,198,162,96]`, respectively exceeding E=79 and E=87 in every branch.
All previously deleted words remain deleted. The exact deleted indices are
`1,2,3,4,5,6,7,8,9,10,11,13,14`; the two survivors are

```
(1,9,16,18,24,32,60), E=116, retained tree floor 114;
(5,8,9,30,32,36,48), E=103, retained tree floor 102.
```

These are necessary-condition survivors only. This audit does not assert
an actual ratio assignment, common phase, or unsafe row for either word.
The proof succeeds by retaining the joint valuation of each actual path,
then the actual graph's exhaustive branch owner. A projected arm or a
selected ratio realization alone does not supply those coordinates.

## 4. Reproduction and frozen dependencies

Run the adjacent `continuing10_20260907_lrc_third_wedge_audit.py`, or after
filing:

```
python 04-computation/continuing10_20260907_lrc_third_wedge_audit.py
python -O 04-computation/continuing10_20260907_lrc_third_wedge_audit.py
```

Both modes pass **2605 always-active gates**, with byte-identical raw LF
stdout. The source locates adjacent dependencies first, then the repository
results directory, preserving the source/data filing contract.

```
audit source 2f4edd8baf15e6e7c610e6a229f217d034bf30e8785c73d2738c7a5b770af366
audit output e396dcfc4126acce9e68c0970d837841bec0e9b1b9263b959344bbc969ec63ca
primary source a768f216da83a8882cd2ef2d52abe08f344a971035e0618aa74242d0b623e6fe
primary output f8ebfb3f7d302bdf24274934c80d06d8c663df3bcc52c0feae6e400cfca9e044
primary certificate 60a6de57a623dbe64a23a77113f140b312f5496247c2253fbe695872740039e9
first-wedges certificate 890d00d44f0195765d62fe1d40b59ad102311f34a84a42a5fddc229d037209e9
full-profile table 935f3f687b6d7c89cc099e536f536238fd753bcc4c1747906d213cef387ca93f
minimum-tree certificate 580a7c930103aab3bea867ad463a90b0e0208323a90ee95a685ff811a761582d
```

The source pins the exact mathematical source/output/JSON dependencies;
the pre-promotion prose hash is recorded above without preventing a later
honest audit-status promotion. All producer files remained untouched.
