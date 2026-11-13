# The H9 rearming flip has a sharp inflation boundary

**PROVED for arbitrary nonempty homogeneous modules; FINITE-EXACT
controls.** Common substitution preserves both the same H and repaired
unit readiness only when all three substituted modules are singleton.
The repaired side nevertheless belongs to
an explicit all-parameter unit-ready family, with an actual primitive NR
exit. This is a marked use of the inherited two-width cyclic family,
not a new unmarked tournament family or an all-target repair theorem.

## 1. Physical source and the preserved root star

Use the history lane's two actual H9 sources

    T=(22,8,2,5,14),      U=(22,8,2,21,6),

with physical root0. They differ only on the arc43: T has4→3, U has3→4.
Both have the same unrooted slot polynomial. Their rooted profiles are
`(9,3,3,9)` and `(9,4,5,10)` respectively, so the flip changes the
indegree-one readiness defect q−s from0 to1 while retaining the root's
unique incoming neighbor3. The birth-return defect `F−h−q+s` is zero
on both sides; these are different notions of readiness.

Inflate vertices1,2,4 to homogeneous transitive blocks of sizes a,b,c≥1.
Keep vertices0 and3 singleton. Let the blocks be A,B,C respectively,
and write

    A0=2^a, B0=2^b, C0=2^c.

The root's unique incoming neighbor is still the same physical singleton3.
The changed cut is the entire set of c arcs between3 and block C; all
other arcs agree between the two inflated sources. Every source is strong,
because its original strong quotient is retained and every block is
nonempty. The full graph has a+b+c+2 vertices.

## 2. Exact rooted responses, without interpolation

Delete root0 and call vertex3 u. On the T side the deletion graph is

    TT_c → V,       V=C3[TT_a,{u},TT_b],
    TT_a → u → TT_b → TT_a.

On the U side it is

    V'=C3[TT_a,{u},TT_(b+c)].

Here the old C and B blocks combine into one transitive module ordered
C then B, because they now have identical external incidence.

For any positive a,e, consider P=`C3[TT_a,{u},TT_e]`, with A0=2^a,
E0=2^e. Its exact rooted state at u is

    ((A0−1)(E0−1)+2, 1, 1, 1+E0(A0−1)).            (14)

The full count is the inherited two-width cyclic formula, also proved in
[Session298 TWO_WIDTH](../../session298_exact_target_geometry_20261113/TWO_WIDTH.md).
A path starting at u must use all of TT_e and then all of TT_a, so there
is exactly one. Deleting u leaves a transitive ordinal sum, also with one
Hamilton path. For F, the path ending at u is either singleton u, giving
one cover, or has a nonempty subset of TT_a immediately before u and an
arbitrary subset of TT_e before that. The second path contains the other
vertices in their forced TT_e-then-TT_a order. Hence F=`1+E0(A0−1)`.

For the T deletion graph, the nonempty TT_c prefix leaves H unchanged
and forces the u-start count to zero. Each prefix vertex independently
chooses which of the two ordered cover slots to precede, so it multiplies
the rooted F by C0. Thus its state at u has

    h_parent=(A0−1)(B0−1)+2,
    s_parent=0,
    F_parent=C0[1+B0(A0−1)].

Reinsert physical root0 as a newborn beaten only by u, and mark0. This is
exactly the inherited N map, not an arbitrary cut fit. It gives

    q_T=(A0−1)(B0−1)+2,
    h_T=(A0−1)[B0(C0+1)−1]+C0+2,
    z_T=(h_T,q_T,q_T,h_T).                           (15)

For the U deletion graph, apply (14) with E0=B0C0 and the same N map:

    q_U=(A0−1)(B0C0−1)+2,
    h_U=(A0−1)(2B0C0−1)+2,
    z_U=(h_U,q_U−1,q_U,h_U+1).                       (16)

Consequently every inflated U is unit-ready, and every inflated T has
q−s=0. The full coordinate F is retained; neither claim follows merely
from the indegree-one root star.

## 3. Exact obstruction to common same-target inflation

Subtract the two counts:

    h_U−h_T = (A0−1)B0(C0−1)−C0
             = B0(A0−2)(C0−1)
               +(B0−2)(C0−1)+(C0−2).               (17)

All three summands are nonnegative because A0,B0,C0≥2. Equality first
forces C0=2, then B0=2 and A0=2. Therefore

    h_U=h_T  iff  a=b=c=1.

This is an all-size obstruction in precisely the common-inflation model
specified above. It is not a finite-search inference and does not exclude
other modules, different inflation parameters on the two sides, or other
target-preserving changes. The structural reason is visible in §2: the
flip merges the C and B modules on one side, whereas the other side has a
nonempty ordinal prefix in its deletion graph.

For example inflating only block C gives `h_T=3·2^c+3` and
`h_U=2^(c+2)+1`, with gap `2^c−2`. Inflating only block A gives
`h_T=5·2^a−1` and `h_U=7·2^a−5`, with gap `2^(a+1)−4`.

**All homogeneous modules, with unit readiness retained.** The obstruction
is not limited to a choice of transitive replacements. Replace vertices
1,2,4 on both sides by arbitrary nonempty tournaments A,B,C, leaving0,3
singleton and using the same uniform exterior incidences. On the U side,
root deletion is still `C3[A,{3},C→B]`. A Hamilton path starting at3 must
traverse the entire C→B module and then A. Hence its number is
`H(A)H(C)H(B)`, and the newborn-root identity gives

    q_U−s_U = H(A)H(B)H(C).                          (17a)

The T deletion graph still has the nonempty ordinal prefix C, so its
corresponding readiness defect is zero. Thus the repaired side is
unit-ready iff every one of A,B,C has exactly one Hamilton path. A
tournament has H1 iff it is transitive: a nontransitive tournament contains
a directed triangle, whose three paths each extend injectively through
successive added vertices by ordinary path insertion; deleting the added
vertices recovers each old word. Transitive tournaments have exactly their
total-order path. Therefore unit readiness forces all three modules to be
transitive, and (17) applies. **The only common homogeneous substitutions
that are both same-target and unit-ready on the repaired side have all
three modules singleton.** This quantifies over arbitrary nonempty module
tournaments, not just bounded orders or a preselected module bank.

## 4. The positive survivor: a two-width unit-ready source

Unmarked, the inflated U is exactly

    C3[TT_a,{u},TT_r],    r=b+c+1,

with the final block ordered `root0,C,B`. Its root is the **first vertex
of the block after the singleton**. Session298's displayed marked chart
places its root in the other block, before the singleton; the unmarked
family agrees, while the rooted response is different.

The direct construction in §2 works independently for every a≥1,r≥2:
start with `C3[TT_a,{u},TT_(r−1)]` and append a newborn beaten only by u.
The source has

    h=(2^a−1)(2^r−1)+2,
    q=(2^a−1)(2^(r−1)−1)+2,
    z=(h,q−1,q,h+1).                                (18)

No empty block is inserted into the positive-inflation proof: r=2 is a
separate instance of this direct construction. At a=4,r=2 it gives the
same full target47 with marked state `(47,16,17,48)`. At a=1,r=3 it is
the repaired H9 source `(9,4,5,10)`.

The actual word N R at this available root has balanced profile

    (3h+2,h,h,3h+2),       gcd(3h+2,h)=1.             (19)

The original inflated T instead has NR profile `(3h_T,h_T,h_T,3h_T)`.
At the original shared target9 these are the distinct exact exits29 and27.
For larger common inflations the starting full counts already differ by
(17), so this is not a same-target comparison there.

Every source in (18) also has a proved route into the repeatable primitive
controller. Each factor `2^j−1` is0 or1 modulo3, so h is0 or2 modulo3.
For the pair in (19), `b−2t=h+2` is therefore a unit modulo3. One C3 cycle
is primitive and reaches `3|t'`, `3∤b'`, after which the inherited C3/C5
continuations apply. The source remains sparse and specified by its two
widths; this does not supply every desired h.

Finally (18) still has birth-return defect zero. Its successful immediate
NR is the k=0 boundary, and one must not replace it by `R^k N R`, k≥1:
that different clock has the inherited zero-defect factor trap. Unit
readiness and the positive-defect chain-renewal chart solve different
continuation obligations.

## 5. The next H29 repair also has a precise inflation boundary

The coordinator's separate three-arc rearming at the NR output29 has,
after deleting its global sink, the six-vertex quotient

    G=(24,1,11,34,46,7),    root4.

In this quotient's label order the vertices are the original names
`(0,2,3,4,5,6)`. Inflate only quotient vertex1, the original vertex2,
to TT_b, keeping b≥1. The root's unique incoming neighbor remains
quotient vertex0. In the root-deletion graph, a Hamilton path starting
at that predecessor is forced to traverse the outside vertices in original
order `0→4→6→3`, then all of TT_b. Consequently the inflated graph
remains unit-ready for every b.

Its exact full rooted state is

    h=15·2^b−1,
    q=3·2^b+3,
    z=(h,q−1,q,h+1).                                (20)

Here is a run-count proof of the full formulas, rather than fitting a
dyadic sequence. The TT_b module has only one outgoing outside neighbor,
vertex0. A Hamilton path can therefore use at most two module runs: each
nonfinal run must be followed by0, which appears only once. Contracting
each nonempty run gives a finite quotient skeleton. The exact counts of
one-run and two-run skeletons are `(29,15)` for G and `(9,3)` after root
deletion. These finite constants are certified by direct enumeration of
the fixed six- and five-vertex skeletons in the script, with no inflated
graphs needed to obtain them. One run has one internal assignment; two
ordered nonempty runs have `2^b−2` assignments, each with its forced
transitive order. Hence H=`29+15(2^b−2)` and q=`9+3(2^b−2)` as claimed.
The forced predecessor-start path then proves s=q−1 and F=h+1 by the
same N construction used above.

The run proof has an all-module form. Replace TT_b by any nonempty
tournament M. Put h_M=H(M), and let c_M count **ordered spanning covers
by two nonempty paths**. A one-run skeleton has h_M internal assignments;
a two-run skeleton has c_M. The modified source therefore has exact state

    (29h_M+15c_M, 8h_M+3c_M, 9h_M+3c_M, 30h_M+15c_M).   (21)

For the root-start defect, the forced outside path in the deletion graph
ends with an arbitrary Hamilton path of M, giving h_M starts at the
predecessor. The indegree-one N identity then gives the last coordinate
in (21), not an assumed independent cover-run bound.

Before the repair, the original U9 source with M substituted at vertex2
has one-/two-run coefficients `(9,4)`, also obtained by the fixed quotient
skeleton enumeration. Its readiness defect is h_M by the same forced
predecessor-start argument, so its rooted F is `10h_M+4c_M`. Its raw NR
output is consequently

    (29h_M+12c_M, 9h_M+4c_M, 9h_M+4c_M, 29h_M+12c_M).   (22)

The exact full-count drift is **3c_M**. Every tournament of order≥2 has
a Hamilton path; splitting that word into two nonempty consecutive pieces
shows c_M>0. A singleton has c_M=0. Hence this proposed common-module
inflation of the H29 repair preserves the target **iff M is a singleton**,
even without imposing unit readiness. The repaired readiness defect is
h_M, so unit readiness alone is equivalent to M being transitive. For
TT_b, h_M=1 and c_M=`2^b−2`, recovering (20) and its drift below.

Before that repair, inflating the original H9 source in this way and then
using NR has full count `12·2^b+5`. The proposed inflated repair therefore
has exact drift

    (15·2^b−1)−(12·2^b+5)=3·2^b−6.

It preserves the target only at b=1. This is another useful unit-ready
family, but does not extend the same-target H29 rearming along its raw
dyadic ray. The coordinator owns the original finite three-arc witness
and its radius-minimality control; no such search is repeated here.

## 6. Reproduction and finite scope

    python3 -B 05-knowledge/results/h17_exact_repairs_20261112/constructors/h9_inflation.py

The checker constructs all128 literal inflated graphs for1≤a,b,c≤4,
counts their rooted profiles by subset endpoint DP, and compares (15)–(17).
It retains literal rows for selected controls and records the exact universe.
It separately constructs the r=2 source at H47 and its literal NR exit143.
It also certifies the fixed quotient skeleton constants in §5 and checks
three small inflated H29 repairs plus a directed-triangle module hostile,
with their target drift and readiness defect retained.
These controls support the implementation; the all-parameter proof is the
deletion/ordered-cover argument above.
