# An arbitrary-module rearm and a target-neutral primitive repair

**PROVED:** the three inherited 462 charts extend to every nonempty
tournament module, with full rooted responses, a protected physical
interface, and an infinite same-target repair of rooted content7 to1.
**FINITE-EXACT:** the small controls and bounded observer search below.
The repair does not restore unit readiness for a nontransitive module.
It does not prove arbitrary-target coverage or Chan–Pak.

## 1. Inheritance and the retained interface

The source is [FAMILIES, F1–F3](../../multicore_metagraph_20261113/FAMILIES.md).
It supplies three literal ready charts with a protected `TT_k` block and
port pair `x→y`, at counts

\[
462\,2^k-5,\qquad1386\,2^k-13,\qquad4158\,2^k-37.
\]

The current task is to change the core at an exact target without
forgetting the physical port. The hostile is the complete homogeneous
inflation obstruction between ready881 and1691 in the same source. The
least-used sidecar here is the deletion count at the port, not merely the
readiness difference. The live objects are a retained module, its one/two
path-cover counts, the two port vertices, selective core changes, and the
four rooted continuation coordinates.

The older HN boundary mechanism is used precisely, not as a scalar
analogy: [WEDGE_GEOMETRY, Section3](../../atomic_tail_reset_20261029/WEDGE_GEOMETRY.md)
separates geometric branch selection from a prescribed boundary extension
relation. Here count equality, port readiness, and primitiveness are also
different predicates. The full marked factorization in
[the earlier ports proof](../../h17_exact_repairs_20261112/ports/PROOF.md)
shows why an unmarked factor identity cannot determine them.

Let `M` be any nonempty tournament, with

\[
h=H(M),\qquad c=C_2(M),
\]

where `C2` counts ordered covers by two **nonempty** directed paths.
Replace only the protected middle `TT_k` by `M`, keeping all exterior
comparisons. Thus `M→x→y→M`, and all its internal arcs are retained.
The fixed cores and tail words are

|chart|rooted core, root0|left tail|right tail|
|---|---|---|---|
|1|`(2,4,1)`|`0`|`1011`|
|2|`(2,36,57,51,35,1)`|empty|`101`|
|3|`(2,4,105,115,103,67,3)`|empty|`011`|

The core dominates the whole tail. The old funnel root `x` beats core
vertex0, loses to all other core vertices, and loses to the entire module;
at fixed tail bits it loses at1 and wins at0. Finally `x→y` and `y` beats
every other vertex. These declarations specify actual tournaments `T_i(M)`.
Their orders are `|M|+10, |M|+11, |M|+12`. Each is strong: replacing a vertex
of the strong singleton quotient by a nonempty tournament preserves strong
connectivity, including between two vertices inside a nonstrong module via
an exterior directed cycle.

## 2. Exact two-run theorem and full rooted response

Write a rooted state as `(H,s,q,F)`, where `s` counts Hamilton paths
starting at `y`, `q=H(T−y)`, and `F` counts spanning two-slot covers with
the first slot nonempty and ending at `y`; the second slot may be empty.

**Theorem.** For every nonempty `M`, the states are

\[
\begin{array}{c|rrrr}
 &H&s&q&F\\\hline
T_1&919h+462c&262h+126c&263h+126c&920h+462c\\
T_2&2759h+1386c&224h+90c&225h+90c&2760h+1386c\\
T_3&8279h+4158c&554h+210c&555h+210c&8280h+4158c.
\end{array}
\tag{1}
\]

In particular `F−H=q−s=h` in every chart.

*Run bound.* Delete `x`. Every out-neighbour of the module is then a fixed
tail vertex to its right. Such a vertex cannot reach the module in the
remaining graph: the core, left tail, module, and `y` all dominate the
right tail, which is internally transitive. Consequently every gap between
two module runs of a Hamilton path uses `x`. Since a path uses `x` once,
there are at most two module runs. The same argument applies after deleting
`y`. Contracting the runs gives coefficients independent of the module:
one run contributes `h`, and two ordered runs contribute `c`. This is an
actual positive path decomposition, not interpolation in module order.

For the full count, the inherited all-k identity `A·2^k+B`, together with
`H(TT_k)=1`, `C2(TT_k)=2^k−2`, determines the two coefficients as
`(2A+B,A)`. For deletion, the finite reduced-word coefficients are respectively
`(263,126)`, `(225,90)`, `(555,210)`. They are independently counted by
`run_coeff` in the accompanying script: each exterior vertex occurs once,
the module placeholder occurs once or twice, consecutive placeholder visits
are forbidden, and every consecutive pair follows its declared arc.
This finite coefficient certificate has no arbitrary module-size cutoff.

The parent `T_i−y` has exactly `h` Hamilton paths starting at `x`. Its first
step must enter core vertex0: entering a tail vertex first would make the
core unreachable. It then takes the unique rooted Hamilton path through
the core, and traverses the ordinal tail; the only freedom is a Hamilton
path of `M`. For an indegree-one root `y`, concatenating a root-ending
two-slot cover proves

\[
F-H=q-s=\#\{\text{parent Hamilton paths starting at }x\}.
\]

Thus `s=q−h` and `F=H+h`, completing (1). This also proves the readiness
assertion independently of the coefficient counts.

## 3. A literal same-target rearm after NR

A birth at root `r` adjoins `v` with `r→v→all other old vertices`. `N`
marks the newborn; `R` keeps the current root. For an indegree-one rooted
source of count `H` and defect `h`, the two births `NR` have state

\[
(3H+2h,H,H,3H+2h).
\tag{2}
\]

Equation(1) gives, for `i=1,2`,

\[
H(T_{i+1}(M))=3H(T_i(M))+2h.
\tag{3}
\]

After those two births, discard their two new vertices and the old fixed
exterior. Insert the next chart's fixed exterior, retaining the actual
module, `x`, and `y`. All comparisons among retained vertices agree exactly:
they are `M→x→y→M` plus the unchanged internal graph of `M`. The completed
macro has the exact same count as (2) and once again has defect `h` at the
same physical `y`, whose sole predecessor is the same physical `x`.

This is a bounded exterior replacement, with size independent of `|M|`.
Its intermediate deletion/insertion cards need not preserve count or
strongness. Both completed graphs are explicit strong tournaments.
For `M=TT_k`, it is the inherited unit-ready chain. For general `M`, it
is a new full-response lift and a rearm of the same nonunit defect.

There is an exact unit boundary: `h=1` iff `M` is transitive. The forward
direction follows because any nontrivial strong component has at least
three Hamilton paths, while SCC counts multiply. Thus arbitrary-module
freedom by itself cannot create new unit-ready sources.

## 4. An infinite exact content7-to1 repair

The rooted content means the gcd of all four coordinates in (1). Since
`q−s=h` and `F−H=h`, elementary gcd reduction gives

\[
g_1=\gcd(h,42c),\quad g_2=\gcd(h,18c),\quad g_3=\gcd(h,42c).
\tag{4}
\]

The actual seven-vertex tournament

    D=(8,5,9,18,71,31,47)

has `H(D)=35`, `C2(D)=606`, and `P_D(2)=676`; the small subset-DP receipt
checks these directly. Let `M_t=D→TT_t` for every integer `t≥0`, including
the empty suffix at `t=0`. Ordinal slot factorization gives

\[
h=35,\qquad c=676\,2^t-70.
\]

This `c` is coprime to35: modulo5 it is `2^t`, and modulo7 it is `4·2^t`.
Consequently the three chart counts and rooted contents are

\[
\begin{array}{c|c|c}
 &H&g\\\hline
T_1(M_t)&312312\,2^t-175&7\\
T_2(M_t)&936936\,2^t-455&1\\
T_3(M_t)&2810808\,2^t-1295&7.
\end{array}
\tag{5}
\]

In the first transition, the raw `NR` output has gcd
`gcd(H(T1),2·35)=7`: its source count is divisible by7 and is coprime
to5. The completed replacement `T2` has **the same H and rooted gcd1**.
Thus (5) is an infinite exact-target primitive repair, retaining a module
whose size grows with `t`, and not merely an equality of two unrelated
numbers. A second independently counted module of H49 gives the same
7-to1 phenomenon; its literal certificate is in the receipt.

The important limit is visible in the last column. Readiness remains35,
and the useful deletion/start residue is lost by the next ordinary `NR`.
In fact `T3` has content7 again. Neither current primitiveness nor a
same-count primitive replacement is automatically a renewable unit-ready
controller. The physical marks and the whole four-coordinate state must
remain in the next inverse problem.

There is now a positive escape from this exact boundary:
[PORT_ESCAPE](PORT_ESCAPE.md) changes the active mark from `y` to the
already retained `x`. A literal RNR then enters a primitive balanced
controller for every clock in both the original middle462 chart and
the H35 family above. Its general necessary-and-sufficient module gate
is given there. This does not restore the specific ternary ready ladder;
it changes the operation language while preserving the actual graph.

## 5. Selective funnel update: a precise additional coordinate

Homogeneous modules cannot solve the unit-readiness problem above. An
actual alternative is a core with a distinguished unique-start port,
contacted only through that port. The inherited six-state funnel compiler
uses

    z=(H,E(1),S(2),P(2),E(2),S(3)).

For an outer word `w`, its final ready count is `rho_w z−1`, where

    rho_w=(B1(w),−D1(w),L2(w),L2(w),−1,0).

If the old core is first changed by one funnel word `u`, the new observer
is `rho_v M_u`. Its `S(3)` coefficient is

\[
2L_2(v)L_3(u)>0.
\tag{6}
\]

Hence there is **no unrestricted six-state linear identity**
`rho_v M_u=3 rho_w`, for any such words, empty words included. This does
not exclude an identity on reachable core states satisfying an additional
relation, or changing both sides, or a more selective multiport interface.
The shortest explicit comparison is

\[
Q(F_1(C),1)=3Q(C,1)+2+12S_C(3).
\tag{7}
\]

The extra marked three-slot response is the exact missing coordinate.
It cannot be discarded on the grounds that the ordinary H agrees.

`intertwiner_probe.py` also tests every two-funnel observer whose two
words are each empty or end1 and each have length≤6: 64 words and4096
distinct rows. None is three times another in that finite universe.
This is a small algebra search and a stopping gate, not an all-language
nonexistence claim. Eighteen rows have gcd3, so row content alone does not
explain the failure.

The independent [periodic-extension theorem](../history/PERIODIC_CORE_OBSTRUCTION.md)
now proves the stronger all-length append-only boundary, and
[REWRITE_DEBT](../history/REWRITE_DEBT.md) gives a nonnegative offset
decomposition for the most immediate common-inner-core graft. The small
observer search above is retained as provenance, not the current limit.

## 6. Reproduction and scope

Run from repository root:

```bash
python3 -B 05-knowledge/results/target_directed_cores_20261112/ports/module_lift.py
python3 -B 05-knowledge/results/target_directed_cores_20261112/ports/intertwiner_probe.py
```

The first script keeps all recipes and literal graphs, counts reduced
one/two/three-run words, checks full rooted profiles independently on
`K1,TT2,TT3,C3`, directly replays `NR` where order≤15, and verifies the
retained induced interface across both replacements. Its two larger
primitive examples use small, independently counted module invariants and
the proved universal coefficient identities; they are not represented as
a new large-graph census. The second script uses exact integers throughout.
The receipts are `module_lift_receipt.json` and `intertwiner_receipt.json`.

The live next question is to replace the factor `h` in the parent-start
channel by a genuine selected unique-start port while controlling the
positive `S(3)` surplus in (6). Signed synchronized selected-cut defects
from [birth repairs B1–B2](../../birth_grammar_repairs_20261112/PROOFS.md)
do not supply that conversion: parent-start readiness is itself a
nonnegative path count, not their signed difference.
