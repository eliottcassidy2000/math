# The two H9 cores exchange an actual unit-ready interface

2026-11-12. **PROVED:** the one-arc marked rebase and its exact two-birth
exit below. The five cores, their two H9 members, and the numerical
attainability of9,15,29 are inherited. The new object is the physical
same-root-star repair and its continuation interface. This is one ready
entrance, not an indefinitely renewable ready-return grammar.

## 1. Same count, same root star, different readiness

Use outgoing-bit rows, with the same physical root0 throughout preparation:

    A=(22,8,2,5,14),       B=(22,8,2,21,6).

The sole edit reverses **4→3 to3→4**. It touches neither root0 nor any
incident arc. The root has the same unique predecessor3 in both graphs;
in particular3→0 survives. Both contain the directed spanning cycle
0→4→2→1→3→0, so both are strong. Their complete rooted profiles are

    z_A=(h,s,q,F)=(9,3,3,9),
    z_B=(h,s,q,F)=(9,4,5,10).                       (1)

Here q=H(T−0) and s counts Hamilton paths starting at0. The inherited
[unit-ready definition](../../session296_multicore_metagraph_20261113/cores/PROOF.md)
is an indegree-one root with **defect q−s=1**. Consequently A has defect0,
whereas B is unit-ready, without changing H, the root name, its predecessor,
or its entire incident star.

For completeness the nine Hamilton paths on each side are as follows;
strings denote vertex sequences, including initial zeroes.

| A | B |
|---|---|
|04132|01342|
|04213|02134|
|04321|04132|
|13042|04213|
|21304|13042|
|30421|21304|
|41302|30421|
|42130|41302|
|43021|42130|

Deleting0 leaves respectively

    A: 4132,4213,4321;
    B: 1342,2134,3421,4132,4213.

In a graph with sole predecessor p of its root, root-start paths are
precisely the deleted paths whose first vertex is not p. Thus the defect
counts deleted paths starting at p. The edit creates exactly one such
path, **3→4→2→1**. This is the missing marked coordinate, not just an
improvement in an unpointed statistic.

## 2. A general indegree-one identity and the exact target15 repair

For every nonempty tournament with an indegree-one root r, one has

    F=h+q−s.                                       (2)

Indeed, concatenate the ordered cover(P,Q) counted by F, splitting after r.
Successful concatenations give each Hamilton path exactly once. Failure
requires Q to start at the sole predecessor p. Then P must be the singleton
r: otherwise its last arc requires p to belong to P. The failures are
therefore exactly the deleted paths starting at p, counted by q−s. This
proves(2), and completes the F-coordinate calculation in(1).

The inherited R/N maps, with the conventions of
[CONTINUATION_FIBRES](CONTINUATION_FIBRES.md), give

    N z_A=(15,6,9,18),       N z_B=(15,5,9,19),
    NRz_A=(27,9,9,27),      NRz_B=(29,9,9,29).       (3)

The first birth adds vertex5 with0→5 and5 beating1,2,3,4; N marks5.
The same named pair(0,5), and all its incident arcs, exist on both sides.
Thus the exact first target **15** is preserved while the root-child pair's
continuation divisor changes3→1, by the inherited
[root-child gcd theorem](../../two_port_return_induction_20261103/ROOTED_RESETS.md).
The individual rooted profile gcds in the first line of(3) are also3 and1.

The second birth adds6 at root5 and retains5. Its new arcs are5→6 and
6→0,1,2,3,4. For every indegree-one source, the maps and(2) imply

    H(NR(T,r))=h+2F=3h+2(q−s).                    (4)

The ready source B therefore gives29=3·9+2, while A gives27=3·9.
This is the **k=0** boundary NR of the often used R^kNR schedule. At this
source the positive-k schedules have zero return defect and are imprimitive;
discarding k=0 would discard the useful repair.

The repaired output is balanced with b=29,t=9 and gcd(b,t)=1. Since
3 divides t and not b, it satisfies the inherited C3-entry guard. Its root
still has indegree one but now q−s=0. Thus(3) provides an accepted primitive
exit and a ready input at9; it does **not** provide another unit-ready
source at29 or prove that the same edit can be repeated there.

There is nevertheless a positive renewable continuation after changing the
interface. Apply the [retained-chain macro L_k](../constructors/CHAIN_REUSE.md),
k≥2, directly to the repaired first N state(15,5,9,19). Its proved formula
gives, with v=2^(k−1),

    s*=10v,       q*=30v+27,
    z*=(q*+3s*,s*,q*,4q*+2s*).

Here gcd(q*,s*)=gcd(27,10v)=1 and3 does not divide s*. At k=2 this is the
ten-vertex state(147,20,87,388). It enters the constructor's closed primitive
chart and accepted exit rules. The active root is now an already retained
chain vertex of indegree three; this preserves the physical two-port
budget but changes the readiness interface. The direct source B→L_k would
instead have content5, so the common target15 birth is essential to this
particular entry. The constructor's checker records the literal entry.

## 3. Why the triangle-overlap distinction matters here

The cyclic triangles are

    A: 013,034,123;       B: 013,123,134.

In A, the first triangle meets each other one in a directed arc, and the
last two share only vertex3. The arc-overlap graph is P3 (A3). In B all
three triangles share1→3, so the arc-overlap graph is K3. These are precisely
the inherited two H9 types in the
[five-core incidence note](../../five_core_web_20261111/history/INCIDENCE_AND_HISTORY.md).
Each has one directed spanning5-cycle; their whole unpointed slot polynomial
is the same:

    P_A(x)=P_B(x)=x^5+6x^3+2x.

The coefficients follow from the three directed triangles and one directed
pentagon in the inherited odd-cycle expansion. The two graphs are therefore
indistinguishable by every uniform context using only that polynomial.
Their different overlap geometry becomes actionable after retaining root0,
predecessor3, and the deleted path3421. The transfer is:

    five-core overlap type → one off-root arc edit → unit-ready marked source
    → exact common first target15 → primitive balanced exit29.

The unpointed polynomial and current count survive. The deletion and start
responses do not. Their difference(0,1,2,1) lies in the exact one-birth fibre
Delta=(0,v,w,v), with v=1,w=2. It changes the two following branches in
opposite directions:

    H_RR:27→25,       H_NR:27→29.

This is a physical branch selector at its specified small baseline, not a
free numerical ±2 actuator at arbitrary targets.

## 4. Provenance, tests, and retained boundary

The older [native reset induction](../../two_port_return_induction_20261103/NATIVE_RESET_INDUCTION.md)
already has a different H9 neutral repair: its N2 example changes the named
birth menu(23,17) to(17,31), with a different retained pair. The older
[reset junction](../../coverage_machine_induction_20261103/RESET_JUNCTION.md)
also contains H15 repair mechanisms. Neither their counts nor the existence
of the five cores is new here. The additional interface in(1)–(4) preserves
the root's whole star, makes the exact H15 newborn target primitive, and
uses the previously excluded NR boundary to realize a ready-source exit.

`hidden_probe.py` checks the displayed literal adjacency rows by independent
permutation/cut counting and subset DP, records all source/deleted paths,
compares every root-incident arc, and recounts both actual N and NR children.
It records the triangle sets, their overlap edges, and the unique directed
pentagons. The fixed five-vertex source is an exact finite object; equation(2)
and its consequence(4) apply at every order with the stated root guard.

The subsequent [all-size inflation theorem](../constructors/H9_INFLATION.md)
identifies a sharp extension boundary: common transitive inflation of the
three other vertices preserves this same-H repair only at the original
five-vertex pair. The repaired side nevertheless has an explicit two-width
unit-ready extension, including a differently marked H47 source. This
retains the positive family while distinguishing it from same-target
inflation of both sides.
