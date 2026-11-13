# Core changes with proved entry and a complete outer-tail obstruction

**PROVED scoped:** the family construction, protected-interface replacement,
and finite bounds for a prescribed dyadic family. **FINITE-EXACT:** the
atlas, its displayed chains, and the resulting exclusion. Chan–Pak and
seeds for every admissible odd coefficient remain **OPEN**.

Here **unit-ready** means E_y(2)−H=1, where E_y(2) counts ordered spanning
two-bin path covers with the first nonempty path ending at y. For an
indegree-one root this equals H(T−y)−s_T(y). The constructions below all
have indegree one; a difference-one deletion/start statistic at an arbitrary
higher-indegree root is not silently treated as the same condition.

Inheritance: [MODULE_RETURN](../cycle_opening_six_return_20261112/compiler/MODULE_RETURN.md)
repairs 45·2^k+7 into a ready representative of 135·2^k+23. Its mod-45
non-entry obstruction is a property of those source families. The full
unique-start grammar and its inverse are inherited from
[FUNNELS F1–F3](../polarized_path_calculus_20261103/FUNNELS.md), and the exact
six-coordinate compiler from
[TARGET_RETURN section 3](../berggren_target_return_20261111/reset/TARGET_RETURN.md).
The new comparison below searches that language with its physical module
interface retained. These are not newly discovered unique-start tournaments.

## A chart is a graph constructor with a specified interface

Choose a strong rooted tournament (C,a) with exactly one Hamiltonian path
starting at a, or the rooted singleton. Write its inherited state

    (h,e,s,p,f,t)=(P_C(1),E_a(1),S_a(2),P_C(2),E_a(2),S_a(3)).

Choose two fixed binary words L,R and b in {0,1}, so L b^k R ends in 1
for every k≥1. Adjoin that transitive tail, dominated by C. Adjoin x with
x→a, C−a→x, and tail vertex v→x exactly at a 1 bit. Finally adjoin y
beaten only by x. Denote this actual tournament Q(C,L,b,R;k).

The marked interface is the ordered transitive block M=TT_k, together with
x→y. The root y beats every vertex of M. The block beats x when b=1;
x beats the block when b=0. The parent at x has a unique starting path,
by the inherited grammar. Therefore Q is strong and

    delta_y = H(Q−y) − number of Hamiltonian paths starting at y = 1.   (F1)

This is **root readiness**. It is not the selected-cut condition H−W(O)=1,
and no unproved conversion between the two meanings of defect is used.

For w=L b^k R, the inherited tail-word coefficients give the full count

    H(Q) = L_2(w)(s+p) + B_1(w)h − D_1(w)e − f − 1.                  (F2)

No path sector is omitted. Write G=s+p. For k≥1 each coefficient is an
affine combination of 2^k and 1, so H(Q)=A·2^k+B. The constructor, not
only (A,B), is retained in every atlas entry.

## Exact family changes keep the same physical endpoints

Suppose two charts have the same b and satisfy A'=3A, B'=3B+d, where
d is 2 or 4. Map their module vertices in transitive order and retain the
same physical x and y. All comparisons among these k+2 vertices agree:
they are determined solely by b and x→y. Delete the other fixed vertices
of the source and insert those of the target with its declared comparisons.
The number of changed vertices is bounded independently of k.

The completed target is literally Q(C',L',b,R';k), so it is ready by F1
and belongs to the next chart by construction. The unique-start inverse
also recognizes its parent. This supplies **entry**, not just a matched
number. At the completed macro boundaries the same physical y is ready
and x is its sole predecessor. Intermediate deletion/insertion cards are
not claimed to remain ready or strong.

For d=2 one may first perform the standard two births, whose count is
3H+2 at a unit-ready source, and then make this same-count replacement.
That excursion's two newborns can be discarded while the original x,y and
M survive. For d=4 the atlas entry is a direct bounded replacement; it is
not falsely labelled a two-birth step from a defect-one source.

The complete bounded search uses every unique-start strong core through
eight (plus the singleton), every fixed context of total length at most
12, and both middle bits. After restricting A≤5000 and |B|≤5000, it retains
28,776 distinct (A,B,b) entries, choosing a shortest representative in each.
There are 265 exact same-interface edges: 159 with d=2 and 106 with d=4.
These counts describe this atlas, not ordinary tournament isomorphism types.

A particularly short three-core chain, for every k≥1, is

    462·2^k−5  →  1386·2^k−13  →  4158·2^k−37,                     (F3)

with H' = 3H+2 at each step. Its orders are k+10,k+11,k+12. The recipes are

|Core (root 0)|Left word|Middle|Right word|Full count|
|---|---|---|---|---|
|C3=(2,4,1)|0|1^k|1011|462·2^k−5|
|(2,36,57,51,35,1)|empty|1^k|101|1386·2^k−13|
|(2,4,105,115,103,67,3)|empty|1^k|011|4158·2^k−37|

The latter cores have ranks 11 and43 in the inherited complete unique-start ranker;
`family_audit.json.gz` also supplies literal graph rows for k=1,2. A second
chain is 455·2^k+59 → 1365·2^k+179 → 4095·2^k+539. These chains occupy
residues outside the original 7/23 mod-45 pair. They do **not** repair the
specific missing entry after 135·2^k+23, nor give unlimited iteration.

## Why the 405 target is absent from the entire outer-tail language

This exclusion allows arbitrary fixed core size and arbitrary fixed words
L,R; only the single growing block in the outer tail is prescribed.

For B_1, use the inherited six word states (last,started), ordered
(0,0),(0,1),(1,0),(1,1),(2,0),(2,1). At a bit b, one transition replaces
last by 1 if b=1 and 2 otherwise; a second sets started=1 when started
already equals 1 or b=0. Both have weight one. Its matrices M_b satisfy

    M_b³−3M_b²+2M_b=0,          M_b²−M_b is entrywise nonnegative.   (F4)

Their entries and exact products are in the receipt. It follows that
B_1(L b^k R)=alpha_B·2^k+beta_B for k≥1, with alpha_B≥0: its coefficient
is a nonnegative row times (M_b²−M_b)/2 times a nonnegative column.
D_1 is simply the indicator that the first tail bit is zero, so is constant.

If l=|L| and r=|R|, the 2^k coefficient of L_2 is

    alpha_L = 2^l (b + sum_j R_j 2^j).

The final 1-bit ensures alpha_L≥1 and alpha_L≥2^(l+r−1). Since
P_C(2)=sum_S H(C[S])H(C[V−S])≥2^m and S_a(2)≥1,

    A = alpha_L G + alpha_B h ≥ G ≥ 2^m+1,
    A ≥ 2^(l+r−1) G.                                           (F5)

For A=405, F5 forces m≤8 and l+r≤8. Both ranges are completely covered
by the atlas. The possible offsets are exactly

    −265, −199, −77, −25, −13, 5, 53, 61, 83, 109, 341, 395.

Thus **405·2^k+71 has no representation in this outer-tail language for
all k≥1**. This is a finite exact exclusion with a proved all-size reduction.
The atlas cutoff cannot hide another A405 offset: at k=1 the word has
l+r+1 letters, so L_2,B_1≤2^(l+r+1); F2 and h≤G give H(Q)≤8A, whence
1−2A≤B≤6A, inside |B|≤5000.

The first value 881 itself is attainable and ready. An inner clock, multiple
growing modules, a changed interface, or an eventual ray beginning after
k=1 lies outside the outer-tail exclusion. The separate inflation theorem
below rules out a much broader operation from its actual ready witnesses.

The empty-prefix, b=1 sublanguage also has a fast target-driven entry test:
A=alpha G forces G to divide A; alpha fixes R uniquely as the least-
significant-bit-first expansion of alpha−1. One then checks
B=h B_1(1R)−f−1. For 37 specified target pairs and all cores through15,
only 147 divisor candidates need their words tested. This gives a complete
recognizer for that sublanguage, not a selector for all odd targets.

## No homogeneous inflation of any ready 881 source reaches 1691

**PROVED using the cited order bound and a complete finite head.** Let T
be any tournament with H(T)=881 and a unit-ready root. Replace any or all
vertices by arbitrary nonempty tournament modules, with uniform comparisons
between modules inherited from T. The result never has H=1691. This even
allows replacing the original port vertices and does not require the result
to remain ready.

First classify the source. For a strong unit-ready881 source its unique-start
parent has H≤(881+1)/2=441, by
[READY](../two_entry_funnels_20261103/READY.md) and
[OBSTRUCTION O2](../ready_return_pentagon_20261103/OBSTRUCTION.md).
[Busch's published sharp bound, Corollary 1](https://www.combinatorics.org/ojs/index.php/eljc/article/download/v13i1n3/pdf/)
gives minimum strong H625 at order13, with increasing minima thereafter.
The parent therefore has order at most12. The complete rank census through18
contains exactly two such parents: (order,rank)=(9,168),(10,174). This is
an all-order classification of strong rooted sources at881, not a cutoff
guess. Every weak unit-ready881 source is one of these followed by a
transitive tail, by the inherited OBSTRUCTION O1.

Substitution is monotone in each module: expand a Hamiltonian path through
its old vertex by a contiguous Hamiltonian path of the inserted module.
In particular a nontransitive module, which has at least three Hamiltonian
paths, gives at least3·881=2643 paths and already overshoots1691. Thus all
potential modules must be transitive. Inflating any further vertex can only
increase the full count.

Of the 21 single-vertex TT2 replacements in the two strong sources, exactly
one has count at most1691: vertex8 of the source with parent (10,174), giving
1559. Every other TT2 replacement overshoots. In the one remaining case,
the positive run expansion is

    H(T[TT_m]) = sum_j a_j j! S(m,j),       a_j≥0.

Here a_1=881 and 2a_2=1559−881=678. For m≥3 the first two channels alone
give at least881+6·339=2915>1691. Hence sizes1,2,3-or-more all miss1691;
adding any second nontrivial module also overshoots. Tail vertices of a
weak source either remain transitive, leaving H unchanged, or acquire a
nontransitive module and multiply H by at least three. This completes the
proof, including arbitrary simultaneous and iterated homogeneous substitution.

Both881 and1691 have unit-ready representatives. The obstruction concerns
the operation connecting them. A realization of the 405 ray therefore
needs an actual change of quotient core already between its first two
values; inflating a fixed ready881 graph, even in several modules, cannot
do it. It is not a forbidden Hamiltonian count theorem.

## A finite inflation budget handles arbitrary module sizes

The preceding argument generalizes. Fix T, h=H(T), and a target h<N<3h.
Let c_v=(H(T[v←TT2])−h)/2. Every nontransitive module already overshoots
N, so only transitive module sizes m_v≥1 need consideration. Contracting
each maximal module run of a Hamiltonian path gives a positive run expansion;
the all-one run vector contributes h, and the vector with two visits only
at v contributes c_v. Consequently every realization obeys

    sum_v c_v (2^m_v−2) ≤ N−h.                              (F6)

This is a necessary budget, not an equality: interactions between split
modules add further nonnegative terms. For positive c_v it bounds m_v
explicitly. Enumerate the remaining integer vectors and count the resulting
graphs exactly. This is a complete target test for the inflation operation.

There is no hidden infinite search at c_v=0. Such vertices are exactly
singleton strong components. If v lies in a nontrivial SCC it lies on a
directed triangle. In T[v←TT2], the two copies can begin and end the path
copy1→a→b→copy2. Repeated vertex insertion extends this to a Hamiltonian
path while retaining the separated copies, proving c_v>0. Conversely,
two separated visits to v imply a directed cycle through v, so are impossible
for a singleton SCC. Replacing that singleton by any transitive chain has
no effect on H and can be normalized back to size one.

Applied to h_k=405·2^k+71 for k=1,...,6, the complete unit-ready source
counts are2,3,3,4,2,2. Busch's bound puts their parent orders at most
12,13,14,16,17,18, respectively. Across all16 strong sources, F6 leaves
only40 nontrivial size vectors; exact counting rejects every one for
N=h_(k+1). Weak sources again add only count-neutral transitive tails.
Thus each of these six adjacent steps requires a change beyond homogeneous
substitution, regardless of endpoint preservation. The small-order controls
also verify the zero-cost criterion on every ordinary tournament through7.

**OPEN hypothesis — adjacent 405-ray inflation separation:** this adjacent
inflation obstruction persists for every k≥1. This is deliberately stronger
than the evidence. Its precise source is the
complete first-six ready classification; its target is a structural lower
bound on how often a constructor must change quotient cores. A positive
inflation witness at the next count pair is a decisive falsification. Even
if true, this hypothesis would constrain constructions, not disprove the
existence of ready graphs or prove Chan–Pak.

## Which endpoint constraints remain real

If arbitrary core replacement is allowed, the bare unit-ready port pair
itself causes no compatibility obstruction. Every strong unit-ready graph
has a root y with sole predecessor x; strongness gives both an in-neighbour
and an out-neighbour of x in the parent. Together with x,y these give both
possible one-vertex module interfaces: M→x→y→M, and x→y→M with x→M.
Thus any two such graphs can retain x,y and a chosen one-vertex interface
while all other vertices are replaced. This does not preserve their separate
deletion and start counts, only their difference1 and the stated comparisons.

The genuine extra requirements are (i) existence of a ready graph at the
target count, (ii) presence of a whole protected TT_k interface, or (iii)
retention of more endpoint-response coordinates than subsequent operations
actually use. These requirements must not be conflated. For example the
census finds ready representatives at881,2645,7937, so a finite +2 path
through those counts can retain its original port pair. The unbounded family
entry problem is stronger than those three scalar instances.

## What a coverage proof still needs

The family address, actual port, and full target now survive each displayed
edge. What is absent is a theorem ensuring that every next target enters
some chart, with a decreasing inverse measure. A larger finite atlas alone
does not supply that theorem. In particular F3 is a uniform finite-depth
chain, and the old 473 repair families are not thereby closed.
