# A hereditary two-type menu inside every unit-ready target

**PROVED:** the reachability extension of inherited start monotonicity and
its balanced/unit-ready induced-root consequence. **FINITE-EXACT:** all
labelled rooted tournaments through order5 in `rooted_heredity.json`.
This is a necessary menu for induced growth, not a source construction.

## Inheritance and the retained root

The insertion proof under a strong ambient hypothesis is already in
[five-core bridges, lemma before Proposition5](../../chanpak_five_core_bridges_20261110/alternatives/NOTE.md)
and [READY_REPAIR RR1](../../lambda_emergent_20261103/READY_REPAIR.md).
The latter also proves an all-additions obstruction using minimum start
multiplicity. We retain the start vertex and weaken the ambient hypothesis
to reachability from that vertex; no new claim is made for the insertion
mechanism itself. Readiness and the indegree-one inverse are inherited
from [READY](../../two_entry_funnels_20261103/READY.md).

## H1. Positive start multiplicity controls all induced restrictions

Let x reach every vertex of a tournament U. Every directed path beginning
at x extends, preserving its first vertex and the relative order of its
old vertices, to a Hamilton path of U. Indeed, if its current vertex set
W is proper, some edge leaves W: otherwise x cannot reach its complement.
Choose an outside vertex z having an incoming edge from W, by a fixed
rule depending only on W. Insert z after the last path vertex beating z;
the following vertex, if present, is beaten by z. Continue until spanning.
Deleting the inserted vertices recovers the original path.

Consequently, for every induced U' containing x,

    starts_(U')(x) <= starts_U(x).                      (1)

It suffices to assume starts_U(x)>0, since a Hamilton path from x proves
the required reachability. The positivity guard cannot be dropped:
in TT2 the last vertex has zero spanning starts, but its singleton
restriction has one. In particular a unique-start parent has at most
one start at its marked vertex in every induced restriction, including
weak restrictions and the singleton boundary.

## H2. Every inherited root is balanced or unit-ready

Let(T,r) be unit-ready: F_T(r)-H(T)=1. The general weighted-cut law gives
delta=F-H-q+s>=d(d-1) and q>=s, where d is the root indegree. A source
has F=H, so d cannot be zero; d>=2 is also impossible. Thus d=1. Let x
be its sole predecessor and U=T-r. The generic newborn identity gives

    F_T(r)-H(T)=starts_U(x)=1.

For an induced T' containing r, there are two cases. If x is absent, r
is a source in T', and its state is(h,h,h,h). If x is present, put
U'=T'-r. By H1, b=starts_(U')(x) belongs to{0,1}. The newborn identities
give

    F_(T')(r)-H(T')=b,
    H(T'-r)-starts_(T')(r)=b.                           (2)

If b=0 the inherited root has balanced state(h,t,t,h). If b=1 it is
unit-ready. This proves the complete two-type menu for every root-containing
induced subtournament. No strongness assumption on T or T' is needed.

## H3. A necessary source menu for subdoubling ready returns

Suppose an h-count tournament V embeds induced in a unit-ready target T
of count N<2h-1. For the target root r, its deletion parent U satisfies

    N=H(U)+F_U(x)-1 >=2H(U)-1,
    H(U)<h.

If V omitted r, induced count monotonicity would give H(V)<=H(U), a
contradiction. If V contained r but omitted x, then r would be a source
of V and H(V)=H(V-r)<=H(U), again impossible. Thus V contains both r
and x. By H2 its inherited root is balanced or unit-ready, and its
deletion count t satisfies t<=H(U)<h.

For N=2h-71 the sharper bound is t<=h-35. Therefore any pure-birth
candidate return, even with an arbitrarily chosen final root, must admit
one of these two rooted source interfaces at an old vertex. This does
not prove that either interface embeds, or that either exists at an
arbitrary h. The complete finite-head exclusions in the containment lane
remain stronger on their five specified pairs.

## H4. A whole forbidden return interval for two actual881 sources

For the literal sources

    A2=(518,388,456,1011,487,455,387,257,1,502),
    B9=(354,228,425,179,231,192,140,1,250),

the independent endpoint subset calculation gives H881 and no vertex
with F-H in{0,1}. Every root state is retained in the receipt. H3 therefore
excludes **every unit-ready target count N with881<=N<1761** from every
pure-birth descendant of either specified graph, at arbitrary order and
with any final root. This conclusion needs no enumeration of target
classes or bound on their orders. It also covers transitive padding,
which still contains the specified source induced.

The fixed initial graph matters. A same-H change of old arcs can alter
its root menu, and is outside this particular interval obstruction. The
separate five-head containment theorem allows arbitrary initial same-H
preparation because its missing count is tested over the complete target
fibre. The two obstructions thus have different, complementary scopes.

## H5. The recovered all-additions mechanism strengthens this to all targets

There is a stronger obstruction for these two sources. The minimum-start
lemma under strong extension and the two root cases below are the actual
mechanism of the older five-core bridges Proposition5, cited above. Its
application to A2 and B9 gives an all-count conclusion.

**Criterion.** Let S be a nontrivial strong tournament such that

    min_x starts_S(x)>=2,    min_r(F_S(r)-H(S))>=2.      (3)

Then no induced extension of S has a unit-ready root, at any count or
order. First suppose the target is strong, with ready root r and strong
unique-start deletion parent U, marked at its sole predecessor x. If
r is outside S, U contains strong S. The inherited minimum-start lemma
forces starts_U(x)>=min starts_S>=2, contradicting unique start. If r
lies in S, it must have exactly one in-neighbour in S: zero contradicts
strongness of S, and two contradict target indegree1. Its sole predecessor
is the same x. The generic newborn identity gives

    starts_(S-r)(x)=F_S(r)-H(S)>=2.

The inherited fixed-x insertion lemma into strong U contradicts
starts_U(x)=1. This exhausts all target-root locations.

A weak ready target of count>1 is a strong ready core followed by a
transitive tail, as independently proved in `CONTAINMENT_AUDIT.md`.
A nontrivial strong S cannot straddle its SCCs or lie inside that tail,
so it embeds inside the strong ready core, already excluded. This proves
the criterion at every order.

The direct literal calculation verifies strongness and gives

| source | minimum starts | minimum readiness defect |
|---|---:|---:|
| A2 |14|10|
| B9 |40|96|

Therefore neither actual source has **any** unit-ready pure-birth
descendant, regardless of target count, number of births, or final root.
The interval conclusion H4 remains a weaker application with different
hypotheses. This all-count result concerns these fixed sources; old-arc
repairs can escape it, exactly as the new three-arc A2 repair does.
It neither forbids unit-ready targets themselves nor allows arbitrary
same-count initial replacements in the hypothesis.
The earlier balanced native-entry constructions remain valid: their
target defect is0, whereas this obstruction requires target defect1.

## H6. Five vertices do not exhaust the hereditary obstruction basis

The regular cyclic tournament R5 has rows(6,12,24,17,3), three starts
at every vertex, and root indegree2. Thus H5 excludes it from every
unit-ready target; this exclusion was already proved, even for2-ready
targets, in the inherited READY_REPAIR RR2. There is a distinct six-vertex
minimal obstruction:

    B=C3[TT2,TT2,TT2],    rows(14,12,56,48,35,3),
    H(B)=45,
    starts=(8,7,8,7,8,7),
    F-H=(13,26,13,26,13,26).

It satisfies H5 and contains no induced R5. Indeed any five selected
vertices contain two in one transitive module. Their comparisons outside
that module agree, while their induced outdegrees inside it differ, so
the induced tournament cannot be regular.

Both R5 and B are **minimal forbidden induced configurations for
embeddability into a unit-ready tournament**, not merely minimal for the
sufficient numerical criterion. The checker exhausts every nonempty
proper induced subset:30 for R5 and62 for B. Each is either transitive
or strong with a vertex having exactly one Hamilton start. A transitive
tournament of order>=2 is unit-ready at its second vertex in transitive
order; a singleton embeds in TT2. A strong unique-start tournament has
an actual unit-ready canonical newborn extension. Those literal newborn
counts and defects are checked, supplying positive embedding witnesses
for every proper induced subset. The full two graphs are excluded by H5.

For B this finite verification also has a simple type explanation: a
strong restriction must meet all three blocks and is C3[TT_a,TT_b,TT_c]
with1<=a,b,c<=2; unless all equal2, one of its checked seven module-size
patterns has a unique-start vertex. Restrictions meeting at most two
blocks are transitive. Every tournament

    C3[TT_a,TT_b,TT_c],    a,b,c>=2,

contains B and therefore forbids unit-ready induced extensions. The same
module-degree argument makes every member R5-free. This is an infinite
family of obstructions invisible to a test for the regular five-core
alone. It does not classify the full forbidden induced basis; it identifies
two actual minimal members and a complete infinite family generated by
one of them. The five/six connection here is an induced-embedding theorem,
not a numerical analogy with unrelated five-element collections.

Run `python -B .../geometry/forbidden_cores.py`, and again with `-O`, for
the whole-core profiles and every proper-subset positive certificate.

Reproduce `python -B .../geometry/rooted_heredity.py`, then with `-O`.
The test includes all weak and strong labelled tournaments through order5,
every root, every induced subset, the positivity-guard hostile, and all
subdoubling cases inside the unit-ready controls. The proofs above are
all-order statements; these finite controls audit their boundaries.
