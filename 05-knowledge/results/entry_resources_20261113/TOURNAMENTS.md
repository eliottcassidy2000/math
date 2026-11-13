# Entry rules and quadratic endpoint budgets

**T1–T2: PROVED within the stated rule classes, using complete finite
counterexample certificates. T3–T4: FINITE-EXACT. T5–T7: elementary scoped
consequences. Chan–Pak and all473 inherited residual repairs remain OPEN.**

For a tournament T let H count directed Hamilton paths. Let q_ij count
ordered pairs of nonempty disjoint paths covering T, with first path
ending at i and second starting at j. Put W_ij=q_ij+q_ji, i≠j. For a
selected set O define

    E_T(O)=sum_{i<j in O} W_ij,   delta(T,O)=H−E_T(O).

The seed condition is delta=1. The separate slope U sums the ordered
two-cover counts whose second path starts in O, allowing the first path
to be empty. Both delta and U are required for the inherited dyadic-ray
constructor and its head threshold. Counts below are these actual covers;
they are not qualified-extension gcds substituted for realizability.

## T1. Every small unconditional selected-cut gadget is padding

A rule adjoins an arbitrary q-vertex tournament G, keeps all old markings,
and chooses which new vertices are selected. Each new vertex has a constant
orientation to O and another constant orientation to its complement.
Require the same rule to preserve delta=1 for every old unit seed.

For 1≤q≤5 the only such rules are unselected transitive padding: a
transitive prefix before every old vertex and a transitive suffix after
every old vertex, with the prefix preceding the suffix. There are q+1
rules per retained transitive-core representative,20 altogether.

The complete universe, modulo isomorphism of G before assigning its
types and markings, has counts8,64,1024,16384,393216:410696 total.
Each of the410676 rejected rules has an actual unit source of order2
or3 whose output defect is not1. The receipt lists every rule and its
witness; a separate two-phase spanning-word automaton verifies every
rejection. No sampling of the rule universe is used.

The survivors have an all-order proof. A Hamilton path must concatenate
the prefix, an old Hamilton path, and the suffix. In a taxed two-cover,
both cut endpoints lie in O: all prefix vertices must be in its first
path and all suffix vertices in its second. Thus both H and E are
unchanged. Adding b suffix vertices multiplies U by2^b; prefix padding
leaves U unchanged. It cannot introduce new odd parts of U.

This theorem says nothing about q≥6, adaptive rules, exceptional old
vertices, or rules on a proper guarded source family. It identifies why
coarse selected/complement interfaces cannot supply the desired entry
within this complete small-gadget class.

## T2. An unconditional rule seeing one stable port also collapses

Let p be an unselected stable port: replacing p by a transitive block of
any positive size preserves delta=1. Adjoin x. On O and on the remaining
old vertices other than p, independently choose x's adjacency to be
constant0, constant1, copy of p, or reversal of p. Choose the arc xp
and whether x is selected. These are64 rules.

Exactly four preserve every such input: an unselected source, an
unselected sink, and the two orders of an unselected transitive clone of p.
The60 other rules have counterexamples in the inherited634 stable-port
presentations through7. Their literal failures are independently checked.
The four survivors are valid by T1 or by the definition of stability.
Their old and new ports remain stable because these operations commute
with further transitive inflation. Only cloning preserves strong input
cores; it retains their existing contracted core.

The input bank is not asserted to contain every stable port at every size.
It is sufficient to disprove each proposed universal rule; the positive
rules are justified separately for all sizes.

## T3. A genuine core-changing move is not yet renewable

The inherited twelve-vertex source is

    [32,5,33,7,15,986,415,287,543,3295,2559,511], O={1,3}, p=0.

Adjoin x whose sole outgoing neighbor is p, then rename x as the root.
The first move is the known physical13/5 transition. Iterating the same
operation gives the following exact (H,delta,U) states:

|Order|H|delta|U|
|---:|---:|---:|---:|
|12|2663|1|11292|
|13|4439|1|29244|
|14|8023|−31|74616|
|15|15125|−29|201060|
|16|27587|−59|535224|
|17|50735|−119|1420740|
|18|93447|−207|3778824|
|19|171769|−385|10048836|

All selected subsets were tested. At12 and13 the only unit masks are
10,18,24, namely {1,3},{1,4},{3,4}. At14,...,19 there are none. Thus
the failed steps cannot be repaired merely by choosing a different marking
on those same cores. This is a finite obstruction, not an all-length theorem.

## T4. The budget complex and first higher obstruction

Define B(T)={O:E_T(O)≤H−1}. Nonnegative pair weights make this a
simplicial complex. Its saturated faces E=H−1 are the unit markings.
The complex is generally not determined by its allowed pairs: a directed
triangle has H=3 and all three pair weights2, so every pair is allowed
but the triple is forbidden. Nor do triple tests suffice at every size.

Complete canonical augmentation through8, with A000568 counts and
independent labelled orbit-mass checks, gives:

|n|Ordinary types|Distinct budget f-vectors|Nonflag complexes|Unit marked orbits|
|---:|---:|---:|---:|---:|
|1|1|1|0|1|
|2|1|1|0|2|
|3|2|2|1|4|
|4|4|3|3|9|
|5|12|6|11|22|
|6|56|12|55|59|
|7|456|25|455|194|
|8|6880|77|6879|1300|

The largest minimal nonface has size3 at every3≤n≤8. The first possible
size for a four-vertex minimal nonface is therefore9, and this actual
order-nine tournament realizes one:

    rows=[454,156,448,405,293,79,410,48,162], H=2151,
    O={1,4,5,6}, delta(O)=−1967.

Deleting one selected vertex gives defects209,11,7,141 respectively.
All four triples are within budget, yet their union is not. These are
defects, not claims of Hamilton counts7 or11. A separate signed-cut
automaton verifies all five values. The source was found on the eighth
fixed-seed random order-nine trial; that search is not a complete order-nine
census. Minimality of its order uses the separate complete census through8.

Paley11 gives a symmetric explanation of the same phenomenon: H=95095,
every pair weight is28842, so every triple is allowed and every quadruple
is forbidden. Unordered-pair transitivity makes these weights constant.
Paley3,7,11,19 have first forbidden sizes3,3,4,4. This finite pattern
does not establish a bound on obstruction size.

Useful intrinsic sequences now have precise definitions: numbers of distinct
budget f-vectors, nonflag budget complexes, and the maximum minimal-nonface
size at each tournament order. Full per-type receipts retain all minimal
nonfaces and saturated-face counts. Face totals sum vertex subsets on one
representative per ordinary class; marked orbit counts divide only saturated
faces by |Aut(T)|. This division uses the inherited marked-unit rigidity
proof in [Session288](../session288_twin_lengths_20261112/audit/REVIEW.md):
the marked automorphism group acts freely on both Hamilton paths and taxed
covers, so its order divides delta=1. No such free action is assumed on
ordinary budget faces.

## T4a. Symmetry turns the full budget complex into a simplex skeleton

**PROVED structural reduction; FINITE-EXACT values below.** If Aut(T) is
transitive on unordered vertex pairs and n≥3, all W_ij equal one integer w. Then

    delta(T,O)=H−binom(|O|,2)w.

If r is the first cardinality with binom(r,2)w>H−1, the entire budget
complex is the complete (r−2)-skeleton of the (n−1)-simplex. Its reduced
integral homology is concentrated in dimension r−2, free of rank
binom(n−1,r−1). To see this, cone every lower-dimensional cycle to one
vertex inside the available skeleton; the top rank then follows from
the alternating binomial identity for its Euler characteristic. The top
cycle group is a subgroup of a free integral chain group, so has no torsion.

This applies to every Paley tournament over a field of order q=3 mod4:
translations and nonzero square multipliers act transitively on unordered
pairs. Choosing the sign of a difference permits the needed square ratio.
Thus an infinite symmetry type has an exact budget description in terms
of only H and w, although their general numerical values remain to be found.

|q|H|w|gcd(H,w)|Budget dimension|Reduced homology rank|
|---:|---:|---:|---:|---:|---:|
|3|3|2|1|1|1|
|7|189|78|3|1|15|
|11|95095|28842|209|2|120|
|19|1172695746915|223117808922|9|2|816|

Only q=3 has a unit marking in this control set. For q=7,11,19 the common
content already excludes delta=1. This is stronger than counting many
allowed faces or observing their topology. In particular

    H(Paley11)=5·7·11·13·19,
    w(Paley11)=2·3·11·19·23,
    delta(Paley11,O)=209[455−138 binom(|O|,2)].

These factors have a precise role:11·19 divides every selected defect.
The per-vertex slope is u=H/q+(q−1)w/2 by vertex transitivity. For a
new vertex whose outgoing neighborhood is O, k=|O|, the exact extension is

    H_child = H + (u−H/q)k − w binom(k,2)
            = H + (w/2)k(q−k).

Indeed the general formula is H_child=H−sum_(i in O)end_i+U(O)−E(O):
split a child Hamilton path at the new vertex, allowing either end to be
empty. The end subtraction is essential. Vertex transitivity makes every
end_i=H/q and every U_i=u, proving the displayed specialization.

Thus this entire one-vertex extension layer has exactly(q+1)/2 count
values. The value indexed by0≤k≤(q−1)/2 has labelled attachment
multiplicity2 binom(q,k), pairing complementary new neighborhoods.
This counts attachments, not child isomorphism classes. Finite differences
give C_1=gcd(H,(q−1)w/2,w)=gcd(H,w), namely1,3,209,9 in the controls.
For Paley11 the six values are

    95095,239305,354673,441199,498883,527725.

Literal children through source order11 are independently counted. The
order19 profile follows from its exact endpoint tables and the proved
formula; no general int64 order20 computation is used. No exclusion for
every larger Paley order is inferred from the four content values.

## T5. Exact search after endpoint tables are available

Unit entry is a nonnegative quadratic subset-sum problem E_T(O)=H−1.
At a partial selection A retain its cost and candidate vertices. A candidate
v adds the marginal sum_i∈A W_iv. Reject over-budget marginals and pairs
whose weight alone exceeds the target. The gcd of all remaining marginals
and affordable candidate-pair weights divides every possible remaining
cost; use it as a necessary congruence guard. Branch over vertices in a
fixed order. Continue after equality to retain any zero-cost additions.

These rules lose no solution: every extension adds exactly its marginals
and internal pair weights, all nonnegative. `weighted_entry_search.py`
agrees with full subset enumeration on114 retained cores; the full
budget census supplies7,412 further per-class controls through8.
For the order19 core in T3 it visits80 subsets rather than524288.
This comparison measures the marking-search stage, not the cost of
constructing W, nor generation of ordinary tournament classes.

Equivalently retain the exact partition polynomial

    Z_T(z,t,x)=sum_O z^|O| t^E_T(O) x^(sum_i∈O U_i).

The coefficient of t^(H−1) is the complete unit-marking/slope supply.
This is a binary quadratic energy model with an additional linear field.
Discarding that field loses ray coefficients; evaluating at one temperature
loses exact energy levels. No general-purpose statistical model is assumed
to solve this coefficient-extraction problem efficiently.

## T6. Large feasible markings exist, without proving unit entry

For a uniformly random labelled n-vertex tournament and any fixed k-set O,

    E[delta(T,O)] = n!/2^(n−1) · (1−2k(k−1)/n).

Each permutation is a Hamilton path with probability2^(−n+1). An ordered
two-path cover is a permutation with one of n−1 cuts. Its n−2 required
arcs have probability2^(−n+2), and the fraction with both cut endpoints in
O is k(k−1)/(n(n−1)). Subtracting gives the identity (with the trivial
n=1 case read separately). Thus for n>2k(k−1) some actual T has
delta(T,O)≥1. The budget-complex dimensions are unbounded over tournaments.

Positive expected defect is not a unit-defect theorem: it neither selects
delta=1 nor supplies the desired U. The strict inequality is a sufficient
bound, not an optimal threshold. This distinguishes plentiful feasible
resources from exact target entry.

## T7. Stable inflation is not a matroid exchange system

The inherited physical J1 seed has two clone lengths with

    delta(a,b)=1−(2^a−2)(2^b−2),  a,b≥1.

Keep one mandatory vertex in each block and give the extra vertices two
disjoint pools A,B. Feasible optional subsets are exactly2^A∪2^B.
This need not be a delta-matroid: for |A|=2, |B|=1 choose feasible sets
X=A,Y=B and e∈B. No f∈X△Y makes X△{e,f} feasible. The two-coordinate
symmetric-exchange axiom fails. This is a direct consequence of the
[physical J1 formula](../joint_stability_forbidden_20261112/README.md#2-j1-the-first-physical-mixed-obstruction),
not a newly claimed J1 obstruction.

If both pools have r vertices and total optional size is fixed at r, the
only feasible states are A and B, at Hamming distance2r. No bounded-size
edit can connect these fixed-order template fibers for arbitrary r.
Without the fixed-size condition both branches connect through the empty
set. This says nothing about leaving this template or changing cores.
That distinction is exactly why the graceful buffer theorem is useful
inspiration: it enlarges the legal state space rather than assuming an
exchange axiom that the original fiber lacks.
