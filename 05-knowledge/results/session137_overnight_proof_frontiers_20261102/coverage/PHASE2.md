# Phase 2 — repair at six, then a deletion-deck barrier at seven

**Status:** FINITE-EXACT, with structural proofs of the quotient and closed-orbit certificates. No all-order interval closure. Source provenance and inherited semantics are as in `RESULT.md`; no further repository reads or external class files were used.

## 1. Three different move graphs

All vertices below are strong tournaments of one fixed order. A move must preserve strongness and change H by `0` or `+/-2`.

1. **Inherited strong-parent atlas:** a cut flip or one-in/one-out exchange over a strong deletion parent, with reparenting allowed.
2. **Unrestricted-parent two-star atlas:** flip one arc, or any two arcs incident with the same vertex; the deletion parent may be nonstrong and the two cut bits may change in the same direction.
3. **Full one-vertex-replacement atlas:** reorient an arbitrary subset of the arcs incident with one vertex. This is equivalently replacement by any other strongly completing ear over its deletion parent.

The second and third graphs are genuine enlargements. Reparenting alone changes only the description of a child, not the child. Relabellings are quotiented explicitly when we work on isomorphism classes.

## 2. The proposed repair works at order six

`RESULT.md` finds five labelled components in graph 1. The explicit diagonal over the nonstrong parent `TT2 => C3` connects its H-images `{23,25}` and `{27,29,31,33}` in graph 2. No unit move can cross a genuine numerical gap in the global spectrum. Consequently graph 2 has exactly four **labelled** components:

| H image | labelled vertices | isomorphism classes |
|---|---:|---:|
| 15:2:19 | 4320 | 6 |
| 23:2:33 | 11760 | 19 |
| {37} | 3600 | 5 |
| 41:2:45 | 2640 | 5 |

These are exactly the maximal odd intervals in the strong order-six spectrum. Graph 3 has the same four numerical component images at order six.

This statement follows from the first packet's complete labelled atlas plus one verified diagonal. The new isomorphism-quotient enumeration independently returns the same four class components.

## 3. Source-preserving class enumeration, not a labelled 2^21 scan

`phase2.py` constructs all tournament isomorphism classes internally by one-vertex augmentation. Given all order-(n-1) classes, append each of the `2^(n-1)` cut patterns and canonicalize. Every tournament has a deletion card isomorphic to a retained parent, so this is a complete induction, not sampling.

Canonicalization first sorts vertices by outdegree and then exhausts every permutation within each equal-outdegree cell. Any isomorphism preserves outdegree, so the minimum LSB-first tournament code among these orders is a complete canonical label. The class tower independently gives

    orders 1..7: 1,1,2,4,12,56,456,
    strong classes at 5,6,7: 6,35,353.

For each strong order-seven representative, graph 2 requires only

    C(7,2)+7*C(6,2)=126

move presentations. Hence all `353*126=44478` presentations are examined, not all `2^21` labelled tournaments. Every retained spanning-forest edge includes a source class code, a literal labelled target code, the flipped pairs, and independently recomputed target H/strongness.

Graph 2 has **23 isomorphism-class components at order seven**, whereas the spectrum has eight maximal odd intervals. The complete records are in `phase2_order7.json`. In particular the unique H=61 class is isolated in graph 2 even though 59 is attained. It first has a unit star move of size three. Allowing larger stars therefore repairs some, but not all, failures.

## 4. Full replacement is exactly deletion-card incidence

Two strong isomorphism classes are adjacent under arbitrary one-vertex replacement if and only if they have an isomorphic vertex-deletion card.

* Forward: deleting the reoriented vertex leaves the same tournament.
* Reverse: choose an isomorphism of the two deletion cards, identify them, and orient the observer to match the two children. Their difference is entirely incident with that observer.

Thus graph 3 can be built exactly from the `353*7=2471` canonical deletion cards. Inside each card's class list, retain pairs whose H values differ by at most two. This gives a compact source-preserving certificate independent of explicitly enumerating every full star.

The order-seven graph 3 has eleven isomorphism-class components with images

    25:2:61,
    65:2:85,
    87:2:105,
    109:2:117,
    121:2:145,
    {147},
    151:2:159,
    {159},
    {171}, {175}, {189}.

Two essential scalar seams remain missing: `85 -> 87` and `145 -> 147`. The extra isolated H=159 class is a different phenomenon: the other H=159 class lies in a component already covering the entire `151:2:159` interval. Component count and numerical interval count must therefore not be conflated.

## 5. A small independent closed-orbit certificate

The cleanest hostile is the strong tournament

    T=(6,112,74,35,13,21,57),
    code(T)=85763,       H(T)=147.                     (1)

It is rigid. All seven of its deletion parents are **strong** and have saturated ear-response lattice `d=2`; the respective parent H values are

    43,29,29,41,43,43,43.                               (2)

Thus neither missing nonstrong parents nor a nontrivial response gcd explains this barrier.

### Finite theorem — the H147 orbit is unit-closed

For each of the seven vertices, reorient all six incident arcs in all 64 possible ways. Among these **448 presentations**:

* no strong replacement has H=145 or H=149;
* exactly eleven have H=147, and all eleven are isomorphic to T;
* the closest distinct H values of strong replacements are 143 and 151.

Therefore the isomorphism orbit of T is closed under every allowed graph-3 move. No sequence of reparentings, arbitrary one-vertex replacements, zero-response moves, and relabellings can take T to a different H value using only changes of magnitude at most two.

This conclusion does not require completeness of the 456-class census: `deck_leaf_audit.py` checks all 448 literal replacements directly by enumerating all `7!` vertex orders, and checks each H147 successor by an explicit full-permutation isomorphism. Its argument is closed under iteration because every zero-response successor is in the original orbit.

Meanwhile the strong tournament

    T'=(38,116,72,35,77,20,41),
    code(T')=610131,       H(T')=145                    (3)

exists. Equations (1)--(3) refute the all-order assertion that every solid global odd interval has a connected physical one-vertex-replacement unit witness, even after admitting every deletion parent and every cut move.

### A compact deck-disjointness certificate

With the canonical degree-cell convention above, the seven cards of T have codes

    1298,1313,8324,280,1298,8280,8280.

There are two H145 isomorphism classes, with codes `147857,17185`. Their combined card support is

    {209,385,2312,2433,8834,8881,12356,16524}.

The two sets are disjoint. This is another exact explanation for the absent 145/147 physical seam. The closed-orbit leaf audit is stronger as an independent certificate because it does not import the classification of the H145 fibre.

## 6. The cheapest repair changes more than one vertex's incident arcs

No set of one or two arbitrary arc reversals takes T to a tournament with H145. Exactly four three-arc reversal sets do so. One is the disjoint matching

    {(0,5),(1,2),(4,6)},

which gives T' in (3). The other three sets are retained in `deck_barrier.json`. None is a star, as required by the full-star obstruction.

Thus the missing information is not merely a larger cut move or an additional one-vertex parent. A physical unit route needs at least a multi-vertex edit at this boundary; scalar adjacency, as already licensed by THM-4118, remains sufficient for an interval certificate without inventing such a route.

The failed implication and survivor are now precise:

    all strongly completing one-vertex presentations, with all cut jumps
       -/-> a physical unit component spanning every scalar interval;

    overlapping/touching scalar component intervals
        -> a valid scalar interval certificate.

The new proof obligation is to organize these **scalar** joins uniformly, or to use genuinely multi-vertex operations with a controlled response. There is no remaining all-order one-vertex-unit-connectivity conjecture to pursue in its former form.

## 7. Replay and artifact contract

Run from this directory:

    python3 -B phase2.py > phase2.out
    python3 -B -O phase2.py > phase2_optimized.out
    cmp phase2.out phase2_optimized.out
    python3 -B fullstar.py > fullstar.out
    python3 -B -O fullstar.py > fullstar_optimized.out
    cmp fullstar.out fullstar_optimized.out
    python3 -B deck_leaf_audit.py > deck_leaf_audit.out
    python3 -B -O deck_leaf_audit.py > deck_leaf_audit_optimized.out
    cmp deck_leaf_audit.out deck_leaf_audit_optimized.out

`phase2_order*.json` retains component/class membership and labelled spanning-forest certificates. `fullstar_order*.json` retains every deletion deck and the common-card spanning forest. `deck_barrier.json` is the compact consequence certificate. `deck_leaf_audit.py` is independent of both engines and of all inherited repository code.

**Scope:** exact orders five through seven, plus the structural equivalence of full vertex replacement and common deletion cards. No all-order global spectrum theorem, global obstruction, or permanent missing odd integer is asserted.
