# LRC14 reflected levels: low components and one common-cell Hunter tree

**codex, 2026-08-01.** Reflection with a proved abstract reduction, exact
finite evidence, and an open analytic invoice.  It is not a proof of the
remaining reflected branch or of LRC(14).

## Inheritance pass

- **Closest proved mechanism:** THM-2941's reflected pair calculus.  Reduced
  level ratios with numerator-plus-denominator at least eight have the
  projective floor `1/105`; the complementary low-ratio graph has clique
  number five.  Its threshold-block referees now close every body having at
  least eight robust edges, namely `2,442/3,003` bodies.
- **Canonical hostile:** `E=(1,2,3,4,6,12)`.  A single high edge can have
  overlap below the whole six-label debt, even at the canonical upper-median
  body-safe cell.  Therefore `there is a high edge` is not a certificate.
- **Corrected near miss:** a deletion/compression induction need not preserve
  the certificate.  On the same hostile body, deleting the lowest level can
  leave zero avoiding gain even though the restored six-word closes.
- **Least-used relevant sidecar:** the Hunter--Kounias spanning-tree union
  bound from THM-1217/2911.  Its pair credits must all be evaluated on one
  common body-safe cell; independently optimized cells cannot be added.

The live concept board is: projective low shape, low component scale, high
complement, common-cell Hunter tree, upper-median boundary label, exact
projected debt, and toothpick/ray drift.

## 1. The intrinsic graph

Let `q=(q_1,...,q_6)` be six distinct positive levels.  For a pair, divide by
its gcd and write `(q_i,q_j)=g(P,Q)`, `(P,Q)=1`.  Define the labelled low graph

```text
ij in L(q)  iff  P+Q <= 7,
```

and let `H(q)=K_6\L(q)`.  The low observable is intrinsic to the level word:
it does not mention a body, a chosen cell, or a transport estimate.  Its
finite symmetric ratio alphabet is

```text
R_low={4/3,3/2,2,5/2,3,4,5,6}^{±1}.
```

Repeated levels have ratio one and belong to the separate same-level branch
already closed on the active body bank.  Distinctness is therefore part of
the reduction below, not a silently discarded tie.

## 2. A proved low-connected / high-connected dichotomy

**Lemma 1 (finite projective atlas).** If `L(q)` is connected, then `q` lies on
one of finitely many labelled primitive integer rays.

**Proof.** Choose a rooted spanning tree of `L(q)` and normalize the root to
one.  Across each oriented tree edge, child divided by parent lies in
`R_low`.  Thus every normalized coordinate is a product of at most five
members of a finite alphabet.  There are only finitely many normalized
labelled rational words.  Clearing denominators and dividing by the common
gcd gives a primitive integer word `P`; the original word is uniquely `gP`
for an integer `g>=1`.  QED.

This is the conceptual completeness argument used by the new edge-eight
referee: normalize a maximum low clique, then grow through a low spanning
tree.  Connectedness always supplies a frontier vertex, and its true branch
survives every exact adjacency/nonadjacency test against the assigned set.

**Lemma 2 (complement escape).** If `L(q)` is disconnected, then `H(q)` is
connected.

**Proof.** Vertices in different low components are adjacent in `H(q)`.  Two
vertices in the same low component both join any vertex in a different
component, so their high distance is at most two.  QED.

**Lemma 3 (common-cell Hunter certificate).** Fix one body-safe cell and let
`A_i` be the six reflected danger events there.  For every spanning tree `T`
on the six labels,

```text
mu(union_i A_i)
 <= sum_i mu(A_i) - sum_(ij in T) mu(A_i intersect A_j).
```

Indeed, root `T` and insert children after their parents.  The new part of a
child is contained in `A_child\A_parent`.  In THM-2941 notation the singleton
sum is `6/7+epsilon(E,q)`, so the cell closes whenever the tree credit is
strictly larger than `epsilon(E,q)`.

Together these lemmas give an unconditional structural split:

```text
L(q) connected
  -> finite labelled primitive shape P times one scale g;

L(q) disconnected
  -> H(q) connected
  -> a high spanning tree exists at one common cell.
```

More recursively, each low component has a finite internal primitive shape
times one component scale, while every cross-component pair is high.  A word
with `c` low components is therefore a finite shape chart over `c` coupled
integer scales, not an unstructured six-dimensional assignment.

## 3. What the reduction preserves and destroys

The source is a labelled six-level word.  The map records reduced pair ratios
and low components.  It preserves the existence of a low spanning tree, the
projective ray in a connected component, and connectivity of the high
complement.  It destroys:

- the ordered label-to-level placement;
- reflected transport orientation and target slope;
- the exact body-safe cell phase and its subdivision unit;
- edge gcds and their cycle consistency inside the single six-vector;
- the exact four-vertex contribution to the debt after an edge is fixed;
- scale residues and boundary/equality flags.

Those are not cosmetic fields.  In particular, quotienting an ordered ray by
the anharmonic `S_3` action loses stabilizer and reflection-sign data needed by
the located overlap.  The faithful atlas must retain the ordered placement.

For Tournament Analysis, the vertices are the six *clauses* at a fixed cell,
the pairwise observable is exact intersection mass, and a directed orientation
may rank the two deletion choices.  Ties must remain ties.  The preserved
target is the Hunter credit; the raw runner tournament is not intrinsic here.

## 4. Exact evidence and the sharp obstruction

The robust-edge-eight computation checks `652,688` exact certificate rows and
closes all `21` bodies in that block for arbitrary positive levels.  This moves
the THM-2941 arbitrary-level count from `2,421` to `2,442`, leaving `561` bodies
with at most seven robust edges.  The proof includes the first connected low
graph missed by the older maximum-clique chart and the genuine disconnected
free-center cylinders.

On the earlier `649`-body bank, the canonical upper-median two-star has two
useful finite signals:

- every injective assignment from levels `1,...,9` closes by its nine-edge
  maximum spanning tree (`55,606,320` exact assignments over spreads six,
  seven, and eight);
- when all nine two-star edges are low, there are exactly `11,856` ordered
  primitive rays, using only `27` integer levels and maximum primitive height
  `180`.  At scale one, all `649*11,856=7,694,544` exact located checks close.

These are **FINITE-EXACT**, not all-scale theorems.

The theorem replay and the two exploratory probes are reproducible by

```bash
python3 04-computation/lrc14_j7_reflected_upper_median_two_star_d6_d8_exact_thm2941.py
python3 04-computation/lrc14_j7_reflected_median_star_chord_scout_20260801.py --ray-audit --body-limit 649 --scale 1
python3 04-computation/lrc14_j7_reflected_median_star_chord_scout_20260801.py --m 1 --spread 9 --sample 128
```

The reason the high branch remains nontrivial is visible on the hostile word
`q=(7,5,4,3,2,1)` for `E=(1,2,3,4,6,12)`.  Its high graph is connected.  The
available cross-determinant transport assigns total maximum-tree floor zero,
but the actual common-cell high-tree credit exceeds the exact debt by

```text
51817283064/549491317375 > 0.
```

Thus the physical aggregate is healthy exactly where the edgewise lower
envelope forgets it.  At the opposite scale, gcd-one words such as
`(101,103,107,109,113,127)` show that a `1/g` tail alone does not control
unbounded primitive height.  A constant per-edge debt cap is also hopeless:
on the hostile `(1,12)` edge it leaves more than `145` million weak primitive
directions by a modest box.  The exact projected debt collapses the observed
failure set dramatically, but no absolute cutoff is yet proved.

## 5. A conditional finite-bank theorem

Suppose, for every labelled stalk edge at the chosen cell, that its exact
failure relation is a finite set of bounded points together with finitely many
homogeneous rays (or finitely many arithmetic subrays after recording scale
residue).  Then on any connected graph of such edges, the simultaneous
all-edge failure set is a finite absolute bank plus finitely many full
projective six-rays.

To see this, a bounded edge bounds one endpoint relative to the other and
propagates boundedness across the connected graph.  If every selected edge is
on a ray, root normalization and cycle consistency leave only finitely many
global ratios.  The conclusion is elementary; its antecedent is the open
mathematics.  Edge scales cannot be chosen independently because all edge gcds
come from one six-vector.

## 6. The next high-leverage target

The missing statement is no longer “find a good pair.”  It is one of the two
following located aggregate lemmas.

1. **Connected-low ray lemma.** For every ordered primitive shape in the
   finite low atlas, prove that a fixed common-cell Hunter tree closes for all
   integer dilations after finitely many explicitly controlled scale residues.
   The exact toothpick/ray recurrence in THM-2941 `(25l)--(25o)` suggests
   looking for a pair-overlap or tree-sum analogue, not sampling scales.
2. **Disconnected-low component lemma.** Normalize the finite shape inside
   each low component and prove that the high complete multipartite graph has
   a located spanning-tree credit exceeding the exact component-coupled debt
   for all component scales.

The cheapest decisive experiment is correspondingly structured: retain one
ordered six-vector, exact cell phase, scale residue, and projected debt;
compute common-cell edge weights; join failure relations as a cycle-consistent
CSP; and separate bounded primitive points from genuine arithmetic rays.  A
positive pattern should be proved through its functional form.  A negative
pattern should be stored as the first failed implication with its exact tree,
cell, debt, and boundary flag.

This reframing turns the remaining reflected wedge from “561 bodies times all
six-tuples” into a finite atlas of low-component shapes coupled by fewer scale
variables.  That is real compression, but the located Hunter inequality over
those scales remains OPEN.
