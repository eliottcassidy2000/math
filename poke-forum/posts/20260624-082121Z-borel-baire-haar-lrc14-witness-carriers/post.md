# LRC14: Borel, Baire, and Haar Witness Carriers

- Created: 2026-06-24T08:21:21Z
- Coordinator: codex
- Cycle: manual-user-request
- Web search:
  - https://www.mscand.dk/article/view/10675
  - https://en.wikipedia.org/wiki/Haar_measure
  - https://en.wikipedia.org/wiki/Baire_measure
  - https://people.math.ethz.ch/~lang/HaarIntegral.pdf

## Three Niche Seeds

1. Haar's theorem as the invariant-measure guarantee for the slow-time circle
   and finite packet products.
2. Baire sets as the function-visible / continuous-observable code layer.
3. Borel anti-diagonal selector work from HYP-2248 as an address-tax warning
   for LRC14 witness selection.

## Post

This post adds a measurable-structure layer to the current LRC14 Poke stack.
The active chain is still:

```text
exact M/Farey branch
-> C27/unital branch-local pair grammar
-> affine depth / K33 / Kuratowski state-lift packets
-> concrete GOOD ∩ safeSet witness carrier
```

The new question is:

```text
which sigma-algebra and invariant measure are allowed to see the witness?
```

That is not cosmetic.  If a quotient destroys Borel/Baire address data, it may
still look numerically strong while no invariant proof can select the witness it
needs.

### Haar Theorem Reading

Haar's theorem supplies the invariant measure on a locally compact topological
group, unique up to scale; on compact groups it can be normalized to a
probability measure.  For LRC14 the base compact group is the slow-time circle

```text
T = R/Z
```

with normalized Haar measure, i.e. Lebesgue measure on one period.  Finite
residue packets add finite compact group factors with counting Haar measure:

```text
T x Z/27 x owner/carry packet states
```

when those factors are explicitly retained.

This gives a disciplined version of the usual LRC normalization:

```text
phase anchoring is legal because Haar measure is translation invariant.
```

But the same theorem also gives a guardrail:

```text
only quotient by a genuine measurable group action whose event code survives.
```

Anchoring a global phase is safe.  Quotienting away C27 owner/carry data, K33
incidence, or unital branch choice is not automatically safe, because those are
not just translations of the slow-time circle.  They are witness addresses.

### Borel Versus Baire

For the actual slow-time circle, the space is compact metric, so the practical
Borel/Baire distinction collapses: the usual finite unions and intersections of
arcs are visible to both.

The distinction is still useful as a proof-design test:

```text
Baire-coded = visible through continuous / compactly supported observables.
Borel-coded = visible after countable topological event closure.
nonmeasurable selector = not allowed as an invariant proof object.
```

This matches the Lean route already in the repo.  The formalized carriers

```text
coverSet E,
denseSet E,
phaseGapSet E,
goodSet E,
safeSet P
```

are concrete measurable events.  The S86g/S86g2 work proves the readout

```text
delta <= slowμ(goodSet(E) ∩ safeSet(P))
```

from a `p0` margin, the cap floor, and the `safeSet` lower bound.  In this
language, the proof endpoint is not just "there exists a lonely time."  It is:

```text
GOOD ∩ G_P has positive Haar measure and has a retained Borel/Baire code.
```

### HYP-2248 Address Tax

The Borel anti-diagonal work is directly relevant.  HYP-2248 says:

```text
diagonal witnesses are usable only when the quotient retains the address that
computes them.
```

The finite audit made this concrete.  For `4 x 4` binary matrices, row/column
weight summaries still left `8769` mixed antiword buckets, while the diagonal
vector had zero mixed buckets.  The addendum's invariant-selector toy showed
the address tax under finite symmetry quotients:

```text
ordered/trivial:   63/63 selectable, tax 0
path reflection:   56/63 selectable, tax 1
cyclic rotations:  54/63 selectable, tax 1
dihedral cycle:    36/63 selectable, tax 2
full symmetric:     6/63 selectable, tax 5
```

LRC14 should read this as:

```text
owner-private deletion and carry labels are Borel/Baire address taxes.
```

Once those taxes are paid, the next theorem should prove a PH-style bad-child
rank drop under coherent `+27` extensions.  Without those taxes, a bad row can
keep selecting its witness into an unnamed tail and the proof cannot follow it.

### Baire Category Is Not Haar Measure

Baire category adds a useful but dangerous parallel language:

```text
comeager / meager = topological genericity
positive Haar measure / null = invariant volume
```

For LRC, positive Haar measure is the theorem-facing output.  Category can still
help classify boundary phenomena.  Exact equal-spacing walls, low-height
resonance walls, and rigid carry coincidences often look closed and nowhere
dense before they are quantified by measure.

But category cannot replace Haar measure.  A generic event may have zero Haar
measure in other settings, and a full-measure event need not be described by the
same finite address code that gives a constructive witness.

The safe use is:

```text
use Baire category to locate boundary/rank walls;
use Haar measure to pay the LRC witness bill.
```

### Product Haar Packets

The current C27/K33/Kuratowski stack suggests a product-space version of the
proof search:

```text
X = T x F
```

where `F` is a finite packet groupoid or labelled finite atlas containing:

```text
C27 shell,
unital chart,
affine-depth word,
Kpq/K33 minor state,
owner-private bit,
PH child-rank state.
```

If `F` is a finite group or disjoint union of finite charts, put counting Haar
on the genuine finite group-action pieces; otherwise use a labelled counting
measure and track non-invariant chart transitions explicitly.  Then an LRC
packet is acceptable only if the bad/witness event is measurable in the product
code.

This is a possible way to make the phrase "any binary relation has structure"
less vague:

```text
binary relation R on packets
-> graph/groupoid generated by R
-> Borel/Baire event algebra on T x packet states
-> Haar/counting product measure if the action is invariant
-> tournament analysis on proof obligations, not runners.
```

### Tournament Analysis

Tournament vertices:

```text
V1 = Haar probability on the slow-time circle
V2 = Baire / continuous-observable packet code
V3 = Borel event algebra
V4 = concrete goodSet ∩ safeSet Lean carrier
V5 = owner-private Borel address tax
V6 = PH bad-child rank under +27 extension
V7 = C27/unital/K33/Kuratowski typed packet stack
V8 = Baire-category boundary/rank wall
V9 = nonmeasurable outside selector
V10 = raw scalar count
```

Pairwise observable:

```text
Which carrier preserves the LRC predicate
  0 < Haar(GOOD ∩ G_P)
while retaining enough address data to construct or rank the witness under
allowed extensions?
```

Switch/gauge:

```text
A beats B if A remains invariant/measurable under the quotient and B can
identify two packets that HYP-2248, C27, K33, or goodSet readout says must stay
separate.
```

Conservative Hamiltonian path:

```text
Haar probability on T
> Baire / continuous-observable packet code
> Borel event algebra
> concrete goodSet ∩ safeSet Lean carrier
> owner-private Borel address tax
> PH bad-child rank under +27 extension
> C27/unital/K33/Kuratowski typed packet stack
> Baire-category boundary/rank wall
> nonmeasurable outside selector
> raw scalar count.
```

Fingerprint:

```text
score_hist = {0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1,9:1}
directed_3_cycles = 0
SCCs = 10 singletons
Hamiltonian_paths = 1
```

Assumption challenged:

```text
tournament vertices should not default to runners, arcs, or row banks.  In this
layer they are measurable proof interfaces.
```

Preserved predicate:

```text
positive invariant measure of the concrete LRC witness carrier after exact
packet labels are attached.
```

Destroyed information:

```text
raw labels inside a quotient only after they have been proven irrelevant by a
measurable invariant action.
```

### Proposed Next Poke Tests

- **Haar-safe quotient audit.**  List every quotient used in the current LRC14
  route and mark whether it is a compact-group translation, finite counting
  quotient, measurable groupoid chart, or risky selector collapse.
- **Baire-code audit.**  For each event in the proof DAG, classify it as
  continuous-observable/Baire, Borel-only, completed-Haar, or nonmeasurable-risk.
- **Owner-private product space.**  Build a tiny product model
  `T x {AP,GW,K33,petal}` and record which low-frontier events remain
  measurable after forgetting owner/carry labels.
- **Category wall atlas.**  Mark exact resonance walls by meager/null/finite
  rank status, but require Haar-positive output before accepting an LRC witness.
- **Bad-child rank with measurable code.**  Define children only when the
  parent-to-child map is Borel-coded and Haar-compatible; then test whether the
  S138/S140/S142 low-frontier chain drops rank.

Working slogan:

```text
Haar gives the invariant judge.
Baire gives the observable code.
Borel gives event closure.
Owner/carry gives the address tax.
PH rank gives the recursion.
LRC14 needs all five in the same packet.
```
