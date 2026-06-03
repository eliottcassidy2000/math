---
source: codex-2026-06-03-S584
status: exact bounded probe + proof-route synthesis
tags: [LRC, additive-circuits, four-term, three-term, hidden-folds, summand-graph, HYP-2115]
---

# Four-term structure as hidden virtual folds

The useful move is to stop asking whether 4-term structure is a second kind of
fold.  It is the same fold seen after completing the summand graph.

After rebasing over the HYP-2114 additive-circuit package, the slogan is
sharper:

```text
3-term folds carry hardness;
4-term energy is their translation-invariant shadow;
S584 measures the hidden virtual node that turns the shadow back into folds.
```

A visible 3-term relation

```text
v_c = v_a + v_b
```

is a literal phase fold:

```text
t v_c = t v_a + t v_b.
```

S577 verified the immediate shield picture: at the pinch clock for `v_c`, the
runner `c` is at an integer.  That is structure, not randomness.

A 4-term relation with no visible triple,

```text
a + b = c + d,
```

is different only because the sum node is missing.  Add the hidden node
`s=a+b=c+d` and the relation becomes

```text
(a,b,s) and (c,d,s).
```

So 3 encodes 4 by adding one virtual apex.

## What the computation says

S584 runs two exact checks.

First, it fixes a hidden sum `s` and deforms rows by choosing several disjoint
pairs `x+(s-x)`.  Visible 3-term folds are filtered out.  The pair-sum flavor
is held fixed, but exact `M` varies.  For `r=2` pair edges, the tested fixed
families all have at least `+0.1333` margin above `delta`; for `r=3`, at least
`+0.1071`.  This is the first important negative result: the 4-term label is
not itself a certificate, because deformation inside the same label changes
the geometry.

Second, it samples no-3-term rows and then adds the highest-multiplicity hidden
sum as a virtual runner.  The original rows remain safely above
`1/(k+1)`, including 4-term-rich rows at `k=9,10`.  But the virtual runner often
lands exactly at distance `0` at an original exact witness.  Example:

```text
V=(6,11,14,15,16,18,19,23,28)
s=34
pairs=(6,28),(11,23),(15,19),(16,18)
M(V)=0.2353
M(V+s)=0.1600
```

The information was present, but not as a real runner.  It was hidden in the
summand completion.

This complements the shifted-AP experiment in HYP-2114.  There the whole
4-term energy profile stays fixed while translation kills the visible folds and
`G` rises.  Here, instead of translating a whole AP, we add back the missing
sum node `s` for arbitrary no-3 rows and see the fold pressure reappear in
`M(V+s)`.

## Back to 3, then back to 4

The 3-term branch is certificate-shaped: a fold creates a shield clock, and
the proof should plug into divisibility, endpoint owners, or Phi.

The 4-term branch is diagnostic-shaped: it says "there is a virtual shield
clock here", but the shield is absent from the row.  That is why 4-term-rich
circuit-free rows are safe in both S584 and the HYP-2114 measure audits.  The
danger has been projected out.

The 4-term data becomes useful only if it is kept as a label:

```text
speed row -> pair-sum fibers -> hidden sum nodes -> virtual folds
```

Then the HYP-2113 stack can decide whether the hidden node becomes an actual
certificate gate.

## Hidden in deformation

The fixed-sum families are the deformation laboratory.  Keeping `s` fixed and
moving from `(1,s-1),(4,s-4)` to another pair selection preserves the same
4-term flavor.  But the exact LRC margin moves.

That means the deformation coordinate is not noise.  It is the missing payload
that raw additive energy forgets.  A proof of Lemma A should not try to show
"4-term energy is harmless" in a scalar way.  HYP-2114 already refutes that
route with shifted APs: 3-free rows can have maximal 4-energy.  The proof
should show:

1. if no hidden node becomes a real fold/shield, the remaining deformation is
   equidistributed enough to give margin;
2. if a hidden node does become real, route it to the visible 3-term machinery.

## Tournament Analysis Note

S584 uses explanatory lenses as tournament vertices:

```text
Phi_gap_gate
visible_3_fold_shield
hidden_4_virtual_sum
circuit_free_randomness
fixed_sum_deformation_fiber
raw_additive_energy
```

The observable is `(certificate_safety, proof_payoff, hidden_exposure, -cost,
maturity)`.  The tournament is transitive with one Hamiltonian path.  This is
not decorative: it ranks which pieces can certify and which merely expose
hidden structure.

The resulting order is exactly the proof instinct:

```text
Phi > visible 3-fold > hidden virtual 4-fold > randomness >
fixed deformation audit > raw additive energy.
```

Raw additive energy is the weakest object.  The hidden virtual fold is the
first object that remembers how 3 encodes 4.
