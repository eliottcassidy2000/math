---
id: HYP-2723
title: LRC14 depth-law convex order and generated-word compatibility are two different layers
status: SYNTHESIS; resolves concurrent HYP-2722 namespace collision
source: codex-2026-06-21-S71
depends_on:
  - HYP-2722
  - HYP-2721
  - HYP-2702
  - HYP-2698
  - THM-558
  - THM-556
related:
  - HYP-2719
  - HYP-2718
  - THM-557
  - OPEN-Q-108
---

# HYP-2723 - Depth Law / Generated Word Two-Layer Split

## Claim

The concurrent S71 work produced two different `HYP-2722` threads.  They are
not contradictory; they live at different layers:

```text
full actual row
  -> depth-law / occupancy convex order
  -> even Krawtchouk bands B2,B4,B6.

decorrelated generated context
  -> fixed-x miss-zeta product word
  -> compatibility filter before q0 scalarization.
```

The mac-mini S10 thread says the global `measS7` crux is a full-row
convex-order occupancy extremality, not a character/Gauss-sum problem.  The
codex S71 thread says the HYP-2721 abstract atom cone is too loose unless an
atom move comes from an HYP-2698 generated miss-zeta word.

Both are needed, but they should not be merged into a single scalar theorem.

## Full-Row Layer

Incoming mac-mini S10 reframes the actual row problem:

```text
pi_E(h) = meas{x : exactly h of 7 sectors are hit},
measS7(E) = pi_E(7) = E[1_{N=0}],
N = number of empty sectors.
```

Bonferroni/Krawtchouk truncations are reads against the depth law:

```text
B_J = E[g_J(N)].
```

Evidence in the pulled commit:

```text
even B2,B4,B6 are consec-max in the tested k=8 box;
odd B3 is dirty;
the crux is convex-order occupancy extremality, not a Gauss-sum or
support-3 character extremality.
```

This is the global actual-row target.

## Generated-Context Layer

Codex S71 tests the decorrelated product-word compression layer:

```text
Z_context,x(A) = product_i Z_i,x(A).
```

On the HYP-2702 sparse-tail frontier, after OR-convolving generated contexts
with candidate clusters, exact atom deltas satisfy:

```text
318 tests
q0_failures = 0
min W1+W2 leakage = 144154/63487
U4/q0 nonpositive = 0
tail45/q0 negative = 0
```

But signed `B2/q0` is not always positive:

```text
B2/q0 nonpositive = 59/318,
min B2/q0 = -1135/282.
```

Therefore full-row even `B2` cleanliness does not transfer directly to
decorrelated generated-context compression.  The generated-word proof should
use:

```text
low-factorial leakage,
B4=U4 / THM-558 transfer tax,
and relation-support handoff,
```

not a blanket `B2` positivity statement.

## Proof Split

The safe proof architecture is:

1. **Full-row convex order.**  Prove actual-row `measS7` extremality through
   the depth law and even Krawtchouk/convex-order bands.
2. **Generated-word compatibility.**  Prove HYP-2698/HYP-2702 context
   compression through fixed-`x` miss-zeta product words.  Use singleton
   death-chain exclusion first, then coherent context merging.
3. **Atom-boundary handoff.**  After product-word compatibility, evaluate the
   HYP-2721 `Q_0` boundary atom.  Cheap abstract LP atom moves are admissible
   only if they survive generated-word compatibility.
4. **Relation support.**  Use HYP-2719 to route support-3/additive-triangle
   packets and high-support tails; do not import OCF reverse-pair cancellation.

## Assumption Challenge

Rejected assumptions:

```text
one HYP-2722 scalar can cover both layers;
full-row B2 positivity implies generated-context B2 positivity;
support-3/Gauss-sum arithmetic is the cover crux;
abstract six-sector atom-cone feasibility implies LRC generated-word feasibility.
```

Productive vertex sets differ by layer:

```text
full row: depth-law states h or empty-count states N;
generated context: miss-zeta product words and coherent context partitions;
atom boundary: missed-count atoms Q_t and low-factorial leakage;
relation layer: support-size carrier packets.
```

The quotient must preserve the predicate appropriate to its layer before
scalarizing.

## Status

HYP-2723 is a synthesis and namespace repair.  It does not claim a new proof,
but it prevents two useful S71 results from being conflated.  The next concrete
step for the codex branch is the singleton generated-word exclusion lemma for
HYP-2722, with `B4=U4` and `W1+W2` leakage as the measured witnesses.
