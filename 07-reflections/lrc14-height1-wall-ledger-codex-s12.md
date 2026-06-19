# LRC14 height-1 wall ledger

Codex 2026-06-19 S12.

The direct support-six Minkowski count did not finish LRC(14), but it taught us
where the false simplification was.  "Large span forces large coefficients" is
wrong because a large offset can sit exactly on a small signed subset-sum wall:

```text
21 = 1 + 2 + 3 + 4 + 5 + 6.
```

That obstruction looked dangerous only while it was being treated as part of
an infinite tail.  Once isolated as a finite wall ledger, it becomes tame.

## What the computation did

For each binding row `k=8,9,10`, I enumerated every primitive one-large row

```text
E = {0} union C union {M},
C subset {1,...,B(k)}, |C|=k-2, M>B(k),
```

with a height-1 type-II support-six relation touching `M`:

```text
+/- M + sum +/- e_i = 0
```

using at least five bounded-core terms.  Then I computed the exact rational
measure of `S7(E)` and compared it with the cap.

The full exact ledger cleared:

```text
k=8:  226046 rows, 0 over cap, worst 0.220151 < 0.381463
k=9:  250264 rows, 0 over cap, worst 0.372109 < 0.494256
k=10:  54173 rows, 0 over cap, worst 0.479762 < 0.604396
```

So the first anti-coset walls are resonant but harmless.

## The proof-shape change

Before this session, the support-six residual split was:

```text
finite bounded-spread check + support-six Minkowski tail
```

After HYP-2612/HYP-2613/HYP-2614/HYP-2616, the honest split is:

```text
finite bounded-spread check
+ finite low-height wall ledgers
+ relative signed support-six theta tail
+ no-scale cluster-collapse quotient
```

That is a better object.  It no longer asks a successive-minima theorem to do
work that is actually finite additive bookkeeping.

## What remains

The result does not prove LRC(14).  It covers height `1`, one large offset, and
bounded cores only.  The next finite target is height `2` one-large walls; after
that, multi-large low-height walls.  The analytic target is the HYP-2613/HYP-2614
relative signed permanent/theta estimate, with HYP-2615's signed-mass sequence
spine as the bookkeeping layer after the finite walls are removed.

## Tournament note

The useful vertices here were not runners.  They were offsets as carriers of
finite wall obligations, with incidence count as the observable.  The incidence
tournaments were transitive in all three binding rows.

This quotient preserves the support-six proof predicate "which offsets carry a
low-height wall" and destroys the phase-time geometry.  That loss is acceptable
for the finite ledger, but it is exactly why the final tail still needs signed
theta information rather than a pure incidence count.
