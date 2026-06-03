---
id: HYP-2115
status: OPEN proof route; S584 supplies exact virtual-node evidence supporting HYP-2114
source: codex-2026-06-03-S584
related: [HYP-2114, HYP-2113, HYP-2112, HYP-2109, HYP-2108, HYP-2107]
---

# HYP-2115: four-term LRC structure is a hidden virtual-fold layer

## Claim

In the additive-circuit branch of LRC, a 4-term relation

```text
a + b = c + d
```

with no visible 3-term fold should not be treated as a hard obstruction by
itself.  It should be completed by a hidden virtual sum node
`s=a+b=c+d`.  Then the 4-term relation becomes two 3-term folds,
`(a,b,s)` and `(c,d,s)`, sharing the same virtual apex.

Thus 3-term structure encodes 4-term structure in the summand completion.  The
information is hidden only because `s` is not a runner in the original speed
row.

## Evidence

S577 showed that visible triples `v_c=v_a+v_b` are literal folds, while sampled
circuit-free rows stay safely above `delta=1/(k+1)`.

The HYP-2114 additive-circuit package sharpened the global split: shifted APs
keep 4-term energy fixed while 3-term folds vanish and `G` rises, and the
additive/Cprime synthesis finds safe measure depressed only by 3-term folds.
Thus 4-term energy is the translation-invariant shadow; visible folds carry the
hardness.

S584 adds the missing intermediate audit.

- Fixed hidden-sum deformation families, after filtering visible 3-term folds,
  have large positive margins.  For two pair edges with hidden sums
  `10,12,14,16,18,20`, the minimum margins over `delta=1/5` are at least
  `+0.1333`; for three pair edges, the minimum margins over `delta=1/7` are at
  least `+0.1071` in the tested families.
- In random no-3-term rows through `k=10`, 4-term-rich rows remain safe with
  margins such as `+0.1045` at `k=9` and `+0.1222` at `k=10`.
- The hidden fold is real: many no-3 rows have a virtual sum `s` that is at
  distance `0` at an original exact witness.  Example:
  `V=(6,11,14,15,16,18,19,23,28)` has hidden `s=34` with pair edges
  `(6,28),(11,23),(15,19),(16,18)`.  The original row has
  `M=0.2353`, margin `+0.1353`; adding the virtual runner gives
  `M(V+s)=0.1600`, augmented margin `+0.0691`, and the hidden runner is at
  distance `0` at the original witness.

So 4-term structure can carry virtual fold pressure, but it is not dangerous
until the virtual node becomes visible to the certificate system.  This is the
exact-node version of HYP-2114's shifted-AP principle.

## Proof Route

The right object is not raw additive energy.  It is a labelled summand
incidence layer:

```text
original speeds -> pair edges -> hidden sum nodes -> virtual 3-folds
```

Rows split into three proof branches.

1. Visible 3-fold: use the S577 fold/shield mechanism and the divisibility
   machinery.
2. Hidden 4-fold only: keep the hidden sum node as a middle label, but do not
   certify from additive energy alone.  Use HYP-2114's gap-structural Lemma A
   route, with hidden virtual pressure as a diagnostic.
3. Dense hidden-sum fibers: push the pair-sum incidence graph into the HYP-2113
   certificate stack.  If the hidden node becomes a real cover/fold gate, pass
   it to fold/shield, Phi, CRT, or endpoint-owner gates.

## Assumption Challenge

Tournament or automaton vertices need not be runners.  For this branch the
candidate vertices are original speeds, pair edges, hidden sum nodes, fixed
hidden-sum deformation fibers, virtual fold clauses, and proof obligations.

The quotient preserves: which hidden virtual sums encode which 4-term
relations, and whether those hidden sums exert pressure at exact witnesses.

The quotient destroys: the individual deformation coordinate inside a fixed
pair-sum fiber.  S584 shows this matters, because rows with the same hidden-sum
flavor can have different exact `M`.

The challenged assumption is that high 4-term additive energy is intrinsically
hard.  In this model, high energy without visible triples is usually a shadow
of absent virtual folds, not a counterexample-shaped obstruction.
