---
id: THM-3041
title: "Projected k3 z243 first-below-floor cardinality descent"
status: >
  CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE AUDIT.
  The candidate reuses THM-3033's hash-pinned exact evaluator on the disjoint
  bank of 151 first-below-high-floor z1=243 rows.  No projected cap or ledger
  change is active at candidate status.
source: modular-farey-quartic-bridge-2026-08-01
depends_on:
  - THM-3033-projected-k3-z246-to-z244-descent-and-z243-high-floor-addendum
  - THM-2984-projected-k3-signed-ray-attainment-and-unit-phase-gate
  - THM-2941-critical-seven-slot-scalar-wall-and-balanced-boundary
script: 04-computation/lrc14_j7_k3_z243_first_below_floor_cardinality_descent_thm3041.py
output: 05-knowledge/results/lrc14_j7_k3_z243_first_below_floor_cardinality_descent_thm3041.out
script_sha256: 91f0b1182c295412706a63c1a92c665be67ca7217b3be28c753090d087569cca
output_sha256: c7f9779633d241d7bb55f29125763144f08cce943014535c2ed2527fe4e0e3fa
semantic_sha256: 9eaffa1eff6004dd7c2a7d12cdede4a7f53d0294cb8aefee2cde7bc3faa35dee
hash_basis: LF-normalized bytes
---

# THM-3041 -- projected k3 z243 first-below-floor cardinality descent

**CANDIDATE + VERIFIED-EXACT; AWAITING INDEPENDENT HOSTILE AUDIT.**

## Candidate statement

In the lossless projected `k=3` scalar atlas inherited from THM-3033, every
one of the `151` `z_1=243` rows whose first label lies strictly below its
integer high floor is empty.

This bank is disjoint from the three `z_1=243` high-floor addendum rows already
proved in THM-3033:

```text
(1,2,3,8,10,12),
(1,3,4,8,10,12),
(2,4,6,8,12,14).                                      (1)
```

Their exact `(states,crude,status,residual)` profiles are respectively

```text
(1,1,0,0), (13,13,0,0), (1,0,1,0).                   (2)
```

No projected-cap or necessary-ledger update is active while this theorem is a
candidate.  If the proof and evidence below pass independent hostile audit,
the prospective arithmetic is the disjoint decrement

```text
375251 - 197 - 3 - 151 = 374900,                      (3)
```

and the projected cap would descend from `z_1<=243` to `z_1<=242`.  Equation
`(3)` is an audit target, not current canon at candidate status.  The result is
only a projected necessary-sector statement and is not LRC(14).

## 1. Frozen universe

The companion pins the LF-normalized source, output, and semantic hashes of
promoted THM-3033.  It parses the same canonical atlas and selects exactly the
rows satisfying

```text
z_1=243,
z_1 < floor(13L/132)+1.                               (4)
```

The exact universe audit is

```text
all z1=243 rows                         154
already-proved high-floor rows           3
first-below-floor rows                  151
distinct bodies                         151
L range                         2520..5045040.         (5)
```

The row selector also verifies the high-floor formula and rejects duplicate
body keys.  It changes no state generator or exclusion criterion inherited
from THM-3033.

## 2. Exact common-status screen

For every selected body, the unmodified THM-3033 evaluator regenerates the
attained quotient-state universe and partitions it into crude capacity,
common-status, and residual states.  The exact totals are

```text
states     79389
crude      30596
status     44159
residual    4634.                                      (6)
```

For each crude state the companion reconstructs the exact target-minus-
capacity gap and checks that it is positive.  For every common-status state it
rebuilds the exact marginals, permitted capacities, and residue-load
histogram, then verifies the returned Farkas identity and strict direction
over `QQ`.  All

```text
44159 / 44159                                           (7)
```

checks pass.  The screen therefore closes `88` bodies outright and leaves
`63` bodies, carrying `4634` residual states, for the terminal proof.

The screen is an upper-relaxation exclusion.  It does not infer literal
carrier realizability from a surviving relaxed state.

## 3. Terminal directions

Fix one of the `63` residual bodies.  THM-2941's projected-safe-arc wall first
forces at least one subsequent label to reach the high threshold.  The exact
duplicate-permitting two-high upper bound then lies strictly below the scalar
floor on every body:

```text
two-high positive gaps = 63 / 63.                       (8)
```

Consequently a genuine residual packet must have exactly one high label.  The
companion enumerates both signs for the two low labels and the height-free unit
ray above the high threshold.  It obtains

```text
zero-high scalar stress passes = 4405
one-high cases                 = 5734.                  (9)
```

The `4405` zero-high cases in `(9)` are not terminal evidence.  They are
deliberate hostile controls excluded only by the inherited THM-2941 wall.
Counting them as scalar exclusions would reverse the proof direction.

For a one-high case of exact denominator `d`, let `S` be the residues of
whole closed body cells safe for the body, the fixed first label, and the two
low labels.  THM-2984 says the strict high-danger band can occupy at most

```text
ceil(d/7)                                               (10)
```

residue classes after any primitive unit scaling and translation.  Every one
of the `5734` cases satisfies

```text
|S| > ceil(d/7).                                        (11)
```

Thus no translated high-danger band contains the safe image, and every
one-high case closes.  The exact totals are

```text
cardinality cases             5734
shape/max-gap cases              0
failed one-high cases             0.                   (12)
```

The stronger affine maximum-gap fallback is retained as a hostile control.
It is compared against literal affine cyclic-block containment for every
subset and every unit modulo `d=1,...,12`, giving `8190` exact checks,
including empty, singleton, and wraparound boundaries.  It is not used in
the candidate closure.

Combining the `88` screen closures with the `63` terminal closures accounts
for all `151` rows, with no survivor.

## 4. Literal infeasibility, not a zero-mass inference

No step argues that `mu(Safe)=0` makes `Safe` empty.  Here an "empty row"
means that no literal pointwise danger cover, hence no LRC counterexample, can
have that projected row.  A nonempty zero-mass safe set would already defeat a
pointwise cover and therefore would close rather than survive the row.

The two proof branches certify literal infeasibility for different reasons.

1. The projected state atlas is lossless for a hypothetical literal cover.
   Such a cover necessarily produces one attained state and, after forgetting
   correlations, a nonnegative common-status table with the recorded exact
   marginals.  A positive crude capacity gap or a strict rational Farkas dual
   proves that even this enlarged table is infeasible.  Since pointwise cover
   implies the relaxed a.e. marginal constraints, failure of the relaxation
   rules out the literal cover; it is not a conclusion drawn from the measure
   of its safe complement.
2. In the terminal branch, `(11)` is pointwise on every complete-cell
   coordinate `y`.  For each `y`, THM-2984 supplies an actual closed-cell
   residue whose point is safe for the body, fixed low labels, and high label.
   Hence the projected safe section is literally the full circle.  A literal
   aligned completion would force that full section into the proper open
   danger union, which is impossible.  This is the pointwise section argument
   of THM-2984 `(23c)-(23d)`, not a positive-Haar-mass surrogate.

This distinction is load-bearing at tight boundaries where a nonempty safe
set can have Haar mass zero.

## 5. Evidence and scope

Reproduction from the repository root:

```text
python 04-computation/lrc14_j7_k3_z243_first_below_floor_cardinality_descent_thm3041.py --processes 8 --checkpoint-dir FRESH_NORMAL_DIR
python -O 04-computation/lrc14_j7_k3_z243_first_below_floor_cardinality_descent_thm3041.py --processes 8 --checkpoint-dir FRESH_OPT_DIR --output FRESH_OPT_OUTPUT
```

The evaluator writes atomic per-row checkpoints whose provenance includes the
current script hash and all three frozen THM-3033 dependency hashes.  A script
or dependency edit therefore cannot silently reuse an old row.  Independent
fresh normal and optimized checkpoint directories are required for the final
evidence claim; checkpoint-only rendering is not counted as optimized
verification.

Fresh normal and optimized runs used disjoint checkpoint namespaces.  Their
`129986`-byte, `229`-line transcripts are byte-identical to each other and to
the stored output.  The frozen hashes are

```text
script LF SHA-256
  91f0b1182c295412706a63c1a92c665be67ca7217b3be28c753090d087569cca
output LF SHA-256
  c7f9779633d241d7bb55f29125763144f08cce943014535c2ed2527fe4e0e3fa
semantic SHA-256
  9eaffa1eff6004dd7c2a7d12cdede4a7f53d0294cb8aefee2cde7bc3faa35dee
```

The output freezes every row record, all residual states, both terminal
directions, failure payloads, the `8190`-case affine control, and a semantic
digest.  No solver-selected basis is part of that semantic digest.

The promoted dependency pins are

```text
THM-3033 source LF SHA-256
  c063e1d94e1f30c4ccd63acfabb172988d23d96f680a961d9042b3f11862fa79
THM-3033 output LF SHA-256
  d012f84d8e380019b4da8efbe9afd6695079c585f47e4aa2c474828ca62014e4
THM-3033 semantic SHA-256
  3f9945d333b13adef8aa8e7b162960bfde97f1e5a77945edcc735013d8ce4231
```

This theorem neither proves a literal LRC row exists nor closes `k<=1`, the
final rung, or LRC(14).  It may not be cited as proved before an independent
hostile audit promotes its status.
