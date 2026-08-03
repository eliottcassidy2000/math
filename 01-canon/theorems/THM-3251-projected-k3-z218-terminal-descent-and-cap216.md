---
id: THM-3251
title: "Projected-k3 z218 terminal descent and composed cap216"
status: "PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED"
source: root/z218-terminal/2026-08-03
depends_on:
  - THM-3139-projected-k3-z225-terminal-and-z224-screen-double-layer-descent
  - THM-3230-projected-k3-z219-common-gcd-three-terminal-descent-and-cap218
  - THM-3242-projected-k3-z217-exact-status-annihilation
  - THM-2941-critical-seven-slot-scalar-wall-and-balanced-boundary
related:
  - THM-3245-pointed-divisor-median-cubes-saturation-band-no-go-and-z219-supplier-support
  - THM-2596-modular-free-factor-farey-gram-owner-cocycle
script: 04-computation/lrc14_j7_k3_z218_terminal_descent_composed_cap216_thm3251.py
output: 05-knowledge/results/lrc14_j7_k3_z218_terminal_descent_composed_cap216_thm3251.out
hash_basis: LF-normalized bytes
script_sha256: d92a6b825268d8fa7147ebcd5229bf8fa1d4cd2ffa3e9c1ec140024bd926832f
output_sha256: 66f7f6bb81992659d5cedf01aa2f22b12872810888a022460eb6a5ce1c05d674
semantic_sha256: fe1f0e0df1c652be79b1e9d1e964bc4db64f1930e25b8fb8a6f7d3c5f47ccedb
audit: >
  An independent hostile audit checked the final theorem against the exact
  companion and transcript, replayed normal and optimized metadata, verified
  all six pinned dependency hashes and semantic sentinels, confirmed that the
  script contains no optimization-sensitive Python assert, and reproduced the
  residual-index and overall semantic digests.  Normal, optimized, generated,
  and stored transcripts have the same LF-normalized hash with empty stderr.
  Every implication direction, quotient boundary, and ledger/cap transition
  was accepted without repair.
---

# THM-3251 -- projected-k3 z218 terminal descent and composed cap216

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## 1. Exact layer and the pointed divisor interval

In the pinned THM-2941 projected `k=3` necessary atlas, the complete layer
following THM-3230 is

```text
z1=218: 119 rows = 117 wall + 2 order,                    (1)
row SHA256 d416b484148e008e9f58273f8533cc9d05a171e5ba4963214d3fc73d7fc20bc7.
```

For every row body `E`, with `L=14 lcm(E)`, put

```text
g=gcd(218,L),                    first_d=L/g.             (2)
```

All `119` rows have

```text
g=2,                            first_d=L/2.              (3)
```

Thus their ambient pointed divisor interval is the two-point chain

```text
first_d | D | L          <-->          D/first_d in {1,2}. (4)
```

This is divisor metadata only.  In particular, `(4)` is not a physical
`C2` action, a carrier transition, or one of the modular free factors of
THM-2596.

## 2. Complete exact screen

The ray/status screen gives

```text
12,833 states = 5,375 crude + 6,522 exact-status + 936 residual. (5)
```

The wall and order pieces are separately

```text
wall:  12,827 = 5,370 + 6,521 + 936,
order:      6 =     5 +     1 +   0.                      (6)
```

All `6,522` status exclusions are verified through the inherited exact
rational Farkas path; none uses the direct shortcut.  The canonical screen
digest is

```text
d89494bb6dca6e5328765782965c0f61d402f8eb05b580dbcc65a1ea665630c3. (7)
```

Exactly `24` wall bodies survive, carrying `936` denominator masks.  Their
body/mask packet and ordered residual-index packet have SHA256 values

```text
e53f23cdce4d26addd82ab9a8321a031a70693c26f3b941ca6a0c86ce033790e,
fa8cfed8aebd06d4f9382a29d8c026d4333465457f770c5ae081bfe6144b0cd7. (8)
```

For every residual mask, if `D=lcm(ds)`, then

```text
D/first_d=2.                                               (9)
```

Hence the complete residual quotient census is

```text
(g,q,count)=(2,2,936).                                    (10)
```

Equation `(10)` says only that every residual occupies the top point of its
rowwise pointed two-chain, abstractly isomorphic to `{1,2}` (indeed `D=L`).
It creates no involution, common carrier, owner map, or phase transport.

## 3. Exact one-high terminal descent

Every residual is a wall row, so an actual residual assignment must have at
least one high slot.  For each of the `24` residual bodies, the
duplicate-permitting two-high union bound has strictly positive gap.  Thus at
most one high slot is possible, and every actual residual assignment has
exactly one.

The complete terminal census is

```text
24 bodies, 936 masks,
867 zero-high hostiles,
1,107 one-high cases over 63 low-label sets,
1,031 coarse closures + 76 exact finite-cell closures,
0 max-gap cases, 0 failures,
minimum finite-cell slack 1.                              (11)
```

Here the `867` zero-high rows are hostile relaxations discarded by the wall
gate, and the `63` low-label-set count is the sum of body-local types rather
than a global identification of labels.

The smallest two-high gap occurs at

```text
E=(1,5,6,9,12,14),
gap=730589/522413892 > 0.                                 (12)
```

For each one-high case, the high ray is replaced by its exact supremum and the
low slots range over a finite complete cell bank.  This enlarges the set of
physical assignments.  Every coarse certificate or exact complete-cell
certificate violates the required projected support cardinality, so all
`24` residual bodies close.  No max-gap fallback is used.

The terminal invariants are

```text
semantic SHA256 18646b406cbe94aa45a1a3bd1abaad338c089c6bf1bf6459bee4df60b03325b7,
tagged closure SHA256 6e09995ebdd03326120531c8d522efec6eb0483242303e6ba4d14cecdfde4272,
case-vector SHA256 689e2709c2de478ffe8b137b7fb1c9f1b90eac9bb304dc0106a18f097ec76da8.
```

## 4. Deferred composition and exact consequence

THM-3230 leaves the projected `k=3` ledger at `373,411`.  Removing the
complete `119`-row layer gives the intermediate consequence

```text
373411-119=373292,                    projected z1 <= 217. (13)
```

Only after `(13)` closes the intervening `z1=218` layer may the already-proved
THM-3242 annihilation of the disjoint eight-row `z1=217` layer be composed.
That gives

```text
373292-8=373284,                      projected z1 <= 216. (14)
```

The next occupied layer is

```text
z1=216: 480 rows = 447 wall + 33 order,
row SHA256 53db9e1d3df2cf2b0398847682d909da81705e43a53ae2553d102fd152337649. (15)
```

Equation `(15)` is a census only; no `z1=216` row is closed here.

## 5. Exact evidence and boundaries

The companion freshly derives and runs all `119` screen tasks from the pinned
atlas, verifies every exact Farkas witness, binds each of the `24` residual
mask banks separately, and then recomputes every terminal body directly
through the THM-3139 complete-cell engine.  Its semantic output is independent
of worker completion order.  The earlier normal checkpoint containing two
interleaved copies of each record is neither an input nor an ordering witness.

The theorem proves only a contiguous decrement in the projected `k=3`
necessary atlas.  In particular:

1. all screen and terminal state spaces are upper relaxations, so excluding
   them is lawful while reversing any implication would not be;
2. the common `g=q=2` quotient is not a physical `C2` action or modular
   free-factor carrier;
3. THM-3242 contributes its eight-row decrement only after `(13)`, never
   before the intervening layer closes; and
4. there is no projected `k<=1`, final-rung, or LRC(14) conclusion.

Normal and optimized replays LF-normalize byte-for-byte to the stored
transcript with empty standard error.  The independent hostile audit also
reproduced every dependency, metadata, quotient, and ledger check and found
no defect.
