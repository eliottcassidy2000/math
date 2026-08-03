---
id: THM-3230
title: "Projected-k3 z219 common-gcd-three terminal descent and cap218"
status: "PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED"
source: root/projected-k3-z219/2026-08-02
depends_on:
  - THM-3218-projected-k3-z220-valuation-product-terminal-descent-and-cap219
  - THM-3207-projected-k3-z221-coprime-terminal-descent-and-cap220
  - THM-3179-projected-k3-z222-composite-divisor-square-terminal-descent-and-cap221
  - THM-2941-critical-seven-slot-scalar-wall-and-balanced-boundary
related:
  - THM-2197-scalar-chord-coverage-has-a-boolean-deficiency-quotient
  - THM-2596-modular-free-factor-farey-gram-owner-cocycle
  - THM-2597-six-vertex-bicycle-modular-abelianization-cycle
script: 04-computation/lrc14_j7_k3_z219_common_gcd_three_terminal_descent_cap218_thm3230.py
output: 05-knowledge/results/lrc14_j7_k3_z219_common_gcd_three_terminal_descent_cap218_thm3230.out
hash_basis: LF-normalized bytes
script_sha256: 0b14b1a86e2a08aae20433ee12d601128090f04658c0934fc6b263fb31d9e723
output_sha256: d8a4e87382da55879ee5f517f79277e4c3913823d5f968931ded22f84ace3bcf
semantic_sha256: 999e0186c706441074e50ca1a2e0d689a9f19c6fa2c13137abc730eb51f92a49
audit: >
  An independent hostile audit rebuilt the complete screen under normal and
  optimized interpreters from the pre-candidate checkpoint and recomputed all
  ten terminal bodies directly through THM-3139. It reproduced every census,
  Farkas split, quotient, digest, minimum gap, ledger decrement, and next-layer
  hash. The canonical normal and optimized replays byte-match the stored output.
---

# THM-3230 -- projected-k3 z219 common-gcd-three terminal descent and cap218

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## 1. Exact layer and the common prime-three fibre

In the pinned THM-2941 projected `k=3` necessary atlas, the complete layer
after THM-3218 is

```text
z1=219: 16 rows = 15 wall + 1 order,                       (1)
row SHA256 6a39c83d47b8080fd52e64479b886f1b6a887150f32cfbf30a676fbb0aaa54ba.
```

For a row with body `E` and ruler `L=14 lcm(E)`, put

```text
g=gcd(219,L),                    first_d=L/g.              (2)
```

Unlike the product-of-valuation-chains layer in THM-3218, all sixteen rows
have

```text
g=3.                                                        (3)
```

Thus the ambient pointed divisor interval is only

```text
first_d | D | L          <-->          q=D/first_d in {1,3}. (4)
```

This prime-three collapse is arithmetically useful, but it is not yet a
`C3` action: `{1,3}` is a pointed two-vertex divisor chain, not a
three-cycle or a torsor.

## 2. Complete exact screen

The ray/status screen gives

```text
2,107 states = 1,022 crude + 661 exact-status + 424 residual. (5)
```

The wall and order pieces are separately

```text
wall:  2,106 = 1,021 + 661 + 424,
order:     1 =     1 +   0 +   0.                           (6)
```

All `661` status exclusions are verified by the inherited exact rational
Farkas path; none uses the direct shortcut.  The canonical digest of the
screen rows through field 18 is

```text
868d5cce16d5f704f783987b4ecba4e55fda009bce1d48611bb8f3dec647e248. (7)
```

Exactly ten wall bodies survive, carrying `424` denominator masks.  Their
body/mask packet has SHA256

```text
69918f14f67d818749581af98f8671d7567da951bf0be52e3274fa558a186881. (8)
```

For every one of those `424` masks, if `D=lcm(ds)`, then

```text
D/first_d=3.                                                (9)
```

Thus the complete quotient census is the singleton

```text
(g,q,count)=(3,3,424).                                     (10)
```

Equation (10) is stronger than merely saying that the top divisor is
available: there is no `q=1` residual at this layer.  It still does not
construct a transition, common carrier, order-three motion, or phase map.

## 3. Exact one-high terminal descent

Every residual is a wall row, so the wall gate requires at least one high
slot.  For each of the ten residual bodies, the duplicate-permitting
two-high union bound has strictly positive gap.  Hence at most one high slot
is possible and every actual residual assignment has exactly one.

The complete terminal census is

```text
10 bodies, 424 masks,
413 zero-high hostiles,
424 one-high cases over 13 low-label sets,
411 coarse closures + 13 exact finite-cell closures,
0 max-gap cases, 0 failures,
minimum finite-cell slack 1.                                (11)
```

The smallest two-high gap occurs at

```text
E=(1,5,7,11,12,13),
gap=176513570591453/67539934271289480 > 0.                 (12)
```

For each remaining one-high case, the high ray is replaced by its exact
supremum and the low slots range over a finite complete cell bank.  This
enlarges the set of physical assignments.  Every complete-cell certificate
violates the required projected support cardinality, so all cases close.  No
max-gap fallback is used.

The terminal invariants are

```text
semantic SHA256 ae07780990ec507a6a3ff883363dd223af4f11cafdc48543d8e5f616c9896049,
tagged closure SHA256 14b8e92a82cfde0daddde1859237f0ec724858a9b993e51208fae1d1b18a0641,
case-vector SHA256 127d5e241f7c6aa4493b158a316890b9d8455205f793dafbdd56725de2c516f5.
```

The closure digest deliberately omits the chosen duplicate-gap maximizing
witness; the semantic digest likewise omits that noncanonical witness.

## 4. Consequence

THM-3218 leaves the projected `k=3` ledger at `373,427`.  Removing the
complete sixteen-row layer gives

```text
373427-16=373411,                    projected z1 <= 218.   (13)
```

The next occupied layer is

```text
z1=218: 119 rows = 117 wall + 2 order,
row SHA256 d416b484148e008e9f58273f8533cc9d05a171e5ba4963214d3fc73d7fc20bc7. (14)
```

Equation (14) is a census only; no `z1=218` row is closed here.

## 5. Exact evidence and boundaries

The companion recomputes all sixteen screens from the pinned atlas, verifies
every raw Farkas certificate, binds each of the ten residual mask banks,
checks `(2)--(10)`, and recomputes all ten terminal bodies directly through
the THM-3139 complete-cell solver.  Normal and optimized modes must agree
byte-for-byte with the stored transcript.

The theorem proves only a projected-`k=3` necessary-atlas decrement.  In
particular:

1. the common `g=3,q=3` state is not a physical `C3` orbit;
2. it supplies neither the invertible order-three map nor the common carrier
   required by THM-2596's modular/Farey picture;
3. it gives no projected `k<=1`, final-rung, or LRC(14) conclusion; and
4. the one-high terminal bank is an upper relaxation, so its successful
   exclusion is lawful, while the converse interpretation would not be.

The independent hostile audit rebuilt the screen from its pre-candidate
checkpoint and the terminal bank directly through THM-3139.  Its normal and
optimized replays reproduce the stored transcript byte-for-byte.  Equations
`(1)--(14)` therefore prove the stated cap.
