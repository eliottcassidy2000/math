---
id: THM-3264
title: "Projected-k3 z216 low-cost gcd8 seventeen-row terminal descent"
status: "PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED"
source: root/z216-gcd8-low-cost/2026-08-03
depends_on:
  - THM-3139-projected-k3-z225-terminal-and-z224-screen-double-layer-descent
  - THM-3261-projected-k3-z216-unique-gcd18-terminal-descent
related:
  - THM-3245-pointed-divisor-median-cubes-saturation-band-no-go-and-z219-supplier-support
  - THM-3259-charge-paired-modular-free-factor-lift-and-idempotent-localization-obstruction
script: 04-computation/lrc14_j7_k3_z216_low_cost_gcd8_terminal_descent_thm3264.py
output: 05-knowledge/results/lrc14_j7_k3_z216_low_cost_gcd8_terminal_descent_thm3264.out
script_sha256: d9841a86850c6609e62d0522b50ea38722cd5850f7173374a3b143ef577e0e3e
output_sha256: 2eebe5f7acb2d1c02fa126ee03166cd274ed209cc52d4b8a729a8fbc5f0a9782
semantic_sha256: 64604458a454f7da3468b07ec5697a6dd62ac4413625a92dbd8e8ffff15e1a7e
hash_basis: LF-normalized bytes
audit: >
  An independent hostile audit reparsed the complete 6,060-row atlas and
  rederived the z216 census and gcd strata, the intrinsic cost selection,
  all 1,746 screened states and 895 exact Farkas exclusions, the three-body
  19-mask top-q8 residual, all three strict two-high gaps, and all 28
  one-high terminal cardinality certificates. It separately checked the
  dependency and semantic pins, ledger arithmetic, stopping boundary, and
  assertion-independent normal/-O/stored byte parity with empty stderr.
---

# THM-3264 -- projected-k3 z216 low-cost gcd8 seventeen-row terminal descent

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## 1. A complete intrinsic sublayer

In THM-3261's first occupied `z1=216` layer, exactly `19` rows have

```text
gcd(216,L)=8.                                               (1)
```

Let `r` be the component count printed in the canonical projected atlas and
use the inherited screen-work invoice

```text
C(E)=L(E) r(E).                                             (2)
```

Exactly `17` of the `19` rows satisfy `C(E)<=2,000,000`.  Their atlas indices
are

```text
8,23,39,57,64,66,115,142,150,152,197,272,277,300,304,306,338, (3)
```

and their total invoice is `12,050,080`.  The two excluded rows are

```text
238: E=(1,8,10,11,13,14),
370: E=(2,8,10,11,13,14),
C(E)=19,059,040 in each case.                                (4)
```

Thus `(3)` is a fully quantified finite sublayer, not a claim about the two
large rows in `(4)` or any other gcd stratum.

## 2. Complete exact screen

Fresh exact screening of all `17` rows gives

```text
1,746 states = 832 crude + 895 exact-status + 19 residual,   (5)
direct Farkas=0, inherited exact Farkas=895.
```

The ordered screen digest is

```text
7ad0a1a8036fefa37a9035d2d2b319206b170dcbea2f7f86b7cef96cd7c1b3f1. (6)
```

Only three bodies survive:

```text
index 64:  E=(1,2,8,10,11,14),  8 masks,
index 152: E=(1,4,8,10,13,14),  9 masks,
index 197: E=(1,5,8,10,13,14),  2 masks.                   (7)
```

The ordered residual digest is

```text
ae75d7b0f6947dbc51dfacaa8b34dea9d2c985b04014eb09ac923ae720a0452d. (8)
```

For every residual mask, with `first_d=L/8` and `D=lcm(ds)`, one has

```text
D/first_d=8, hence D=L.                                     (9)
```

The pointed divisor interval is therefore the four-point divisor chain of
`8`, and all `19` residuals occupy its top.  This is not an iterated `C2`
action, parity tower, common carrier, owner map, or phase transport.

## 3. Complete terminal exclusion

All rows in `(3)` are wall rows, so an actual residual contains a high slot.
For the three bodies in `(7)`, the duplicate-permitting two-high gaps are

```text
64:  945034906/306786044565,
152: 35479656293/19991919811500,
197: 368677303/146893907160,                                (10)
```

and are strictly positive.  Hence every actual residual has exactly one high
slot.  The complete enlarged terminal bank has

```text
3 bodies, 19 masks,
19 zero-high hostile relaxations,
28 one-high cases over 12 body-local low-label sets,
28 coarse-cardinality closures,
0 exact-cell fallbacks, 0 max-gap fallbacks, 0 failures,
minimum cardinality slack 5.                                (11)
```

The terminal semantic digest is

```text
9328d556544c230503b6860ea27958d71bfa3f8af7816d02d350da360e92acac. (12)
```

The wall condition discards the zero-high relaxations.  Each remaining
one-high state lies in a complete upper relaxation and fails the necessary
projected support cardinality, so all three residual bodies and all `17`
selected rows close.

## 4. Consequence and stopping boundary

Composing only after promoted THM-3261 gives

```text
373283-17=373266.                                           (13)
```

The projected cap remains `z1<=216`.  The occupied layer becomes

```text
462 rows = 429 wall + 33 order.                             (14)
```

The two high-cost gcd-eight rows `(4)` remain open, as do every row in the
gcd `24,36,72` strata.  The theorem proves no physical-cover classification,
`k<=1` result, final rung, or LRC(14).

## 5. Exact evidence

The companion pins THM-3139 and promoted THM-3261 by source, output, and
semantic hashes; freshly derives all `480` rows and component counts; checks
the intrinsic selection `(2)--(4)`; reruns all `17` screens; binds the screen,
residual, quotient, gap, and terminal packets; and verifies `(13)--(14)`.
Normal and optimized runs must LF-normalize byte-for-byte to the stored
transcript with empty standard error.  All truth-bearing checks use explicit
exceptions rather than Python assertions.
