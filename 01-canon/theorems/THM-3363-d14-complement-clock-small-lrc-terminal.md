---
id: THM-3363
title: "The D=14 complement-clock small-LRC terminal"
status: >
  PROVED + VERIFIED-EXACT + CITED-INPUT.  In the one-aligned/six-drift
  support-transfer branch, the canonical row E={1,2,3,4,6,12}, L=168, D=14
  is impossible for every aligned multiplier and every one of its 26
  denominator shapes, with no numerator-height cutoff.  Its projected carrier
  plus one complement clock would make eight integer danger combs cover the
  circle, contradicting the cited lonely-runner theorem for eight speeds.
source: codex-major-frontiers-2026-08-14
external_input: >
  Tanupat Trakulthongchai, "Nine and Ten Lonely Runners," Electronic Journal
  of Combinatorics 33(2) (2026), P2.46, DOI 10.37236/14972.  The proof uses
  only its nine-runner case: eight nonzero integer speeds have a common lonely
  time at distance at least 1/9.
depends_on:
  - THM-2928-critical-seven-comb-grid-tensorization-and-drift-tariff
related:
  - THM-2174-endpoint-phase-scale-obstruction
  - THM-3360-uniform-reflected-high-pair-floor-by-admissible-affine-tails
script: 04-computation/lrc14_d14_complement_clock_terminal_thm3363.py
output: 05-knowledge/results/lrc14_d14_complement_clock_terminal_thm3363.out
script_sha256: b3069105298192bf0ebae590ae16ca2cc73d5f238e4fc6079730f8df32419fef
output_sha256: 78df70987b3495d103fe42709cf7fcb8794afe43540088f929bc6be802cf8f01
hash_basis: LF-normalized bytes
---

# THM-3363 -- the D=14 complement-clock small-LRC terminal

**PROVED + VERIFIED-EXACT + CITED-INPUT.**

## 1. Statement

Write

```text
D_s={x in T: ||sx||<1/14},              T=R/Z.             (1)
```

Consider the one-aligned/six-drift support-transfer row with

```text
E={1,2,3,4,6,12},       L=168,       D=14.                (2)
```

Let `alpha>=1` be the aligned multiplier and let `b_1,...,b_6` be the six
positive integer quotient speeds supplied by THM-2928.  Then the exact
pointwise containment

```text
Y_14(alpha) subset D_(b_1) union ... union D_(b_6)         (3)
```

is impossible.  This holds uniformly over all quotient numerators and all
`26` denominator multisets of size six with entries in `{2,7,14}` and lcm
`14`.  No height bound or asymptotic ray decomposition is required.

Consequently this canonical body/divisor row is deleted from the arbitrary-
residue `k=1` necessary ledger.  The ledger changes from `27,240` to `27,239`
body/divisor rows, and its denominator-shape occurrence count changes from
`38,954,725,590,760` to `38,954,725,590,734`.  The remaining `27,239` rows,
projected sectors, physical entry, the rung, and LRC(14) remain open.

## 2. Exact projected carrier

The body-safe cell calculation in THM-2928 gives the twelve half-open integer
ranges

```text
[12,13), [15,26), [30,39), [45,52), [60,69), [71,78),
[90,97), [99,108), [116,123), [129,138), [142,153), [155,156). (4)
```

They contain `88` addresses, and their support modulo `14` is exactly

```text
S_14={1,2,...,12}.                                        (5)
```

Let `R_alpha=T\D_alpha`.  THM-2928's support projection is

```text
Y_14(alpha)=union_(r=1)^12 (r+R_alpha)/14.                (6)
```

For `x=(r+u)/14`, integer `r` disappears from the aligned clock:

```text
||14 alpha x||=||alpha(r+u)||=||alpha u||.                (7)
```

The selected cells in `(6)` fill the closed middle arc from `1/14` to
`13/14`; their cell boundaries lie in `D_(14 alpha)`.  Equations `(5)`--`(7)`
therefore give the pointwise identity

```text
Y_14(alpha)=[1/14,13/14] \ D_(14 alpha).                 (8)
```

The companion exact referee checks `(4)`--`(8)` at every arrangement endpoint
and open atom for `1<=alpha<=64`, totaling `118,336` rational probes.  This is
a hostile endpoint audit of the symbolic identity, not the reason it holds.

## 3. Add the missing complement clock

The complement of the open unit-speed danger arc is

```text
T\D_1=[1/14,13/14].                                      (9)
```

If `(3)` held, `(8)` and `(9)` would imply

```text
T = D_1 union D_(14 alpha) union
    D_(b_1) union ... union D_(b_6).                     (10)
```

After deduplicating, `(10)` is a cover by at most eight positive integer
speeds.  (In the physical row they are in fact distinct: THM-2928 gives
distinct positive `b_i`; exact denominator `d_i>1` excludes `14|b_i`; and
`b_i=1` would duplicate the body speed `12`.)  Distinctness is not
load-bearing because collisions only reduce the number of combs.

Trakulthongchai's published nine-runner theorem supplies `x in T` for these
at most eight nonzero speeds with

```text
min(||x||,||14 alpha x||,||b_1x||,...,||b_6x||)>=1/9.    (11)
```

Since `1/9>1/14`, this `x` lies in none of the eight open sets in `(10)`, a
contradiction.  Thus `(3)` is impossible.

The cited source is [*Nine and Ten Lonely Runners*](https://www.combinatorics.org/ojs/index.php/eljc/article/view/v33i2p46),
Electronic Journal of Combinatorics **33** (2026), P2.46.

## 4. Mechanism and boundary

The useful operation is **complement-clock completion**:

```text
projected carrier + clocks covering its complement
                    -> a lower-dimensional global LRC cover.              (12)
```

It preserves the pointwise open-set semantics and avoids every numerator-
height coordinate.  In general, if the complement of a projected carrier is
covered by `r` positive integer danger combs and the carrier were covered by
`p` more, deduplication would produce a global cover by at most `p+r` combs.
Any already-known LRC case at radius at least `1/14` then rules it out.

This theorem uses the exceptional contiguous support word `(5)`.  It does not
assert that every one-aligned support complement has a five-clock completion,
nor does it bound the quotient speeds on other rows.  Finding such small
complement covers, or proving that none exist, is the next finite support-level
test.

## Reproduction

```text
python 04-computation/lrc14_d14_complement_clock_terminal_thm3363.py
python -O 04-computation/lrc14_d14_complement_clock_terminal_thm3363.py
```
