---
id: HYP-9029
title: "The strong-interval tiling law: H-spectrum completeness reduces to overlapping strong intervals"
status: >
  OPEN ALL-ORDER TILING LAW. THM-4097 promotes the complete order-nine strong
  spectrum and the exact coverage prefix: strong(7)>=[65,105],
  strong(8)>=[69,609], strong(9)>=[85,2881], so every odd in [65,2881]
  has a strong witness and every allowed odd through 2885 is globally
  realized. The next interval, order ten, and global continuation remain
  open.
source: kind-pasteur-2026-07-26-S134
related:
  - THM-1370-h-spectrum-omits-7-21-all-n
  - THM-1862-order-join-reduction-principle
  - THM-1936-signed-redei-join-multiplicative
  - HYP-9028-circulant-hamiltonian-excess-tends-to-e
  - THM-4097-order-nine-strong-ear-spectrum-solid-interval-and-lane-extension
script: 04-computation/strong_h_spectrum_intervals_kps_S134.py
output: 05-knowledge/results/strong_h_spectrum_intervals_kps_S134.out
data: 05-knowledge/results/strong_H_spectrum_m9_values_kps_S134.out
---

# HYP-9029 -- the strong spectra tile the odd line

The H-spectrum's gap half is PROVED (`{7,21}` forbidden); the completeness
half ("every other odd occurs") is THM-1370's standing conjecture. THM-4097
now proves exact global coverage through `2885`, with `2887` the first target
not forced by the known finite strong values and multiplication. Canon already
reduces completeness to strong tournaments
(spectrum = multiplicative closure of strong H-values, THM-1862/
THM-1936, machine-checked). This hypothesis records the missing
structural mechanism, sitting unremarked in the exhaustive censuses:

**The strong spectra contain giant contiguous odd intervals that
overlap consecutively.**

```text
strong(7) >= [65, 105]     (islands below: [25,61], hole only at 63)
strong(8) >= [69, 609]     (islands: 45, [49,53], [57,65])
strong(9) >= [85, 2881]    (islands: 75, 81; 1,399 consecutive odds)

junctions: 69 <= 105, 85 <= 609  -- every odd in [65, 2881] has a
STRONG witness; products of m <= 8 strong values give all odd <= 400
except {7, 21} (machine-checked in the monad s5 out).
```

**Conjecture (tiling law).** For every `n >= 7`, `strong(n)`
contains a contiguous odd interval `[c_n, d_n]` with
`c_{n+1} <= d_n`. Consequently (with the canon semigroup law and
the finite base) `spectrum(infinity) = N_odd \ {7, 21}` -- the
completeness conjecture in full.

Why the tiling should hold with room to spare: the left edge `c_n`
tracks the Busch floor `f(n) = min{3^a 5^b : 2a+3b = n-1}` (growth
`~ 5/3` per step: 25, 45, 75, 125, ...), while the right edge `d_n`
tracks near-maximal H, i.e. `~ e n!/2^{n-1}` up to the HYP-9028
excess constant -- factorial against geometric. The observed slack
grows `36 -> 524` in one step. The two S134 threads meet here: the
excess constant governs `d_n`, the Busch floor governs `c_n`.

Proof route (recorded): a surgery/insertion move within strong
tournaments that changes `H` by exactly `+-2` inside a band would
give contiguity directly; THM-1370's insertion lemma and the
THM-1900 insertion-response calculus are the natural instruments.
The single-arc-core family provably cannot do it (its gap density
is 1/2, OPEN-Q-055), so the surgery must use a richer family.

THM-4097 identifies the richer exact carrier. Every strong tournament is a
nonconstant ear of a strong parent, and the response is
`H(parent)+cut_w(S)+sum_S h`, where `w` is an integral symmetric cut weight
and `h` is an integral zero-sum orientation field. At `8->9`, cut weights
`{3,4}` generate the entire 1,482-value strong spectrum, while `w=4` alone
misses exactly `89,93,105,125`. Any all-order surgery must control the field,
not merely the cut energy.

Independent-path note (same session): a LABELED brute-force
enumeration (2^28 tournaments, different algorithm from monad's
iso-class generation; strong_h_spectrum_m8_labeled_kps_S134.cpp)
reproduces the m = 8 strong spectrum EXACTLY -- 297 values, min 45,
max 661, identical sets -- independently re-verifying both the
MISTAKE-055 correction (min 45, not the Leonardo 41) and this
hypothesis's central interval [69, 609].

Cheapest decisive tests: (i) m = 10 strong interval via canonical
augmentation (the monad s6 generator extends; expect
`c_10 ~ 100-130`, `d_10 ~ 15,000`); (ii) attempt the +-2 surgery
lemma on the doubly-regular-minus-arc family; (iii) prove
`d_n >= f(n+2)` for one explicit construction family, which already
chains two floors ahead.
