---
id: HYP-9029
title: "The strong-interval tiling law: H-spectrum completeness reduces to overlapping strong intervals"
status: >
  OPEN ALL-ORDER TILING LAW. THM-4097 promotes the complete order-nine strong
  spectrum and the exact coverage prefix: strong(7)>=[65,105],
  strong(8)>=[69,609], strong(9)>=[85,2881], so every odd in [65,2881]
  has a strong witness. THM-4102 adds a selected order-ten interval
  [249,14649]; THM-4104 adds selected order-eleven [429,80265], proving
  every allowed odd through 80405 globally realized. Complete order-ten and
  order-eleven images and global continuation remain open. THM-4111 proves
  that every full-cut recursively selected bank has unbounded maxima;
  THM-4115 strengthens the recurrence using exact Walsh variance. Neither
  proves interval overlap or unbounded solid-interval endpoints.
source: kind-pasteur-2026-07-26-S134
related:
  - THM-1370-h-spectrum-omits-7-21-all-n
  - THM-1862-order-join-reduction-principle
  - THM-1936-signed-redei-join-multiplicative
  - HYP-9028-circulant-hamiltonian-excess-tends-to-e
  - THM-4097-order-nine-strong-ear-spectrum-solid-interval-and-lane-extension
  - THM-4102-selected-order-ten-strong-ear-solid-interval
  - THM-4104-selected-order-eleven-strong-ear-solid-interval
  - THM-4111-uniform-ear-average-and-recursive-selected-bank-growth
  - THM-4115-uniform-ear-cut-walsh-variance-and-sharp-growth-refinement
script: 04-computation/strong_h_spectrum_intervals_kps_S134.py
output: 05-knowledge/results/strong_h_spectrum_intervals_kps_S134.out
data: 05-knowledge/results/strong_H_spectrum_m9_values_kps_S134.out
---

# HYP-9029 -- the strong spectra tile the odd line

The H-spectrum's gap half is PROVED (`{7,21}` forbidden); the completeness
half ("every other odd occurs") is THM-1370's standing conjecture. THM-4097
proves exact global coverage through `2885`; THM-4102/4104's selected banks
extend it through `80405`, with `80407` and `7*11527=80689` the first two lane
targets not forced by current finite strong values and multiplication. Canon reduces
completeness to strong tournaments
(spectrum = multiplicative closure of strong H-values, THM-1862/
THM-1936, machine-checked). This hypothesis records the missing
structural mechanism, sitting unremarked in the exhaustive censuses:

**The strong spectra contain giant contiguous odd intervals that
overlap consecutively.**

```text
strong(7) >= [65, 105]     (islands below: [25,61], hole only at 63)
strong(8) >= [69, 609]     (islands: 45, [49,53], [57,65])
strong(9) >= [85, 2881]    (islands: 75, 81; 1,399 consecutive odds)
strong(10) >= [249,14649]  (selected bank; 7,201 consecutive odds)
strong(11) >= [429,80265]  (selected bank; 39,919 consecutive odds)

junctions: 69<=105, 85<=609, 249<=2881, 429<=14649.
Thus every odd in [65,80265] has a STRONG witness; multiplication bridges
the remaining allowed values through 80405.
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

THM-4102 confirms that a deterministic one-witness-per-value order-nine bank
already tiles the much longer odd interval `[249,14649]` at order ten. This is
a selected image rather than a complete order-ten census, so it supports but
does not prove the all-order tiling law.

THM-4104 repeats the selection at order ten and obtains `[429,80265]` at order
eleven. Its directed quadratic `C+sum_S L-sum_(SxS)Q` is exactly THM-4097's
symmetric-cut-plus-zero-sum-field carrier in unsymmetrized coordinates. This
suggests a weaker sufficient theorem than full spectral induction: choose one
labelled representative per attained value so that the next bank's cut
quadratics contain an overlapping interval. Equal-`H` parents can have
different quadratics, so the labelled `(Q,L)` sidecar is indispensable.

THM-4111 removes growth of the scalar maximum from the list of conjectural
steps. Its uniform cut identity implies `M_(n+1)>=(n+3)M_n/4` for every
one-representative-per-value bank that expands all nonconstant cuts, hence
factorial-over-exponential and therefore unbounded maxima. This does not yet
advance the solid right endpoint `d_n`: the averaging operation destroys the
individual-cut distribution, and its maximum can sit outside every long
interval. The open structural coordinate is now explicitly **dispersion and
overlap**, not raw growth.

THM-4115 sharpens this diagnosis. Its exact degree-two Walsh expansion gives

```text
Var=1/4(sum_i h_i^2+sum_(i<j)w_ij^2)
```

and, using the support bound `H(T+x_S)>=H(T)`, improves the universal bank
recurrence to

```text
M_(n+1) >= ((n+1)(n+2)/(4n))M_n.
```

The directed triangle attains equality, so this invariant-free order-three
coefficient is sharp. Variance separates the known equal-mean hostile, but it
still collapses the labelled coefficient arrangement and hence does not force
a small-ball estimate, a `+/-2` chain, or an interval. In fact, `400/544`
labelled strong order-five parents have no `M-2` one-bit neighbor at any
maximizing cut. The cheapest next test is therefore a **multi-cut or
multi-parent overlap** theorem conditional on the full `(w,h)` geometry, not
single-bit descent from the scalar maximum.

Independent-path note (same session): a LABELED brute-force
enumeration (2^28 tournaments, different algorithm from monad's
iso-class generation; strong_h_spectrum_m8_labeled_kps_S134.cpp)
reproduces the m = 8 strong spectrum EXACTLY -- 297 values, min 45,
max 661, identical sets -- independently re-verifying both the
MISTAKE-055 correction (min 45, not the Leonardo 41) and this
hypothesis's central interval [69, 609].

Cheapest decisive tests: (i) compile the selected order-twelve bank
(`43,251*(2^11-2)=88,491,546` ears) with the quadratic recurrence; (ii) prove
the bank-selection overlap lemma by adding a variance/small-ball or local
`+-2` control to THM-4111's exact mean, without enumerating every equal-`H` fibre;
(iii) attempt the `+-2` surgery on the doubly-regular-minus-arc family; and
(iv) prove `d_n>=f(n+2)` for one explicit construction family.
