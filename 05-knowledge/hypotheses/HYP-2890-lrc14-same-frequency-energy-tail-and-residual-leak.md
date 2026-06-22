---
id: HYP-2890
status: COMPUTATIONAL SIGNAL / proof-target; exact bounded scans for k=8,9
source: codex-2026-06-22-S104
tags: [lrc14, additive-energy, fourier-tail, fejer, p0, residual, tournament-analysis]
related:
  - HYP-2885
  - HYP-2889
  - HYP-+2888
  - HYP-2887
  - HYP-2636
  - THM-534
  - OPEN-Q-108
results:
  - 04-computation/lrc14_gamma_frequency_tail_codex_s104.py
  - 05-knowledge/results/lrc14_gamma_frequency_tail_codex_s104.out
---

# HYP-2890: the additive-energy frequency tail is positive, but the closing constant is a residual-leak inequality

KPS S31k proved that the `m=1` additive-energy coefficient

```text
Gamma_k(1) = sum_A (-1)^|A| |hat f_A(1)|^4 (1-|A|/7)^(k-4)
```

is positive for `k=8..12`.  S104 extends this to the whole same-frequency
packet

```text
Gamma_k^sf = sum_{m>=1, 7 not divide m} Gamma_k(m),
Gamma_k(m) = sum_A (-1)^|A| |hat f_A(m)|^4 (1-|A|/7)^(k-4).
```

For fixed residue `r mod 7`, `m^4 Gamma_k(m)` is constant.  Therefore this is
an explicitly convergent `1/m^4` tail, not one of the conditionally convergent
support-level tails.

## Confirmed packet facts

The residue constants are all positive for `k=8..13`; the `H=12` partial sum
already certifies positivity with a tiny absolute tail:

```text
k=8:  Gamma_sf=+1.869969145e-03, H=12 tail <= 1.084e-06
k=9:  Gamma_sf=+1.869969145e-03, H=12 tail <= 1.084e-06
k=10: Gamma_sf=+1.631180083e-03, H=12 tail <= 9.458e-07
k=11: Gamma_sf=+1.288806723e-03, H=12 tail <= 7.473e-07
k=12: Gamma_sf=+9.354693542e-04, H=12 tail <= 5.424e-07
k=13: Gamma_sf=+6.212037851e-04, H=12 tail <= 3.602e-07
```

The full same-frequency packet is a fixed multiple of the KPS `m=1` coefficient
in this range:

```text
Gamma_sf / Gamma(1) = 2.098724...
```

Thus the leading additive-energy monotonicity is not a one-frequency accident.
It has a clean positive same-frequency tail.

## Guardrail: the packet overpredicts AP

Let

```text
R_sf(E) = p0(E) - p0_decorr(k) - Gamma_k^sf A*(E).
```

For AP rows, the same-frequency packet overpredicts the actual deviation:

```text
k=8 AP:  actual dev +0.269509, samefreq +0.418873, R_sf=-0.149364
k=9 AP:  actual dev +0.311244, samefreq +0.628310, R_sf=-0.317066
k=10 AP: actual dev +0.341382, samefreq +0.782966, R_sf=-0.441584
```

So additive energy alone is still not the proof.  The missing term is a
negative hidden-fold / support-cycle correction that grows with AP-likeness.
This agrees with HYP-2889: scalar additive energy is an AP-facing carrier, not
a monotone sufficient statistic.

## The sharper proof target

For a fixed `k`, AP extremality is equivalent to the residual-leak inequality

```text
R_sf(E) - R_sf(AP_k)
  <= Gamma_k^sf (A*(AP_k) - A*(E)).
```

The left side is the amount of negative AP correction that leaks back when the
row is deformed away from AP.  The right side is the positive same-frequency
energy advantage that AP has.

S104 exact bounded scans found no violations:

```text
k=8, anchored max<=14: rows=3432, violations=0
  corr(A*,p0)=+0.5603, corr(A*,R_sf)=-0.4502
  worst leak ratio = 0.469476
  row (0,2,4,6,7,8,10,12)

k=9, anchored max<=14: rows=3003, violations=0
  corr(A*,p0)=+0.5433, corr(A*,R_sf)=-0.3518
  worst leak ratio = 0.932737
  row (0,2,4,6,7,8,10,12,14)
```

The k=9 pressure row is especially informative: it is the even AP
`0,2,4,6,8,10,12,14` plus the midpoint bridge `7`.  This is exactly the
scaling-invariant tiling signal from HYP-+2888, not a random high-energy row.

## Proof shape

The next theorem should not be scalar `p0 <= G(A*)`.  It should be:

```text
positive same-frequency additive-energy packet
  + AP-facing Fejer/difference-profile majorization
  + residual-leak bound for labelled hidden-fold/support-cycle terms
  => p0(E) <= p0(AP_k).
```

Candidate routes for the residual-leak bound:

1. Low relation-depth / near-scale-reducible rows:
   finite AP/Freiman atlas, with the worst k=9 bridge row as the first target.
2. High relation-depth rows:
   HYP-2636 channelwise Abel/L2 cancellation after exact-period or coimage
   packet projection.
3. Repeated-packet rows:
   HYP-2887 octahedral current/Hodge decomposition, where the residual is a
   face-curl / wall-ledger correction rather than scalar energy.

## Tournament Analysis

Vertices are proof lenses, not runners:

```text
m=1_lead,
same_frequency_tail,
hidden_fold_shape,
support_cycle_tail.
```

The tested observable is `(explained_AP_deviation, monotonicity_payoff,
remaining_residual, formal_tail_cost)`.  The useful path is

```text
m=1_lead -> same_frequency_tail -> hidden_fold_shape -> support_cycle_tail.
```

Assumption challenged: once the leading additive-energy coefficient is
positive, summing its frequency tail should close the proof.  It does not.  The
frequency tail is clean and positive, but it exposes the real constant as a
residual-leak inequality over labelled Fourier/sector packets.

## Post-pull synthesis with KPS S31l

The incoming KPS S31l moment-coefficient audit is exactly the missing warning
label for `R_sf`.  S104 proves the same-frequency `s=2` additive-energy packet
is positive and has a clean `1/m^4` tail.  S31l checks higher additive moments
at frequency `1`:

```text
k=9:  s=2 +8.910e-04, s=3 -7.138e-05, s=4 -4.448e-05
k=12: s=2 +4.457e-04, s=3 +3.689e-05, s=4 -9.843e-06, ...
```

So the residual after the positive same-frequency energy packet is genuinely
signed.  This rules out a Savchenko-style "every term is maximized by the
regular/AP object" proof.  The surviving analogy is the H-max/Jensen template:
keep the signs, prove the AP-facing convexity/majorization inequality for the
whole labelled functional, and use the residual-leak inequality above as the
finite quantitative form of that convexity.

## S105 carrier refinement

HYP-2891 supplies finite quotient objects for this residual-leak inequality.
After choosing a tangent sector, the 64 residual masks fold to 16 Clebsch
classes with exact pair design balance, but every class mixes missed depth.
Thus the residual should be decomposed by signed covariance class, not by scalar
`q_t`.  Local AP-facing compression should be tracked on the Bruhat `S4`
carrier: commutation squares are the likely design-balanced components, while
braid hexagons are the finite low-depth AP/Freiman packets to inspect first.
