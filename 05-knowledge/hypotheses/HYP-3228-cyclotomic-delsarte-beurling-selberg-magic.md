---
id: HYP-3228
title: Cyclotomic Delsarte/Beurling-Selberg magic function
status: EVIDENCE / exact bounded-bank synthesis; not an LRC14 proof
source: codex-2026-06-28
tangent: T1326
technique: LTI-326
tournament_technique: LTT-226
script: 04-computation/lrc_cyclotomic_delsarte_bs_magic_codex_20260628.py
result: 05-knowledge/results/lrc_cyclotomic_delsarte_bs_magic_codex_20260628.out
reflection: 07-reflections/cyclotomic-delsarte-beurling-selberg-magic-codex-20260628.md
related:
  - HYP-3245
  - HYP-3242
  - HYP-3241
  - HYP-3240
  - HYP-3239
  - HYP-3238
  - HYP-3237
  - HYP-3236
  - HYP-3235
  - HYP-3234
  - HYP-3233
  - HYP-3232
  - HYP-3230
  - HYP-3229
  - HYP-3227
  - HYP-3226
  - HYP-3225
  - HYP-3224
  - HYP-3223
  - HYP-3222
  - HYP-3221
  - HYP-3215
  - HYP-3214
  - HYP-3213
  - HYP-3212
  - HYP-3210
  - HYP-3205
  - HYP-3204
  - HYP-3203
  - HYP-3202
  - HYP-3200
  - HYP-3153
  - HYP-3138
  - HYP-3132
  - OPEN-Q-108
---

# HYP-3228: Cyclotomic Delsarte/Beurling-Selberg Magic Function

## Claim

The k=8 `L_y` dual has a concrete finite magic-function carrier:

```text
f(n) = ((n-1)(n-2)(n-4)(n-5))/4,  n=0,...,6.
```

Its shell values are

```text
[10, 0, 0, 1, 0, 0, 10],
```

so `E[f(N)] = 10q0 + q3 + 10q6 = 10 L_y`.  This is the finite
Delsarte/Beurling-Selberg packet behind the current k=8 coefficient face:
zeros at shells `1,2,4,5`, positive mass at endpoint shells `0,6`, and one
central Worpitzky repair at shell `3`.

The same object has two other useful forms:

```text
f(n) = 10 - 10 C(n,1) + 10 C(n,2) - 9 C(n,3) + 6 C(n,4),
P(z)=10+z^3+10z^6,
z^-3 P(z)=10(u^3-3u)+1,  u=z+z^-1.
```

Modulo the real seventh-cyclotomic cubic `u^3+u^2-2u-1`, the Joukowski form
reduces to `-10u^2-10u+11`.  This pins the "magic function" language to the
same level-7 Chebyshev/de Moivre arithmetic as HYP-3212/HYP-3213.

## Exact Scout

The script recomputes the exact anchored bounded bank

```text
E = {0} union A, A subset {1,...,14}, |A|=7
rows_all = 3432
rows_primitive = 3431
```

For consecutive speeds,

```text
q = [481/1470, 359/1470, 25/147, 26/245, 17/210, 5/98, 1/49]
magic = 2633/735
L_y = 2633/7350
```

Extremality matches HYP-3204 exactly:

```text
primitive magic beaters = 0
primitive magic ties = 1
all-bank magic beaters = 0
all-bank magic ties = 2
ties = consecutive and doubled AP
```

The all-bank doubled AP tie is the known nonprimitive dilation artifact.

## AP Support Alignment

HYP-3203's AP-polarized support normal is still the useful root-locus normal.
The script verifies:

```text
primitive AP-support beaters = 0
primitive AP-support ties = 1
all-bank AP-support beaters = 0
all-bank AP-support ties = 2
max support = 39766/540225
```

For every primitive non-AP row in the bounded bank, the AP-support deficit
controls the magic deficit by a positive exact ratio.  The observed range is:

```text
min magic_deficit/support_deficit = 7390327/146434  ~= 50.468654821
max magic_deficit/support_deficit = 6857550/28337  ~= 241.999858842
```

This is strong alignment, but not identity.  The cosine between the centered
magic vector `[7,-3,-3,-2,-3,-3,7]` and the AP-support direction is only
`0.182704091`.  AP support is a coercive normal; the full magic vector still
needs the ordered-tail/central-exchange sidecar from HYP-3204.

## Cyclotomic Guardrail

A naive reading would demand that the Laurent polynomial

```text
10 + rho z^3 + 10z^6
```

be nonnegative on the seventh roots.  The exact scout shows why that is not
the right terminal certificate: cyclic nonnegative completion would require

```text
rho >= 18.019377358048,
```

while the actual Delsarte/`L_y` packet has `rho=1`.  Literal cyclic
PSD/Fejer positivity therefore overprices the central `q3` repair by a large
factor.  This is compatible with HYP-3221: config-blind algebraic positivity
is not the final proof.  The finite magic function should instead be glued to
AP support, ordered-tail exchange, Toeplitz trap discharge, and the
Joukowski/Hermite-Biehler sidecar.

Post-rebase HYP-3214 sharpens the language: the positive sector-side
cyclotomic magic function is the Fejer kernel `F_7`, while HYP-3228 is the
shell `L_y` Delsarte dual.  The two objects should be treated as distinct
faces of the same level-7 packet rather than identified coefficient-by-
coefficient.

The later HYP-3215 literature-route note adds a proof-routing warning: the
Fejer/Cohn-Elkies LP face may need a separate induction-base or
polyhedron-flatness handoff before it can close LRC(14).  HYP-3228 should be
used as a finite shell-dual coordinate inside that larger route, not as a
standalone replacement for the base-case check.

The later mac-mini S75 HYP-3227 comb-overlap Gram kernel adds the spatial
dual of HYP-3214's Fejer kernel and proves a single-arc peeling lemma for
speed `1`.  HYP-3228 should be tested against that Gram/peeling coordinate:
the shell `L_y` dual is a row-level coefficient witness, while the Gram kernel
is a measure-domain PSD witness.

Post-rebase companion HYP-3245 supplies the lag-space projection of the same
split.  Its finite trap table says non-AP local maxima move ordinary
speed-support autocorrelation mass from short lags to outward lags.  This is a
natural sidecar candidate for decomposing `magic_deficit`: HYP-3228 names the
shell-contact equation, while HYP-3245 records the transport that may connect
that equation to AP support, Toeplitz slack, HYP-3236 Green resistance slack,
Vitali/Brouwer core-wall sign, Brouwer trace-sign/SOS factorization, or
ordered-tail pricing.  Incoming HYP-3232 is the separate apex-fold/modulus-covariance
companion.
Incoming HYP-3238 and HYP-3239 sharpen the sign bookkeeping for this handoff:
shell magic cannot be compressed through the positive/even channel unless the
odd/negative payload is zero, restored, dual-annihilated, or explicitly kept,
and the `n=14` reflection packet should be tested in the dihedral sign irrep /
Borsuk-Ulam chart, not only the Brouwer fixed-point chart.
Incoming HYP-3240 and HYP-3241 sharpen the core witness bookkeeping: the shell
deficit should distinguish primitive `Phi_14` core witnesses shared by AP and
Goddyn-Wong from promoted `Phi_{14d}` dilation witnesses, and it should not
collapse the dip to a single imaginary-quadratic norm scalar.
Incoming HYP-3242 adds the topological reading: shell slack and lag transport
should also be compared with danger-cover Euler characteristic, lonely-hole
survival, and Cech/Betti sidecars.

## Proof-Frontier Use

The next theorem-facing statement should not be "find a positive Fourier
kernel."  It should be a multi-chart certificate:

```text
quartic shell contact:
  f(n)>=0 and f has contact zeros at 1,2,4,5

Delsarte/Newton face:
  10S0 -10S1 +10S2 -9S3 +6S4

AP-support face:
  <q-u, q_AP-u> maximized by AP/doubled AP

ordered-tail face:
  central q3 gain is priced by q0+q6 loss

moment/spectral face:
  Toeplitz and covariance layers discharge local traps

Joukowski/HB face:
  10(u^3-3u)+1 carries the level-7 Chebyshev arithmetic
```

This narrows HYP-3224's requested dual: the visible coefficient vector is no
longer unknown.  The remaining work is to prove that every legal primitive row
has nonnegative slack against this packet, with AP as the only primitive
equality row, without replacing the analytic/equidistribution obligation by a
false cyclic PSD claim.

## Tournament Analysis

The script uses proof currencies, not runners or raw arcs, as vertices.
Alternate vertex sets considered were miss-count shells, Laurent modes,
AP-support normal, exchange traps, and sector-boundary events.  The chosen
quotient preserves the `L_y` shell functional and AP-support direction but
forgets sector identity, covariance-pair labels, and exchange-neighbor
witnesses.  The challenged assumption is explicit: cyclotomic magic need not
mean nonnegative cyclic Fourier kernel.

The proof-carrier tournament is transitive:

```text
quartic_shell_magic_contact
-> delsarte_newton_binomial_dual
-> joukowski_chebyshev_cubic
-> ap_polarized_support_normal
-> ordered_tail_exchange_price
-> toeplitz_trap_discharge_sidecar
-> cyclic_psd_regularization
-> raw_cyclotomic_energy
-> raw_q3_scalar
```

-> HYP-3245, HYP-3243, HYP-3242, HYP-3241, HYP-3240, HYP-3239, HYP-3238, HYP-3237, HYP-3236, HYP-3235, HYP-3234, HYP-3233, HYP-3232, HYP-3230, HYP-3229, HYP-3228, HYP-3215, HYP-3214, HYP-3226, HYP-3224, HYP-3223, HYP-3222, HYP-3221, HYP-3219, HYP-3213,
HYP-3212, HYP-3210, HYP-3205, HYP-3204, HYP-3203, HYP-3202, HYP-3200,
HYP-3153, HYP-3138, HYP-3132, T1326, LTI-326, LTT-226, OPEN-Q-108.
