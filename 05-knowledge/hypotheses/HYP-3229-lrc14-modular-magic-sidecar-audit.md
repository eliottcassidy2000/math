---
id: HYP-3229
title: Modular magic sidecar audit for the LRC14 Fejer certificate route
status: EVIDENCE / exact algebraic scout plus proof-target triage; not an LRC14 proof
source: codex-2026-06-28
tangent: T1329
technique: LTI-329
tournament_technique: LTT-229
script: 04-computation/lrc14_modular_magic_sidecar_audit_codex_20260628.py
result: 05-knowledge/results/lrc14_modular_magic_sidecar_audit_codex_20260628.out
reflection: 07-reflections/lrc14-modular-magic-sidecar-audit-codex-20260628.md
related:
  - HYP-3227
  - HYP-3226
  - HYP-3225
  - HYP-3224
  - HYP-3215
  - HYP-3214
  - HYP-3213
  - HYP-3212
  - HYP-3210
  - HYP-3205
  - HYP-3203
  - HYP-3201
  - HYP-3162
  - HYP-3161
  - HYP-3160
  - THM-577
  - OPEN-Q-108
---

# HYP-3229: Modular Magic Sidecar Audit

## Claim

The current LRC14 proof route should keep two kernels at the core and four
arithmetic analogies as labelled sidecars.

Core:

```text
F_7 sector kernel
  = Fejer / de-Moivre / Chebyshev / positive-definite slack

comb-overlap Gram kernel
  = S75 spatial kernel K(p,q)=meas(D_p cap D_q)=<1_Dp,1_Dq>

Johnson J(14,2) cap kernel
  = binomial cap mass on the 14-clock
```

Proof-facing sidecar:

```text
Gamma0(7) Eisenstein candidate
  E_7(tau) = (7 E_2(7 tau) - E_2(tau)) / 6
```

Guarded inspiration sidecars:

```text
Dirichlet-L / Stark conductor-7 arithmetic
Beraha/Tutte chromatic-root constants
Mahler/Lehmer height of the de-Moivre unit
subshift / transfer-operator autocorrelation
```

This refines HYP-3214 rather than replacing it.  HYP-3214 already showed that
the cyclotomic Delsarte magic function is the Fejer kernel.  HYP-3229 says
the next useful modular move is not a raw closed-form cap guess; it is a
coefficient-generating sidecar for an explicit finite LP/Toeplitz certificate.

Post-fetch integration with HYP-3215: the literature-facing proof route names
Fejer/Cohn-Elkies LP construction as Gap A and separately flags the need to
verify the LRC induction base.  HYP-3229 supplies a concrete candidate
coefficient engine for that Gap A finite LP.  It does not settle the HYP-3215
base flag, and it should not be used to blur the bounded certificate problem
with the global descent/base problem.

Second fetch-time integration with mac-mini S75: the spatial build of the
magic function is the comb-overlap Gram kernel
`K(p,q)=meas(D_p cap D_q)=<1_Dp,1_Dq>`, PSD by construction.  S75 proves
`K(1,q)=1/(7q)` for `q<=13`, `K(7,q)=1/49`, and the single-arc peeling
recursion

```text
cap(P) = cap(P\{1}) - (1/7)(1 - 1/min(P\{1}))
```

for `1 in P`, speeds `<=13`.  This is the finite spatial dual of HYP-3214's
Fejer Fourier kernel.  HYP-3229 should now match the Gamma0(7) coefficient
engine to this Gram basis, with the remaining debt isolated in order-3
triple-overlap constants.

## Exact Readout

The scout verifies the Chebyshev/de-Moivre identity exactly:

```text
V_7(u) = u^7 - 7u^5 + 14u^3 - 7u
m(u) = u^3 + u^2 - 2u - 1
(V_7(u)-2)/(u-2)
  = u^6 + 2u^5 - 3u^4 - 6u^3 + 2u^2 + 4u + 1
  = m(u)^2
disc(m) = 49 = 7^2
```

The Fejer Fourier coefficients are the triangular positive sequence
`(7-|n|)_+` for `|n|<=6`.  This is still the strongest certificate payload:
nonnegative trigonometric function, nonnegative Fourier coefficients, and
double zeros at the de-Moivre angles.

The explicit level-7 modular sidecar is:

```text
E_7(tau) := (7 E_2(7 tau) - E_2(tau))/6
a_n = 4 sigma_1(n) - 28 sigma_1(n/7)  if 7 | n
a_n = 4 sigma_1(n)                    otherwise
```

First coefficients:

```text
[4, 12, 16, 28, 24, 48, 4, 60, 52, 72, 48, 112, 56, 12,
 96, 124, 72, 156, 80, 168, 16, 144, 96, 240, 124, 168,
 160, 28]
```

This is proof-relevant because it is a divisor-fibered coefficient source at
level 7.  It can be used to construct finite certificate rows with retained
prime-power and conductor labels instead of a scalar smoothing guess.

The spatial finite certificate side is the S75 Gram row:

```text
K(1,q) = 1/(7q), q=1,...,13
K(7,q) = 1/49
```

So the next certificate should be a three-way compatibility object:

```text
Fejer Fourier weights
Gamma0(7) divisor coefficients
comb-overlap Gram entries
```

## Guardrails

The Dirichlet-L/Stark closed-form route is not yet a proof route.  The raw
Bernoulli samples `B_2(a/7)` do carry a conductor-7 grid, and the de-Moivre
field has discriminant `7^2`, but the normalized `L(-1, chi)` values for the
even primitive mod-7 characters contract to denominator `7`:

```text
principal:          L(-1,chi) = 1/2
even primitive chi2: L(-1,chi) = -4/7 + 2/7 omega
even primitive chi4: L(-1,chi) = -6/7 - 2/7 omega
```

So the honest conclusion is:

```text
7^2 is real as discriminant/double-root arithmetic.
It is not yet a direct Stark/Dirichlet-L cap formula.
```

The Beraha and Mahler sidecars are exact but not terminal:

```text
B_7 = 2 + 2 cos(2pi/7) has polynomial B^3 - 5B^2 + 6B - 1
Mahler(m) = B_7 - 1
```

This gives a height monitor for perturbations of the de-Moivre cubic, not a
certificate inequality.

The subshift sidecar is similarly useful only after labelling the retained
predicate:

```text
P(z) = 1 + z + ... + z^6
P(z)P(z^-1) has autocorrelation coefficients 7-|n|
```

Thus the AP row is the rank-one transfer state whose autocorrelation is the
same Fejer kernel.  Perturbed rows should be studied by transfer-matrix
Perron defects, but that is a route to a certificate comparison, not a proof
by analogy.

## Proof Target

The sharper LRC14 proof target is:

```text
consec is the unique configuration simultaneously sharp for
  (1) Fejer/de-Moivre sector slack,
  (2) S75 comb-overlap Gram positivity and single-arc peeling,
  (3) Johnson J(14,2) cap mass,
  (4) HYP-3205/HYP-3224 Toeplitz/covariance coordinates,
  (5) HYP-3227 trap-discharge conductance graph.
```

The immediate theorem-facing tasks are:

1. Match `E_7(tau)` q-coefficients to the S75 comb-overlap Gram basis and
   turn the match into explicit finite LP/Toeplitz certificate rows.
2. Test whether Fejer/Gamma0(7)/Gram slack dominates the HYP-3227 Green-only
   trap conductance weights and the precision-defect island.
3. Isolate the remaining S75 order-3 triple-overlap constants as named debt.
4. Interface the resulting finite certificate rows with HYP-3215's
   Fejer/Cohn-Elkies/polyhedron-flatness Gap A, while keeping the induction
   base verification as a separate global-route obligation.
5. Treat the Stark/Dirichlet-L, Beraha/Tutte, Mahler/Lehmer, and subshift
   observations as sidecars until each preserves a named LRC predicate.

## Tournament Analysis

The tournament vertices are proof carriers, not runners, roots, arcs, or raw
constants:

```text
fejer_demoivre_kernel
comb_overlap_gram_kernel
johnson_14_cap_kernel
gamma0_7_eisenstein_sidecar
toeplitz_green_conductance_bridge
subshift_transfer_operator
dirichlet_l_stark_sidecar
beraha_tutte_sidecar
mahler_lehmer_height_sidecar
raw_arithmetic_numerology_guardrail
```

Pairwise observable:

```text
certificate payload retained,
then formal checkability,
then negative risk.
```

Fingerprint:

```text
score_hist = {0:1, 1:1, 2:1, 3:1, 4:1, 5:1, 6:1, 7:1, 8:1, 9:1}
directed_3_cycles = 0
Hamiltonian path =
  fejer_demoivre_kernel
  -> comb_overlap_gram_kernel
  -> johnson_14_cap_kernel
  -> toeplitz_green_conductance_bridge
  -> gamma0_7_eisenstein_sidecar
  -> subshift_transfer_operator
  -> beraha_tutte_sidecar
  -> dirichlet_l_stark_sidecar
  -> mahler_lehmer_height_sidecar
  -> raw_arithmetic_numerology_guardrail
```

Assumption challenge:

```text
Rejected vertices:
  runners, arcs, roots, q-coefficients, and famous constants alone.

Chosen vertices:
  proof carriers, because the preserved LRC predicate is certificate payload.

Preserves:
  Fejer positivity/PSD, double-zero sharpness, cap-vs-sector split,
  comb-overlap Gram positivity, conductor-7 labels, Toeplitz/Green trap
  discharge.

Destroys if scalarized:
  literal runner identity, height/numerology, raw closed forms,
  and unlabelled extremality.
```

## Sources Checked

- DLMF section 1.15 for Fejer kernels and summability context.
- Szmidt--Urbanowicz--Zagier, "Congruences among generalized Bernoulli
  numbers", for the generalized Bernoulli / Dirichlet-L relation.
- Royle's Beraha-number chromatic-root note for Beraha/Tutte context.

-> HYP-3229, HYP-3227, HYP-3215, HYP-3214, HYP-3213, HYP-3212, HYP-3205,
HYP-3203, HYP-3201, HYP-3162, HYP-3161, HYP-3160, THM-577, T1329,
LTI-329, LTT-229, OPEN-Q-108.
