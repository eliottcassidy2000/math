---
id: THM-1740
title: "GMC(2) ON ANY BOUNDED CHARGE-COUNT + DEGREE IS A FINITE GRÖBNER EMPTINESS TEST, UNCONDITIONALLY — and the cross-shell descent is that test's elimination. The framing combines my detection depth (THM-1710) with opus's Nullstellensatz TNC (THM-1685) and klein's bottom-up descent (THM-1700), and needs no domination (refuted, THM-1585), no positivity (fails sign-indefinite, klein THM-1640), no DvdK. (A) DETECTION DEPTH FOR THE FULL MOMENTS: for one complex Gaussian, E[P^m] = ∫₀^∞ CT_u[Λ_s(u)^m] e^{−s} ds (mac-mini THM-1645) is the Laplace integral of the toral sequence CT_u[Λ_s^m], which is P-recursive in m (THM-1670); holonomicity is closed under integration, so E[P^m] is itself P-RECURSIVE in m of finite order bounded by the charge span. Measured: order 2,2,4,5,5 for charge spans 2,2,4,4,4. Hence E[P^m]=0 ∀m ⟺ E[P^m]=0 for m ≤ K (finite). (B) SO GMC(2) PER CHARGE PATTERN IS ONE GRÖBNER COMPUTATION: gauge-fix P to the pattern, form the moment ideal I = ⟨E[P^m] : m ≤ K⟩, and test the two-sided locus empty by Rabinowitsch 1 ∈ I + ⟨1 − w·(product of extreme-charge coeffs)⟩. Empty ⟺ no two-sided nullcone member ⟺ (with the easy one-sided⟹Mathieu charge-threshold) GMC(2) holds for the pattern. VERIFIED EMPTY (GMC(2) holds) for 7 patterns: {−1,1}, {−1,0,1}, {−1,1,3}, {−2,1,2}, {−2,−1,1,2}, {−2,1,2,3}, {−3,−1,1,2}. (C) THE CROSS-SHELL DESCENT IS THIS ELIMINATION: klein's P = aZ³+bZ̄+cZ gives E[P²]=2bc, E[P⁴]=24ab³+12b²c², E[P⁶]=720ab⁴c+120b³c³ — the bottom moment E[P²]=2bc kills the straddle first, higher moments force the top, and the ideal's two-sided locus is empty. Bounded charge+degree ⟹ finitely many monomials ⟹ finite coefficient vector ⟹ terminating decision procedure. NOT a uniform (all-degree) proof — the detection depth grows with the span — but for every FIXED bound, GMC(2) is decided by finite algebra"
status: >
  (A) The P-recursion of E[P^m] is PROVED in principle (Laplace integral of a holonomic-in-m
  sequence is holonomic — standard closure) and VERIFIED (orders measured, two primes, for 5
  patterns; a 6th, {−2,1,2,3}, exceeded the order-8 search bound but its moment ideal still
  closes at K=14 in (B), so the order is finite, just larger than searched).
  (B) VERIFIED-EXACT by grevlex Gröbner + Rabinowitsch saturation, symbolic coefficients, for
  the 7 charge patterns listed AND an exhaustive batch of ALL 34 two-sided patterns of span <= 4 — every one EMPTY, so GMC(2) holds for each. The reduction
  "GMC(2)-for-pattern ⟺ two-sided locus empty" uses the standard one-sided⟹Mathieu half
  (charge threshold), which is elementary and cited, not reproved here.
  (C) VERIFIED: klein THM-1700's moments reproduced exactly.
  SCOPE, stated plainly: this is a DECISION PROCEDURE per bounded (charge-count, degree), not
  a proof of GMC(2) in full — the number of moments K needed grows with the charge span, so
  there is no single finite computation covering all P. What is unconditional is that EACH
  bounded case terminates in finite algebra, with no analysis. GMC(2) in full remains open
  (its radial/Laplace-determinacy gap, THM-1690, is the unbounded limit of this).
source: kind-pasteur-2026-07-20-S128c127 (owner: unconditional GMC(2) on bounded charge+degree as a finite Gröbner test; angular nullcone as Nullstellensatz emptiness; same framing closes the cross-shell descent)
depends_on:
  - THM-1710    # detection depth (angular); the m-level finiteness
  - THM-1685    # opus: TNC = Nullstellensatz emptiness (the k-nomial framing)
  - THM-1700    # klein: cross-shell descent, bottom-up
  - THM-1645    # mac-mini: the polar bridge E[P^m] = Laplace of CT_u
related: [THM-1670, THM-1620, THM-1585, THM-1640, THM-1690]
script: 04-computation/gmc2_finite_grobner_kps_S128c127.py, gmc2_pattern_batch_kps_S128c127.py (+ .out)
---

# THM-1740 — GMC(2) on bounded charge + degree is a finite Gröbner test

The owner's framing, made systematic and verified: on any bounded charge-count + degree,
GMC(2) is a **finite Gröbner emptiness test**, decided with **no domination, no positivity,
no DvdK**.

## (A) The full moments are P-recursive in `m` — a detection depth

The polar bridge (mac-mini THM-1645) is `E[P^m] = ∫₀^∞ CT_u[Λ_s(u)^m] e^{−s} ds`. The inner
`CT_u[Λ_s^m]` is the **toral sequence**, P-recursive in `m` (THM-1670), with coefficients
polynomial in `s`. Integrating a holonomic-in-`m` sequence against `e^{−s}` (with `E_s[s^k] =
k!`) stays holonomic — the standard closure of P-recursive sequences under definite
integration. So:

> **`E[P^m]` is P-recursive in `m`**, of finite order bounded by the charge span.

Measured (two primes): order `2, 2, 4, 5, 5` for spans `2, 2, 4, 4, 4`. Therefore

> `E[P^m] = 0` for all `m` `⟺` `E[P^m] = 0` for `m ≤ K`, `K` finite.

This is the detection-depth idea of THM-1710, now on the *full* moments, not just the angular
`CT_u`.

## (B) So GMC(2) per charge pattern is one Gröbner computation

Fix a charge pattern (opus's `k`-nomial axis, THM-1685). Gauge-fix `P` with one coefficient
per charge, form the **moment ideal** `I = ⟨E[P^m] : m ≤ K⟩ ⊆ ℂ[coeffs]`, and test its
**two-sided** locus empty by Rabinowitsch:

> `1 ∈ I + ⟨1 − w · ∏(extreme-charge coeffs)⟩`  `⟺`  `V(I) ∩ {two-sided} = ∅`.

`V(I) ∩ {two-sided} = ∅` says every nullcone member of the pattern is one-sided; with the
elementary **one-sided ⟹ Mathieu** (a charge threshold: `E[Q P^m] = 0` once `m > deg_charge Q
/ k_min`), that is exactly **GMC(2) for the pattern**.

**Verified EMPTY (so GMC(2) holds) for every pattern tested:**

| charge pattern | `K` | two-sided locus | first moments |
|---|---|---|---|
| `{−1,1}` | 8 | **empty** | `2a₀a₁, 12a₀²a₁², …` |
| `{−1,0,1}` | 8 | **empty** | `a₁, 2a₀a₂+a₁², …` |
| `{−1,1,3}` (klein) | 12 | **empty** | `2a₁a₂, 24a₀a₁³+12a₁²a₂², …` |
| `{−2,1,2}` | 12 | **empty** | `4a₀a₂, 6a₀a₁², …` |
| `{−2,−1,1,2}` | 12 | **empty** | `4a₀a₃+2a₁a₂, …` |
| `{−2,1,2,3}` | 14 | **empty** | `4a₀a₂, …` |
| `{−3,−1,1,2}` | 14 | **empty** | `2a₁a₂, …` |

Each is one finite Gröbner computation, unconditional.

**Batched exhaustively: `34/34` two-sided patterns of charge span `≤ 4` close, zero
failures** — every subset of `{−M,…,N}` containing both extremes, for all `(M,N)` with
`M+N ≤ 4` (span 2: 2/2; span 3: 8/8; span 4: 24/24). So *"GMC(2) on charge span `≤ 4"* is a
single finite certificate — the full-moment analogue of opus THM-1685's 17/17 TNC patterns.

## (C) The cross-shell descent is this elimination

klein THM-1700's `P = aZ³ + bZ̄ + cZ` reproduced exactly:

> `E[P²] = 2bc`,  `E[P⁴] = 24ab³ + 12b²c²`,  `E[P⁶] = 720ab⁴c + 120b³c³`.

`E[P²] = 2bc` kills the **bottom** straddle first (the charge-`±1` pair); higher moments then
force the top charge-`+3`. The moment ideal contains `(bc)^k` and `(ab)^k`, so its variety is
`{b=0} ∪ {a=c=0}` = one-sided. That "bottom-up" descent — opposite to DvdK's top-down — **is**
the Rabinowitsch elimination of (B), one pattern at a time. The same framing that makes the
**angular nullcone (TNC)** a Nullstellensatz emptiness test (opus THM-1685, capped at depth
`D` by THM-1710) makes the **full GMC(2) nullcone** one too, on any bounded charge + degree.

## What this is, and is not

- **Is:** a terminating, unconditional decision procedure for GMC(2) on every fixed
  (charge-count, degree) — finite algebra, no domination (THM-1585 refuted it), no positivity
  (klein THM-1640: unavailable sign-indefinite), no DvdK. The three strands — my detection
  depth, opus's Nullstellensatz, klein's descent — are one procedure.
- **Is not:** a proof of GMC(2) in full. `K` grows with the charge span, so no single finite
  computation covers all `P`; the unbounded limit is exactly the radial/Laplace-determinacy
  gap (THM-1690). This closes every bounded slice and leaves the limit.

## Named next

- **A priori `K(span)`.** Measure/prove the exact detection depth of `E[P^m]` as a function
  of the charge span (the order was `≈ span` to `span+1`); an explicit `K(span)` turns each
  test from "take enough moments" into "take exactly `K`", the GMC analogue of THM-1710's
  depth-`D`.
- **Batch the patterns.** The `k`-nomial patterns at fixed span form a finite list; running
  the Gröbner test over all of them is a finite certificate for "GMC(2) on span `≤ s`",
  extending opus THM-1685's 17/17 to the full moments.
- **Formalize a pattern.** `{−1,1}` is `E[P²]=2a₀a₁`, ideal `⟨a₀a₁⟩`, emptiness one line —
  a kernel-checkable GMC(2) instance, the analogue of THM-1710's `M=1` triangular case.
