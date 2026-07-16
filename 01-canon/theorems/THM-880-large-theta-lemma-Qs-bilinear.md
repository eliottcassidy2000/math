---
id: THM-880
title: THE LARGE-θ LEMMA AND THE EXACT Q_s BILINEAR FORM — (L1) Σ_{ℓ≥1} sin²(πℓθ)/ℓ² = π²θ(1−θ)/2; (L2) the Möbius-sinc closed form V(e) = Σ_{j mod e} tent(j/e) = e/49 + {e/7}(1−{e/7})/e (each aliasing sum is ONE large-θ term; E(e) ≤ 1/(4e), zero iff 7|e); (L3) Σ_{ℓ≠0} e(ℓt)/ℓ² = 2π²B₂({t}); hence THM-729's density second moment is an EXACT FINITE RATIONAL BILINEAR FORM: Q_s = −2π² Σ_{k,k'} ε_k ε_{k'} {wΔ_{kk'}}(1−{wΔ_{kk'}}) — no ℓ-sum, no truncation; the diagonal recovers THM-729's rigorous backbone and the "off-diagonal cancellation" is sign cancellation in ε_kε_{k'}
status: (L1)-(L3) PROVED (Fourier of B₂; L2 verified EXACT e ≤ 2000); the Q_s identity PROVED (expand |U|², apply L3, use Σε = 0; verified against ℓ-sums on random arc systems and on THM-729-style 7-section clusters, now with EXACT rational Q_s); Q_s = O(diam) still OPEN as a uniform theorem — now reduced to sign-equidistribution of the large-θ form, with the O(1)·diam law verified EXACTLY (ratios 2.9–8.0, no LMAX artifacts); dipole-transport bound proved but weak in the observed regime (w·maxwidth ≈ 6–29)
source: klein-2026-07-16-S314; upgrades THM-729 (its named open: "a clean bound must exhibit the arc-differencing self-cancellation" — the bilinear form IS that exhibition)
depends_on: [THM-727/728/729]
verification: 04-computation/clock_moduli_largetheta_Qs_klein_S314.py (8/8)
---

# THM-880 — the large-θ lemma; Q_s as a finite bilinear form

**(L1)** Σ_{ℓ≥1} sin²(πℓθ)/ℓ² = π²θ(1−θ)/2. *Proof:* sin² = (1−cos)/2 and
Σcos(2πℓθ)/ℓ² = π²B₂({θ}) (the Fourier series of the Bernoulli B₂); subtract. ∎
**(L2)** V(e) = e·Σ_{e|n} t̂_n with t̂_n = (sin(πn/7)/(πn))² (tent = 1_I ⋆ 1_I, Fejér):
V(e) = e/49 + (2/e)Σ_m sin²(πme/7)/(π²m²) = e/49 + {e/7}(1−{e/7})/e by (L1). So the whole
Möbius/aliasing error is ONE large-θ term: 0 ≤ E(e) ≤ 1/(4e), = 0 iff 7 | e. (Exact for
e ≤ 2000.) This powers THM-878's divisor arithmetic and LEM-020's (RAM).
**(Q)** U(N) = Σ_k ε_k e(−N p_k), Σε_k = 0 (arc endpoints, signed). Then
Q_s = Σ_{ℓ≠0}|U(ℓw)|²/ℓ² = Σ_{k,k'}ε_kε_{k'}·2π²B₂({w(p_k−p_{k'})}) by (L3); the constant
1/6 dies against (Σε)² = 0, leaving
> **Q_s = −2π² Σ_{k,k'} ε_k ε_{k'} {w(p_k−p_{k'})}(1 − {w(p_k−p_{k'})})** — exact, finite,
> rational for rational data. Diagonal (k, k' in the same arc) = THM-729's rigorous
> 2π²Σ{w·u}(1−{w·u}); everything else is the off-diagonal WITH ITS SIGNS — the
> self-cancellation THM-729 asked to exhibit, exhibited.

**Consequences.** (i) All Q_s numerics are now exact-rational and decide-shaped (Lean
batch candidate); the O(1)·diam law re-verified with zero truncation error. (ii) The
dipole-transport bound Q_s ≤ 2π²Σ_{i,j}min(1/2, 2w·min(u_i,u_j)) is rigorous but weak in
the observed regime (w·u = O(10)): sharp O(diam) is genuinely about sign-equidistribution
of {wΔ} over the endpoint-difference multiset — the honest remaining content, now stated
as a finite exponential-sum-free question about one bilinear form.
