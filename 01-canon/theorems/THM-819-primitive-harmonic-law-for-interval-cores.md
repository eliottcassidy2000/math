---
id: THM-819
title: THE PRIMITIVE HARMONIC LAW — for every k ≥ 1, the interval core's good measure at the sub-tight threshold is EXACTLY the primitive harmonic sum — m({1,…,k}; λ = 1/(k+2)) = (2/((k+1)(k+2))) · Σ_{1≤u≤k, gcd(u,k+1)=1} 1/u. COROLLARY: it equals H_k/C(k+2,2) (THM-805's harmonic-measure form) ⟺ k+1 is PRIME — the "mod-6 law" of HYP-6885 was the shadow of primality (primes > 3 are ≡ ±1 mod 6, so n = k+2 = p+1 ≡ 0,2 mod 6; the mod-6 exceptions n ∈ {10,16,22} are n−1 ∈ {9,15,21} composite). The deep well is the p = 13 instance: m({1..12}; 1/14) = 2H₁₂/(13·14) = H₁₂/91 (mac-mini THM-736's base measure — now an instance of a totient-structured law)
status: PROVED (three-step proof below: Dirichlet confinement + per-witness two-sided extents + unit-inversion bijection; every step elementary and self-contained) + REFEREED (exact ℚ, k = 1..30, zero exceptions; also re-derives every EQUAL/differ entry of the cont.8 ladder including the exact failure ratios)
source: kind-pasteur-2026-07-15-S128 (cont.9; owner: prove the mod-6 law)
depends_on: []
related:
  - THM-805 (the Tutte/harmonic bridge whose part (ii) this completes and corrects: the harmonic-number form is the prime case)
  - HYP-6885 (a) — RESOLVED by this theorem (the mod-6 characterization was correct data, wrong lens)
  - mac-mini THM-736 (the p=13 deep-well instance |G'| = H₁₂/91)
  - HYP-6865 (the staircase's electrical resistance H_{n−2}: matches the good measure exactly iff n−1 prime — the staircase network "knows" primality through the LRC)
verification:
  - 04-computation/thm819_primitive_harmonic_referee_kps_S128c9.py
  - 05-knowledge/results/thm819_primitive_harmonic_referee_kps_S128c9.out
---

# THM-819 — the primitive harmonic law

**Theorem.** For k ≥ 1, set q = k+1, δ = 1/(q(k+2)). Then

>  m({1,…,k}; 1/(k+2)) = 2δ · Σ_{1≤u≤k, gcd(u,q)=1} 1/u.

**Step 1 (Dirichlet confinement).** Let t be good (‖jt‖ ≥ 1/(k+2) for all j ≤ k). By Dirichlet's
theorem with Q = k+1 there is 1 ≤ j ≤ q with ‖jt‖ < 1/(k+2). Goodness forces j = q: t lies within
δ' := ‖qt‖/q < δ of some a/q. If d = gcd(a, q) > 1, then (q/d) ≤ k and ‖(q/d)t‖ = ‖a/d + ε/d‖ =
|ε|/d < 1/(k+2) (a/d ∈ ℤ), contradicting goodness. So **every good t lies within δ of a PRIMITIVE
a/q** — the good set is contained in the δ-neighborhoods of the φ-witnesses.

**Step 2 (per-witness extents).** Fix a primitive a/q and write t = a/q + s, |s| ≤ δ. For j ≤ k:
‖j·a/q‖ = (ja mod q)/q or its complement; the residues ja mod q run over nonzero residues. The
UNIQUE binding speeds are j⁺ ≡ a^{−1} (mod q) (residue +1: ‖j⁺t‖ = 1/q + j⁺s) and j⁻ = q − j⁺
(residue −1: ‖j⁻t‖ = 1/q − j⁻s). All other speeds have ‖ja/q‖ ≥ 2/q, and over |s| ≤ δ they move
at most kδ = k/(q(k+2)), with 2/q − 1/(k+2) = (k+3)/(q(k+2)) > kδ — never binding. Speeds sharing
a factor with q sit at ≥ 2/q automatically (their residues are nonzero multiples of the gcd).
Hence the good part of the neighborhood is exactly a/q − δ/j⁺ ≤ t ≤ a/q + δ/j⁻ (each side: the
binding constraint reaches 1/(k+2) precisely at excursion (1/q − 1/(k+2))/j± = δ·(k+2−q)... =
δ/j± after simplifying with 1/q − 1/(k+2) = δ). Neighborhoods are disjoint (2δ < 1/q for k ≥ 2;
k=1 checked directly), so each witness contributes length δ(1/j⁺ + 1/j⁻).

**Step 3 (unit inversion).** As a runs over the primitive residues mod q, j⁺ = a^{−1} mod q runs
over the primitive residues bijectively, and j⁻ = q − j⁺ likewise. Therefore
Σ_a (1/j⁺ + 1/j⁻) = 2·Σ_{u primitive mod q, u ≤ k} 1/u, and m = 2δ·Σ_{gcd(u,q)=1} 1/u. ∎

**Corollary (the prime law).** Σ_{u≤k, gcd(u,k+1)=1} 1/u = H_k ⟺ every u ≤ k is coprime to k+1
⟺ k+1 prime. So THM-805(ii)'s harmonic form m = H_k/C(k+2,2) holds iff n−1 = k+1 is prime —
{4,6,8,12,14,18,20} are exactly the n ≤ 22 with n−1 prime, and the deep well (n=14, p=13) is one
of them. The failure ratios of the cont.8 referee are reproduced exactly as H^prim/H_k.

**LRC reading.** The interval-core good measure — the base quantity of every peel bound on
interval bodies (THM-733/735/738/741 ladder) — is a TOTIENT-weighted harmonic object: only the
speeds coprime to k+1 contribute measure at the sub-tight threshold, each in proportion 1/u, and
each witness interval is "powered" by the inverse pair (a, a^{−1}). Via HYP-6865, the staircase's
electrical resistance equals this measure precisely when n−1 is prime: the Smith network of the
staircase detects primality of n−1 through the Lonely Runner.

## Evidence log

- [x] referee (tutte-harmonic ladder, cont.8 + this session k=1..30): exact ℚ equality of the
      primitive formula at every k; prime/composite split exactly as the corollary states
- [x] j=4 run co-harvest at time of writing: 49/2002 bodies complete, 49/49 clean (fails=[])
