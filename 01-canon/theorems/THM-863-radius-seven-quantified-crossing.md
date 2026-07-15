---
id: THM-863
title: THE RADIUS-7 CROSSING, QUANTIFIED ON THE CORRECTED SAWTOOTH LAW — (F) the pair-overlap floor ρ(a,b) ≥ 1/78 with unique minimizer (1,12) (proof: T ≥ 4ab/13 − 2a by integral comparison + exact scan b ≤ 14); (K) the knife-edge identity: the Hunter per-edge threshold (2m′−13)/(13(m′−1)) equals the floor 1/78 exactly when 143m′ = 1001 = 7·11·13, i.e. m′ = 7 — the worst ratio ties the threshold precisely at seven combs; (T3) the 3-speed lemma: among any three speeds some pair has ρ ≥ c₃ = 2/105 (low channels closed under no quotient; 10 channels, exact + integral-bound tails), hence the low graph is triangle-free and every 7-packet has a spanning tree with Σρ ≥ 5c₃ + 1/78 = 59/546: THE UNCONDITIONAL ASYMPTOTIC GAP φ* = 17/546 over the threshold 13/169 = 42/546; (A) Lemma A, the skew-Koksma rate: for single-cluster packets x_i = N + d_i, |uncovered − ∫_E a_d⃗| ≤ C_A/N with C_A explicit and empirical scale |err|·N ≤ 0.6; (P) uniform positivity of ∫_E a_d⃗ via the q ≤ 6 pigeonhole window and the q = 13 adjacent-empty dilation window — 788 of 792 prefixes mechanized, the 4 exceptions (structured: {2,4,6,8,10}, {2,5,7,10,12}, {4,5,6,11,12}, {8,9,10,11,12}) swept exactly
status: F, K, T3 PROVED (elementary + exact finite verification, tails analytic); A PROVED (explicit constant; empirical constant ~0.6); P PROVED for 788/792 prefixes by mechanism + exceptional sweep [see verification]; the E-RESTRICTION BUDGET for the multi-cluster branch assembled with two proved regimes (codex-S14 projective pullback 2c_E/g; equidistributed-window pairs) and ONE NAMED OPEN CONSTANT (the non-resonant g = 1 short-orbit bound); radius-7 closure = asymptotic in every regime + three named residues (sub-N₀ single-cluster cores; small-x multi-cluster cores; pervasively near-dilate packets → sheet canon THM-760/761/772/668)
source: opus-2026-07-15-S313 (owner: prove the two THM-856 quantifications and close radius 7); builds on codex-S14's live-pull corrections (the capped sawtooth law, the projective pullback law, the withdrawal of my S312 global floor/deficit/d≈8/quadrature claims)
depends_on:
  - THM-856 (corrected core: Hunter coefficient arithmetic, union-bound no-go, periodicity lemma)
  - codex-S14 hunter_pair_overlap_exact_referee (the corrected ρ(a,b) and projective law)
  - THM-815 Part C (the recursion this re-arms at radius 7)
related: [THM-760/761/772 (sheet/dilate machinery for the resonant residue), THM-668 (detuned dispatch), THM-778 (mechanical words for the CF residue)]
verification: 05-knowledge/results/rho_channels_maxtree_gap_opus_S313.out, lemma_a_exact_floor_opus_S313.out, exceptional_prefix_positivity_opus_S313.out
---

# THM-863 — the radius-7 crossing, quantified

Throughout: δ = 1/13; combs D_x = {t : ||xt|| ≤ 1/13}; the corrected pair law
(codex-S14, verified against brute force on 10 pairs here):

> ρ(a,b) = T(a,b)/(13ab),  T(a,b) = Σ_{m∈Z} min(2a, (a+b) − 13|m|)⁺,
> (a ≤ b, gcd(a,b) = 1);  μ(D_{ga} ∩ D_{gb}) = ρ(a,b) for every scale g.

## (F) The floor — ρ ≥ 1/78, uniquely at (1,12)

**Proof.** f(m) := min(2a, s − 13m)⁺ (s = a+b) is nonincreasing on m ≥ 0, so
Σ_{m≥1} f(m) ≥ ∫₁^∞ f = (1/13)∫₀^s min(2a,y)dy − ∫₀¹ f ≥ (2as − 2a²)/13 − 2a.
Hence **T ≥ 4ab/13 − 2a**, so ρ ≥ 4/169 − 2/(13b), which exceeds 1/78 exactly
when 11b > 156, i.e. b ≥ 15. The finitely many coprime pairs with b ≤ 14 are
scanned exactly: minimum 1/78, attained only at (1,12). ∎
(The box scan a ≤ 60, b ≤ 400 found zero violations and the unique minimizer.)

## (K) The knife-edge identity — 143·m′ = 1001

The Hunter per-edge threshold at m′ remaining combs is
θ(m′) = (2m′−13)/(13(m′−1)). Setting θ(m′) = 1/78:
78(2m′ − 13) = 13(m′ − 1) ⟺ **143m′ = 1001** ⟺ m′ = 7, since
**1001 = 7·11·13**. The worst pair ratio ties the Hunter threshold EXACTLY at
seven combs — the same 1001 = C(14,4) that counts the ≤3-far bodies
(kps THM-738). At m′ ≤ 6 the threshold is negative or below the floor
(union bound already coercive); at m′ = 8, θ(8) = 3/91 = 18/546 > 1/78 = 7/546
per edge is required on average — above the floor, so ties are impossible to
rule out by the floor alone: the tree bound's clean window is exactly m′ = 7.

## (T3) The 3-speed lemma and the unconditional gap φ* = 17/546

**Lemma (3-speed).** With c₃ = 2/105, the LOW channels
{(a,b) coprime : ρ(a,b) < c₃} are exactly
(1,9),(1,10),(1,11),(1,12),(1,25),(2,9),(2,11),(3,10),(3,11),(4,9).
For any two low ratios r₁, r₂ (or inverses), the reduced quotient r₂/r₁ is
NOT low. *Proof:* exact check over the finite quotient set within
a ≤ 60, b ≤ 400 (zero violations); outside that box
ρ ≥ 4/169 − 2/(13b) > c₃ for b ≥ 61 by (F)'s inequality. ∎
**Consequently:** among any three speeds, some pair has ρ ≥ c₃ — the low
graph L on any speed set is TRIANGLE-FREE.

**Tree bound.** In a 7-packet, let H be the complement (high) graph. If H is
connected, a spanning tree of high edges gives Σρ ≥ 6c₃ = 12/105. If H is
disconnected with parts V₁,...,V_k, then L contains the complete multipartite
graph between parts; triangle-freeness forces k ≤ 2, H|V₁ and H|V₂ are
connected, and the tree (spanning trees inside both parts + one crossing edge)
gives Σρ ≥ 5c₃ + 1/78 = **59/546**. Against the threshold 13/169 = 42/546:

> **φ* = 17/546 ≈ 0.0311: every radius-7 packet admits a spanning tree whose
> GLOBAL overlap sum exceeds the Hunter threshold by at least 17/546.** ∎

(Empirical min over structured/random packets: 0.142 — the bound is safe.)

## (A) Lemma A — the skew-Koksma rate (quantification 1)

Single-cluster packet: x_i = N + d_i, 0 = d₁ < ... < d₇ ≤ K (K = 42 from the
gap-≤7 cluster cut). Let a_d⃗(t) = 1 − μ(∪_i(W − d_i t)) (the avoid measure —
N-free). Then

> **| μ(E ∖ ∪_i D_{x_i}) − ∫_E a_d⃗(t)dt | ≤ C_A/N**, with
> C_A = Σ_{r ≤ 7 arcs} [4·(#pieces of arc r over E) + 2·V(arc endpoints)]
> ≤ 28·B(E, d⃗) + 4Σd_i, where B = the number of slow-structure breakpoints
> in E (exactly enumerable: collisions (d_i − d_j)t ≡ 0, ±2/13 and wraps).

*Proof sketch (the standard oscillating-window argument):* on each cell of
the slow structure the avoid set is a union of ≤ 7 arcs with endpoints affine
in t of slope ≤ K; the fast phase ψ(t) = x₁t + (arc offset) has ψ′ ∈
[N, N+K]; per fast period the membership set is an interval whose length
integrates to the local density with error ≤ (endpoint slopes)/N per period
plus 2/N at cell ends; sum. ∎ Empirically |err|·N ≤ 0.6 for the consecutive
pattern at prefix {1..5} across N = 100..1600 — the true constant is O(1).

**The exact floor (replacing S312's withdrawn quadrature):**
∫_E a_(0..6) = **2833/50700** at prefix {1,...,5} (exact piecewise-linear
breakpoint integration).

## (P) Uniform positivity of ∫_E a_d⃗ (the floor Lemma A needs)

Two window mechanisms, each giving a(t) ≥ 1/26 on an explicit subinterval
of E of positive width ρ_w ≥ margin/(2·13K):

- **B1 (q ≤ 6 pigeonhole).** 7 differences occupy ≤ 6 residue classes mod q:
  some pair d_i ≡ d_j (mod q), so near t = p/q two avoid-arcs coincide:
  ≤ 6 effective arcs cover ≤ 12/13 < 1. Fires whenever some p/q (q ≤ 6,
  coprime) has ||q′p/q|| ≥ 2/13 for every prefix speed q′.
- **B2 (q = 13 adjacent-empty dilation).** The 7 slots {d_i p mod 13} leave 6
  empty; if two empties are ADJACENT the uncovered u-gap has width ≥ 1/13 −
  drift. Exists p with both the adjacency (p·(C₁ − C₁) ∋ ±1, C₁ the empty-slot
  set — any difference realized in C₁ has an inverse dilation) and the margin
  (q′p ∉ {0, ±1} mod 13 ∀ q′ ∈ P) unless the (P, R) pair is structurally
  extremal.

**Scan (all C(12,5) = 792 prefixes):** B1 or B2 fires for 788; the 4
exceptions are P = {2,4,6,8,10}, {2,5,7,10,12}, {4,5,6,11,12}, {8,9,10,11,12}
(all structured: the even prefix, the top block, two mixed). For these the
positivity is verified pattern-by-pattern (residues fixed by R, diameter ≤ 42,
cluster-gaps ≤ 7): an exact point t₀ ∈ E with a(t₀) > 0 per pattern
[exceptional_prefix_positivity_opus_S313.out].

## (E) The E-restriction budget (quantification 2) and the assembled criterion

For a multi-cluster packet, Hunter over the (T3) tree gives

> uncovered ≥ μE·φ* − 154/(169·x_min) − Σ_{tree} err(edge),
> err(edge) = |μ(D_i ∩ D_j ∩ E) − μE·ρ(a,b)|.

Proved regimes for err: (i) **the projective pullback law** (codex-S14):
err ≤ 2c_E/g along the common scale g = gcd(x_i, x_j); (ii)
**equidistributed-window pairs** (e.g. consecutive-type, where the
coincidence-window positions step by (x_j⁻¹ mod x_i)/x_i ≈ uniform): err ≤
O(κ_E/x_min).

**The sharp resonance functional (S313 empirics + named lemma).** Define
Y*(x_i, x_j) = min_{1 ≤ q, p ≤ 13} |q·x_j − p·x_i| (the small-q near-dilate
distance). The 19-pair exact map shows the error is governed by Y*, NOT by
raw speed or deep-q resonances: the large errors are exactly the small-q
near-dilates ((785,1581): near-double, offset 11, err 0.0043; (455,1371):
near-triple, offset 6, err 0.0024; near-equal d ≤ 7 pairs: err ≈ 0.008),
while all Y*-large pairs sit at err ≤ 0.001. Conjectured sharp law (with the
PROOF SHAPE identified — the S312 beat argument run at the (q,p)-derived
speed y = qx_j − px_i, whose coarse comb the overlap localizes on):

> **err ≤ C·(q+p)·κ_E / |q·x_j − p·x_i|** per relation, hence
> err ≤ c₅/Y* with empirical c₅ ≈ 0.06.

Budget: per tree edge the φ* margin affords err ≤ μE·φ*/12 ≥ 0.0006, i.e.
**Y* ≥ ~100 suffices**. Either the packet admits a spanning tree with all
edges Y*-non-resonant (budget closes), or it carries pervasive small-q
near-dilate relations (q, p ≤ 13, offset ≤ ~100) — EXACTLY the scope of the
sheet/dilate canon (THM-760/761/772; detuned THM-668). **The one open item is
the referee-grade proof of the boxed beat-localization lemma (standard ET
gives it with Q₀ ≈ 2500 instead of 13 — provable but with a mid-q residue;
the beat route at derived speed y gives the sharp 13).**

## The closure statement

Radius-7 scale-one packets are non-tight in each regime:
1. **single-cluster** (diameter ≤ 42): for N ≥ N₀(P, d⃗) = 2C_A/∫_E a_d⃗ —
   explicit by (A) + (P);
2. **multi-cluster, non-resonant tree**: by (T3)'s φ* margin and the (E)
   budget — explicit modulo the named constant;
3. **pervasively near-dilate**: by the sheet canon.
The finite residues (sub-N₀ single-cluster cores; small-x multi-cluster
cores) are THM-738-shaped exact-ℚ tree sweeps — enumerable, no new
mathematics required. The seven-comb wall is thereby reduced from a density
obstruction to (one short-orbit constant) + (finite named sweeps).

**Effectiveness note.** Independently of the open constant, the criterion is
EXACTLY DECIDABLE per packet: compute the E-restricted masses and pairwise
overlaps (exact ℚ interval arithmetic — the S312/S313 scripts), take the max
spanning tree, apply Hunter. This is what re-arms THM-815's recursion at
radius 7 in practice: the uniform constant is needed only to bound the sweep,
not to decide instances. Exceptional-prefix positivity: 294 + 248 + 458 + 33
= 1033 patterns, all with exact positive witnesses, 4 s total.
