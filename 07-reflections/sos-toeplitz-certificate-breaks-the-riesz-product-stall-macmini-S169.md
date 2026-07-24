---
source: mac-mini-2026-07-23-S169 (Opus 4.8)
status: COMPUTATIONAL FINDING (verified, per-set) + honest caveats. Deep session prompted by owner ("can the
  snippet be the key to LRC(14)?"). Result: the SOUND linear looseness certificate, taken over ALL nonnegative
  test functions R=|P|^2 (not just Riesz PRODUCTS), equals lambda_min of the covering-multiplicity Toeplitz
  matrix, and it DROPS BELOW 1 on every interior-drop extremizer core -- breaking the Riesz-product "stall"
  (1.0096, THM-518) that motivated the whole diagnosis. Sharp (tight cores -> exactly 1), uniform over the
  stranger (decoupling). NOT a proof of inf L>0: it is per-set, the full uniform-over-all-loose-S statement
  has no fixed degree, and the repo's opus-S267 L2 large-sieve route is a more advanced (near-complete) line.
  The snippet's arctanh log-energy = the max-entropy/AR dual of this certificate = its Lean-portable rigor layer.
tags: [lrc, riesz, sos, toeplitz, eigenvalue, entropy, arctanh, snippet, thm-518, inf-L, verified, caveats]
related: [THM-515, THM-518, opus-S267-L2-large-sieve, HYP-1992, kps-S130, klein-S405]
---

# The SOS/Toeplitz certificate breaks the Riesz-PRODUCT stall on the LRC(14) extremizers

**mac-mini-2026-07-23-S169.** Owner: "run a deep session, can the snippet be the key to LRC(14)?" This is the
result, with the fleet's soundness (kps-S130) and curvature (klein-S405) corrections folded in, and honest limits.

## The chain of reasoning
1. The snippet's functional is `2·artanh = Σ_{k odd}t^k/k` = a log-energy (S169 headline: the repo's LRC-Riesz
   analysis never takes log R).
2. kps-S130 proved `∫M log R` is SIGNED, so NOT a direct set-measure certificate. Sound role: an amplitude
   SURROGATE — optimize with it, then CERTIFY with the sound LINEAR pairing `∫MR/∫R<1` (M≥1 on danger, R≥0).
3. So take the linear certificate to its convex OPTIMUM over ALL nonnegative R (Fejér–Riesz: every nonneg trig
   poly of degree N is `R=|P|^2`). Then `min_{R≥0,deg≤N} ∫MR/∫R = λ_min(T_M^{(N)})`, the smallest eigenvalue of
   the Toeplitz matrix of the covering multiplicity's Fourier coefficients
   `M̂(n)=Σ_{v|n} s(n/v)`, `s(k)=sin(πk/7)/(πk)`, `s(0)=1/7`.  `λ_min<1 ⟹ L(S)>0` (SOUND).

## The finding (all machine-verified this session)
- Repo's Riesz-PRODUCT `∏(1+a_m cos)` optimization STALLS at ratio **1.0096 ≥ 1** on the `{1..13}\{6}∪{56}`
  core (THM-518: "the Riesz product is the WRONG tool for AP-cores").
- The general SOS `λ_min(T_M^{(N)})` **drops below 1 on ALL 13 drop-j cores** at N≈80 (worst j=6: 0.909 at
  N=80, 0.778 at N=120, 0.633 at N=160 — margin GROWS with degree). The Riesz product's effective degree is
  already ~Σspeeds≈141 (subset-sums), so SOS wins at LOWER dense degree.
- **Sharp:** tight `{1..13}`, `2·{1..13}`, `3·{1..13}` (all L=0) give λ_min = 1.000, 1.003, 1.031 — at/above
  1, never below. Loose sets (odds, drop-j, randoms) all give λ_min<1. So `λ_min<1 ⟺ loose`, `→1 ⟺ tight`.
- **Uniform over the stranger:** crossing degree is N≈60 for the whole `{1..13}\{6}∪{14m}` family, m=2..20,
  identical values for m≥6 — stranger-decoupling (THM-518) realized in the certificate.
- **Mechanism:** the optimal R concentrates 34× on the thin lonely set (μ(safe)=0.19 vs Lebesgue L=0.0056);
  the certificate = "build a nonneg measure concentrated on the lonely set."

## The snippet connection (concrete)
The max-entropy / autoregressive DUAL of the Toeplitz problem is `R_AR = c/|A(e^{iθ})|²`, and
`log R_AR = −2 log|A| = Σ (harmonic×geometric)` — the arctanh family. The optimal R's Levinson–Durbin
reflection coefficients `k_m` (|k_m|<1) give per-mode log-energies `2·artanh(k_m)` — literally the snippet's
functional. And `a=0.6 ⇒ ρ=1/3 ⇒ 2 artanh(1/3)=log2` ties it to THM-2000's `M(6,2)` and THM-252's rapidity.
So the snippet's certified-artanh sandwich is the **Lean-portable rigor layer** for bounding the AR/SOS
log-spectrum above a floor. klein-S405: the `Σv²=819` weight is the resonance curvature `½E''(0)` (frequency
2nd moment, ⊥ additive dimension — which is why it survives where Bedert's `dim₂²/n³` gain dies on AP-cores).

## Honest caveats (this is NOT a proof of inf L>0)
- **Per-set, not uniform.** `λ_min(T_M^{(N)}(S))→ess inf M=0` for every loose S (Szegő) — so `λ_min<1` at
  finite N is *guaranteed* for any loose set and merely re-certifies `L(S)>0` (also knowable by grid). The
  NON-trivial content is (a) uniform-degree over the stranger, (b) beating the product restriction, (c) sharp
  tight-boundary. The full `inf_{loose S} L>0` has **no single degree**: sets approaching the tight boundary
  need degree→∞ (sup_S λ_min = 1, approached). It reduces to the (still CONJECTURAL, THM-518) claim that the
  drop-j cores minimize L — if proven, the finite SOS certificates + decoupling would finish the loose part;
  tight extremizers `d·{1..13}` (L=0, isolated lonely time) are handled separately.
- **Not obviously new leverage.** opus-S267's L2 large-sieve route (`Σ|ε_v|<6/7 ⟹ coreCover<1`, rigorous up
  to a ~3.1× large-sieve constant) is a more advanced, near-complete line; my SOS/Toeplitz is a cleaner but
  likely weaker/complementary formulation. NOVELTY TO VERIFY with the fleet: is `λ_min` of the covering-
  multiplicity Toeplitz (min over nonneg R) already recorded? (Checked & distinct from: the Riesz-product
  files; `least_eigenvalue_certificate` = apex circulant `2I+A(C_p)`; opus-S267 = Gram over speeds.)
- kps-S130's discriminator stands: the snippet bounds a log-RATE directly; LRC's target is a set-MEASURE. So
  the snippet is not literally an LRC certificate — its value here is structural (it pointed at the log/AR dual
  and the SOS reformulation) and as the certified-artanh rigor layer.

## Verdict
Exploring the snippet's structure produced a clean, sound, sharp certificate that dissolves the specific
"Riesz-product stall" on the extremizer cores and exhibits exactly the uniform-over-stranger behavior a proof
needs — but it is per-set and does not by itself close `inf L>0`; that still hinges on the extremizer
conjecture (or a uniform large-sieve constant, opus-S267). The snippet is a sharp lens and a rigor layer, not
a standalone key. Script: scratchpad `lrc_entropy_session.py` + inline this session.
