---
id: THM-517
title: The explicit Mayer cluster expansion of the OCF hard-core gas — the cumulants of log I(Ω,z), the excluded-volume integral |E(Ω)|, and the forbidden values {7,21} as the unique excluded 3-cluster (7) and a multiplicative gap (21); answers OPEN-Q-100(a)
status: PROVED (the cluster identities; the forbidden-value characterization) — VERIFIED exhaustively n=4,5,6 (33,856 tournaments, 0 failures) + n=7 sample, with an independent adversarial re-derivation (C4 counterexample pinning the c_k vs graphical-b_k distinction). The {7,21}-non-achievability-for-all-n is CITED from the already-closed Busch/A005517 result (HYP-2271, monad-compute-2026-06-06), here re-confirmed, not newly closed.
source: kind-pasteur-2026-06-15-S7 (workflow: develop→adversarial-verify→synthesize)
depends_on:
  - THM-002   # OCF: H = I(Ω,2), the hard-core lattice gas (T006)
  - THM-029   # H=7 impossible: Ω=K₃ non-realizable (the unique excluded cluster)
  - THM-510   # the free-gas reframe (H_free=3^{α₁}); OPEN-Q-100 posed the explicit expansion
related:
  - HYP-2271  # {7,21} = genus-2 semigroup gaps; strong-min=A005517 (Busch) — the non-achievability half, already closed
  - HYP-2544  # the LRC-14 singular series as a parallel (conditional) Mayer gas
  - HYP-2545  # the apex-prime-7 dual gate (OCF realizability vs LRC weight-zero)
  - HYP-2546  # E₇ as a Fano-substrate structural cousin (+ HYP-2530 retraction)
  - OPEN-Q-100
  - reflection: the-ocf-is-a-mayer-gas-lrc-is-its-conditional-twin-e7-is-its-fano-cousin-kps
# NOTE: renumbered from THM-515 (collision with mac-mini-S4's THM-515 = LRC lonely-measure θ-form)
---

# THM-517 — the explicit Mayer cluster expansion of the OCF gas

**One line.** `H(T)=I(Ω,2)=Σ_k α_k 2^k` (THM-002) is the grand partition function of a
single-species **hard-core lattice gas** on the odd-cycle conflict graph `Ω` at fugacity
`z=2`; its cluster (Ursell/Mayer) expansion is explicit, the forbidden value `7` is the
unique excluded cluster `Ω=K₃`, and `21` is a multiplicative gap.

## A. The gas and its cluster expansion (PROVED, verified 0 failures)

- **Ideal gas.** `Ω` edgeless ⟹ all odd cycles independent ⟹ `H_free=(1+2)^{α₁}=3^{α₁}`;
  `H ≤ 3^{α₁}` termwise, equality iff `Ω` edgeless. (Verified all 33,856 tournaments `n≤6`.)
- **Cumulant expansion.** `log I(Ω,z)=Σ_{k≥1} c_k z^k` with the **analytic Taylor cumulants**
  `c_k`, and `exp(Σ_k c_k z^k)=I(Ω,z)` is an **exact formal power-series identity** (verified
  0 failures). The order-2 cumulant is `c_2 = α_2 − α₁²/2 = −|E(Ω)| − α₁/2`.
- **The excluded-volume integral (the clean object).** `|E(Ω)| = C(α₁,2) − α_2` = the number
  of conflicting odd-cycle pairs = the **Ursell pair integral** `−b_2`. The repo defect `p33`
  is the 3-cycle/3-cycle block of `|E(Ω)|` (at Paley `T₇`, `|E|=3153` splits as
  `(3,3):84,(3,5):588,(3,7):336,(5,5):861,(5,7):1008,(7,7):276`; `p55` already dominates `p33`).
  The 3rd Ursell integral is `b_3 = P₂(Ω) − t(Ω)` (cherries minus triangles).

> **CORRECTION the adversarial pass forced (and I had independently caught):** the clean
> `−|E(Ω)|` is the **graphical/Ursell excluded-volume integral** `α_2 − C(α₁,2)`, **not** the
> order-2 coefficient of `log I` (that is `c_2 = −|E| − α₁/2`). The graphical `b_k` do **not**
> exponentiate to `I` — decisive `C₄` counterexample: graphical `b=[4,−4]` gives `z²`-coeff `4`
> but `α_2(C₄)=2`. Only the `c_k` satisfy `exp(Σc_k z^k)=I`.

> **Convergence (formal, not numeric at z=2).** The radius of convergence is
> `R = min|root of I(Ω,z)|` — for Paley `T₇`, `R≈0.0125`, so `z=2` is far outside and the
> numeric series `Σ c_k 2^k` **diverges** (partial sums `~10^174`). "`H` = the cluster series"
> holds only as a **resummed/formal** identity — exactly the project fact that `λ=2` lies
> outside hard-core uniqueness for `Ω` of max-degree `≥4`.

**Not this gas (separate generating function):** the Witt cumulants `W_k` and `TQ`
(`tr A⁷=7(c₇+TQ)`) belong to the spectral/Bowen–Lanford zeta `det(I−uA)^{−1}=∏(1−u^k)^{−W_k}`
of the adjacency matrix `A`, **not** to `log I(Ω,z)`. The only clean bridge is
`α₁ = c₃+c₅+c₇+… =` the sum of odd primitive-necklace counts.

## B. The forbidden values {7,21} in cluster language (PROVED / cited)

- **`H=7` = the unique excluded cluster.** `7 = I(K₃,2) = Φ₃(2) = 2³−1`. The **only**
  graph-realizable cluster profile with `Σα_k2^k=7` is `α=(1,3,0)`, i.e. `Ω=K₃` (three
  pairwise-conflicting odd cycles, zero independent pairs — the maximal-excluded-volume
  3-cluster), and THM-029 proves exactly that `Ω` is **non-realizable**. So `7` is forbidden by
  **single-cluster rigidity**: its lone realizing `Ω` lies outside the realizable-conflict-graph
  cone. (Verified: exhaustively no `T`, `n≤6`, has `Ω=K₃`; 200k random `n=7` gave none.)
- **`H=21` is different — a multiplicative gap.** `21` has **four** graph-realizable profiles
  (`[1,4,3],[1,6,2],[1,8,1],[1,10]`, each built explicitly), so it is **not** single-cluster
  rigid. It is forbidden **multiplicatively**: `H` factors over strong components, `21=3·7`,
  and every odd factorization routes through the forbidden strong-`7`.
- **Non-achievability for all `n` (CITED, re-confirmed).** Already closed in canon (HYP-2271,
  monad-compute-2026-06-06): `strong-min(m)=A005517=3,5,9,15,25,45,75,…` (Busch 2006, EJC 13
  #N3 = Moon's upper bound), strictly increasing with `a(m+3)=5a(m)` for `m≥3`, so `{7,21}` sit
  below the strong-min floor for `m≥7` and are absent at `m≤6`. The cluster picture **cleanly
  re-explains** the `7` half (unique excluded cluster) and **re-describes** the `21` half
  (multiplicative); the surjectivity half of HYP-2271 (every odd `≥23` achievable) stays open
  and is **provably outside** cluster/multiplicative reach (a prime `H` needs a strong
  tournament of that `H`, which closure cannot supply).

## Scope / honesty

PROVED & VERIFIED: A (the cluster identities, with the c_k/|E| correction and the formal-only
convergence) and B (the 7=unique-excluded-cluster / 21=multiplicative characterization). CITED,
not new: the {7,21}-for-all-n non-achievability (Busch/A005517, HYP-2271) — the workflow
independently re-derived it (good cross-validation) but it was already canon. The LRC-parallel
(HYP-2544), the apex-7 dual gate (HYP-2545), and the E₇-Fano cousin (HYP-2546) are separate
claims of decreasing rigor (parallel-by-re-reading → resonance → flagged analogy). Answers
OPEN-Q-100(a). Scripts: `04-computation/{mayer_ocf_cluster,forbidden_cluster_characterization,
forbidden_21_profiles,forbidden_semigroup_closure,forbidden_surjectivity_route}_kps.py`
(+ adversarial re-derivations `adversarial_mayer_check_kps.py`).
