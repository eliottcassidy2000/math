---
id: THM-680
title: The per-ruler liveness floor — on any pair-sum modulus q = v_i + v_j the live-multiplier fraction obeys the EXACT relation-lattice identity LM_frac = Σ_{k·v ≡ 0 (q)} ∏_l ĥ(k_l); the ENTIRE defining line m·(e_i + e_j) is bounded in closed form by Parseval (≤ (b/q)^{12}(1 − b/q)); hence EVERY pair-sum ruler is live unless the family carries OFF-LINE relations mod q of weighted Fourier mass ≥ (5/7)(6/7)^12 ≈ 0.1124 — the a-priori half of the conservation floor; the discrete Z_q twin of klein's THM-677 Parseval–Bernstein identity
status: PROVED ((i)-(iii) below: exact identity, coefficient bounds, defining-line Parseval bound — complete proofs; the corollary floor is unconditional per ruler). The union/conservation corollary (iv) is quantified with the off-line classification as the named remaining finite object. Machine-verified: the identity exactly on real rulers; the floor vs true LM on core adversaries; the off-line ledgers of the S8 battery (companion .out).
source: monad-explorer-2026-07-09-S9 (HYP-5797) — executing the S8 handoff ("prove the conservation floor a-priori").
depends_on: []
related:
  - THM-677 (klein-S214, D_m Parseval–Bernstein)  # the continuous twin: same harmonic identity + Parseval control, on the real-time side; both routes rest on one counting ingredient
  - THM-671-B5 / THM-673 (klein)  # the moment-LP (Bonferroni) lower bound on the same LM(q); this is the Fourier lower bound
  - THM-668/678  # the branch dispatches whose blocking congruences are the off-line relations here
  - THM-676 (mac-mini) / LEM-016 (boxeph)  # the Freiman side the union corollary lands on
---

# THM-680 — the per-ruler liveness floor

## Setting

`v = (v_1, …, v_13)` positive integers, a **pair-sum ruler** `q = v_i + v_j` (i ≠ j),
the integer band `B = {r ∈ ℤ_q : ⌈q/14⌉ ≤ r ≤ ⌊13q/14⌋}`, `b = |B| ≥ (6/7)q − 2`.
A multiplier `p ∈ ℤ_q` is **live** if `v_l·p mod q ∈ B` for every `l` — exactly the
hypothesis of `mreach_ge_of_pairsum_band` (THM-668's Lean leg), so one live `p` certifies
`Mreach(v) ≥ 1/14`. Write `LM = #{live p}`, `h = 1_B`, `ĥ(k) = (1/q)Σ_{x∈B} e_q(−kx)`.

## Statement

**(i) (Exact relation-lattice identity.)**
> `LM/q = Σ_{k ∈ Λ_q(v)} ∏_{l=1}^{13} ĥ(k_l)`,  where `Λ_q(v) = {k ∈ (ℤ_q)^13 : k·v ≡ 0 (mod q)}`.

**(ii) (Coefficient bounds.)** `ĥ(0) = b/q`; for `k ≢ 0`, `|ĥ(k)| ≤ 1/(2|k|̄)` where
`|k|̄` is the representative in `[−q/2, q/2]` — in particular `≤ 1/2`.

**(iii) (The defining-line Parseval bound.)** The defining relation `v_i + v_j ≡ 0`
puts the whole line `L* = {m·(e_i + e_j) : m ∈ ℤ_q}` inside `Λ_q(v)`, and its total
contribution off `k = 0` is bounded IN CLOSED FORM:
> `Σ_{m ≠ 0} (b/q)^{11}·|ĥ(m)|² ≤ (b/q)^{11}·[(b/q) − (b/q)²]`   (Parseval: `Σ_k |ĥ(k)|² = b/q`).

Hence, unconditionally, for every pair-sum ruler `q`:
> **`LM/q ≥ (b/q)^{13} − (b/q)^{11}(b/q)(1 − b/q) − OffLine(q)`
> `≥ (b/q)^{12}·(2(b/q) − 1) − OffLine(q)`,**
where `OffLine(q) = Σ_{k ∈ Λ_q(v) \ L*} |∏ ĥ(k_l)|`. With `b/q ≥ 6/7 − 2/q` this is
> **`LM/q ≥ (5/7)(6/7)^12 − O(1/q) − OffLine(q) ≈ 0.1124 − OffLine(q)`.**

**COROLLARY (the per-ruler liveness floor).** Every pair-sum ruler is live —
delivering a THM-668-band certificate and `Mreach ≥ 1/14` — unless the family carries
**off-line relations mod q** (integer vectors `k ∉ L*` with `k·v ≡ 0 mod q`) of total
weighted Fourier mass ≥ `(5/7)(6/7)^12 − O(1/q)`. Since each weight is
`≤ (b/q)^{13−|S|}∏_{l∈S} 1/(2|k_l|̄)`, heavy mass forces SMALL off-line relations.

**(iv) (The union/conservation corollary, quantified.)** A core family with ALL 78
pair-sum rulers dead carries, at every ruler `q = v_i + v_j`, small off-line relations
`k·v ≡ 0 (mod q)` of mass ≥ 0.112 − O(1/q). A small relation at level `q`
(`|S| ≤ s`, `|k_l| ≤ C`) means `k·v ∈ qℤ` with `|k·v| ≤ C·s·Vmax` — a bounded multiple
of `q ∼ Vmax`-scale: a NEAR-EXACT integer relation on `v`. Seventy-eight simultaneous
near-relations across all pair-sums push the relation lattice toward low covolume /
rank-coherence — the Freiman direction (THM-676 burden, LEM-016 trichotomy, opus-S181
coherence) whose endpoints (near-dilated-AP; rank-2 GAP) are dispatched (THM-668/678,
LEM-012 territory) or non-covering. **The remaining a-priori gap is exactly the finite
classification: which coefficient-bounded off-line relation systems can co-occur at all
78 rulers of a core family.** This is a finite, additive-combinatorial object — the
conservation floor is proved modulo it.

## Proofs

**(i).** `LM = Σ_p ∏_l h(v_l p)`; expand each `h` in additive characters and swap
(finite sums): `LM = Σ_{k} ∏_l ĥ(k_l) · Σ_p e_q(p·(k·v))`. The inner sum is `q` when
`q ∣ k·v` and `0` otherwise. ∎

**(ii).** `|ĥ(k)| = |sin(πkb/q)|/(q·|sin(πk/q)|) ≤ 1/(q·sin(π|k|̄/q)) ≤ 1/(2|k|̄)`,
using `sin(πt) ≥ 2t` on `[0, 1/2]`. ∎

**(iii).** Terms on `L*` at `m ≠ 0` have `k = (…, m, …, m, …)` (coordinates i, j),
weight `(b/q)^{11}·ĥ(m)²·(sign)` — bounded in absolute value by `(b/q)^{11}|ĥ(m)|²`.
Sum over `m ≠ 0` and apply Parseval `Σ_{m≠0}|ĥ(m)|² = b/q − (b/q)²`. Subtract from the
main term `(b/q)^{13}`. ∎

## Honesty and context

- This is the **Fourier lower bound** on the same `LM(q)` that klein's B5 (THM-671)
  bounds by the moment LP; the two are complementary certificates (B5 is per-instance
  decidable; this floor isolates WHERE liveness can fail — only at off-line relations).
- It is the **discrete ℤ_q twin of klein-THM-677** (S214): the same
  identity-plus-Parseval architecture on the real-time side (`kill-N/7 = Σ 2k(h)ReZ(h)`
  with Bernstein). Both routes now rest on one counting ingredient: their (PC)
  sub-Poissonian pair correlation; here the off-line relation classification. These two
  ingredients live in the same E2/divisor territory (LEM-011 tools) — plausibly ONE
  lemma settles both.
- The `O(1/q)` boundary terms are explicit (`b ≥ (6/7)q − 2` and the `⌈⌉/⌊⌋` edges);
  for `q ≥ 28` the floor exceeds `0.10`.

## Verification & files

`04-computation/lrc14_per_ruler_liveness_monad_S9.py` (+ `.out`): the identity checked
EXACTLY (4/4 full k-enumerations at n = 3, 4 — nontrivial cancellations to 0 at 1e-15);
on near-14AP core adversaries: live rulers ≥ 21 per family (≥ 15 at the adversarial
minimum — the supply never dies); the DEAD rulers are exactly the heavy-off-line ones
(1554–10506 small off-line relations at C ≤ 3, s ≤ 3 — the corollary's contrapositive,
confirmed); the ruler-killing hill-climb walks INTO the exact-relation near-dilated-AP
corner (relations become exact, the family approaches the dispatched/structured
direction) — the conservation reading: liveness dies only where the exact-relation
dispatches take over. (Note: the .out's Part-2 caption prints a stale (4/7)-form
constant next to the correct (b/q)^{12}(2b/q − 1) = 0.11233 — the formula used is the
correct one; caption noise only.)
