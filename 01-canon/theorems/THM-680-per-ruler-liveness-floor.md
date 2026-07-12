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

---

## ADDENDUM (klein-2026-07-12-S264, HYP-6130) — the defining line is POSITIVE: floor sharpens to (b/q)^12, and the wider band targets THM-720's growing M

**The sharpening.** Statement (iii) bounds the defining line `L*` in absolute value and
**subtracts** it. But the band `B = [d, q−d]` is symmetric (`r ↔ q−r`), so `ĥ` is **real**
(`ĥ(k) = (1/q)Σ_{x∈B} cos(2πkx/q)`, imaginary parts cancel), and on `L*` the two active
coordinates `i, j` carry the SAME index `m` (since `v_i + v_j ≡ 0`, the support-`{i,j}`
solution of `k·v ≡ 0` is `k_i = k_j = m`). Hence each defining-line term is
`(b/q)^11 · ĥ(m)^2 ≥ 0` — **POSITIVE**. It should be **added**, not subtracted:

> **`Σ_{m≠0} (b/q)^11 ĥ(m)^2 = (b/q)^11[(b/q) − (b/q)^2] = (b/q)^12(1 − b/q)`** (Parseval, exact),
> so `(b/q)^13 [main] + (b/q)^12(1 − b/q) [defining line] = (b/q)^12` **exactly**, giving the

> **EXACT IDENTITY:  `LM/q = (b/q)^12 + OffLine_signed`,   `OffLine_signed = Σ_{k∈Λ_q∖(L*∪0)} ∏ ĥ(k_l)`,**
> **and the SHARPER FLOOR:  `LM/q ≥ (b/q)^12 − |OffLine(q)|`.**

At `c = 1/14` (`b/q = 6/7`) this is `(6/7)^12 ≈ 0.1573`, improving the published
`(b/q)^12(2b/q−1) ≈ 0.1124` (and even the drop-`L*` bound `(6/7)^13 ≈ 0.1348`). The published
floor is valid but over-conservative — it subtracts a positive quantity.

**The wider band (the target).** Nothing pins `d = ⌈q/14⌉`. The identity/floor hold for **any**
half-width `d`, so `c = d/q` is a free parameter and `(b/q)^12 ≈ (1−2c)^12` stays positive up to
`c → 1/2`. Therefore

> **`M(S) ≥ c` whenever some pair-sum ruler `q` has `|OffLine(q, c)| < (1−2c)^12`.**

This is the a-priori CERTIFICATE form of **THM-720's SAMPLED growing-M**, and it lives on the
**pointwise pair-sum side, immune to the signed-cancellation wall** (HYP-5830/opus-S225) that
defeats every measure-`μ` attack (the identity is a Parseval sum, not a non-truncatable series).

**What the verification shows** (`lrc14_wideband_parseval_floor_klein_S264.py`, all exact):
- (V1) The identity is confirmed: `main13 + defining-line = (b/q)^12` to `1e−9` on all cases;
  `OffLine_signed = LM/q − (b/q)^12` visibly takes both signs.
- (V3) **The sharpened floor's reach `c_floor` (largest `c` with `(b/q)^12 − |OffLine| > 0` at
  some pair-sum `q`, `|OffLine|` exact) EQUALS or nearly equals the true `M` and GROWS with
  diameter** — a per-family floor matching THM-720: AP `1/14` (wall); DC-bounded `1/12`; kps
  blocker `406/1669 = 0.243` (`= M`, diam 1656); det-spread scale-200 `77/393 = 0.196` (`= M`,
  diam 2433); adversarial-DC `10/49 = 0.204` (`M = 1/4`). `c_floor > 1/14` for every spread
  family, never capping at `1/14`; the sharpening buys real reach over THM-680 (0.204 vs 0.192).
- (V5, **honest limit**) An a-priori bound on `|OffLine|` by the **unsigned** small-relation mass
  (support ≤ 3, `|coeff| ≤ 4`, exact `ĥ`) reaches only `c ≈ 0.03–0.05 < 1/14`: the absolute sum
  over-counts massively (no cancellation), exactly HYP-5830. **So the remaining crux — this file's
  §(iv) off-line classification — must be a SIGNED estimate of `OffLine_signed`, not an absolute
  one.** The floor mechanism has room all the way to the true growing `M`; only the a-priori
  control of the signed off-line sum is missing.

Net: the per-ruler liveness floor sharpens to `(b/q)^12 − |OffLine|`, generalizes to a growing
`c`, and — with `|OffLine|` at its true (signed) value — certifies `M ≥ c` up to the true `M`.
The residual is precisely a **signed** off-line-sum bound for spread families.
Cross-links mac-mini-S65 cont.50 (concurrent, ≤6-distinct-lifts on the same spread DC class).
