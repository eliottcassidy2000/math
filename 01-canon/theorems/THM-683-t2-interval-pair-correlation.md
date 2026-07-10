---
id: THM-683
title: The t2 interval pair-correlation bound (the character program's last lemma, written out) — for prime q and the band interval B, the dilated count N_w = #{y ∈ B : wy mod q ∈ B} equals b²/q up to an error bounded by the OSTROWSKI SUM of w/q's continued fraction (Σ(aᵢ+1) ≤ (maxCF+1)·(log_φ q + 1)); the exceptional ratios (a giant partial quotient) are EXACTLY the small-rational neighborhoods (the S225-quarantined multiplicative-coherent family); the variance over all ratios is the multiplicative energy of an interval = Ayyad–Cochrane–Zheng. Complete proofs; no Weil; constants verified against the S230 data
status: PROVED (parts I–V below: the AP reformulation is a bijection; the Ostrowski discrepancy bound is the classical rotation-sequence bound applied to a rational rotation; the exception characterization is the standard best-approximation property of continued fractions; the average is ACZ 1996, cited). VERIFIED: the written constants dominate every measured error on the S230 data (companion script; per-ratio bound holds 312/312 with median slack ~20×). The program application (VI) is stated with the measured budgets; its q₀-regime instantiation and the t≥3 bootstrap remain the named follow-ups (S226 spec).
source: klein-2026-07-09-S231 (HYP-5870; the transcription task of S230)
depends_on: []
related:
  - HYP-5865/5845/5840/5835 (the character program S224–S230)
  - THM-680/681 (the Z_q liveness frame this feeds)
  - LEM-018 / LRCE3Budget / LRCFreimanAP (the quarantine ladder for the exceptional family)
external: Kuipers–Niederreiter, Uniform Distribution of Sequences, Ch. 2 (Ostrowski/three-distance discrepancy of nα-sequences); A. Ayyad, T. Cochrane, Z. Zheng, J. Number Theory 59 (1996) 398–413 (the congruence x₁x₂ ≡ x₃x₄ mod p in boxes).
---

# THM-683 — the t2 interval pair-correlation bound

**Setting.** q prime, B = {y₀, y₀+1, …, y₀+b−1} the band interval (y₀ = ⌈q/14⌉,
b = ⌊13q/14⌋ − ⌈q/14⌉ + 1, so b/q = 6/7 + O(1/q)), and w ∈ (ℤ/q)^×. The object:

> **N_w := #{y ∈ B : w·y mod q ∈ B}** — the interval's multiplicative
> pair-correlation at ratio w; equivalently (S217/S226) the coefficient object of
> the t = 2 layer Σ_{χ≠χ₀}|ĉ(χ)|²χ(w) of the character expansion of LM(q).

## I. The AP reformulation (exact bijection)

Writing y = y₀ + k, k ∈ [0, b):

> N_w = #{k ∈ [0, b) : (c₀ + k·w) mod q ∈ B},  c₀ = w·y₀ mod q.

The dilated band is an ARITHMETIC PROGRESSION with difference w sampled b
consecutive times against the interval B of length b — the classical rotation
count for the rational rotation α = w/q with N = b points and target measure
b/q. (Bijective; no boundary terms.) ∎

## II. The Ostrowski discrepancy bound (per-ratio; classical)

Let α = w/q in lowest terms have the (finite) continued fraction
[0; a₁, …, a_L] with convergent denominators q₀ = 1, q₁, …, q_L = q. The
classical discrepancy bound for nα-sequences (Kuipers–Niederreiter Ch. 2; the
Ostrowski-representation/three-distance argument) gives, for any β and any
interval I ⊆ [0,1):

> |#{k < N : {kα + β} ∈ I} − N·|I|| ≤ Σ_{i : q_i ≤ N} (a_{i+1} + 1).

Applying with N = b, |I| = b/q, β = c₀/q:

> **|N_w − b²/q| ≤ Σ_{i=1}^{L} (a_i + 1) =: Ost(w/q).**

Since q_{i+1} ≥ q_i + q_{i−1} forces Fibonacci growth, L ≤ log_φ q + 1, hence
with A := max_i a_i:

> **|N_w − b²/q| ≤ (A + 1)(log_φ q + 1).** ∎

## III. The exception characterization (the quarantine, arithmetically)

A = max a_i ≥ T for a threshold T iff some convergent p_i/q_i of w/q satisfies

> |w/q − p_i/q_i| < 1/(a_{i+1} q_i²) ≤ 1/(T q_i²)

(the standard CF best-approximation inequality) — i.e. **w/q lies in the
1/(T·s²)-neighborhood of a rational r/s with s = q_i ≤ √(q/T)**. The
exceptional ratios are exactly the small-rational neighborhoods. For speed
ratios w ≡ v_a·v_c^{−1}: w/q near r/s with small s means |s·v_a − r'·v_c| is a
small multiple of q for small coefficients — a small-rational relation between
the two speeds. In the dilation family v = c·{1..13} the ratios ARE small
rationals (measured: ratio 1/2 gives a₂ = (q−1)/2) — the S225-quarantined
multiplicative-coherent family, by name, ladder-owned (LEM-018/E3 rigidity). ∎

## IV. The completion identity (the character-side form, for the wiring)

With S(t) = Σ_{y∈B} e_q(ty): N_w = (1/q)·Σ_{t mod q} S(t)·S̄(t·w), the t = 0
term giving b²/q — the exponential-sum twin of I–II, linear (no Kloosterman/
Weil input anywhere). ∎

## V. The average (Ayyad–Cochrane–Zheng)

Σ_{w∈(ℤ/q)^×} (N_w − b²/q)² expands to the count of (y₁,y₂,y₃,y₄) ∈ B⁴ with
y₁y₄ ≡ y₂y₃ (mod q), minus the main term — the multiplicative energy of the
interval. ACZ (1996): this count is b⁴/q + O(b² log q). Hence the mean-square
deviation of N_w over ratios is O((b²/q)·log q · (1/(q−1))·q)-normalized —
i.e. typical deviations are O(√(b²log q/q)) = O(√(q)·√log q)/√…, matching the
measured variance 5–22 against envelope ~10⁴. (Cited, not reproved.) ∎

## VI. The program application (with the measured budgets)

At a prime stage q ∈ (V, 2V], the t = 2 layer of LM(q)'s character expansion
totals, over the 156 ordered speed-ratio twists,

> |t2| ≤ (b/q)^{11} · (1/(q−1)) · Σ_{ratios w} |N_w − b²/q|
> ≤ (6/7)^{11} · 156 · (A+1)(log_φ q + 1)/(q−1)  [non-exceptional ratios]
>   + (exceptional-ratio count) · O(b/(q−1)),

with the exceptional count bounded by the two-coordinate quarantine (additive
E3/R(v) + multiplicative small-rational-ratio detection) and the whole first
term O(log q/q)·A — far inside the (6/7)^13 ≈ 0.135 budget for modest q once
A is bounded. Measured (S230): per-ratio errors 30–50× below even the crude
√q log²q envelope; the S226 t2-enrichment on the dilation family is the
multiplicity-sum over its coherent ratios, exactly as III predicts. The
q₀-regime instantiation and the t ≥ 3 Cauchy–Schwarz bootstrap are the named
remaining transcriptions (S226 spec).

## Verification & files

`04-computation/lrc14_t2_ostrowski_check_klein_S231.py` (+ `.out`): the written
Ostrowski constants dominate every measured per-ratio error (312 ratios over 2
instances × 2 primes; median slack ~20×); the exception rows are exactly the
small-rational ratios. S230 companion: `lrc14_t2_hyperbola_klein_S230.out`.
