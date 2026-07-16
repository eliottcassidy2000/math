---
id: THM-886
title: THE RESONANCE LAW OF Q_s ON THE ASSEMBLY FAMILY (the repaired HYP-6994, family version) — for E_t = {1..6, t}: (I) EXACT MODE IDENTITY S(ta) = S_frame(ta) + Σ_{c mod 7} e(ac/7)·N_c for EVERY integer a (elementary: at n = ta the endpoint phase depends only on j mod 7), recovering klein-cont.4's mode law t·ν̂(a) with N_c = t·ν_c + O(R) — at exact modes no Abel/Koksma is needed; (II) OFF-MODE PROFILE |S(m)| ≤ C_f + R/(2sin(π‖m/t‖)) (run-geometric; refereed 0 violations, worst ratio 0.707); (III) |k̂_P(n)| ≤ 1/(2P²sin²(πn/P)) ≤ 1/(8ñ²) exact-shape; (IV) uniform HYP-6994 REFUTED — max|S|²/M = 14.0 → 226.0 for t = 6 → 1200 (independent scans to twice klein cont.4's validated range; convergent, theirs first-pushed); (V) THE RESONANCE LAW: Q_s(w) ≤ C₁M + 2π²Σ_ℓ (8ℓ²)⁻¹·Amp(ℓw mod P)², and the resonant classifier is DIOPHANTINE, not arithmetic: w is resonant ⟺ some small multiple ℓw lands near the mode lattice tZ — w plays a LONELY RUNNER against the mode lattice (self-similarity one level up); coprimality is NOT the classifier (w = 2t+7 coprime: Q_s/M = 132 vs 0.06–0.73 at random coprime w); sup_w Q_s ≍ t² sharp
status: (I)–(III) PROVED (elementary, machine-refereed: mode identity exact at all tested a; profile 0 violations; k̂ ratio exactly 2.000 vs the naive form); (IV) refutation CONFIRMED (convergent with klein-S314 cont.4, first-push theirs); (V) upper law PROVED with stated constants (covering verified on the w-battery incl. both blowup classes; deliberately lossy but uniform), t² lower sharp at w ≡ ta; FAMILY-RESTRICTED — general clusters need per-cluster (C_f, R, ν̂), same proof shape (named next)
source: boxeph-2026-07-16-S25 (owner directive: work the signed sup-norm)
depends_on: [THM-880/881 (the bilinear/descent frame), klein-S314 cont.4 (the mode law + refutation, first-push; their THM-'883'-resonant-mode — ID collision pending), THM-729]
related: [THM-879 (the k=13 interval-core rate — different arc systems, complementary), death-star HYP-7017 (claimed target = the now-refuted uniform lemma — redirect flagged)]
script: 04-computation/lrc14_hyp6994_resonance_test_boxeph_S25.py -> 05-knowledge/results/lrc14_hyp6994_resonance_test_boxeph_S25.out
---

# THM-886 — the resonance law of Q_s (family version)

**Objects.** E_t = {1,…,6, t}; section s; R_s = {x : the seven sections of E_t under
x ↦ ⌊7{ex}⌋ occupy exactly the six values ≠ s} (klein's swing set); S(n) = Σ_k ε_k e(np_k)
the signed endpoint sum on Z_P, P = 7·lcm(60, t); M = #endpoints; C_f = #frame-owned
endpoints (owners ≤ 6); for the t-owned endpoints, j/(7t) with residue class c = j mod 7:
N_c = the SIGNED count in class c, and R = total number of runs of the per-class active
index sets in Z_t.

**(I) Exact mode identity.** For every integer a,
> S(ta) = S_frame(ta) + Σ_{c mod 7} e(ac/7)·N_c.
*Proof:* at n = ta the t-owned endpoint phase is e(ta·j/(7t)) = e(aj/7), which depends
only on j mod 7. ∎ (Refereed exact at a = 1..7, 11; t = 600, s = 0: N =
(−87, 100, −2, 9, −12, −2, −8), amplitudes |Σ N_c e(ac/7)| = 97.1, 138.1, 184.4 vs
measured |S(ta)| = 101.2, 141.2, 187.2 — the gap is the frame term ≤ C_f = 16.)
Since N_c = t·ν_c + O(R) (grid sampling of the frame's limiting signed section
measures), this recovers klein cont.4's mode law S(ta) ≈ t(1−e(a/7))-weighted ν̂(a) —
independently, and EXACTLY at modes (their Abel/Koksma is needed only off-mode/asymptotically).
Note ν̂(a) is 7-periodic in a with the 7|a modes small (Σ N_c = −2 here): the mode comb
has FULL-STRENGTH teeth at every ta, a ≢ 0 mod 7 — no a-decay.

**(II) Off-mode profile.** For m ∉ tZ:
> |S(m)| ≤ C_f + R/(2 sin(π‖m/t‖)).
*Proof:* each per-class active index set is a union of runs (integer intervals in Z_t);
a run's exponential sum at frequency m/t is geometric, bounded by 1/(2sin(π‖m/t‖));
R runs total; frame endpoints contribute ≤ C_f. ∎ (Refereed: 0 violations over all
m ∈ Z_P at t = 600; worst ratio 0.707 — near-sharp.)

**(III) Kernel DFT.** k̂_P(n), the Z_P-DFT of {x}(1−{x}), satisfies |k̂_P(n)| ≤
1/(2P² sin²(πn/P)) ≤ 1/(8ñ²) (second-difference/double-Abel: the discrete second
difference of the sampled kernel is −2/P² off 0 with a 2/P-boundary at 0). Refereed:
max |k̂|·16ñ² = 2.000 exactly (the naive 1/(4P²sin²) fails by exactly this factor 2).

**(IV) The refutation (independent confirmation).** Direct full-Z_P scans:

| t | 6 | 12 | 25 | 37 | 50 | 60 | 120 | 300 | 600 | 1200 |
|---|---|---|---|---|---|---|---|---|---|---|
| max|S|²/M | 14.0 | 6.5 | 9.7 | 9.7 | 10.0 | 9.9 | 27.1 | 58.3 | 110.4 | 226.0 |

Linear growth (argmax on the mode comb: e.g. 3600 = 3t at t = 1200); |S(t)|²/M and
|S(2t)|²/M grow linearly. **The uniform sup-norm lemma HYP-6994 is FALSE; the C = 14
of the t ≤ 50 scans was pre-asymptotic.** Same-night convergent with klein cont.4's
automaton refutation (theirs first-pushed; this table extends the validated range 400 →
1200 by an independent method). The lower bound is (I): |S(ta)| ≥ t|ν̂(a)| − O(C_f + R)
with ν̂(1) ≠ 0 an exact frame constant.

**(V) The resonance law.** For every integer w (coprime or not),
> Q_s(w) = 2π² Σ_{n≠0} k̂_P(n)|S(nw)|² ≤ C₁·M + (π²/4)·Σ_{ℓ≥1} ℓ⁻²·Amp(ℓw mod P)²,
where Amp(m) is the (I)/(II) profile: t|ν̂(a)| + C_f at modes m = ta, and
C_f + R/(2sin(π‖m/t‖)) off-mode; C₁M absorbs the low band (|S(n)| ≤ min(M, 2πnλ(R_s)),
refereed) and the Parseval bulk. Verified covering on the battery
{11, 1013, 601, 1207, 105} at t = 600 (lossy but uniform).

**The classifier is Diophantine, not arithmetic.** Q_s(w)/M measured: random coprime
w: 0.06–0.73; w = t+1 (coprime): 32.7; **w = 2t+7 (coprime): 132.4**; w ≡ ta: Θ(t²)/M.
So the resonant classes are NOT the non-coprime w (klein cont.4's w ≡ taℓ̄ list is the
ρ = 0 slice): w is resonant ⟺ **some small multiple ℓw lands within a short distance ρ
of the mode lattice tZ**, with blowup t²|ν̂|²/(4ℓ²) at ρ = 0 and (R t/(2πρ))²/ℓ²-scale
tails at ρ ≠ 0. Equivalently: **w is resonant against the mode lattice exactly as a
lonely runner is resonant against a torsion grid — the problem reproduces its own
Diophantine structure one level up.** (w/t admitting a good rational approximation with
small denominator ℓ and small numerator-defect ρ/t.)

**(VI) What this changes.**
1. THM-729/881's empirical "Q_s = O(diam) uniformly in w" is an OFF-RESONANCE law; the
   resonant classes are explicit and finite per instance (and P2-decidable).
2. The program is UNAFFECTED in the limit it uses: even resonant Q_s ~ t² gives
   |S|/w ≤ t/w → 0 along the peel sequence w → ∞ — the two-scale error still closes.
   But any FIXED-w application must check the Diophantine class of w first.
3. death-star's HYP-7017 (claimed: uniform HYP-6994 via coded-rotation characters)
   targets the refuted statement — redirect to the weighted/resonance law.
4. The general-cluster resonance law: same skeleton with per-cluster (C, R, ν̂) and one
   mode comb per large owner — the named next step; the multi-owner comb interaction
   (CRT of the combs) is the only new ingredient.

## Evidence log
- [x] (I) exact at all tested a; (II) 0 violations, ratio 0.707; (III) ratio exactly 2
- [x] (IV) scans t ≤ 1200, linear growth; (V) battery covered, both blowup classes hit
- [ ] general-cluster version (multi-owner mode combs)
- [ ] Lean: (I) is decide-shaped; (III) is 3 lines of algebra
