---
id: THM-877
title: THE RAMANUJAN TRUNCATION SECOND MOMENT — the truncation error of the exactly-cancelling series Σ_l c_l(h)/l = 0 has closed form T_L(h) = Σ_{d|h} M(⌊L/d⌋) (M = Möbius partial sum), and its limit mean square over h is EXACTLY Σ_{g≤L} φ(g)/g² = (6/π²)·log L + O(1) — LOG-GROWING, NOT BOUNDED: the "bounded mean square = textbook territory" hope of the sharp-rate program is corrected to a (6/π²)·log L law, proved by Gauss's totient identity + the unit identity Σ_{a≤y} M(⌊y/a⌋)/a = 1 applied twice
status: PROVED (three short steps below, fully classical) + machine-exact (closed form vs direct sums; S(L) = totient sum for 8 values of L to 800; finite-H mean squares match; kps cont.26's concentration numbers reproduced to 4 decimals)
source: mac-mini-2026-07-15-S112 (owner: "work the ramanujan truncation mean square classical estimate"); completes the residual named in kind-pasteur cont.26 ("a divisor-sum second-moment problem, textbook territory")
depends_on:
  - THM-873 (kps: Ramanujan-Fourier expansion of interval-core good sets; cont.26's braid — this file quantifies its truncation)
related:
  - THM-874 (the Möbius-log² grammar — same coprime filtration, profile face)
  - klein-S280 sharp rate Q_s = O(r) — SHARPENED TARGET: the raw second moment carries a (6/π²)log L; linearity needs the v-grid restriction to absorb it (or the rate is r·log)
script: 04-computation/ramanujan_truncation_wallwords_macmini_S112.py -> 05-knowledge/results/ramanujan_truncation_wallwords_macmini_S112.out
---

# THM-877 — the Ramanujan truncation second moment

Let c_l(h) be the Ramanujan sum, M(y) = Σ_{m≤y} μ(m)/m, and
T_L(h) = Σ_{l≤L} c_l(h)/l — the depth-L truncation of Ramanujan's identity
Σ_{l≥1} c_l(h)/l = 0 (h ≥ 1), i.e. the exact error of the interval-core spectrum
truncation at l = L (kps cont.26 (B)).

**(1) Closed form.** c_l(h) = Σ_{d|(l,h)} d·μ(l/d); swapping sums with l = dm:

> T_L(h) = Σ_{d|h, d≤L} M(⌊L/d⌋).

So the truncation error is a divisor transform of Möbius partial sums; since M(1) = 1
and M decays, T_L(h) ≈ #{divisors of h in (L/2, L]} + lower dyadic layers — the exact
mechanism of cont.26's "concentrates on divisor-rich h" (reproduced: T_12(12) = 1.9662,
T_12(30) = 2.2996, T_12(1) = −0.0004).

**(2) The second moment, exactly.** As H → ∞, (1/H)Σ_{h≤H} T_L(h)² → S(L) :=
Σ_{d,e≤L} M(⌊L/d⌋)M(⌊L/e⌋)/lcm(d,e) (the O(1)-per-pair boundary is O(L²/H)). Then:

> **S(L) = Σ_{g≤L} φ(g)/g².**

*Proof.* 1/lcm = (d,e)/(de) and Gauss's (d,e) = Σ_{g|d, g|e} φ(g) give
S = Σ_g (φ(g)/g²)·[Σ_{a≤L/g} M(⌊(L/g)/a⌋)/a]². The bracket is 1 for every integer
y = ⌊L/g⌋ ≥ 1 by the unit identity Σ_{a≤y} M(⌊y/a⌋)/a = Σ_{n≤y}(1/n)Σ_{m|n}μ(m) = 1
(group the pairs am = n; only n = 1 survives; floors compose: ⌊L/(ga)⌋ = ⌊⌊L/g⌋/a⌋). ∎
(Machine-exact for L = 6, 12, 25, 50, 100, 200, 400, 800.)

**(3) Asymptotics and the corrected target.** Classically Σ_{g≤L} φ(g)/g² =
(6/π²)(log L + γ − Σ_p log p/(p²−1)) + o(1): the mean square is **Θ(log L)** — numerics:
increments ≈ 0.42 per doubling = (6/π²)·log 2 ✓ (S = 2.22, 3.08, 3.50, 3.92, 4.34, 4.76
at L = 12, 50, 100, 200, 400, 800). Hence: the truncation error's variance on the full
integer grid is NOT bounded; it grows logarithmically with ζ(2)⁻¹ slope. Consequence for
the sharp-rate program (klein-S280 Q_s = O(r), kps cont.26): the raw second moment
contributes a factor (6/π²)log L, so either the v-grid restriction (h confined to the
resonance grid, kps cont.26 (A)) absorbs the log — the remaining thing to prove — or the
true rate is r·log r. kps's prime-window empirics (disc·v² ∈ [1.15, 3.2] over v = 17..199,
a 2.8× spread across a log-ratio-1.87 window) are CONSISTENT with the log law, not with
boundedness. The open analytic core, in final classical form: **does restricting h to the
v-grid remove the (6/π²)log L?**

## Appendix — the wall-word answer (owner's question 3)

For the k-truncated three-distance gap word W(a, q, k) of {i·a/q : 0 ≤ i ≤ k}
(exhaustive q ≤ 60, k = 5, 8, 12): **negation q−a acts by word reversal, always**
(one-line proof: t ↦ −t reflects the point set). **Inversion a ↦ a⁻¹ does NOT act
dihedrally** (first failures q = 9/11/15 at k = 5/8/12): the classical CF-reversal
duality is a full-orbit statement, and the k-truncation breaks it — inversion maps the
index-window to the value-window (Ostrowski transpose), a genuinely different word.
Generic units are lawless (dihedral coincidence rate ≈ 0.001). So: numerator
multiplication permutes wall-words as labels, but the MECHANICAL symmetry group of the
truncated wall-word is exactly {±1}; the inverse-pair witness structure (THM-819) is a
threshold/full-window phenomenon, not a truncated-word symmetry. Refined transpose-duality
statement: backlogged.
