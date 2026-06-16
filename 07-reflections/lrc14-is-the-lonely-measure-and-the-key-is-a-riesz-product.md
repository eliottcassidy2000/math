# LRC(14) is the lonely measure, and the key is a Riesz product

**Source:** mac-mini-2026-06-15-S4. Dispatch: spend a long session applying recent
innovations to the LRC(14) proof. Canon: THM-515, T826, HYP-2540..2543, OPEN-Q-104.
Builds on THM-501 (singular series), THM-503 (almost-Sidon class), THM-504 (the
|T|≥3 wall), and the session-stack innovations (free-gas/Mayer, FKN/Walsh, theta/
positive-definite, det/permanent wall).

## What the singular series really is

The LRC(14) singular series, stripped of its circle-method clothing, is one clean
object:
```
L(S) = Σ_{t ∈ Λ} ∏_i h(t_i)        (theta-form, Λ = relation lattice, h = −sinc band)
     = ∫_0^1 ∏_i 1_safe(v_i τ) dτ    (the Lebesgue MEASURE of the lonely set).
```
Because each `1_safe ∈ {0,1}`, `L ≥ 0` is structural (it is a measure) — `L` is just *the fraction of time all 13 runners are simultaneously ≥ 1/14 from the
origin*. C'(14) (hence LRC(14)) is the single statement `inf_S L(S) > 0` over primitive
multiple-of-14 speed sets: **the lonely set never has zero measure unless the speeds are
exactly tight.** Everything else — the (6/7)^13 main term, the (−1)^{|T|} relation
corrections, the sinc weights — is the harmonic expansion of this one measure.

## Why every clean reframe hits the same wall

I brought the session's innovations to bear, and they converge on the same diagnosis:

- **Theta / positive-definite:** gives `L ≥ 0` for free and the lonely-measure identity,
  but the naive lattice sum diverges (the box truncation reproduces the |T|≥3 blow-up).
  Positivity is structural; the *uniform lower bound* is not.
- **Free-gas / Mayer cluster expansion:** `L = (6/7)^13 + Σ_T (−7/6)^{|T|}R_T` is a
  *signed polymer gas* (the runners are particles, the additive relations are polymers).
  But standard Mayer/Kotecký–Preiss needs absolutely summable activities, and `A_3 = ∞`
  (THM-504). The free gas works for the *bosonic* OCF (`H = I(Ω,2)`, T006) and breaks on
  the *fermionic* LRC sign.
- **Pfaffian / determinant (the BSD/Hodge alternating→square spine):** the `(−1)^{|T|}`
  *looks* fermionic, but the within-level signs are positively correlated (half-period
  positivity, THM-504-A), so they do not assemble into an antisymmetric kernel; and the
  project's Pfaffian objects are A-affine/spectral — the *wrong side of the Valiant
  det/permanent wall*. The LRC content lives on the speed-relation lattice, not the
  tournament adjacency matrix. Closed door.
- **FKN / second moment:** Paley–Zygmund needs a nonnegative integrand; the raw signed
  sum is not one. It is a *finishing* tool, not an opener.

The unanimous verdict: the saving that keeps `L > 0` is **cross-level**, in the
`(−1)^{|T|}` alternation over the support filtration, with each level conditionally
convergent. No absolute method, and no reframe that merely *recovers* positivity, can
reach it. This is not a failure of imagination — it is the precise reason LRC(14) is
hard, the same "improve the union bound by a factor of two" that Tao identifies.

## The one innovation that is the actual key: a Riesz product

There is a published, working precedent on the *cousin* problem — the lonely-runner
*gap* `ML(V) ≥ 1/(2n) + 1/n^{5/3}` (arXiv:2511.16636, 2025), which beat Tao using
**exactly** a signed-sinc-lattice-sum tested against a **positive** measure. The trick is
a **Riesz product**:
```
R(τ) = ∏_{m∈D} (1 + a_m cos 2π m τ) ≥ 0   over a dissociated set D,
```
a probability density whose Fourier coefficients are **signed** `(a/2)^{level}`. Pair the
covering `M = Σ_v 1_danger` against `R`:
```
∫ M · R = Σ_v Σ_k s(k) R̂(v k).
```
The main term (`k=0`) is `Σ_v s(0) = 13/7 ≈ 1.857 > 1`. The *signed* Riesz coefficients,
sitting on the relation frequencies, must subtract `≥ 0.857` to force `∫MR < 1 = ∫R` —
and since `M ≥ 1` on any cover, that **certifies the lonely set has positive measure**.
The dissociated set `D` is tuned to the speed-relation lattice — precisely the
*additive-energy* structure that governs `L` (the cores have dense relation lattices;
that density is what the Riesz product must dissociate). Bonami's hypercontractive
inequality controls the higher levels.

This is the right key because it **keeps a positive test measure while extracting
cancellation from a signed character sum** — exactly the cross-level cancellation
THM-504 says is essential and no absolute method captures. My feasibility probe (signed
cosines pull `∫MR` from `1.857` to `1.41` on the `j=6` extremizer) shows the pairing
machinery works; reaching `< 1` needs the optimized dissociated construction — the
concrete open step (HYP-2540, OPEN-Q-104).

## The shape of the remaining problem

Three facts now organize the endgame:

1. **It suffices to handle the extremizers.** `L` is governed by the relation-lattice
   additive energy; the infimum sits at the maximal-energy interior-drop cores
   `{1..13}\{j}∪{14m}`, `j=6` worst (`L ≈ 0.0053`). Dominance is orthogonal to `L`
   (HYP-2526): the hard core is bounded and balanced. So `inf L>0` reduces to a Riesz
   product for a *bounded family* of dense-relation cores.
2. **Tao's third moment locates the difficulty.** `∫M³` (triple danger-overlaps) is large
   exactly at the cores — the arithmetic-progression structure Tao's method flags as the
   obstruction. The Riesz product must dissociate these APs (Bedert's `E_k` level-sets are precisely sum-sets of APs, bounded by Bonami).
3. **Positivity is no longer the issue — uniformity is.** `L ≥ 0` is structural; the
   content is the constant-factor gap between `13/7` (union bound) and `1` (the truth),
   the same factor-of-two everywhere in lonely-runner theory.

The session did not prove `inf L>0` — it is genuinely as hard as the cousin gap problem.
But it relocated the problem onto a **tool with a 2025 success on the neighbor**, set up
the exact pairing, reduced the target to a bounded family of cores, and shut three
plausible-looking doors (polymer-gas, Pfaffian, raw second-moment) so the next attempt
spends its effort on the Riesz product, where the open problem actually lives.
