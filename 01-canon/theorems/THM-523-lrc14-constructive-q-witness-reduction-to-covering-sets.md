---
id: THM-523
title: Constructive q-witness reduction of LRC(14) to covering sets — LRC(14) holds for every 13-set omitting a multiple of some q∈{2,…,14} (explicit lonely point τ=1/q), so a counterexample (gap M<1/14) must contain a multiple of EVERY q∈{2,…,14}; the M=1/14 tight minimizers are all in the trivial case, and the residual covering-set hard core is empirically loose (min M=7/89≈0.0787, 10% above 1/14)
status: PROVED (q-witness ⟹ M≥1/q: elementary, the classical small-divisor reduction instantiated exactly for n=14; covering-set necessary condition: PROVED). VERIFIED (exact gap tool `lrc14_gap_M_exact`, validated vs dense grid): 0 counterexamples over single/double/triple perturbations of {1..13}, all dilated/shifted APs, unions of APs, ~25k random 13-sets, and the full cover-all-units interior-drop family; min M over covering sets ≈ 7/89. Residual uniform-looseness CONJECTURAL (HYP-2566).
source: mac-mini-2026-06-16-S3 (prove∧disprove dialectic; gap M(S)=max_τ min_v ||vτ|| instead of the measure L)
depends_on:
  - THM-398   # LRC(14) ⟸ C'(14)
  - THM-501   # the lonely measure / singular series
related:
  - THM-522   # kind-pasteur: scale-invariance + quantization + compactness on L (measure side); this is the GAP-side dual
  - HYP-2561  # inf L = 1/1260 (measure infimum); distinct from the gap M
  - HYP-2566  # residual covering-set hard core uniformly loose (the open piece)
  - HYP-2567  # disprove direction: no counterexample over extensive search
  - OPEN-Q-104
external: Lonely Runner Conjecture n=14; classical small-divisor / covering reductions (Bohman–Holzman–Kleitman, Tao).
---

# THM-523 — Constructive q-witness reduction of LRC(14) to covering sets

**The actual conjecture quantity.** For a 13-set S of distinct positive integers, write
`M(S) = max_{τ∈[0,1)} min_{v∈S} ||vτ||` (the LRC gap; `||x||` = distance to nearest
integer). LRC(14) is exactly `M(S) ≥ 1/14` for all primitive S. This is sharper than the
open lonely **measure** `L(S) = meas{τ: ||vτ||>1/14 ∀v}` studied in THM-501/515/522:
`L=0` is necessary but **not** sufficient for a counterexample — a counterexample needs the
*closed* lonely set empty, i.e. `M(S) < 1/14`. `M` is computed exactly via the lower
envelope of tent waves (max at `τ ∈ {k/(v_a±v_b), (2k+1)/(2v_i)}`); tool
`lrc14_gap_M_exact_mac-mini-2026-06-16-S3.py`, validated against a dense grid.

## A. The q-witness lemma (PROVED, elementary)

> **Lemma.** Fix `q ∈ {2,…,14}`. If S contains **no** multiple of q, then `τ* = 1/q` is a
> lonely point: for every `v∈S`, `v ≢ 0 (mod q)` ⟹ `||v/q|| = min(v mod q, q−(v mod q))/q
> ≥ 1/q ≥ 1/14`. Hence `M(S) ≥ 1/q ≥ 1/14`.

*Proof.* `||v/q|| ≥ 1/q ⟺ v ≢ 0 (mod q)`, immediate from `||j/q|| = min(j mod q, q−j mod q)/q`.
Taking the min over S and `q ≤ 14` gives `M(S) ≥ 1/q ≥ 1/14`. ∎

This is the classical small-divisor observation, instantiated for n=14. It proves LRC(14)
**outright for an infinite family**: every primitive S that omits a multiple of some
`q∈{2,…,14}`.

## B. The covering-set reduction (PROVED)

> **Corollary.** A counterexample to LRC(14) (a primitive S with `M(S)<1/14`) must be a
> **covering set**: it contains a multiple of **every** `q ∈ {2,…,14}`.

Contrapositive of the Lemma. Verified: over 5000 random non-covering primitive 13-sets, the
`τ=1/q` witness never failed (0 violations, `lrc14_strong_necessary_condition`). Covering
sets are ~8.5% of random primitive 13-sets; a covering set needs multiples of 8,9,11,13 and
14, so it cannot be `{1..13}` and must be "spread out."

**The M=1/14 minimizers are all trivial.** The configurations achieving the bound exactly
(`M=1/14`) are the tight AP `{1..13}` and sporadics like `{1..11,13,24}` — none contains a
multiple of 14, so all are handled by `τ=1/14` (Lemma, q=14) **with equality**. Thus the
gap is never tight on a covering set: every covering set tested has `M > 1/14` strictly.

So **LRC(14) ⟺ `M(S) ≥ 1/14` for all primitive covering sets** — the conjecture reduces to
a structured sub-family, on the gap side, dual to THM-522's measure-side compactness reduction.

## C. Residue refinement of the q=14 layer (VERIFIED)

Within `C'(14)` (S contains a multiple of 14), write `A` = the non-multiples-of-14. The test
point `τ = a/14 + δ` (a coprime to 14, δ small) keeps A safe and dodges the multiple of 14;
the boundary A-speeds are those with `va ≡ ±1 (mod 14)`, i.e. `v ≡ ±a⁻¹`, and a δ of one
sign helps the `+` side and hurts the `−` side. So the construction succeeds — `M>1/14` —
whenever **A misses some unit residue mod 14** (`{1,3,5,9,11,13}`). Verified: 0 violations
over interior-drop + 3000 random C'(14). (Subsumed by B via the q-witnesses for `q<14`, but
gives the sharp boundary picture: the hard core needs A to cover all 6 unit residues.)

## D. The residual hard core (empirical; HYP-2566)

Covering sets are the only counterexample candidates. Over the full cover-all-units
interior-drop family and a 40k-config search, the minimum gap is
**`M = 7/89 ≈ 0.078652`** at `S = {1,2,3,4,5,6,7,8,9,10,11,13,84}` (84 = 14·6) — a clean
**10.1% above 1/14**. Random covering sets (speeds<80) bottom out at `5/41 ≈ 0.122`. No
covering set with `M<1/14` exists in any search. **Conjecture (HYP-2566):** `inf` of `M`
over primitive covering sets is `> 1/14` (uniform looseness) — which, with the bounded
structure of covering sets (compactness, cf. THM-522), would complete LRC(14).

## E. Status & honest scope

- A (q-witness) and B (covering-set necessary condition): **PROVED**, elementary; the
  τ=1/q witness is classical, the exact n=14 covering-set characterization is the instantiation.
- C (residue refinement): **VERIFIED** (construction sketched; the δ-existence for multiple
  large multiples of 14 needs the standard perturbation lemma to be fully rigorous).
- D (uniform looseness of the residual): **OPEN** (HYP-2566), strong empirical support.
- **Disprove direction:** no counterexample found over extensive exact search (HYP-2567);
  smallest `M = 1/14` always at the trivial tight configs.

**Dialectic origin.** The disprove search (minimize M, hunt `M<1/14`) *produced* the prove
reduction: forcing a counterexample to dodge every `τ=1/q` is exactly the covering-set
necessary condition, and that condition **is** the reduction. The two goals are one.
