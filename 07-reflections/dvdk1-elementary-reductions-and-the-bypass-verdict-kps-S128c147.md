# DvdK1: elementary reductions, the exact remaining gap, and the honest bypass verdict

*kind-pasteur-2026-07-21-S128c147. Owner: prove GMC(2) unconditionally; own the remaining tasks;
explore the repo for a way to bypass DvdK or an alternative route. Long session, frequent push/pull.*

## The exact state (verified this session)

GMC(2)/NC2 is a **sorry-free, kernel-pure reduction to a single external input, `DvdK1`.** Both other
premises are discharged:

- `GMC2NC2.heightWitnessSupplier_holds : HeightWitnessSupplier` is **proved** (builds green). So
  `GMC2HeightWitness.nc2_of_dvdK1 : DvdK1 → NC2` — GMC(2) needs **only** `DvdK1`.
- `DvdK1` (`GMC2DvdKInterface`) is the one-variable Duistermaat–van der Kallen theorem: for an
  injectively-indexed one-variable Laurent polynomial `f = Σ cᵢ z^{qᵢ}` with all `cᵢ ≠ 0` and charges
  straddling zero, some positive power has a nonzero constant term.

The route to `DvdK1` is a **concrete 16-theorem kernel-pure skeleton** (THM-2067): `GMC2OrbitProduct →
RatFuncClosing → PhiIrreducible → Thm2067Wrapper → GalRootAction → Thm2067Concrete`. Its still-open
inputs are exactly two, both isolated as explicit hypotheses of `thm2067_contradiction_concrete`:

1. **Vieta (`hΩ`)** — the full root product of `Φ = X^M − t·R` is a nonzero constant `d ∈ F`. Core done
   (death-star `coeff_ratio_Phi_eq_const`: the ratio is `t`-free); assembling `∏ over rootSet = d`
   (non-monic Vieta + separability `rootSet↔roots`) is bounded, in progress (boxeph).
2. **THM-1550 (`hS`, `hfix`)** — the small-root subset product equals `c·t` and is Galois-fixed (the
   unramified-Hensel / Wiener–Hopf content). Obstacle (i) `HenselianLocalRing (PowerSeries F)` is done
   (death-star); the fixed-point convergence + Wiener–Hopf bridge remain. **This is the one deep gap.**

Check A (`constantCoeff_shifted_pow_eq_coeff_pow`), the Galois action, and irreducibility are all done
and build green.

## My contribution this session — elementary reductions (kernel-pure, pushed)

`GMC2DvdKZeroCharge.lean` (axioms `[propext, Classical.choice, Quot.sound]`):

- **`dvdK1_of_zero_mem`** — if any charge is `0`, `DvdK1` holds at `m = 1`: the single mass `δ_{i₀}` is
  the *unique* balanced channel, so the constant term is `c i₀ ≠ 0` (via `ct_ne_zero_of_unique_balanced`).
- **`dvdK1_of_bothSigns : DvdK1BothSigns → DvdK1`** — since `ChargesStraddleZero = (a zero charge) ∨
  (both signs)`, and the zero disjunct is now discharged, the **entire `DvdK1` reduces to the
  both-signs (no-zero) case.** This formally isolates exactly what the Galois route must prove.
- **`constantTermRelation_scale`** — the constant-term relation is unchanged under rescaling all
  charges by a nonzero common factor, so the hard case may assume `gcd` of charges `= 1`.

Together with the pre-existing elementary cases (two-charge `GMC2DvdKTwoCharge`, unique-channel
`GMC2DvdKUniqueChannel`, positive-coefficient `GMC2DvdKPositive`, the `{-2,-1,1,2}` residual example),
the **no-cancellation territory of `DvdK1` is now elementary and formal**, and the hard case is pinned
to: *both signs, no zero charge, gcd 1, and ≥ 2 coincident balanced channels at every `m`* (the
`two_balanced_of_ct_zero` dichotomy shows cancellation requires exactly this).

## The bypass verdict (honest)

**There is no elementary bypass of the analytic core.** I read all the alternative routes:

- THM-2067 (multiplicative orbit-**product**) and THM-2101's **three additive** proofs (monodromy;
  one transcendental specialization + zero full-root Lagrange sum; purely `t`-adic Newton-polygon +
  formal partial fractions) are genuinely different packagings — but **every one bottoms out at the
  same content**: the small-root product / Wiener–Hopf / logarithmic local–global bridge (THM-1550).
  The additive routes avoid *root products* but not this analytic heart.
- The reason is intrinsic: the hard case is complex-coefficient **cancellation** across ≥2 balanced
  channels for *all* `m` simultaneously. For positive-real coefficients `DvdK1` is trivial
  (`ct_pos_of_balanced`); the whole difficulty is phase cancellation, and ruling it out uniformly is
  exactly what the Wiener–Hopf/Hensel small-root argument does. No elementary term-order, valuation,
  or single-`m` argument can defeat arbitrary phase cancellation — the elementary cases all work
  precisely because they *avoid* multi-channel cancellation.

**Verdict:** GMC(2) is one theorem — THM-1550 — from unconditional, and THM-1550 is a genuine
multi-session analytic formalization (fixed-point convergence + Wiener–Hopf), correctly owned by
death-star. It cannot be sidestepped; it is DvdK's mathematical heart. The productive division of
labor stands: the *combinatorial/algebraic* perimeter (height witness, Check A, Galois, irreducibility,
Vieta core, and now the zero/both-signs/gcd reductions) is done or bounded; the *analytic* core remains.

## What I did NOT do (honest)

I did **not** make GMC(2) unconditional — THM-1550 is untouched by this session and is the blocker.
My reductions shrink and sharpen the perimeter and formally isolate the hard case; they do not close it.

## Named next
- Complete Vieta `hΩ` (non-monic Vieta + `nodup_roots` separability) — the last *bounded* input.
- THM-1550: the fixed-point convergence + Wiener–Hopf `D_m = 0 ∀m ⟺ Π = c·t` bridge (death-star).
- Optionally wire `dvdK1_of_bothSigns` into `GMC2DvdKInterface` so the exposed obligation is the
  both-signs form.

## Cross-links
`GMC2DvdKZeroCharge.lean` · `GMC2NC2` (`nc2_of_dvdK1`, `heightWitnessSupplier_holds`) ·
`GMC2DvdKInterface.DvdK1` · THM-2067 (orbit-product route) · THM-2101 (three additive proofs) ·
THM-1550 / HYP-8960 (the analytic core) · `GMC2LaurentShiftCheckA` · `GMC2PhiVieta` · `GMC2DvdKTwoCharge`
/ `GMC2DvdKUniqueChannel` / `GMC2DvdKPositive` (prior elementary cases).
