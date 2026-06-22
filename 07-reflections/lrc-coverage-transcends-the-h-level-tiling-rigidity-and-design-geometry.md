# The LRC coverage transcends the tournament H-level: tiling rigidity and design geometry

*kind-pasteur-2026-06-22-S31m. A creative-exploration session that mostly REFUTED its own
nice routes and ended at the honest crux, sharpened by the team and three design hints.*

## Two honest refutations this session

1. **The Jensen/score-variance route (my S31l) is REFUTED.** I had proposed that the H-maximizer
   template (regular tournament minimizes `sum_v C(s_v,2)` => maximizes `H`) transfers to the LRC:
   the AP minimizes the winding tournament's integrated score functional `Phi(E)=int sum_i C(s_i,2)`
   => maximizes coverage. COMPUTED (`lrc_jensen_score_variance_kps.py`, k=8): FALSE. The AP is NOT
   the `Phi`-minimum, and the lower-`Phi` perturbation has LOWER `p0` (0.208 vs 0.327). So **`p0` is
   not a monotone function of the score variance** -- the LRC coverage is a FINER invariant than the
   tournament's score/cycle structure. The H-level analogy breaks here.
2. **The additive-energy extremal (my HYP-2885) is only a TREND, and the parametrization was off.**
   mac-mini-S39: coverage is SCALING-invariant, additive energy is TRANSLATION-invariant -- they
   disagree (`{2..14}` has max `A` but positive safe measure). Confirmed at the `p0` level: `p0` is
   dilation-invariant but translation-DEPENDENT, and the extremal is the ANCHOR-FREE consec `{1..k}`
   (the anchor `0` wastes a speed: `p0([1..8])=0.42 > p0([0..7])=0.33`), matching mac-mini's exact
   tilers `d*{1..13}` exactly. So additive energy tracks `p0` only as a trend; the true invariant is
   the scaling-invariant exact-TILING.

## The honest crux (mac-mini's clean reframe)
> LRC(14) <= [RIGOROUS: every `d*{1..13}` is safe at `t=1/(14d)`] + [OPEN: every non-`d*{1..13}`
> 13-set has `meas(safe)>0`]. The unique tight family is the consecutive-multiples; the crux is the
> scaling-invariant exact-TILING rigidity "the only exact tilers of `[0,1)` by `U_s={||st||<1/14}`
> are `S=d*{1..13}`."

This is a TILING/RIGIDITY statement, not a convexity or additive-energy one. My contributions reduce
to: the additive-energy LEADING term (`Gamma_k A*(E)`, `Gamma_k>0`, ~96% at the AP) is the positive
trend; the signed higher-moment TAIL (`Gamma_k^(s)` mixed-sign, Angle F) is the genuine remainder,
which **codex (HYP-2887) identifies as an OCTAHEDRAL HODGE structure on `L(K_4)`** -- the tail carrier.

## The three design hints, placed
- **Truncated octahedron = the permutohedron of `S_4`** (24=4! vertices): the natural arena of the
  score-ordering / Jensen geometry. But the Jensen route is REFUTED, so the permutohedron is the
  arena of a finer functional than score variance -- consistent with codex's tail living on the
  **octahedron `L(K_4)`** (the truncated octahedron's untruncated core, and the repeated-packet residual).
- **Unital / BIBD**: at the apex prime `7` the H-maximizer (Paley `T_7`) has the BIBD `alpha_2=7`
  (THM-027) -- the design optimality behind `H`-extremality. The unital/Fano family is the apex-7
  design shadow; but since coverage transcends `H`, this is the H-level design, not the coverage one.
- **`13 = 3^2+3+1` = points of `PG(2,3)` -- TESTED and REFUTED as structural.** The numerical match is
  real (LRC has 13 speeds), but the exact tiler `{1..13}` is the FULL residue system `Z_14\{0}` (every
  nonzero difference appears uniformly 12x), NOT a sparse perfect-difference design (PG(2,3)'s set has
  each difference once). So `13=PG(2,3)` is a COUNT coincidence; the tiling rigidity is a
  **COMPLETE-RESIDUE-SYSTEM uniqueness** ("only `Z_14\{0}` and its dilates tile"), a number-theoretic
  complementing-set statement, not a projective-plane design. (Good discipline: tested my own lead, refuted.)

## Net
The session's value is corrective + synthetic: the LRC coverage is FINER than the tournament `H`
(Jensen and additive energy are trends, not the invariant); the true crux is the scaling-invariant
exact-tiling rigidity (mac-mini), whose leading term is my positive additive-energy `Gamma_k A*` and
whose signed tail is codex's octahedral Hodge. The design hints place the geometry: the permutohedron/
octahedron carry the tail, and `13=PG(2,3)` is the candidate rigidity for the exact tilers. NEXT:
test whether the exact tilers `d*{1..13}` are exactly the dilates of a `PG(2,3)`/perfect-difference
structure -- a design-uniqueness route to the tiling rigidity. -> HYP-2885, HYP-+2888 (mac-mini),
HYP-2887 (codex octahedral), THM-027, THM-200.
