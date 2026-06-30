---
id: HYP-3726
title: The tighter LRC floor margin = 2/(2n-1) - 1/n = 1/(n(2n-1)) is the reciprocal of a HEXAGONAL number, and that one denominator wears five hats at once (all verified n=2..15): n(2n-1) = H_n (n-th hexagonal) = T_(2n-1) ((2n-1)-th triangular) = C(2n,2) (arcs of K_(2n) = the TOURNAMENT side) = dim so(2n) (even orthogonal Lie algebra / skew 2nx2n) = 2*B(2n-1,2) (a Beta moment 2*int_0^1 x^(2n-2)(1-x)dx on the LRC circle). The margins SUM to a constant: Sum_n 1/(n(2n-1)) = 2 ln 2 = ln 4 (= log-det of the square Cartan Z^2, det 4; the per-level margin lives on the triangle T_(2n-1), the SUM on the square ln4). DOUBLING BRIDGE (a new Mode C, n->2n): Phi_6(2n) = 2*[margin-denominator at level n] + 1 = 2*n(2n-1)+1; so the convergent modulus at 2n = twice the margin denom at n, plus 1. LRC14 HINGE: Phi_6(14)=183=2*T_13+1, T_13=91=margin-denom(7)=(Phi_6(14)-1)/2, 14=2*7 (the apex-7 doubled). LEVERAGE LEADS: (1) tournament embedding LRC_n -> K_(2n) (margin=1/arcs => bridge to H(T)/OCF, the dual mandate); (2) summable safe-measure (Sum=ln4<inf => Borel-Cantelli/union-bound control); (3) Beta-moment LP (margin=2B(2n-1,2) = explicit slack for a Beurling-Selberg/moment lower bound on the floor)
status: VERIFIED identities (exact, n=2..15; sum to 2ln2 numerically to 7 digits; Beta-moment exact). The three LEVERAGE leads are CONJECTURAL/exploratory (not yet executed). The margin is the MEDIANT margin (covering-min = mediant confirmed only n=7,8; the form 1/(n(2n-1)) is the gap mediant-value minus floor regardless).
source: mac-mini-2026-06-30-S48
related:
  - HYP-3725  # the refutation (convergent NOT covering-min); this is what survives -- the hexagon in the margin
  - HYP-3723  # klein: Phi_6 = 2T+1 (the convergent modulus is triangular); margin denom = T_(2n-1) -- triangle at two scales
  - HYP-3615  # lonely-measure (the safe sliver); leverage (2)
  - THM-579   # the floor as a 2nd moment; leverage (3) Beta-moment LP
reflections:
  - 07-reflections/the-hexagons-revenge-floor-margin-is-one-over-a-hexagonal-number.md
results:
  - 04-computation/floor_margin_hexagonal_macmini_20260630.py
  - 05-knowledge/results/floor_margin_hexagonal_macmini_20260630.out
---

# HYP-3726 -- the floor margin is one over a hexagonal number

S47 refuted the hexagonal/Eisenstein **covering-min**. This HYP records what survives: the hexagon is woven
into the **floor margin** itself. The owner asked to leverage the tighter margin creatively and hunt a cheeky
out-of-the-box connection -- here it is.

## The identity chain (verified exact, n=2..15)
`margin(n) = 2/(2n-1) - 1/n = 1/(n(2n-1))`, and the denominator is simultaneously:
- `H_n` -- the n-th **hexagonal** number (the hexagon I "killed" in the covering-min returns here);
- `T_{2n-1}` -- the (2n-1)-th **triangular** number (*everything is the triangle*, at the signed-LRC index `2n-1`);
- `C(2n,2)` -- the **arcs of `K_{2n}`**: the LRC floor margin = 1/(number of arcs in a `2n`-vertex tournament)
  -- the project's TWO MANDATES (runners and tournaments) meeting in one number;
- `dim so(2n)` -- the even **orthogonal Lie algebra** = the skew-symmetric `2n x 2n` space (tournaments ARE
  skew sign-matrices);
- `2 B(2n-1, 2)` -- a **Beta moment** `2 int_0^1 x^{2n-2}(1-x)dx` on `[0,1]` (the LRC circle; the Wallis/Beta
  family that gives the project its `pi`).

## The sum, and the doubling bridge
- `Sum_{n>=1} 1/(n(2n-1)) = 2 int_0^1 dx/(1+x) = 2 ln 2 = ln 4`. Per-level the margin sits on the **triangle**
  `T_{2n-1}`; the **total** is `ln(det Z^2)` -- the **square** Cartan (det 4), the `+1`-off-diagonal cousin of
  the `A2`/Eisenstein lattice (det 3) that ran the covering-min story. Triangle vs square, per-level vs total.
- **Doubling (Mode C, `n -> 2n`):** `Phi_6(2n) = 2*[margin-denom at level n] + 1 = 2 n(2n-1) + 1`. The
  convergent modulus at `2n` is twice the margin denominator at `n`, plus one (klein's `Phi_6 = 2T+1` with
  `T = T_{2n-1}`). A third reduction mode beside Mode A (`n->n-1`) and Mode B (`n->n-2`).
- **LRC14 hinge:** `Phi_6(14) = 2 T_13 + 1`, and `T_13 = 91 = margin-denom(7) = (Phi_6(14)-1)/2`. The apex-7
  (genus-1 boundary, forbidden `H=7`, Fano plane, last solved LRC) is `n=14`'s margin hinge under `14 = 2*7`.

## Three leverage leads (conjectural -- next sessions)
1. **Tournament embedding.** `margin = 1/(arcs of K_{2n})` + the `n->2n` doubling suggest mapping a covering
   `(n-1)`-set to a `2n`-vertex tournament so the floor margin becomes an `H(T)`/OCF quantity -- making the two
   mandates compute the same number. Test: is there a natural covering-set -> tournament map with
   `margin = 1/(arc count)` and `M` tied to `H(T)`?
2. **Summable safe-measure (Borel-Cantelli).** The margin is the measure of the safe sliver at the
   covering-min; `Sum = ln 4 < infinity`. A union-bound / Borel-Cantelli argument with total budget `ln 4`
   could lift per-level positivity to all `n` (ties to HYP-3615 lonely-measure, THM-579 floor-as-2nd-moment).
3. **Beta-moment LP.** `margin = 2 B(2n-1,2)` is an explicit moment on the circle -- a ready test function for
   a Beurling-Selberg / moment-LP lower bound on the floor, with the margin as closed-form slack.

## Caveat
This is the **mediant** margin. The mediant is the covering-min only at n=7,8 (proved); at n>=10 the true
covering-min trajectory is open (HYP-3725). But the *form* `1/(n(2n-1))` is the gap (mediant value − floor)
regardless, and the identity chain is unconditional. The cheeky content -- hexagon = triangle = tournament-arcs
= so(2n) = Beta-moment, summing to `ln 4`, doubling to `Phi_6(2n)` -- stands on its own.
