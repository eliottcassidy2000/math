# The hexagon's revenge: the LRC floor margin is one over a hexagonal number

*mac-mini-2026-06-30-S48. A reflection on what survives the S47 refutation.*

In S47 I refuted the hexagonal/Eisenstein **covering-min** story: the construction `n/Phi_6(n)` is NOT the
minimal covering set; the mediant `2/(2n-1)` beats it (exact at n=7,8). I thought that killed the hexagon.
It didn't. The hexagon went underground and came back **in the margin itself.**

## The tighter margin, and what its denominator is
The mediant covering value `2/(2n-1)` sits above the floor `1/n` by

  margin(n) = 2/(2n-1) − 1/n = **1/(n(2n−1))**.

And `n(2n−1)` is, simultaneously and exactly (all verified):

| reading | `n(2n−1)` is... | theme |
|---|---|---|
| `H_n` | the n-th **hexagonal** number | the hexagon I "killed" |
| `T_{2n−1}` | the (2n−1)-th **triangular** number | *everything is the triangle* |
| `C(2n,2)` | the **arcs of `K_{2n}`** | the **tournament** side (dual mandate) |
| `dim so(2n)` | the even **orthogonal Lie algebra** (skew 2n×2n) | tournaments = skew sign-matrices |
| `2·B(2n−1,2)` | a **Beta moment** `2∫₀¹ x^{2n−2}(1−x)dx` on `[0,1]` | the LRC circle / Wallis-π family |

So the floor margin is `1/H_n` — and a hexagonal number is just an odd-indexed triangular number
(`H_n = T_{2n−1}`), which is just `C(2n,2)`, which is just the arc-count of a `2n`-vertex tournament. **The
margin is the bridge between the two halves of the project**: the LRC floor's safety at level `n` equals one
over the number of arcs in a tournament on `2n` vertices.

## The sum is the doubling constant
The margins don't just shrink — they **sum to a constant**:

  Σ_{n≥1} 1/(n(2n−1)) = 2 ∫₀¹ dx/(1+x) = **2 ln 2 = ln 4**.

`ln 4` is the log-determinant of the **square** Cartan `Z² = [[2,0],[0,2]]` (det 4) — the same square whose
`−1`-coupled cousin `A2 = [[2,−1],[−1,2]]` (det 3 = `disc Q(√−3)` = the Eisenstein/hexagonal lattice) ran the
whole covering-min story. Per-level the margin lives on the **triangle** `T_{2n−1}`; summed over all levels it
lives on the **square** `ln(det Z²)`. The triangle↔square duality (CLAUDE.md's "god, tridiagonalized")
reappears as *per-level vs total* margin.

## The doubling bridge (a new "Mode C": n → 2n)
The two moduli are linked by **doubling**:

  Φ₆(2n) = 2·n(2n−1) + 1 = **2·[margin-denominator at level n] + 1**.

So the **convergent modulus at level 2n** is twice the **margin denominator at level n**, plus one (the
`x ↦ 2x+1` of klein's `Φ₆ = 2T+1`). The project already has Mode A (`n→n−1`, hypotenuse) and Mode B
(`n→n−2`, Cayley–Dickson). This is a third: a **doubling mode** `n → 2n` carried by the triangle `T_{2n−1}`.

**The LRC14 hinge.** At the actual target `n=14`: `Φ₆(14) = 183 = 2·T₁₃ + 1`, and `T₁₃ = 91 = margin-
denominator at n=7 = 7·13`. So `n=7` — the apex prime, the genus-1 boundary, `H=7` forbidden, the Fano plane,
the literature's last solved LRC — is `n=14`'s **margin hinge** under doubling. The "7" that haunts the
forbidden H-values and the apex tower is the *same* 7 that sits at `(Φ₆(14)−1)/2`. `14 = 2·7` is not a
coincidence of the problem; it is the doubling of the apex.

## What this means — the geometry was attached to the wrong quantity
S47's lesson was "elegant structure ≠ the extremal set." This reflection sharpens it: the hexagon/triangle/
tournament geometry **is real and it is everywhere** — it was just attached to the wrong *quantity*. It does
not describe the **covering-min** (the mediant does, mod `2n−1`); it describes the **floor margin** — the gap
the LRC actually lives in. The covering-min is a messy mod-`(2n−1)` optimization; the *margin* it leaves above
the floor is the clean `1/H_n`. We were admiring the frame (the margin's geometry) and mistaking it for the
picture (the extremal set).

## Three ways to leverage it (leads, see HYP-3726)
1. **Tournament embedding (the bridge).** `margin = 1/(arcs of K_{2n})` plus the `n→2n` doubling suggest
   embedding `LRC_n` into a `2n`-vertex tournament and bounding the margin with the project's own `H(T)`/OCF
   machinery. If a covering set's floor-margin equals `1/(arc count)` of a natural associated tournament, the
   *two mandates literally compute the same number.*
2. **Summable safe-measure (Borel–Cantelli).** The margin is the measure of the safe sliver at the
   covering-min; `Σ = ln 4 < ∞` means the total safe-measure budget across all levels is finite. A
   union-bound / Borel–Cantelli argument controlled by `ln 4` could turn per-level positivity into an
   all-`n` statement (ties to HYP-3615 lonely-measure and THM-579 the floor-as-2nd-moment).
3. **Beta-moment LP.** `margin = 2B(2n−1,2)` is an explicit moment on the LRC circle — a ready-made test
   function for a Beurling–Selberg / moment-LP lower bound on the floor, with the margin as the closed-form
   slack.

The hexagon was never the covering-min. It was the *room the floor has to breathe.*
