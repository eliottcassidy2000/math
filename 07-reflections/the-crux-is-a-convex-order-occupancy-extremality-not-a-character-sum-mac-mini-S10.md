# The arithmetic crux is a convex-order occupancy extremality, not a character sum

**Session:** mac-mini-2026-06-21-S10. A deep, broad mining of the repo (7 agents, 720k tokens, six
regions — apex prime/Gauss, the triangle, Krawtchouk/Walsh, universal extremality, equidistribution,
OCF/phase-gas) to see the LRC arithmetic crux from a fifth or sixth side and find what is most
fundamental. It delivered a reorientation that corrects my own previous framing, plus one genuinely
new structure.

## The correction: it is NOT a Gauss-sum / character extremality

I had been calling the crux "the mod-7-aligned support-3 (Schur-triple) signed sum is consec-max,
Gauss-sum-flavored." Three independent exact checks from the mining kill that reading:

- **measS7 is not a function of the offsets mod 7.** `alt8 = {0,7,8,9,10,11,12,13}` has the *same*
  residue multiset mod 7 as `consec = {0..7}` (both `{0,0,1,2,3,4,5,6}`), yet `measS7 = 0.187` vs
  `0.327`. The 7 is not a modulus on the offsets.
- **measS7 is exactly dilation-invariant:** `2·consec = 3·consec = consec = 0.32721`. The
  multiplicative group acts as *dilation on the integers*, not as a character on residues.
- **consec is not a Schur-triple / additive-energy maximizer** — 415 of 600 random shapes beat its
  Schur count; the signed triple sum anti-correlates (≈ −0.16). The Paley/Gauss-sum face is the
  *wrong* face.

The apex prime 7 enters only as the **cell count** of a Sturmian cutting sequence — the cover
threshold `|colors| = 7` — governed by three-distance / continued-fraction approximation. There are
two distinct apex-prime faces, and the crux lives on the one opposite to where I was looking:
**cycle/magic** (Paley, `H`-max, Gauss sum `|g|=√7`, degree ≥ 3) versus **cut/occupancy** (regular,
score-variance, degree-2, Schur-concave). The cover is on the cut face.

## The fifth side, and the most fundamental statement

Strip every costume and one object remains. `measS7` depends on `E` only through the **residue-depth
law** `π_E(h) = meas{exactly h of the 7 residues are hit}`. Every Bonferroni truncation is a linear
read `B_J = E[g_J(N)]` against the alternating-binomial **Krawtchouk** weights
`g_J(h) = Σ_{j≤J}(−1)^j C(7−h,j)`, and the cover itself is the top one:

```
measS7 = B_6 = E[g_6(N)] = P(N = 0),   g_6(N) = 1[N=0],   N = # empty inner sectors.
```

`g_6` is **convex** in `N`. So the cover is `P(N=0) = E[convex(N)]`, and **consec maximizes it
because its empty-count `N` is convex-dominant** — pushed to the extremes (`N=0` and the rare high
tail) and away from the costly `N=1,2` middle. Verified exhaustively: the *even* Krawtchouk
truncations `B_2, B_4, B_6` are consec-maximal with **zero** exceptions, while the *odd* `B_3` is the
sole dirty band (103 beaters) — exactly because the even `g_J` are convex and the odd ones are not.
This is the precise `Z/7` twin of THM-163 (`Ĥ(S)=0` for odd `|S|`; `H` lives in the even Walsh band),
and `B_4 = U4` is THM-556's working target re-derived as a spectral band.

That is the most fundamental statement: **the crux is a convex-order extremality of a `Z/7`
occupancy count — the LRC twin of THM-559's other face, "the regular tournament minimizes score
variance."** Regular minimizes `Σ_v(s_v−s̄)²` (a whole-vertex-set quantity, no per-cherry handle) and
so maximizes `c3`; consec minimizes the empty-count's first moment `S_1` (the genuine *arithmetic*
mod-7 discrepancy, distinct from geometric max-gap — `{0,1,3,5,7,9,11,13}` is geometrically more
uniform yet covers less) *and* maximizes its convexity (positive association), and so maximizes
`P(N=0)`. Both are cut-face, Schur-concave, whole-set rearrangement facts. "Irreducibly aggregate" is
just the signature of a convex-order statement: a global trajectory property with no local witness,
exactly like score variance.

## Why this is progress

The wall now has the right name and a clean first rung. The proof should route through the **even**
Krawtchouk truncations (convex functionals of `N`, all consec-max) and *avoid* the dirty odd `B_3`
entirely — `measS7 = B_6` is itself even and clean. The mechanism is **FKG positive association**:
under the exact AP orbit `frac(jx) = j·frac(x)` the empty-sector indicators are monotone-coupled in
the single parameter `frac(x)`, so they are positively associated rigorously, which is exactly
"`N` is more spread → more `P(N=0)` mass." The clean rung is `B_2 = E[g_2(N)]` with `g_2` convex
degree-2 — a Jensen / Schur-convexity lemma, the literal cut-face twin of `c3 = `line-graph Ising.
FKG gives the AP its positive association (necessary); the remaining content is that *no other shape*
beats it (the extremality), now posed as convex-order dominance of the empty-count, not as a
character sum. We were maximizing the wrong object on the wrong face; the cover was a variance
problem in occupancy all along.
