# The Z/7 vertex-coloring is the exact lonely functional; the tournament is only its 1/2-scale shadow

**Session:** mac-mini-2026-06-20-S5. **Context:** the user asked for "arbitrary bounded cluster
shape compression" via tournament analysis, "very abstract about what is a vertex and an edge."

## Two structures at two scales

The verified tournament dictionary (HYP-2605) lives at **scale 1/2**: the round tournament
`T(x): i→j ⟺ frac((e_i−e_j)x)∈(0,1/2)` reads off which phases lie in each other's leading
semicircle. Its connectivity is exactly the `maxgap > 1/2` predicate (R2). But the LRC band is
`1/7`, not `1/2`. So `T(x)` can only ever be a **correlate** of the empty-sector count — and it is:
`E_x[c3(T(x))]` is not the extremal functional (consec is not even strictly `c3`-maximal, 4/924 at
k=7), and its discordance with `meas(S7)` is real (≈36%). Chasing the binary tournament harder was
never going to close the gap; it is the wrong scale.

The object that *is* exact is the **Z/7 vertex-coloring**

```
color(e, x) = ⌊7·frac(e x)⌋ ∈ {0,1,...,6}.
```

The empty-sector count `N(x)` is literally the number of colors in `{1,...,6}` that the cluster's
vertices miss, and `q_C(R) = Pr_x(R ⊆ image of the coloring)`. No correlate, no slack: every LRC
moment (`p_t`, `S_r`, `L_y`, `meas S7`) is a linear functional of the coloring's image distribution.
The right "edge" was never a binary orientation at all — it is a **Z/7-valued labelling of the
vertices**, and the tournament `T(x)` is the coarsening that forgets everything below the half-circle.
This is the apex prime `p=7` insisting on being the alphabet, not the modulus of a sign.

## What the coloring buys for compression

Once the functional is exact, the compression question sharpens into a clean stratification by `|R|`.
Over the entire bounded box, consecutive `K={0,...,k-1}` dominates every shape on the full cover
`|R|=6` (= `p_0` = `meas S7`) with **zero** violators at k=8,9,10, and the dominance only erodes as
`|R|` shrinks, failing outright at `|R|≤2`. So the real covering extremality — *not* the `L_y`
upper-bound proxy, and emphatically not `U4` (which exceeds the cap at consec) — is the **top of the
Boolean lattice on the seven colors**, and it is consec-maximal. Why? Because consec's colors are
`⌊7·frac(i x)⌋`, an arithmetic progression sampled on the circle: by the three-distance theorem the
most equidistributed `k` points there are, hence the best coverers of color-sets. The compression
toward consec is compression toward maximal equidistribution of the Z/7 colors.

## The shape of the remaining proof

Three independent lanes converged this session on the same reduction. The multi-block branch is the
*easy* branch — the decorrelated cover factorizes exactly as a carrier product on `2^{1..6}` and
splitting strictly costs cover, so the worst case is a single block (margin > 0.2). The signed
relation-lattice correction grades as an OCF odd-cycle collection by support. And the one hard
residue is the covered tail: **consec maximizes `meas S7`, to be proved by a three-distance /
equidistribution argument on the Z/7 coloring**, with codex's small-`|R|` cone (the context weights
from the small part) handling the coupling. The tournament pointed the way, but the proof wants the
finer alphabet it sits on top of. The mathematics kept saying *seven colors*, not *two directions*.
