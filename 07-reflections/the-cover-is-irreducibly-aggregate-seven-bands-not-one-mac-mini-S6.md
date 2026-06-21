# The seven-sector cover is irreducibly aggregate — seven bands, not one

**Session:** mac-mini-2026-06-20-S6. The prompt asked to think in colorings and cellular automata,
and to mine the repo for inspiration. Both paid off — not with a proof, but with the cleanest map
yet of *why* the proof is hard, and the right coordinates to attempt it in.

## The object is a coloring, and the question is arithmetic

Last session pinned the exact functional: the Z/7 vertex-coloring `c(e,x)=⌊7·frac(ex)⌋`, with
`measS7(E)=meas{x : the colors hit all of Z/7}`. This session asked: among `k`-element offset sets,
why does the consecutive block `{0,…,k-1}` cover Z/7 the most? Every geometric instinct says
"because the AP is the most uniform point set." That instinct is **wrong here**: consec does not
maximize `meas{maxgap < t}` (13–178 shapes beat it). The cover is not a geometry of gaps on the
circle — it is the **arithmetic of residues mod 7**. The `e=0` clock pins residue 0; the measure is
scale-invariant but *not* translation-invariant; the seven matters because it is the modulus, not
because it sets a scale. The right tools are Sturmian partial sums (THM-536) and the multiplicative
structure of Z/7, not the three-gap theorem.

## Seven bands, each a twist of one

The one durable structural gain is the **slope-band decomposition** (HYP-2703, proved): splitting
`θ=7x∈[0,7)` into its seven unit slope bands `s=⌊θ⌋` gives

```
measS7(E) = (1/7) Σ_{s=0}^{6} bandcover_s(E),
bandcover_s(E) = meas{ f∈[0,1) : {(s·e + ⌊e·f⌋) mod 7 : e∈E} = Z/7 }.
```

Every band is the **same** pure (slope-0) Sturmian cover problem, twisted only by the multiplicative
map `e ↦ s·e mod 7`. So the intimidating 64-term inclusion–exclusion certificate collapses to a
**seven-term signed sum**, and the seven terms are one problem seen through the seven multipliers of
Z/7. This is the kind of object the project keeps finding — the apex prime insisting that its whole
multiplicative group is the symmetry, the way QR_7/Paley is the extremiser for the tournament count
(THM-126). I hoped the quadratic character would organize the trade. It does not: consec loses band
`s=1` (a residue) and wins band `s=3` (a non-residue). The split is not QR; it is **slow versus
fast**.

## Why it cannot be cut into pieces

The session's real yield is negative, and unusually complete. Consec maximizes `measS7` — verified
exactly, zero violators, k=8,9,10 — but it maximizes it through **no sub-structure whatsoever**.
It is not the per-IE-block extremiser. It is not the per-band extremiser. There is no monotone
compression path to it, in value space or index space. It does not stochastically dominate the
empty-sector count, and its coverage-count distribution does not majorize — consec dominates *only*
at the very top, the event `|colors|=7` that is `measS7` itself, and loses every weaker threshold.
Each of these was a plausible proof route; each is now a dead end with an exact counterexample. The
extremality lives only in the full aggregate.

And the mechanism, finally legible, is a **trade across the seven bands**. Consec reaches less far in
`e`, so in the slow bands `s=0,1` — where the Sturmian walk is still loitering near residue 0 — a
shape with a larger top element covers more. But in the fast bands `s=2..6` the twist `e↦se` flings
consec's regular residues clear across Z/7, and it wins those five bands by more than it loses the
two. Five small wins outweighing two small losses: that is the whole theorem, and it is why no single
band, block, or step ever sees it. The cover is a cellular automaton of `k` clocks on a seven-ring,
and consec is not the configuration that is locally best anywhere — it is the one whose seven global
exposures net out ahead. The proof, when it comes, must add the seven bands up with their signs; the
mathematics has been telling us, in every refuted shortcut, that there is no smaller place to look.
