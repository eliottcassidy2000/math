# The census splits: a rigid unit-skeleton plus one Jacobsthal doubling

*kind-pasteur-2026-06-28-S256. The owner asked to chase the ℚ(√−7) route and a few other LRC(14)
obligations, switching when one stalls. The ℚ(√−7) route gave a partial floor result; switching to the
census produced the cleaner finding — the entire tight locus is one rigid skeleton with a single degree
of freedom. Along the way a breakpoint bug nearly fabricated a false census, which is its own lesson.*

## The setup: what the census is, sharply

A 14-free 13-set has `M(S) ≥ 1/14` (THM-523, recovered as the equioscillation at the six units `a/14`,
`a ∈ (ℤ/14)*`). **Tight** means equality, `M(S) = 1/14`: the safety function never rises above `1/14`,
so the danger arcs `{t : ‖st‖ ≤ 1/14}` cover the circle with gaps exactly the six unit points. The
census conjecture is that the only tight sets are the AP `{1,…,13}` and GW `{1,…,11,13,24}`.

A prior session showed the tight/loose distinction is **not** a congruence condition: `12→26` keeps every
residue mod 14 (`26 ≡ 12`) yet jumps `M` to `1/12`. So the census lives in the *magnitude* layer. This
session asks: how much magnitude freedom is there, exactly?

## The finding: the runners split into a rigid skeleton and one free coordinate

Classify each AP runner by whether it **binds** a unit — i.e. whether `s·a ≡ ±1 (mod 14)` for some unit
`a`. That congruence forces `s` odd and `s ≠ 7`, so:

> **Binding runners = `{1,3,5,9,11,13}` = the units `(ℤ/14)*`.** Each pins one antipodal unit-pair (`s`
> and `14−s` pin the same pair). **Covering runners = `{2,4,6,7,8,10,12}`** — the six evens and the apex
> prime 7 — carry no unit; they exist only to cover the gaps between units and hold `M ≤ 1/14`.

And the rigidity is split cleanly along this line (verified by exhaustive single-replacement search,
`v ≤ 100`, with the *complete* breakpoint set):

> **Every binding runner is rigid. Every covering runner is rigid except 12.** The lone flexible
> coordinate is `12 → 24` — which is exactly GW. GW is itself single-swap-rigid (its only tight neighbor
> is `24 → 12`, back to the AP); there are no double-doublings and no genuinely new 2-swaps.

So the **entire freedom of the tight locus is one runner's Jacobsthal doubling.** The census is a rigid
unit-skeleton with a single hinge.

## Why 12, and only 12

The doubling `12 → 24` survives for a clean reason: `24 = 2·12` kills every point `12` killed
(`24·(j/12) = 2j ∈ ℤ`) *and* `24` is still 14-free (`24 ≡ 10`). It is the unique covering runner whose
double is 14-free and preserves the full cover. The other large covering runners fail — but *not* by
leaving a kill-point uncovered. `8→16`, `10→20`, `11→22`, `13→26` each create a lonely point at a
fractional location (`t = 20/23, 19/27, 9/25, 25/27`) where **no** AP runner vanished; the cover
redistributes globally and a sliver rises above `1/14`. So the failure of the other doublings is a
*global* covering fact, not a local one. That is the precise reason the census is hard: the surviving
hinge is explained locally, but the rigidity of everything else is global.

## What this buys, and the honest boundary

The split reduces the census to two sub-rigidities:

1. **The binding skeleton must be the units.** This is the residue layer — the same `(ℤ/14)*` /
   ℚ(√−7) structure that organizes the floor's resonances (the three binding pairs are the three
   Galois-conjugate `{QR,NQR}` pairs mod 7). Tight ⟹ `S` contains a runner of each unit residue; the
   open part is that their *magnitudes* must be the units themselves, not large lifts (`1`, not `15`).
2. **The covering layer admits only `12→24`.** This is the magnitude layer — the global covering /
   three-gap rigidity (HYP-2917's apex-lock Steinhaus rigidity, the gcd-window `[2,3]`).

Both sub-rigidities are still global covering statements; the split does not dissolve them. But it
**localizes** the census's entire degrees of freedom to one explicit hinge and names the skeleton, which
is a far more structured target than "minimize over all integer 13-sets." The apex prime 7 sits, tellingly,
in the *covering* layer (odd but non-binding, since `7a ≡ 0` or `7`), not the skeleton — the prime that
names the field is the prime that does the covering.

## The lesson that almost wasn't

The first pass used candidate maximizers `t = k/(2s_i)` and `k/(s_i−s_j)` only, and reported 23 "tight"
single-replacements — a false census. The maximizer of `min_i ‖s_i t‖` can sit at an `s_i+s_j`
*crossing*, so the complete breakpoint set is `{k/(2s_i)} ∪ {k/(s_i−s_j)} ∪ {k/(s_i+s_j)}`. With it, the
census collapsed to the single correct hinge. Incomplete breakpoints *underestimate* `M` — a one-sided
error that fabricates tights, the worst direction for a census. (MISTAKE-86.) The anchors `M(AP)=1/14`,
`M(GW)=1/14`, `M({1..11,13,26})=1/12` catch it.

## The pointer beyond

If the tight locus is one rigid skeleton with a single hinge, the natural object is the **deformation
space** of the cover: the skeleton (units) is infinitesimally rigid, and the one nonzero deformation
mode is the `12→24` doubling. A proof of the census would be a rigidity theorem — the cover's
"infinitesimal flex" is one-dimensional, generated by the Jacobsthal hinge — in the spirit of bar-joint
rigidity, but for a circle covering by arithmetic progressions. The three-gap theorem controls a single
runner's gaps; the census is the statement that the *joint* gap structure of thirteen progressions has
exactly one flex. That is the magnitude-layer rigidity the residue layer (ℚ(√−7)) cannot see.

— Related: [[lrc14-thread]], HYP-3258 (this split), HYP-3256 (residue layer / ℚ(√−7) skeleton),
HYP-3254 (floor splitting), HYP-2920/2921 (single-swap census), HYP-2917 (apex-lock Steinhaus rigidity),
MISTAKE-86 (breakpoints), THM-523; `the-tournament-spectrum-is-the-object.md`,
`lonely-runner-as-chebyshev-equioscillation.md`.
