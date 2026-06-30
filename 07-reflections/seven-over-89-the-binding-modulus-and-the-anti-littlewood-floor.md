# 7/89, the binding modulus, and the anti-Littlewood floor

*mac-mini-2026-06-29-S18. Asked to understand where 7/89 comes from, with the Littlewood conjecture as inspiration. Understanding it turned up a tighter covering set — and reframed the whole conjecture as the opposite of Littlewood. New: HYP-3551.*

## What 7/89 actually is

`7/89 = M(\{1,\dots,11,13,84\})` was the project's tightest known covering set. Take it apart and it is completely transparent. The set is the **skip-12 core** `\{1,\dots,11,13\}` — which has `M = 1/12`, because resonance 12 has no killer — completed to a covering set by `84 = \mathrm{lcm}(12,14)`, the smallest speed that kills both 12 and 14 at once. The optimizing time is `t^* = 37/89`, and there the runners `5` and `84` sit at the *same* distance `7/89` from the wall, because `84 \equiv -5 \pmod{89}` — so the **binding modulus** is `89 = 84 + 5`, the sum of the binding pair. The numerator `7` is the **packing radius**: at the best rotation, the thirteen speeds land on residues mod 89 that all avoid the 7-neighbourhood of 0, and 7 is the largest gap you can force. So `M = j/D` with `D` a binding-pair sum and `j` the covering radius of the speeds in `\mathbb{Z}/D` — the general law, made concrete. (`89 = F_{11}` is a Fibonacci red herring; the next set in the family is `14/173`.)

## Understanding it produced a tighter one

The skip-12 reading immediately raises a question: why skip 12? The denser the core, the closer its `M` to `1/14`, and the gentler the killer needs to be. The densest single-skip core is `\{1,\dots,12\}` — skip 13 — with `M = 1/13`, closer to `1/14` than the skip-12 core's `1/12`. Its required killer is `182 = \mathrm{lcm}(13,14)`, the *largest* minimal killer (13 and 14 are coprime), so it equidistributes and barely punctures the `1/13` structure. The result, grid-verified to four-million-point resolution:

> `M(\{1,\dots,12,182\}) = 14/183 \approx 0.07650`, only **+7.1%** above `1/14`,

**tighter than `7/89` (+10.1%)**. The recorded extremal was not extremal. The tightest covering set is the densest coverable core plus the minimal far killer that perturbs it least, and the asymmetry between skipping 12 and skipping 13 — `1/12` versus `1/13` — is the whole difference. The covering infimum is pinned near `1/13`, a fixed margin above `1/14`: the structural reason the conjecture's "uniform looseness" (HYP-2566) should hold.

## The Littlewood inspiration, and where it actually lands

Littlewood's conjecture is that for every pair of reals, `\inf_q\, q\,\lVert q\alpha\rVert\,\lVert q\beta\rVert = 0` — you can simultaneously approximate two numbers so well that even the *product* of the errors, weighted by `q`, vanishes. I went in expecting the tight LRC configs to live at the badly-approximable directions that make Littlewood hard — golden-ratio continued fractions, bounded partial quotients. They do not. The binding ratios `84/5 = [16;1,4]` and their siblings have *large* partial quotients; the tightness has nothing to do with bad approximability and everything to do with core density.

What Littlewood gives instead is the right *frame*. At the binding, `\lVert v_a t^*\rVert = \lVert v_b t^*\rVert = M`: the time `t^*` simultaneously approximates one wall by both runners, and the Littlewood product is exactly `M^2`. Littlewood says such a product can be driven to zero. The Lonely Runner says `M \ge 1/14`, so the product is `\ge 1/14^2` — it **cannot** be driven to zero. **The Lonely Runner is the anti-Littlewood**: where Littlewood asserts a simultaneous-approximation product vanishes, the LRC floor asserts the same kind of product stays positive. The tightest covering set is the LRC's closest approach to a Littlewood collapse — the place where the speeds most nearly trap every time near a wall — and it stops at `14/183`, not at zero.

This is the same theme the floor wears on the other side. Last session the floor itself turned out to be a positive Euler product (`\prod_p(1-p^{-2})`), positive because each factor is. Here the gap stays positive because the binding product cannot collapse. Both are statements that a *product* the rest of analysis expects to vanish is, for the Lonely Runner, bounded away from zero — the multiplicative-positivity signature of the whole problem. Littlewood is the conjecture that products of approximation errors collapse; the Lonely Runner is the theorem (where proved) that a dual product does not. They are two readings of how small `\lVert q\alpha\rVert` can be, pointed in opposite directions.

## The correction, recorded

`7/89` is not the tightest covering set; `14/183` is (so far), at `\{1,\dots,12,182\}`, and the gap-line margin is `+7.1%`, not `+10.1%`. The understanding is more useful than the number: tightest = densest coverable core + minimal far killer, the infimum sits near `1/13`, and the whole thing is the anti-Littlewood — a positive floor on a simultaneous-approximation product.

See [[the-two-razor-thin-lines-of-lrc14]] (HYP-3548, the gap line this sharpens), [[the-cut-prefix-is-a-continued-fraction-times-an-euler-product]] (HYP-3550, the Euler-product floor — the multiplicative-positivity twin), [[the-resonance-killing-game-and-the-zeta-duality-of-the-lonely-runner]] (the killer/survivor picture). New: HYP-3551.
