# The edge is isolated, not thin

*mac-mini-2026-06-30-S41. The owner asked me to push the disproof mode while being precise about the razor-thin edge between proof and disproof. Being precise dissolved the "razor-thin."*

## What I expected, and what is true

We have called it the razor-thin edge for many sessions — the line a counterexample would have to walk, the floor that must stay positive, the measure that sinks toward zero. The picture in my head was a continuum: the floor getting thinner and thinner, a disproof creeping up on the boundary until, at the limit, it crosses. So I went to characterize that boundary exactly.

The boundary is `gap(O) = 0`. By the irreducibility of `Φ_p`, the apex gap vanishes at exactly one core: the full `Z_p`, the complete mod-`p` covering. One point. Every proper, partial covering has a strictly positive gap. And then the surprise: how close can a partial covering get to that point? For `p=7`, the answer is the doublet, `4cos²(3π/7) = 0.198`, and there is *nothing* between `0.198` and `0`. A spectral gap. The disproof boundary is not the end of a continuum the counterexample can slide along — it is an isolated point, fenced off by `0.198`, with empty space between it and everything else.

So for LRC14 the edge is **not razor-thin at all**. It is a clean, isolated point. A would-be disproof cannot approach it; it must land on it exactly — on the complete mod-`7` covering — making a jump of `0.198` with no foothold in between. And the complete covering is the apex cusp, where the lonely measure is zero but the `φ(n)` unit touch-points are still there: a lonely set of measure zero, nonempty, the counting witness. The disproof reaches the one point it could live on and finds the units already occupying it.

## Where the thinness actually is

Then where did "razor-thin" come from? Two places, both real, neither the discrete edge. The gap-`M` line is comfortably far — covering sets sit `10%` above `1/14`. The measure line is genuinely thin: `inf R' = 0` over the infinite family, the floor sinking to zero. But that thinness is a *product* — the total floor is the product of per-level gaps over the 2-adic descent, and a long descent drives a product of numbers each `≥ 0.198` down toward zero by sheer accumulation. It is soft thinness, a vanishing product, not a discrete approach to the boundary. The per-level factor never gets near zero; only their product does, and a vanishing product is closed by existence, not by a positive lower bound. LRC14 is **thin in measure, isolated in gap.**

## The genus is the thinness

The cleanest thing I learned is that the isolation is the genus. I tabulated the minimum positive gap by apex prime and it tracks `genus(X_0(2p))` exactly: `1.0, 0.382, 0.198` for `p=3,5,7` (genus `0,0,1`), the doublet each time; then it falls off a cliff to `0.0078, 0.0049` for `p=11,13` (genus `2`), and the minimizer is no longer a doublet but a large core like `{1,2,4,6,8,9}` mod `11` whose character sum nearly vanishes. At genus `≥ 2` the edge *does* go razor-thin — proper cores crowd up against zero. The genus literally measures how thin the disproof edge is, and `p=7` is the last value where it is still isolated. That is the precise sense in which LRC14 is the frontier: not the hardest case where the floor is smallest, but the **last case where the boundary is a clean isolated point** rather than a thin approach. One step further, to `n=22`, and the edge frays.

So the razor-thin edge, looked at exactly, is three different statements wearing one phrase. The gap-`M` edge: far. The measure edge: thin, but softly, a vanishing product closed by counting. The discrete gap edge: for LRC14, isolated by `0.198` — and that isolation is the genus, which holds at one only through `p=7` and then breaks. The disproof has nowhere thin to walk at `n=14`; it has only the single complete-covering point to stand on, where the units are waiting. The thinness that frightens lives one prime further out.

See [[the-small-measure-regime-is-the-heart]] (HYP-3615, the units close the edge), klein-S11 (HYP-3587, genus = the binding core), klein-S16 (HYP-3597, `inf R'=0` is a product), HYP-3548 (the two razor-thin lines), THM-590 (the apex gap). New: HYP-3700.
