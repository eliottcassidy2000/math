# The r=6 sharp horn is itself a mini Lonely Runner problem (death-star-S58)

A self-similar observation from work toward the r=6 uniform bound (THM-1132,
renumbered after THM-1123 had already been claimed).  The sharp horn and the
one-variable `G(σ)` bands are proved; the all-core, all-scale landing/drift
lemma is still open for arbitrary five-killer shapes.  THM-1134 has since
closed the entire step-two family over all cores and scales, and supplies a
general multiplier-chart cone plus a separated-ratio gate.

The r=6 covering-killer stratum, after the sharp measure horn, asks: within a core-safe arc,
do 5 far killers `{b, b+2, …, b+8}` leave a gap wider than one killer's danger arc? Writing
the killer phases at time `t` as `bt + m·(2t)`, they form an **arithmetic progression of step
σ = 2t**. So the question becomes:

> Can 5 points in arithmetic progression (step σ) on the circle simultaneously avoid a fixed
> arc of width `1/7` around the origin?

That is *exactly* a Lonely Runner / view-obstruction problem in miniature — 5 "runners" whose
positions are a dilation of a single time, all trying to be `≥ 1/14` from a marked point. The
whole apparatus we are using to attack LRC(**14**) reappears, one level down, as a question
about **5 runners against a width-`1/7`(= `2/14`) window**. The relevant quantity `G(σ)` — the
largest gap the 5 arcs leave — is minimized (`2/35`) exactly when the AP is equidistributed
(`σ=1/5`), the same "resonant/equidistributed is extremal" principle that governs the parent
problem, and maximized when the AP *collapses* (`σ=1/3` gives only 3 distinct points, `σ=1/2`
only 2), the same way LRC extremals are the arithmetically degenerate families.

Two things worth keeping:

1. **The 7 is doing double duty.** `1/14 = 1/(2·7)`; the danger window the mini-runners avoid
   has width `1/7`, and the sharp horn constant is `1/(7L)`. The modulus 7 that makes the
   Fourier kernel vanish on multiples of 7 (THM-1061) is the same 7 setting the window here.
   The problem keeps folding `14 = 2·7` back on itself at every scale.

2. **MAX not MEAN, again.** The reduction lands on an *existence* condition — `∃ t: G(2t)>1/7`
   — not an average. The bad set (`G<1/7`, near `σ=1/5,2/5,…`) has positive measure; loneliness
   survives only because the core-safe region can *reach* a good band. This is MISTAKE-129's
   lesson (good-period existence is a supremum) reappearing as the crux of a different route.

The recursion — a runner problem inside the runner problem, with the same extremal principle and
the same `2·7` arithmetic — is the kind of self-reference this project keeps surfacing. Whether
the descent continues (does the 5-runner window-problem reduce to a 3-runner one?) is a natural
question the sharp-horn frame now makes askable.

## Multiplier charts resolve the fixed-chart arity mirage

The mini-runner recursion sharpened after allowing all core-safe charts
`t=u/13`, not only `u=1`.  Every at-most-five-point residue pattern has some
nonzero multiplier with a six-unit cyclic gap; the exact proof has ten affine
orbits.  This produces a fixed Kakeya rectangle and the cone
`B>=17 max(A,80)`.  For the step-two AP pattern, a stronger 792-core rectangle
atlas plus exact finite complement closes every legal scale.  The lesson is
that the faithful vertex is not a killer at one frozen time, but an affine
residue orbit together with the selectable chart.  A single-chart tournament
forgets precisely the maximization that supplies the wide gap.
