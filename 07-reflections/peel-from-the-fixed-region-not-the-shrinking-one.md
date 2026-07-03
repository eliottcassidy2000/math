# Peel from the fixed region, not the shrinking one

*mac-mini-2026-07-02-S18. Reflection on HYP-3875 (the simultaneous union-bound far-peel).*

The intermediate band needed a bound on what a handful of fast runners do to a lonely set, and the program
already had a peel — one runner at a time, each removing its danger, each leaving the region a little
smaller and a little more fragmented. It is the obvious way to compose perturbations: apply the first, look
at what is left, apply the second to that, and so on. And it does not work past the first step. A fast
runner's danger is a comb with as many teeth as its speed; subtracting it from a region shatters that region
into as many pieces as there are teeth. So after one peel the region has a number of pieces of order N, and
the next peel's boundary fee — which is proportional to the number of pieces — is of order N times one-over-N,
which is a constant. The fee stops shrinking. Each further peel costs the same fixed amount no matter how fast
the runner, and a fixed cost per peel, summed over the runners, sinks the whole estimate. The iterated peel
is a good idea that eats itself.

The fix is a change of what each peel looks at. Do not peel the second runner from the shattered remains of
the first. Peel every runner from the *original* region — the one lonely set of the bounded core, fixed once,
never fragmented — and take the union of the damage. The measure of what survives all the dangers is at least
the measure of the region minus the sum of the measures of the individual dangers on it. Each of those is the
single-runner rate applied to the same fixed region with its same small number of pieces, so each fee is the
small original number of pieces times one-over-speed, and the total is that number of pieces times the sum of
the reciprocal speeds. Nothing compounds. The fee is linear in the number of runners and vanishes as the
speeds grow. It is the difference between a quantifier that reads "for each runner, its effect on what the
previous runners left" and one that reads "the sum over runners of each effect on the original" — and the
second is bounded where the first is not.

In Lean the two are almost the same length and utterly different in what they cost. The iterated peel is a
theorem that peels one element and must then be folded, and the fold carries the growing piece-count as an
invariant that has to be controlled at every step. The union peel is a single subtraction of the whole
danger-list from the fixed region, one application of the lemma that subtracting a set loses at most its
intersection, and then the intersection distributes over the concatenated dangers into a plain sum. The hard
part — the growth — simply is not there, because the region never changed. Three short structural lemmas and
the single-runner rate, and the multi-runner bound falls out. The proof is easy because the mathematics was
made easy, and the mathematics was made easy by refusing to let the object evolve under the operations.

There is a price, and it is worth naming honestly. The union bound is crude: it gives a floor that decays
linearly in the number of runners, one minus the runner-count over seven, where a sharper analysis that
tracked the true joint density would give a floor decaying like a product, seven-sixths to the power of the
count. The linear floor dies at seven runners; the product floor would survive further. But the band that
actually binds has at most six fast runners, and for six the linear floor is still positive, and positive by
a finite margin that a finite search can reach. The cruder instrument, applied to the fixed region, closes
the case the sharper instrument was being built to close — and closes it now, sorry-free, on machinery that
already existed, instead of after the long telescoping the sharp bound demands. Sharp is not the same as
necessary.

The pattern that transcends the theorem: **when you compose many perturbations, do not apply them in sequence
to the evolving object — each step would inherit all the damage of the ones before, and the damage compounds;
apply every perturbation to the original, fixed object and sum the effects.** The union bound is not a weaker
tool that you settle for. It is the tool that keeps the object still, and an object that does not move under
the operations is one whose cost you can actually count. Reach for the fixed region before you reach for the
sharp constant; often the case is already closed by the crude estimate that never let anything move.
