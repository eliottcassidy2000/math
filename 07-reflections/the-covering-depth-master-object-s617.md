# The covering-depth master object (S617)

For a long time I treated "loneliness certificate," "tight configuration," "resonance relation," and "apex" as four
different objects that happened to interact. This session collapsed them into one, by asking the dumbest possible
question: what is a lonely time, *concretely*? It is a point of the clock that misses every runner's forbidden arc.
So LRC is a circular-arc covering problem, and there is exactly one object worth studying — the **depth
distribution** `p_k = meas{t : depth(t) = k}`. Lonely times are `{depth = 0}`, and everything else is a functional
of `{p_k}`.

The first thing the master object tells you is why LRC is hard, in one line. Each arc has measure `2δ`, so the depth
integrates to `2nδ = 2n/(n+1)`, which is bigger than 1. The arcs over-cover the circle by a factor of nearly two.
The union bound on `p₀` is *negative*. There is no first-moment proof, ever — the loneliness is entirely a
correlation effect. I find that clarifying: the difficulty isn't hidden, it's the very first inequality, and it's
visibly hopeless. You are forced upward, to the structure of the overlaps.

And the overlaps have a name. The collapse family — the configs where `p₀ = 0`, where the over-covering is so
perfectly arranged that no point escapes — turned out, on an exhaustive small search, to be *exactly* the additive
chains. Not just the arithmetic progression: the sporadic `(1,3,4,7)` and `(1,3,4,5,9)`, each top the sum of two
below. Every collapse set is closed under `a+b=c`; every set without such a relation has `p₀ > 0`. That is the whole
sub-problem the user pointed at, and it has a clean mechanism. The clock distance is subadditive — `‖(a+b)t‖ ≤
‖at‖ + ‖bt‖` — so a sum-relation literally chains three arcs together: pin `a` and `b` at the boundary and `c` is
forced inward. An additive chain is a tower of these pins, and a tall enough tower leaves no slack. I formalized the
subadditivity; it is one triangle inequality, and it is the engine of the entire collapse family.

What makes me trust the picture is that `a+b=c` is a `Σ mᵢ vᵢ = 0` resonance with coefficients `(1,1,−1)`. So the
collapse family *is* the resonance family. The thing the user has kept telling me to sidestep is not an obstruction
that appears in the analysis — it is the definition of the extremal set. `p₀ = 0` and "rich in small resonances" are
the same sentence read twice. And that snaps onto last session's apex sheaf: `p₀ > 0` iff the global certificate
section is nonempty iff `H⁰ ≠ ∅`. The depth distribution is the *measure-theoretic refinement* of the sheaf's `H⁰` —
the sheaf sees empty-or-not, `p₀` sees how empty. The apex (the multiple-of-n, σ-fixed lane) is the resonance that
the whole circle conspires around; the covering-depth `p₀` is what that conspiracy costs.

The singleton-wall test was the satisfying part. The altitude principle says a wall built from one relation should
open linearly — `(loglog)¹`, one averaging level. I perturbed `{1,2}` off its wall and `p₀` opened as `ε^{0.99}`.
Clean. The multi-relation walls bent below linear, the fingerprint of stacked relations. The exponent counts the
relations you break, which is the codimension, which is the altitude. Four sessions of "number of logs = number of
nested averagings" and here it is again, measured to two decimal places in a covering problem.
