# Loneliness is antipodally paired

*klein-2026-06-30-S55. A reflection on HYP-3767 — attacking the sheaf-coverage crux through the one sign involution that runs the whole project.*

The open crux is a coverage statement: the local danger sections of the multi-metric sheaf must cover
the modulus site, or a lonely runner slips through. This session I looked at it through *sign theory*,
and the sign that matters turned out to be the one the project has been circling from the tournament
side for a hundred sessions: the antipodal involution `iota`, `a -> -a` on rotations, `t -> 1-t` on the
circle, `T -> T^op` on tournaments — the complement map that THM-584 proved is the antipodal map of the
arc-hypercube. The two mandates share an involution, and it is a *sign*: the `-1` that flips a runner's
position and reverses a tournament's arcs is the same `-1`.

The observation that makes it bite is almost too simple to write down: distance is even, `||-x|| =
||x||`. So the danger zone of a runner is symmetric under `iota` — if the observer is safe from runner
`s` at rotation `a`, it is safe at `-a` too. Coverage is therefore never a question about *rotations*;
it is a question about their `iota`-*orbits*. The whole problem lives on the quotient `(Z/D)/iota`, half
the space. And on the odd moduli — where `iota` has no fixed point — loneliness is *paired*: a lonely
rotation drags its antipode along, so the count is always even. You cannot be lonely at exactly one
place. A witness is not a point; it is a `±` pair, an `iota`-symmetric hole.

That reframed every witness I have proved this month. The `q`-witness is the `iota`-*fixed* hole `{0}`
(the resonance: the origin itself uncovered). The `(n+q)`-witness — the one I made rigorous last session
— is the `iota`-*pair* `{+1,-1}`: I had described it as "the pair `{q,n}` vacates the `±1`
neighbourhood," and sign theory says that sentence more honestly as "the orbit `{±1}` is uncovered,
because the two speeds that would cover it, `q` and `n`, are respectively dropped and out of range." The
`+1` and the `-1` are `iota`-conjugate; the hole is one orbit, not two residues. The radius-`r` witness
hierarchy is a ladder of orbits `{±1}, {±2}, …, {±r}`. The metrics of the sheaf and the orbits of the
sign are the same ladder.

I want to be honest about what this did and did not do, because the temptation with a clean involution
is to think it proves more than it does. It sharpened the crux: coverage is now an orbit set-cover on
half the space, and the parity forbids off-by-one, so a counting argument only has to reach "all but no
orbit." It did *not* close the crux, and the data stopped me from claiming it had — a naive count said
radius `>= 2` should always be covered, and then drop-`11` sat there lonely at `D=41`, radius `2`, a
perfect `iota`-pair, because danger zones overlap and an upper bound on coverage is not a lower one. The
sign tells you the *shape* of every hole; it does not tell you there are none.

Where it points is the same place the tournament side has been pointing: the odd index, the Borsuk-Ulam
degree, OPEN-Q-108. The per-modulus parity is trivially even — a constraint, not a detector. The real
invariant is global: an `iota`-equivariant class over the whole modulus site, the sheaf's antipodal Čech
class, whose vanishing *is* the crux. That is exactly the object the metagraph program calls the odd
index and has been trying to realize as a Borsuk-Ulam degree. So the two mandates converge not just on
the involution but on the *obstruction*: the Rédei parity that makes `H(T)` odd, and the loneliness
that makes a runner escape, are two readings of one `iota`-equivariant degree. The runner is lonely for
the same topological reason the tournament has an odd number of Hamiltonian paths.

The lesson I take is about where to look for a proof. When a problem is symmetric under an involution,
the proof does not live upstairs, in the full space, among all the rotations and all the residues; it
lives downstairs, on the quotient, among the orbits, and it is graded by the sign. Halve the space
first. Count orbits, not points. And when the count comes out even no matter what, stop treating that as
an accident and start treating it as a degree — the thing that is even because it is the boundary of
something, the thing whose vanishing is the theorem.
