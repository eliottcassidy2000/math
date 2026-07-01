# Five twos are one corner

*klein-2026-07-01-S73. A reflection on HYP-3806 — bringing Chebyshev equioscillation to the covering-min,
and finding every framework saying "2."*

Asked to find more abstract things to care about and apply them to proofs, I reached for the oldest lens in
approximation theory: Chebyshev's. The covering-min `M(S) = max_t min_v ||vt||` is a minimax, and minimax
problems have a signature — the extremizer *equioscillates*: it is held in place by a set of active
constraints, the alternation set, and the size of that set measures how rigid the optimum is. So I computed
the alternation set of the construction — the runners that bind at the optimal time `t*` — and it was, for
every `n`, exactly two runners: the slowest, `1`, and the killer, `n(n-1)`. Two, and no more.

Then I noticed that "two" was already everywhere. The lonely measure has two atoms, at `t*` and `1-t*`. Its
Verblunsky sequence terminates at `|alpha_1| = 1`, the orthogonal-polynomial signature of a two-point
measure. The minimax has exactly two global maximizers. And the two binding runners are the two solutions
of `v ≡ ±1 (mod Phi6)`. Five different frameworks — Chebyshev alternation, atomic support, OPUC
termination, the maximizer count, the residue equation — each independently produced the number two. When
that happens, the number is not a coincidence of one calculation; it is the invariant, and the frameworks
are five photographs of one object. Here the object is a length-two alternation at the hexagonal binding,
and the single fact underneath all five is the killer identity `n(n-1) ≡ -1 (mod Phi6)`: the killer is not
an ad-hoc large speed the construction bolts on, it is precisely the `v ≡ -1` partner of the slowest runner
`v ≡ +1`, the second foot of the alternation.

Two is the smallest a nondegenerate alternation can be, and that is the point. A minimax whose extremizer
equioscillates on the *minimal* number of constraints is at a *corner* — a vertex of the feasible region,
isolated, admitting no continuous deformation that preserves optimality. This is exactly the "deep-well
isolation" the search sessions kept hitting: five thousand random covering sets, none within striking
distance, the construction alone at the bottom of a narrow well. The reason is now structural rather than
empirical — the construction is a corner solution of a Chebyshev problem, and corners are isolated by
definition. The rigidity people measured by sampling is the rigidity Chebyshev theory predicts from the
alternation length.

The second gift of the session was a name for a thing that spans the whole project. Both halves of this
repository — lonely runners and tournaments — turn out to have a *covering excess*: the amount by which a
covering constraint pushes the extremal value above the free floor. For the runners it is exactly
`M_C - 1/n = (n-1)/(n·Phi6)`, and the numerator is the dropped speed, the very partial quotient of the
continued fraction `[0; n-1, n]`. For the tournaments it is `rho(n) - ceil(log2|G_n|)`, the flip-rank's
overshoot past the information floor. In both, "free" is easy and "covering" is dear, and the excess is set
by the arithmetic (the modulus `Phi6`) or the symmetry (the group `S_n`). Naming it makes the two problems
visibly the same shape: an extremal value equal to a floor plus a covering excess, with the excess carried
by the object the construction is forced to add — the killer, the seventh arc.

What this buys the proof is a sharper target than I have had. The covering-min lower bound is now a
Chebyshev / linear-programming duality statement with a *two-point* dual: a certificate supported on the
runners `{1, killer}` at the times `{t*, 1-t*}`. A beater would have to break a length-two alternation that
the residue equation `v ≡ ±1 (mod Phi6)` pins to the hexagonal lattice, and do it while covering every
modulus up to `n` — two obligations that meet on the same finite lattice and, the structure suggests,
cannot both be satisfied. I did not close that; it is OPEN-Q-108 still, now wearing its Chebyshev clothes.
But the clothes fit better than any I have put on it: a finite, two-dimensional dual on a finite phase
group is the kind of object one can actually hope to construct. The lesson is that the most useful abstract
concept is often the one that tells you how *small* the remaining problem is — Chebyshev did not solve the
runner, but it told me the answer hinges on two points, and two points is a place you can stand.
