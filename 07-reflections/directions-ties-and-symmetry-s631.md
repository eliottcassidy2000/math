# Reversible runners, the unit-distance trienement, and where the symmetry hides (S631)

Two reframes were handed to me, and the instruction was to work them back and forth until they
opened. They opened onto the same room.

**Reframe one: let the runners reverse.** The first thing you notice is that you can't — a single
runner flipping direction is invisible, because `‖−v t‖ = ‖v t‖`. The observer cannot tell clockwise
from counterclockwise for any one runner. The user warned that tricks live here, and they do: the
direction is invisible *individually* but decisive *pairwise*. Two runners going the same way separate
at the difference of their speeds; two going opposite ways separate at the sum. So a choice of
directions — every other runner clockwise, say — is not a symmetry to be quotiented away; it is a
*tournament orientation*, and it decides, pair by pair, whether the relevant gap is a difference or a
sum. I computed the reversible mutual-loneliness gap over all orientations of the arithmetic
progression, and the modulus told the whole story: all-same-direction lives at the difference modulus
`n` with the familiar gap `1/(n−1)`, and the moment you reverse runners the witness denominators pick
up exactly `2n−1` — 7, 9, 11 at n = 4, 5, 6. That is the pair-sum-sieve modulus the canon has been
orbiting for a year. The shell-`2n−1` framework, the prime-shell dodge, the `±`-involution that is
complex conjugation — all of it is the *reversible* Lonely Runner. We had been studying the
opposite-direction problem without saying the word "opposite." The `±` was the reversal the whole time.

**Reframe two: read unit distances as a trienement.** Orient each pair of points by whether they are
closer than one (arrow), exactly one (tie), or farther than one (arrow back). The unit-distance graph
is the tie subgraph, and Erdős asks for the most ties. In this language the disproof of the grid stops
being a list of lemmas and becomes a single sentence: maximize the tie-degree. And the tie-degree of a
lattice is exactly its count of norm-one vectors — the units, the roots of unity, the modulus-one
elements. The triangular lattice gives six ties per interior point and the square lattice four, which
is precisely why triangular wins at small scale: `ℚ(√−3)` has six units, `ℚ(i)` has four. The tie-graph
is a Cayley graph, and its tie-degree *is* its symmetry. So "complex-multiplication beats the grid"
becomes "a tie-graph with a larger automorphism group beats one with a smaller one," and the engine of
the disproof — bounded-height modulus-one algebraic numbers, manufactured in quantity by a class field
tower — is just the statement that you can keep the tie-degree growing while the field grows.

Does the trienement simplify the disproof? I want to be honest about this, because it is the user's
direct question. The *statement* simplifies, cleanly and usefully: maximize ties = maximize the
modulus-one supply = maximize the symmetry. That is genuinely the right frame, and it is *our* frame —
the perspective/rigidity key, automorphism as the controlling invariant. But the *proof* does not
simplify. That class field towers actually deliver a superlinear supply of bounded-height
modulus-one elements — Golod–Shafarevich, bounded root discriminant, split primes — is number theory
the combinatorial trienement cannot reach down and touch. The reframe tells you what to maximize and
why it is a symmetry question; it does not hand you the tower. I would rather say that plainly than
claim a simplification that isn't there.

What the two reframes share is the gem. Both problems are tournaments with a tie relation. In both,
the tie is the *resonance* — loneliness pinned at `1/n`, distance pinned at one — and the arrows are
the metric order around it. And in both, the density of ties is set by a *symmetry group*: the
orientation group for the runners, the automorphism group for the points. The runner reframe shows the
orientation choosing the modulus; the point reframe shows the automorphism choosing the tie-count.
Underneath both sits the same arithmetic object — the modulus-one elements of a complex-multiplication
field, appearing literally as norm-one vectors for unit distances and as the `(ℤ/(2n−1))*` multipliers
of the shell for the Lonely Runner. Trienement-plus-symmetry is the skeleton; CM modulus-one elements
are the muscle. The user asked whether reversing runners and tying distances were two ideas. They are
one idea, seen from the two problems we have been carrying, and the place they meet is the perspective
key wearing a new coat: the symmetry that ties resonances together is the symmetry we have been calling
rigidity all along.
