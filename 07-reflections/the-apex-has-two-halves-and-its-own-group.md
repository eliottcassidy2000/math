# The apex has two halves, and its own group

*klein-2026-07-01-S85. A reflection on HYP-3830/3831 — the hard apex splits into a √p certificate half and a
π measure half, its symmetry group is a subgroup of the Hurwitz group, and the game framing names the crux as
a discrepancy.*

This session had the shape of a place rather than a step: the owner pointed at the apex — the `{7,21}`
obstruction, the Paley heptagon, the `√21` residual — and asked me to look at it from every side at once. What
came back is that the apex is not a single thing but a place where several structures meet, and the useful move
is to name each face rather than to reduce them to one.

The cleanest new fact is that `√p` has four faces and they are literally the same number. At `p=7`, the value
`2.6458` is the imaginary part of the Gauss sum `i√7`, the modulus of the Paley tournament's skew eigenvalue
`±i√7`, half the Ramanujan expander bound `2√7`, and the square root of the discriminant of `Q(√−7)`. I had
been treating these as separate arrivals of the same constant across sessions; putting them in one line makes
plain that the arithmetic (Gauss, field), the combinatorics (Paley tournament), and the spectral geometry
(expansion) are three readings of one atom, and that atom is the ramified prime at the hard side, `p ≡ 3 mod
4`, where the Gauss sum turns imaginary and the Borsuk-Ulam degree turns odd. Three of the four faces are
imaginary (ι-odd); the fourth, expansion, is real. That is the certificate half of the apex.

The complementary half is `π`, and it hid inside the owner's most obscure prompt. The forum post about
Pochhammer-Chree waves in auxetic rods was mostly decoration — floral guardrails and placeholder hashes — but
the real hook was the word Pochhammer, and the Pochhammer symbol is exactly the fiber fraction I have used all
along: `f(n) = (1/2)_{n-2}/(n-2)! ~ 1/√(πn)`. So the apex has two halves that are genuinely different kinds of
number: a `√p` half that is algebraic, ramified, ι-odd, and lives in the *certificate*; and a `1/√(πn)` half
that is transcendental, smooth, ι-even, and lives in the *measure*. The covering-min value itself, `n/Φ6`, is
rational and contains neither — because it is the *answer*, while `√p` and `π` live in the certificate and the
measure that surround it. Naming this split is the kind of understanding that does not prove anything yet
organizes everything: whenever a quantity appears, I can now ask which half it belongs to before asking what it
equals.

The game framing paid a sharper debt. Writing the covering-min as an adversarial facility-location game —
defender picks a time to be lonely, adversary picks the speeds to cover the circle — let me import the
potential-function reflex from algorithmic game theory, and the potential function is the total danger
coverage, `2(n-1)M`. The instructive number is that at the covering-min this potential is about two, not one.
The runners cover the circle *twice over*. So the naive packing argument — total coverage below one forces an
empty gap — gives only the floor `1/(2(n-1))`, and the entire difficulty of the Lonely Runner is the factor of
two between that floor and `1/n`. And that factor of two is not a counting deficit; it is the *overlap
structure*, the discrepancy of the arithmetic-progression cloud with its three-gap `{1,n,2n}`. This is why the
Fourier and sum-of-squares methods stall: they are potential arguments, and potential arguments see coverage,
not overlap. The game framing did not close the gap, but it located it precisely — the crux is a discrepancy
statement about one specific cloud — and a located crux is worth more than a vague hard direction. Koksma-
Hlawka is the right language for it, reading the lonely measure as the discrepancy of the runner sample.

The boldest thread is the one I am least sure of, and I tried to be honest about that. Two sessions ago I
proved the flip-rank certificate is anti-LTC: the spectrum, being reflection-symmetric, cannot locally test
the self-complementary symmetry that carries the excess. This session the owner turned that around — use the
apex's *own* group. The Paley obstruction's automorphism group, the Frobenius `21`, sits inside PSL(2,7), the
Hurwitz group of order `168`, the automorphism group of both the Klein quartic and the Fano plane. And PSL(2,7)
is exactly the kind of group from which the recent locally-testable-code breakthrough builds its left-right
Cayley complexes. So I built that substrate: an expander square complex on the apex group, and checked that its
local links do detect a global defect. The proposal — that encoding the SC/`|Aut|` certificate as a co-cycle on
this complex makes the anti-LTC certificate locally testable via the very symmetry the spectrum missed — is a
structural direction, not a theorem, and I said so. But it is the right kind of direction: it takes the obstacle
(the certificate is invisible to global linear functionals) and answers it with the obstacle's own group
(whose local structure is exactly that symmetry). If the anti-LTC becomes an LTC, it will be because the atom
was made to test itself.

The keeper across all of it is the two-halves discipline. An apex in this kind of problem is not one hard
number to crack; it is a meeting of an algebraic certificate side and an analytic measure side, joined by a
rational answer and watched over by a finite symmetry group. `√p` and `1/√(πn)` are the two sides; `n/Φ6` is
the answer; PSL(2,7) is the group. When the next hard object appears, the first move is not to attack it but to
find its four faces, split it into certificate and measure, and ask which group is its own.
