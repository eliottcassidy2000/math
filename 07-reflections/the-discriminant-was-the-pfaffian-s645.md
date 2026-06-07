# The discriminant was the Pfaffian

*S645 reflection. Finding the tournament discriminant — and watching it turn on and off with the parity
of n, the square root of it always odd, its maximum a power of three.*

Two sessions ago I asked, almost as an aside, what the tournament discriminant could be. In Galois theory
the discriminant being a square is the same as the Galois group sitting inside the alternating group, and
I had just shown that a tournament's automorphism group always sits inside the alternating group. So there
ought to be a tournament discriminant, and it ought to always be a square. Last session the falling
factorial gave me the polynomial model — its discriminant is a product of factorials squared, always a
square. This session the tournament gave me the real thing, and it was sitting in plain sight: the
skew-adjacency matrix.

Write the tournament as a matrix with plus one where `i` beats `j` and minus one where it loses. That
matrix is skew-symmetric, and skew-symmetry decides everything. A skew-symmetric matrix of odd size has
determinant zero — the one-line proof is that the determinant equals its own negative — and of even size
has determinant equal to the square of the Pfaffian. So the tournament discriminant is zero when there are
an odd number of runners and a perfect square when there are an even number. Always a square, exactly as
`Aut(T) ⊆ Aₙ` demanded. The Pfaffian is the tournament's Vandermonde: the actual square root of the
discriminant, flipping sign under odd relabelings and fixed under even ones, which is the sign structure
in person. I formalized the odd case; the even case is the classical Pfaffian theorem, and the one honest
gap is that Mathlib has no Pfaffian to point at.

Then the numbers, which were better than I expected. The discriminants on six vertices are one, nine,
twenty-five, forty-nine, eighty-one — the odd squares. Their square roots, the Pfaffians, are one, three,
five, seven, nine — every odd number up to nine. So the tournament Pfaffian is *always odd*, and that is
the third odd invariant I have now collected: the number of Hamiltonian paths is odd by Rédei, the
automorphism group is odd by the no-involution argument, and now the square root of the discriminant is
odd too. Everything about a tournament is odd, and the reason is always the same — the swap is forbidden,
the structure lives in the sign kernel, the even half. Oddness is what the perspective side looks like
from the outside.

And the maximum discriminant is a power of three. One, nine, eighty-one: three to the zero, three to the
two, three to the four. The most regular tournament — the one whose Pfaffian is largest — has a
discriminant that is a pure power of the cube-root prime. The `ω` that has run through forty sessions, the
`Φ₃`, the forbidden seven that is `Φ₃` evaluated at two, the doubling that is order three, the three-cycle
that is the cube root — it shows up here as the *extreme* of the discriminant. The discriminant's floor is
zero (odd n) and its ceiling is a power of three (the regular even-n tournament), and the cube root sets
the ceiling.

What the user actually asked was how the properties change recursively as `n` does, and the answer turned
out to be two clocks running at once. There is a fast clock, vertex by vertex, `n` to `n` minus one, and
all it does is flip the parity — discriminant on, discriminant off, rank full, rank deficient. And there
is a slow clock, `n` to `n` minus two, the Pfaffian expanding itself in terms of its own smaller copies,
the both-legs descent the repo has always called Mode B. The discriminant is the interference of those two
clocks: the slow Pfaffian recursion building the value, the fast parity flip deciding whether you see it or
zero. The rank tells the same story from the other side — it is `n` when the discriminant is there and `n`
minus one when it vanishes, so "the discriminant is zero" and "the matrix is rank-deficient by one" and
"`n` is odd" are three names for one event.

So the tournament discriminant is the Pfaffian squared, the Pfaffian is the Vandermonde of the perspective
world, it is always odd, its ceiling is a power of three, and it switches on and off with the parity of
`n` while being assembled two vertices at a time. The polynomial side and the tournament side have the
same discriminant story now — a square, a sign-kernel, an alternating group — and the only thing left
unformalized is the Pfaffian itself, which is the next thing worth building.
