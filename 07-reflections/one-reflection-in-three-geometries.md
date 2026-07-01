# One reflection, in three geometries

*klein-2026-07-01-S81. A reflection on HYP-3814 — the Cayley transform glues the staircase to the
circle, complement is a single involution seen in three coordinate systems, and the odd/even duality
is that involution's fixed-point parity.*

For most of this project the tournament and the runner have lived on opposite shores. The tournament is
combinatorial: a tiling of a staircase, a point in a hypercube, a thing you canonicalize by trying every
relabeling. The runner is analytic: a point on a circle, a measure with Verblunsky coefficients, a phase
`nv mod Phi6`. We have crossed between them by analogy — "the complement is like the antipode," "the
self-complementary class is like the palindromic measure" — and the analogies have been productive. But an
analogy is a rope bridge; you can carry ideas across, but you always feel the gap underneath. This session,
prompted to take the gluing seriously, I found that the bridge is not a rope but a diffeomorphism, and it has
a name older than any of this work: the Cayley transform.

The move is almost too simple to have waited this long. A tournament has a skew-symmetric sign matrix
`S = A - A^T`, entries `+-1` off the diagonal. Skew matrices are the Lie algebra; the Cayley transform
`U = (I - S)(I + S)^{-1}` is the classical retraction to the group, and it sends the imaginary spectrum of
`S` to the unit circle. So every tournament, with no choices made, becomes a set of points on the circle —
the eigenvalues of `U`, the runners. The staircase and the circle are not two shores; they are one object in
two coordinate systems, and the Cayley transform is the change of coordinates. I checked that `U` is
orthogonal and its eigenvalues lie on `|z| = 1` for every tournament up to `n = 5`, but the check was a
formality: it is forced by `S` being skew. The content is that the dictionary exists at all.

What the dictionary makes rigid is the complement. On the staircase, complement is the grid reflection
`sigma` (opus proved this years of sessions ago). On the matrix, complement is `A -> A^T`, which is
`S -> -S`, a negation. On the circle, `U(-S) = U(S)^{-1}`, which for points on the unit circle is
`theta -> -theta`, a reflection across the real axis. Three descriptions that used to be three facts are now
one fact: `sigma`, `S -> -S`, and `theta -> -theta` are the same involution, transported by Cayley. When I
say "complement is a reflection," I no longer mean it three times loosely; I mean it once, precisely, and the
three appearances are shadows. The rope bridge became a mirror, and a mirror is the same on both sides by
definition.

The payoff I did not expect was the odd/even duality falling out as the mirror's fixed-point structure.
A self-complementary class is exactly a class fixed by the reflection. Rédei's theorem — the number of
Hamiltonian paths `H` is odd — is the statement that the fixed fiber has odd size. The non-self-complementary
classes come in mirror pairs, so their merged fiber is `2H`, even. That is the entire SC-odd / NS-even
parity the owner has been circling for a dozen sessions, and here it is not a computation about merged
metagraph nodes but a sentence about a reflection: *fixed points are odd (Rédei), moved points pair up (even).*
The same sentence is true on the staircase (odd `sigma`-fixed fiber), on the circle (odd palindromic mass),
and on the spectrum (`U` conjugate to `U^{-1}`). One fact, three geometries — which is exactly what the
prompt asked me to see, and I believe it because the geometries are provably the same geometry.

Two corollaries came for free, and one of them is a small confession. The skew spectrum is `+-`-symmetric,
so `spec(-S) = spec(S)` always: the reflection acts *trivially* on the spectrum. That is precisely why the
skew-spectrum was a weak invariant back in S72 — it cannot see the complement, because the complement is a
reflection and the spectrum is reflection-blind. A weakness I had merely observed is now explained by the
same structure that explains everything else. And the `n`-parity: an odd-size skew matrix is singular, so
`U` has a `+1` eigenvalue — a runner pinned at angle zero — exactly when `n` is odd. The observer's fixed
point, in the tournament's own spectrum.

The confession is about the Paley heptagon. I had carried a memory that its Cayley transform gives the
seventh roots of unity, `U^7 = I`. It does not. The Paley skew circulant has eigenvalues `+-i*sqrt(p)` — the
Gauss sum — so the Cayley eigenvalues sit at `cos(theta) = -(p-1)/(p+1)`, which is `-3/4` at `p = 7`, an
irrational angle, no root of unity in sight. The roots of unity were real, but they live on a *different*
circle: the vertex loop, where a circulant tournament's vertices sit at `n`-th roots of unity. Two circles,
and I had glued the wrong pair. Checking HYP-3802 showed the file itself was careful; the error was only in
my recollection. But the correction is the better fact: the Cayley angle of a Paley tournament *is* the Gauss
sum, `cos(theta) = -(p-1)/(p+1)`, exact and clean, and it says the arithmetic of the tournament is written on
the circle in a coordinate I had not been reading. The mistake was worth making because its repair was
sharper than the thing I misremembered — which is the recurring shape of this project, and the reason the
reflections directory exists.

The keeper is methodological, and it is about when to stop reaching for analogy. For a long time the right
move between the tournament and the runner was to notice a resemblance and carry it across by hand. That was
correct while the bridge was a rope. But when two objects keep resembling each other in fact after fact, the
resemblance is usually a functor waiting to be named, and the work shifts from carrying ideas across to
finding the map that makes the carrying automatic. Here the map was the Cayley transform, sitting in plain
sight in the skew matrix the whole time. The lesson for the next accumulation of analogies: before proving
the twentieth parallel by hand, ask whether the nineteen already proved are telling you the two sides are the
same side, and go looking for the mirror.
