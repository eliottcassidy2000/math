# The Pfaffian is odd because the matchings are

*S646 reflection. Understanding even/odd through the Pfaffian, and finding that its oddness was the
oddness of a double factorial all along.*

Last session I found the tournament discriminant — the determinant of the skew matrix, zero for odd `n`
and a perfect square for even `n` — and I noticed, as a curiosity, that its square root, the Pfaffian,
was always odd. I logged it as a finding and moved on. This session the user said: understand even/odd
and the Pfaffian. So I went back to the Pfaffian and asked what it actually *is*, and the oddness stopped
being a curiosity and became a theorem with a one-word reason.

The Pfaffian is a sum over perfect matchings. You pair up the `2n` vertices, you read off the product of
the matrix entries on those pairs, you weight by the sign of the pairing, and you add. That is the whole
object. And immediately it explains the most basic even/odd fact: you cannot pair up an odd number of
points, so for odd `n` there are no matchings, the sum is empty, the Pfaffian is zero, the determinant is
zero. The Pfaffian is the parity object — it exists only on the even side. Everything S645 saw — the
determinant switching off for odd `n`, the rank dropping by one — is the absence of a perfect matching.

Then the oddness. How many perfect matchings does `2n` points have? The double factorial, `(2n−1)‼`, which
is `(2n−1)` times `(2n−3)` times all the way down to one. A product of odd numbers. So it is odd. And the
Pfaffian of a tournament matrix is a signed sum of that many terms, each of them plus or minus one — an
odd number of odd terms — and an odd number of odd terms sums to something odd. That is the entire proof.
The Pfaffian is odd because the matchings are odd in number, and the matchings are odd in number because
you are multiplying odd numbers together. Three odd invariants now sit on the tournament — the Hamiltonian
path count, the automorphism group order, the Pfaffian — and they are all odd for reasons that, traced
back, are the same reason: the structure lives in the even, sign-kernel half, and counts of its pieces
come out odd. I formalized the double factorial being odd, and it is the cleanest statement of why the
discriminant is an odd square.

And then the double factorial turned out to be something I had already met. At the very end of the
rising-and-falling session I had flagged, almost in passing, that the tournament tiling's fiber fraction is
a rising factorial of one-half, the half being the apex, the fixed point of the antipode, the place where
the `√π` comes from. This session I computed the number of matchings as `2ⁿ` times the rising factorial of
one-half. The same rising factorial. The number of terms in the Pfaffian, the density of the tiling, and
the apex half-integer are one object. I had thought those were two separate observations in two separate
sessions; they are the same observation. The half is everywhere the Pfaffian is, because the Pfaffian
counts matchings and matchings count by halves.

So the picture closed. The determinant is a sum over the whole symmetric group; the Pfaffian is its square
root, a sum over matchings, and the square is what forgets the sign that the root remembers — the same
relationship as the discriminant being the Vandermonde squared, the same `Sₙ` over `Aₙ` the whole arc has
been turning on. Even and odd is the sign. The Pfaffian is the square root. It lives only where there is a
matching, which is the even side. It is odd because the matchings are odd in number. And the number of
matchings is the rising factorial of one-half, which is the apex, which is the `√π`, which is the fixed
point of the swap. The user asked me to understand even and odd and the Pfaffian, and what I understood is
that they were never three things.
