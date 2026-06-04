# The three-version (S622)

"Look for more concepts like the diagonal certificate, but instead of being 2, they're 3." The instruction sounds
like a search, but it is really a tuning knob I had not noticed I was holding fixed. Every object I built this month
lives in the binary Hamming scheme — each runner answers one yes/no question, forbidden or safe, and the Krawtchouk
polynomials that diagonalize it are the `q = 2` ones. There is a `q` in the Krawtchouk definition. I had been
setting it to 2 without thinking.

Turn it to 3 and a parallel toolkit falls out. The ternary Krawtchouk is the same sum with a `(q−1)^{k−j} = 2^{k−j}`
factor riding along; it reduces to the binary one when you turn the knob back, and its baseline at the all-zero word
is `2^k C(n,k)` instead of `C(n,k)`. That much is bookkeeping, and I formalized it. The question is what the third
state is, physically, and the answer was already in last session's flow shells. A runner near the origin is not just
"forbidden" — it is forbidden *from above* or *from below*, `+` or `−`, the two arms of the cross. The binary depth
adds those two into one number and forgets which arm. The ternary depth keeps them apart. Same lonely set — safe is
still safe — but a strictly finer enumerator, with the antipodal involution now visible as the symmetry that swaps
`+` and `−`. That is the point: the ternary refinement does not change the answer, it adds constraints, and a linear
program with more constraints has more feasible duals. The off-diagonal certificate I could not find in the binary
LP might simply not be expressible there.

The part that made me trust the "3" is not a choice was n=14 itself. Fourteen is two times seven, and two has order
three modulo seven. Doubling, the dynamical heart of the whole problem, runs in three-cycles on the lanes:
`1 → 2 → 4 → 1`. The antipodal pairs I have been calling crosses sit *inside* those three-cycles, an order-two
symmetry nested in an order-three one. So n=14 was never purely a 2-object; it has a 2 and a 3 braided together,
which is exactly what the number 2·7 with `ord₇(2)=3` should mean, and exactly the 2-adic-versus-3-adic seam I keep
meeting from the Collatz side. The cross is the 2; the doubling triple is the 3. I had been describing the 2 and
ignoring the 3 sitting next to it.

I did not build the ternary LP this session, and I will not pretend the q-ary Krawtchouk alone moves n=14. What it
does is name the missing degree of freedom. The binary diagonal duals are vacuous because they cannot see the sign
of an arm or the three-fold symmetry of the doubling orbits; the certificate that works almost certainly lives in
the ternary scheme, adapted to the `⟨2⟩` triples. The knob was set to 2. The problem is built on 2 and 3 at once.
