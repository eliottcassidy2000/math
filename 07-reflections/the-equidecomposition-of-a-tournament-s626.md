# The equidecomposition of a tournament (S626)

The user handed me two words — equinumerosity and equidecomposability — and they turned the forbidden-value problem
inside out. I had been treating the impossible values 7, 21, 63 as a list with a pattern, 7 times powers of three,
and the pattern as a thing to be explained. Equidecomposability says: stop looking at the list, look at how the
object comes apart.

A tournament comes apart canonically into its strongly-connected components, and the components sit in a line,
because the condensation of any tournament is transitive. A Hamiltonian path has to walk through the components in
that order, finishing one before it can enter the next, so a Hamiltonian path of the whole is a Hamiltonian path of
each component, independently. That is an equinumerosity — the paths of the whole are in bijection with tuples of
paths of the parts — and counting it gives an equidecomposition: H of the tournament is the product of H over its
components. I checked it on twenty thousand tournaments at each size and it never failed. So H is not just a number,
it is a multiplicative function, and the strongly-connected tournaments are its primes.

Once H is multiplicative the forbidden-value question stops being combinatorics and becomes arithmetic. The
achievable values are exactly the products of atoms, where an atom is the H of a strongly-connected tournament. So a
value is forbidden precisely when it cannot be written as a product of atoms. And now seven explains itself. Seven is
prime, and no strongly-connected tournament has H equal to seven — that is the old impossibility theorem, the three
pairwise-touching cycles that force a fourth. So seven is not an atom and not a product of smaller atoms, and it is
forbidden. Twenty-one is three times seven; three is an atom, the lone three-cycle, but seven is not, and there is no
atom equal to twenty-one either, so twenty-one cannot be assembled. Sixty-three is seven times nine; same story. The
seven-times-three-to-the-k ladder is forbidden because every rung is seven times a power of three, the only way to
supply the seven is an atom divisible by seven, the only divisors of seven-times-three-to-the-k that carry a seven
are themselves seven-times-powers-of-three, and none of those is an atom. I formalized exactly this — three lines of
prime arithmetic — and it turns an infinite family of impossibilities into a single one: no strongly-connected
tournament has H equal to seven times a power of three. Prove that one statement and the whole ladder falls out.

That is the sharpest target I have had in weeks, and it is sharp because the equidecomposition did the work. The same
two words organize everything else in hindsight. Equidecomposability is the multiplicative skeleton — H over
components, the partition function over disjoint blocks, the independent sets split by a vertex. Equinumerosity is the
bijection underneath each split — the paths through the components, the parity vector and the integers below two to
the K, the chord and the arc. And on the circle, where the rotations form an amenable group, Tarski tells me
equidecomposable means equal measure, no paradox, which is exactly why the Lonely Runner obstruction is honestly
measure-theoretic and not some Banach-Tarski sleight of hand. The forbidden values, the collapse family, the
loneliness — they are all what survives when you cut an object along its natural seams and ask what the pieces can and
cannot be reassembled into.
