# What we keep track of is the transfer spectrum (S720)

Six sessions arrived at one table. The lonely runner's dual conformal symmetry, the discrete-log twistor,
the unit-distance angle twistor, the Pfaffian even/odd seam, the autocorrelation operator, and the
`A+B+C-D-E-F+G` tournament recursion looked like six investigations. They are one object looked at from
six sides, and the object is the spectrum of a transfer operator. The instruction was to reframe what we
keep track of and what it means; the answer is that we have been keeping track of the elementary symmetric
functions of a transfer spectrum all along, and they split cleanly into what is frozen and what runs.

Start with the recursion the owner keeps returning to. `A+B+C-D-E-F+G` is inclusion-exclusion over the
three corner subtriangles, characteristic polynomial `(x-1)^3`, the third finite difference. Read its
eigenvalues — all three are `1` — and read their symmetric functions: the sum is three, the product is
one, the middle pairwise-sum is three. Now look at the other order-three recursion the project found, the
base-path Hamiltonian-path count of THM-337, characteristic polynomial `x^3 - 3x^2 - x - 1`, dominant root
about three point three eight. Its eigenvalue sum is again three, its product is again one, and only the
middle symmetric function has changed, from three to minus one. Two recursions, identical except in one
coordinate. So they are not two recursions; they are two points of a single one-parameter family,
`x^3 - 3x^2 + a x - 1`, with the sum of eigenvalues frozen at three, the product frozen at one, and the
middle coefficient `a` free to move.

The three frozen and running pieces have meanings, and the meanings are the whole reframe. The sum of
eigenvalues, three, is the number of corners of the triangle — the geometry, the simplex dimension. It is
frozen because the shape does not change; a tetrahedral recursion would have sum four, a `d`-simplex sum
`d+1`, and the additive crystalline point would be `(x-1)^{d+1}` with binomial middle coefficients. The
product of eigenvalues, one, is the determinant of the transfer matrix, and it is `±1` — unimodular —
which is the recursion-level face of the Pfaffian. S713 showed `det(I+2A) = Pf^2` with the Pfaffian forced
odd; here that parity reappears as the frozen unit product of the transfer spectrum. The geometry lives in
the top symmetric function, the parity in the bottom, and both are frozen. What runs is the middle: the
single coordinate `a`, which I will call the temperature.

The temperature is exactly the additive-versus-multiplicative axis the project keeps meeting. At `a = 3`
all three eigenvalues collapse to one — the crystalline point — and growth is polynomial; this is the
additive regime, where the invariant is a sum over tiles and the inclusion-exclusion is perfectly tuned.
This is `A+B+C-D-E-F+G`. Cool it no further; it is already frozen solid. Any `a` below three cracks the
triple root apart: one eigenvalue climbs above one, the others sink (their product with it staying one),
and growth turns exponential. At `a = -1` you are at the Hamiltonian-path count, dominant root three point
three eight, fully molten. The crystalline point is unique — sum three and product one and every
eigenvalue on the unit circle force them all to be exactly one — so additivity is a measure-zero miracle,
the single temperature at which the inclusion-exclusion closes; everything else is hot. Additive
invariants are rare and cold; multiplicative invariants, like the number of Hamiltonian paths, are the
generic hot case.

And now the twistors fall into place as the cooling operation. The discrete log that linearized the lonely
runner's multiplier symmetry, and the angle that linearized the unit distance rotation symmetry, both do
the same thing: they take a multiplicative structure and make it additive. At the level of the spectrum
this is literally lowering the temperature. A geometric sequence `g^n` has a single eigenvalue `g`, hot and
off the unit circle; its discrete log is the arithmetic sequence `n`, characteristic polynomial `(x-1)^2`,
eigenvalue one, cold. The twistor moves a structure from the hot end of the ladder toward the crystalline
end. So the relation the owner asked for — how do the twistors relate to `A+B+C-D-E-F+G` — is that
`A+B+C-D-E-F+G` is the cold crystalline limit and the twistors are the maps that cool toward it. The
inclusion-exclusion is where you arrive when you have linearized everything; the discrete log and the angle
are how you get there.

This rewrites the ledger. We have been recording the number of Hamiltonian paths, the Pfaffian, the lonely
gap, the unit-distance count, the autocorrelation, the recursion coefficients, the edge directions, the
half-systems, as if they were separate quantities. They are coordinates on one thing. The autocorrelation
operator `MM*` is the convolution that the transforms diagonalize. The transforms are two: Mobius
inversion on the boolean lattice, which is `A+B+C-D-E-F+G` and handles the tiling geometry, and Fourier on
the cyclic group, which is the twistor and handles the dual conformal symmetry. The spectrum they expose
has its top symmetric function fixed by the geometry (the corners), its bottom fixed by the parity (the
Pfaffian, unimodular), and its middle free — the temperature, the one real degree of freedom, the
additive-to-multiplicative dial. Hamiltonian paths are that dial turned to hot; the inclusion-exclusion is
it turned to cold; the twistor is the hand that turns it; the Pfaffian is the seal that the determinant
stays one while it turns.

What things mean, then, is this. Geometry is the count of corners and lives in the leading coefficients,
frozen. Parity is the determinant and lives in the constant term, frozen and equal to the Pfaffian. The
content of a problem — whether it grows polynomially or explodes, whether it is solvable additively or
only multiplicatively, whether the lonely runner dodges or the unit distances pack — is the temperature,
the middle of the spectrum, the single number that the inclusion-exclusion freezes and the twistor thaws.
The next thing to track is not another invariant. It is where on the temperature ladder a given problem
sits, and which transform moves it.
