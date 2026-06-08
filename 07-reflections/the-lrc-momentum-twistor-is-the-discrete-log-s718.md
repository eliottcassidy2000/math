# The LRC momentum twistor is the discrete logarithm (S718)

Last session found the LRC's dual conformal symmetry: the multiplier group `(Z/m)*` acting on residues,
with inversion as the special-conformal generator, a symmetry hidden in spacetime and manifest only in the
arithmetic dual. The follow-up asked for a momentum twistor — the change of variables that, in scattering
amplitudes, makes dual conformal symmetry *linear* and trivial to act with. The answer is forced by one
observation: the LRC dual conformal action is multiplicative, `v -> a v`, and the thing that linearizes a
multiplicative action is a logarithm. The LRC momentum twistor is the discrete logarithm.

Pick a primitive root `g` modulo the shell `m` and send each residue to its discrete log
`ell(v) = log_g v` in `Z/phi(m)`. Three checks, all exact across every shell I tried including the
ramified `27`: the multiplier `v -> a v` becomes the translation `ell -> ell + ell(a)`; the inversion
`v -> v^{-1}` becomes negation `ell -> -ell`; and `ell(-1) = phi(m)/2`. The hidden symmetry that scrambled
the worldlines beyond recognition in spacetime is, in the twistor coordinate, just a shift and a flip.
This is exactly what a momentum twistor is supposed to do — Hodges' twistors turn the nonlinear dual
conformal group into a linear `GL(4)` action; the discrete log turns the LRC's multiplicative dual
conformal group into the additive group `Z/phi(m)`.

Once the symmetry is linear the hard object simplifies. A multiplier is bad for runner `i` exactly when it
pushes that runner onto the central band, `a v_i = ±1`, which in logs reads `ell(a) in {-ell(v_i),
c - ell(v_i)}` with `c = phi/2`. So the entire set of bad multipliers is `B = (-L) cup (c - L)`, the union
of two translates of the negated log-set `-L`. A dodge exists if and only if these two translates fail to
cover `Z/phi(m)`. The multiplicative dodge problem — find a rotation that keeps every runner off the band
— has become a problem about whether two translates of a set cover a cyclic group. That is the whole point
of a twistor: the question that was tangled in the original variables is a clean covering question in the
new ones. And it is exact: the criterion `|B| < phi(m)` matched a direct multiplier search six hundred out
of six hundred times at every shell.

The transversal core, the genuinely hard frontier, falls out as a tiling condition. The two translates
cover `Z/phi(m)` with no overlap precisely when `L` contains one element of each pair `{x, x+c}` — a
half-system, a transversal of the order-two subgroup `{0, c}`. In residues this says `L` contains no pair
`{v, -v}`, which is exactly THM-420's `±`-transversal; the twistor recovers the theorem and upgrades it to
the statement that the core is where `{-L, c-L}` *tiles* the group. This also reconnects to THM-415, whose
signed-LRC half-system is the same object seen through the same log coordinate. The twistor is the common
chart in which the dual conformal symmetry (S717), the transversal core (THM-420), and the signed-LRC
half-system (THM-415) are one picture.

The most satisfying thing the twistor does is explain, by a size count, why some `n` are hard and others
ramified. The half-system characterization requires `|L| = n - 1` to equal the number of pairs,
`phi(2n-1)/2`. That is the single equation `phi(2n-1) = 2n - 2`, and it holds if and only if `2n-1` is
prime. At a prime shell the sizes match, the transversal core exists as a unit half-system, and THM-420 is
the clean statement — this is where the open frontier `n = 15, 19, 21, 22` lives. At a ramified shell like
`27` for `n = 14`, `phi(27) = 18` is strictly less than `26`, the sizes cannot match, no unit half-system
exists, and the unit residues *always* dodge. The difficulty has not vanished; it has moved to the
residues divisible by three, which are not units and therefore not on the twistor at all. The momentum
twistor literally separates `n = 14` into its two heads: a unit part it linearizes and dispatches, and a
ramified off-twistor part that is the real obstruction — the same fiber the divisor and 7-clock work
(S643, S710) isolated by other means. Three roads, one place.

It does not close the conjecture. At a prime shell the half-system core is still there, still open, and the
twistor only tells you it is a tiling of `Z/phi(m)` rather than trivializing it. But it is the right
coordinate. The dual conformal symmetry is manifest and linear; the dodge is additive covering; the prime
vs ramified split is a totient size match; and the hardness of `n = 14` is exactly the part of the residue
set that the logarithm cannot reach. The amplitudes borrowed us a word, and the word turned out to name the
discrete log.
