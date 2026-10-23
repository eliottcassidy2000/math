# Independent proof and consumer audit

October 23, 2026. **PROVED:** the elementary reconstruction, transport, and
component-necklace statements with the scope below. **FINITE-EXACT:** the
declared bounded control universes. This review does not promote scalar
bucket data to an all-order metagraph classification.

## Root curl, its exact ambiguity, and the decoder

The difference of two arc arrays with equal signed root curls is a gradient
`d_ij = u_i - u_j`. Its entries belong to `{-2,0,2}`. Consequently the
potential has at most two distinct levels: three levels would force an
extreme difference of magnitude at least 4. A nonconstant potential changes
exactly one uniformly directed cut. Such a cut is a prefix of the unique
ordered strong-component decomposition. Reversing it cyclically rotates
that decomposition and leaves each labelled component internally intact.
Conversely, every such rotation has unchanged curl. These implications give
the claimed exact fibre, rather than merely a necessary condition.

Distinct members of this fibre have distinct scores at every fixed root.
Compare any pair: the root lies on one side of a nonempty proper directed
cut, and reversing it changes its score by the nonzero size of the opposite
side. The argument applies to every pair of fibre members. Thus one root
score restores the lost cut without requiring all vertex scores.

The decoder is complete for arbitrary submitted tensors, including rejection.
Writing `x_i = a_0i`, the equality `a_ij = kappa_ij - x_i + x_j` gives the
four exact constraints stated in the proof. Along an arc of the sign
tournament, a feasible x cannot increase. It is therefore constant on every
strong component and selects a prefix of its component chain. Testing all
those cuts and the literal arc alphabet tests every possible realization.
This also explains why `(3,3,3)` fails: one shared root-star coordinate would
have to be both +1 and -1. It is not enough that each coordinate separately
belongs to the permitted curl alphabet.

The code's original decoder accepted Boolean/float aliases through Python
equality. This audit requested literal-integer guards for the order, curl
coordinates, and root score; those guards and hostile controls were added.
The mathematical signature is indexed and signed. Sorting its entries or
discarding their signs would destroy the proved injectivity.

## Ports, orbits, and the small-order boundary

On a fixed directed path, the displayed recurrence reconstructs the root
star from consecutive root curls, so no score sidecar is needed. The update
`a' = -a + 2p` gives the claimed curl and score transport by linearity. An
actual directed Hamiltonian path is a separate source guard; it cannot be
inferred from a score list.

For two paths on one tournament, equality of their normalized adjacency
arrays is equivalent to an automorphism taking the first ordered path to
the second. One direction follows by applying the automorphism; the other
uses the vertex permutation defined by corresponding path positions. The
automorphism action on ordered spanning paths is free, giving the inherited
`H/|Aut|` port count. Hence the test matches literal path orbits, not only
their scores or unsigned triangle types.

The all-order pairing claim needs `n >= 3`: only then is the free-edge set
nonempty and the endpoint-closing edge free. The audit caught the omitted
boundary; the proof now states it explicitly. At order 2 the partner is the
identity, although the reconstruction and affine transport remain valid.
Blue preservation follows because path reversal is an anti-symmetry of the
path flow itself, so it commutes with this affine update in the required
signed form.

The uniqueness conclusion is correctly conditional on a **complete
authenticated port inventory**. It does not say that margins, node sizes,
parity, electrical quantities, or any bounded motif profile determine that
inventory. It supplies an exact rejection test for a false proposed port
attachment. It supplies neither a smaller enumeration algorithm nor a proof
of black-core connectivity.

## Component cuts and even/odd buckets

An independent permutation calculation on the three labelled rotations of
`C3 -> K1 -> K1` gave `(H, |Aut|, number of anti-automorphisms)` equal to
`(3,3,0), (3,3,0), (3,3,3)`. All three have the identical root-curl vector
`(3,1,1,-1,-1,1)`. Thus the same labelled curl fibre really reaches both a
self-converse class with one tiling and a non-self-converse pair with two
merged tilings. This is a concrete obstruction to identifying curl fibres
with converse-merged nodes.

For the general necklace, unique ordered SCC decomposition identifies
unlabelled classes with the distinct linear cuts of its word of strong
component types. If its least period is d and converse-reversal is rotation
by a, then conjugating rotation by that involution reverses its direction.
A cut t is self-converse exactly when `2t = a (mod d)`. This proves the
one / zero-or-two parity alternatives and the merged-node count `(d+s)/2`.
If converse-reversal lies in a different necklace, all d classes instead
pair with that other necklace. The common product `H/|Aut|` is unchanged;
its oddness uses the inherited odd Hamiltonian-path and automorphism-order
facts. These are statements about an SCC necklace, not about all blue-line
multiplicities.

The initial finite-test description used `K1 -> C3` as a strong-component
letter. It is not strong. The audit flagged this type error: the formal
three-letter involution calculation is valid, but such letters do not give
the asserted SCC realization. A valid concrete replacement supplied during
review is the strong six-tournament with edge-mask 18 in this packet's
convention and its converse. It contains the directed Hamilton cycle
`0 -> 5 -> 4 -> 3 -> 2 -> 1 -> 0`; checking all 720 permutations gives no
anti-automorphism. Its ordered score vector is `(2,1,1,3,4,4)`. This supplies
a genuine pair of distinct converse component types for the formal test.
The final proof and checker use this corrected strong letter and explicitly
verify both strong connectivity and absence of an anti-automorphism.

## Verification and literature scope

The independent optimized-Python replay verifies every labelled tournament
through order 6, every fixed-path port through order 7, the explicit
five-vertex witness, and the finite necklace controls. The proof above does
not extrapolate those finite checks to larger orders; it gives the missing
all-order arguments separately.

The final independent replay passes **388733 explicit checks** under
optimized Python. Audited `check.py` SHA256 (also its LF-normalized hash):
`e008c8b0ad044ac866d26d5d65e0c5c9bff6387701f64999c4ddf88891a544d6`.
The producer's normal run and this optimized run agree on the declared
finite universes, 1374 necklace cases, and the corrected order-six letter.

The formal-series identity `C = log T` also follows directly from the proved
cyclic sequence of labelled strong components. Its equality with the
connected-labelled-graph series is an equinumerosity, not by itself a
canonical graph/tournament bijection. The cited embedding-based Bernardi
construction is appropriately described as prior art and is not implemented
by this checker. This audit does not independently certify that external
implementation or an observer-free/unlabelled version of it.
