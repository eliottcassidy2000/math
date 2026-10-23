# Independent audit: exact graph/tournament quotient coupling

PASS on the sharp all-size theorem and signed-moment equality classification,
subject to the explicit carrier restrictions below. No literature-priority
claim is inferred.

## Probability normalization and operators

The underlying sample space is the entire labelled edge-bit cube with uniform
measure. Graph and tournament relabelings are different measure-preserving
actions, the latter affine. Averaging each group gives its orthogonal conditional
expectation. On graph-invariant functions the operator P_G P_T P_G has transition
matrix K_ik=sum_j W_ij W_kj/(g_i t_j). Thus pi_i=g_i/2^N, not a uniform measure
on unlabelled classes. Direct checks confirm stochasticity, positive diagonal
and detailed balance. The displayed weighted Gram form proves positivity.

## Walsh orbit normalization

An ordinary graph orbit of edge subsets spans a monomial invariant space for
the twisted action. Its ordinary invariant subspace is one-dimensional. A
negative stabilizer sign kills the twisted invariant line; otherwise its
coordinates can be chosen as signs. Projection between the uniform orbit sum
and that signed orbit sum has coefficient mu(F) in both directions. Their
normalizations agree, so the resulting eigenvalue is exactly mu(F)^2, with no
missing orbit-size or n! factor. Relabeling can change mu's sign but not its
square. Signed orbit averages may vanish even when the stabilizer permits an
invariant line.

The independent checker constructs graph and tournament quotient classes
using a separate permutation-image implementation and tests EACH Walsh orbit
vector's exact eigenvector equation, not merely matching traces or floating
spectra, for n=2..5. Both complement/converse-paired and unpaired spaces pass.

## Sharp signed-moment theorem

The first-vertex recurrence has the correct sign: every neighbor smaller than
the first vertex contributes one minus sign. If no vertex covers all edges,
every deletion is nonempty and induction bounds the average by1/3. Otherwise
the graph is a star plus isolates, with alternating rank sum0 for odd degree
and magnitude1/(d+1) for even degree. This closes every case, including a single
edge and arbitrary disconnected graphs.

The equality classification is also sound. Relative orders of disjoint vertex
sets are independent, so moments multiply over components and two nontrivial
components cannot attain1/3. Connected nonstar equality forces every deletion
to be extremal. Their maximum degrees are at most2; nonzero moments force even
edge counts both before and after deletion, hence every original degree is
positive even. Deleting a neighbor then forces that degree to be2. The graph
is a cycle, and its deletion path is extremal only for C4. P3 and C4 attain the
bound, including the stated C4 sign under the chosen labels. Isolates preserve
these values.

Independent direct permutation averages verify the signed recurrence and
extremal shapes for all labelled graphs through n=5. The producer's integer
Walsh–Hadamard controls through n=7 are replayed in normal and optimized modes;
they support but do not replace the structural proof.

## Gap, pairing, and boundaries

Every nonconstant eigenvalue is at most1/9 and the P3 sector attains1/9 for
n>=3. Therefore the gap is exactly8/9. The constant eigenspace is one-dimensional,
so reversibility with positive stationary masses implies irreducibility; positive
diagonal then implies aperiodicity. At n=2 the ordinary chain has nonconstant
eigenvalue0, while the paired chain has a single state and trivial mean-zero
space. A gap of1 in that latter case is a convention, not an extra eigenvalue.

Complement acts by (-1)^|F| on Walsh characters. Pairing keeps even-edge sectors,
and the discarded odd-edge sectors already have mu=0. Consequently nonzero
spectra, including multiplicities, agree. The exact L2/variance contractions
and TV upper bound have the correct normalization. In particular the point-mass
factor sqrt(1/pi_i-1) prevents interpreting a fixed gap as uniform-in-n TV
mixing from arbitrarily rare states.

The result concerns the full labelled cube and its two orbit quotients. It does
not apply without a new proof to the fixed-Hamiltonian-path metagraph, regular
tournaments, a black-core walk, shared-port-preserving moves, or arbitrary
restricted inventories. It gives no faster exact tournament enumeration by
itself. Scalar spectra do not restore discarded native observables.

## Exact lost-observable addendum

The rank formula is correct because the nonzero singular values of the
cross-projection have squared values equal to K's positive eigenvalues.
Consequently tournament invariant dimension minus rank(K) counts invisible
linear observable directions, not missing tournament objects. The recorded
ranks1,2,3,10 and the paired/unpaired nullities were independently checked.
At order4 the cubic star observable is indeed4 on a source over C3, -4 on C3
over a sink, and0 on the other types. Its ordinary graph conditional mean is
zero on each class. A direct proof is that its four star-character coefficients
alternate +1,-1,+1,-1 in reference vertex order, and ordinary relabeling
averages the four star supports equally. It therefore demonstrates genuine
information loss despite the large coupling gap.

## Reproduction

`check_independent.py` is standard-library-only, self-contained and portable.
It verifies exact quotient eigenvectors through n=5 and direct signed moments
through n=5. Producer replays and immutable source pins are retained in this
packet. No source files or Git state were modified by the auditor.

## Attribution-only final delta

The added inheritance paragraph accurately credits THM-4721 section5/TheoremW
and THM-4531/TheoremG for the signed Walsh basis. Both source files were read
at pinned commit5258d52 and their SHA256 hashes verified. Automorphism-even
is correctly distinguished from even vertex degrees. The unsigned/signed
overlap and sharp coupling gap remain the present extension. Checker and
executable receipts are byte-identical to the earlier audited versions.
