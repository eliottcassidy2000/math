# Positive Pfaffian kernels and regular prime/composite composition

**PROVED elementary composition rule and recurrence / FINITE-EXACT
witness / CITED inherited modular decomposition.** This refines the
regular-card sector of tournament enumeration. It does not replace the
existing all-order Burnside formula or produce a new unrestricted A000568
term. Same-agent review; no literature-priority claim.

## 1. Recover the stronger existing tools first

The repository already has the
[prime cycle-index composition formula](../session86_relations_attachment_carriers_20261028/A000568_COMPOSITION.md),
the canonical prime/ordered decomposition and actual decorated
automorphism action in
[Session68](../session68_tournament_composition_classification_20261026/session68-classification/CLASSIFICATION.md),
and the [regular-parent child theorem](../regular_parent_symmetry_20261029/PROOFS.md).
We do not derive a weaker substitute for any of them.

The underused sidecar is the SIGNED Pfaffian cofactor vector. Its exact
root covariance was recovered in
[Session102's parity interfaces](../session102_tournament_classification_20261029/pfaffian/PARITY_INTERFACES.md),
with dependencies THM-174 (`01-canon/theorems/THM-174-pfaffian-determinant-identity.md`)
and THM-3729 (`01-canon/theorems/THM-3729-rooted-pfaffian-response-and-sign-root-deletion-average.md`).
We use only the elementary kernel identity, re-established below; no
unrelated path-count consequences are imported.

## 2. Exact weighted regularity test

Let Q be a tournament with skew adjacency K (K_ij=1 for i->j). Substitute
a nonempty tournament T_i of size w_i at each labelled vertex. Cross
edges follow Q. The resulting tournament is regular if and only if

    each T_i is regular, and K w=0.                           (1)

Indeed, the signed row sum at a vertex in T_i is its internal row sum
plus (Kw)_i. If all total row sums vanish, all internal row sums in T_i
are the same. Their sum is zero by skew symmetry, so they all vanish,
and then Kw=0. The converse is immediate. In particular all w_i are odd.

For even-order Q, K is nonsingular. For odd order m, it has rank m-1 and
kernel vector

    p_i=(-1)^i Pf(K with row and column i deleted),  i=0,...,m-1. (2)

Every even tournament Pfaffian is odd: modulo two it is the number of
perfect matchings, an odd double factorial. Thus it is nonzero. Even
determinants are squared nonzero Pfaffians; odd skew rank is at most
m-1, and an even principal submatrix attains that bound. Expanding a
Pfaffian with a repeated row/column gives Kp=0, proving (2) and its
one-dimensionality. Every coordinate in (2) is nonzero and odd.

Consequently Q admits a regular inflation iff the coordinates of p all
have the same sign. After changing the common sign and dividing by their
odd gcd, obtain the unique primitive positive odd integer vector u.
All possible child-size vectors are then EXACTLY

    w=t u,  t a positive odd integer.                         (3)

Primitivity forces an integral multiplier, and odd sizes force it odd.
Every vector (3) is realizable: every odd size admits a cyclic regular
tournament. The smallest regular inflation has order W_Q=sum u_i.
Its child labels, orientations and actual ports still matter for later
attachments; (3) determines size feasibility, not the whole graph.

Every automorphism of Q preserves u coordinatewise under its vertex
permutation. It sends the positive one-dimensional kernel ray to itself;
the primitive positive integral normalization rules out rescaling.

## 3. A weighted prime-action recurrence

Let R_n count regular tournament isomorphism classes, with R_1=1 and
R_n=0 for even n. Let P^+ be the ordinary isomorphism classes of prime
tournaments Q of odd order at least3 whose kernel has positive generator
u_Q. Retain their ACTUAL automorphism actions on vertices together with
u_Q. For n>1,

    R_n = sum_(Q in P^+, n=t W_Q, t odd)
             (1/|Aut Q|) sum_(sigma in Aut Q)
                 product_(cycles c of sigma) R_(t u_Q(c)).    (4)

Here u_Q(c) is the common coordinate on the cycle. This is a finite
recurrence: every child size is smaller than n and |Q|<=W_Q<=n.

Proof. A regular tournament with more than one vertex is strongly
connected. If it had a source strong component of size a and a sink
component of size b, regular degree (n-1)/2 would imply both a and b
at least(n+1)/2, impossible for disjoint components. Its canonical
maximal strong-module quotient is therefore prime. Its children must
be regular by (1), and their sizes satisfy (3). Conversely, substitution
over a prime quotient has these blocks as its canonical maximal proper
strong modules, so the prime skeleton and child types are recovered.
Isomorphism classes of child choices are the orbits of the actual Aut Q
action. Burnside's lemma gives (4). QED.

No mere multiset of weights or aggregate automorphism GROUP ORDERS
suffices for (4). The weighted action is the required observer. It
refines the older ordinary prime-cycle-index grammar specifically for
regularity. It still needs the positive-kernel prime inventory; this
paper does not eliminate that bottleneck or prove a faster unrestricted
A000568 algorithm.

Converse sends K to -K, leaves the positive kernel ray and size vector
unchanged, and reverses each child. It does not identify a non-self-
converse class with its converse in the ordinary enumeration.

## 4. Concrete prime7 -> composite11 ->505 order12 children

Here is an actual prime seven-vertex skeleton, with zero-based ports:

    K = [ 0 -1 -1  1  1 -1 -1 ]
        [ 1  0  1 -1 -1  1 -1 ]
        [ 1 -1  0 -1  1 -1  1 ]
        [-1  1  1  0 -1  1  1 ]
        [-1  1 -1  1  0  1 -1 ]
        [ 1 -1  1 -1 -1  0  1 ]
        [ 1  1 -1 -1  1 -1  0 ].

Its primitive positive kernel is (3,1,1,3,1,1,1). Every pair-module
closure is the full vertex set, and all5040 permutations have been
tested: Aut Q is trivial. Put a directed triangle at ports0 and3 and
a singleton at the other five ports. The resulting eleven-vertex
tournament is regular, with full automorphism group C3 x C3. Canonical
prime decomposition ensures this is the full group, not merely a
displayed subgroup. The checker directly verifies all eleven row sums.
Multiplying all sizes by three gives another regular inflation at33;
increasing only the first size by two destroys regularity. This is the
cheapest hostile against independent child-size selection.

For the eleven-vertex graph the number of binary subset orbits is
4² 2^5=512, and the number of vertex orbits is7. The incoming regular-
parent theorem therefore gives exactly

    512-7=505                                                   (5)

ordinary isomorphism classes of one-vertex extensions. These505 classes
have order12; they are one parent's exact contribution, not A000568(12).
Different regular parents' extension families are disjoint by that same
theorem. The checker constructs the nine actual child rotations and
independently counts all2048 subset orbits and all eleven point orbits.

This connects a signed linear invariant to a canonical composite graph,
then to an unmarked enumeration count. It preserves ports until the
specific justified quotient by automorphisms. No metagraph path data is
silently discarded: (4)--(5) count objects, not native edge multiplicities.

## 5. Hostiles, switching and the spectral lane

Every odd tournament can be switched by the signs of its cofactor vector
to give a positive kernel. This follows from diagonal covariance of K
and its one-dimensional kernel. Switching was used only to FIND the
displayed witness. It is not a legal native complement move, need not
preserve primality, and changes named orientations. We test primality
again AFTER switching, and then keep the resulting literal graph.

The balanced ramified primes in [LIFTS.md](LIFTS.md) are already regular,
so their positive kernel vector is the all-ones vector. They supply exact
prime skeletons for uniform-size regular composition. The seven-vertex
witness shows why the larger weighted grammar is useful: its uninflated
skeleton is not regular, but its composite is. A common positive kernel
does not fix spectrum; and a small spectrum does not supply an admissible
Hamilton-path relocation. Those remain separate predicates.

The exact test universe consists of all1024 labelled five-vertex
tournaments,2000 seeded seven-vertex samples and300 seeded nine-vertex
samples:3324 kernel checks. Only the first universe is exhaustive.
The recurrence (4) is a proof from canonical decomposition, not an
extrapolation from these samples.
