# Three black moves around every rotation-symmetric endcap presentation

Research continuation: October 20, 2026.

**PROVED:** every marked circulant endcap of odd order at least five has
an explicit three-black replacement for its original marked path-reversal
edge, with an arbitrary retained middle tournament. Equivalently, every
blue endcap presentation whose full cyclic rotation is an automorphism,
at endcap order at least five, has such a replacement. The final labelled
tournament is exactly the original path-reversal target.

**Scope:** this explicit macro treats marked translation cycles. The
companion [regular-endcap theorem](regular.md) covers every marked blue
presentation of every regular endcap of order at least five, including
noncirculant regular endcaps and other Hamiltonian paths of circulants,
using a cited Hamiltonian-connectivity theorem and exact small boundaries.
General large nonregular endcaps and global black connectivity remain
open; no global black-connectivity theorem is inferred here.

## 1. Inheritance, marked objects, and scope

The current sources were read before deriving:

* [scc_basins_20261020d/README.md](../scc_basins_20261020d/README.md) and
  [moves.md](../scc_basins_20261020d/moves.md), Sections 3--4, already
  cover every paired outer strong SCC of order four or five by finite
  boundary identities. Their singleton-padding theorem is also all-order.
* [The general endcap criterion](../deficiency_three_20261020/ear_probe.md),
  Section 3, gives a three-move macro with closing/chord guards. Its
  first move fails to be black precisely when the full cyclic rotation
  is an automorphism. A regular pentagon is the explicit hostile.
* [The all-order descent and basin reduction](../black_descent_20261020c/scc.md)
  reduce black connectivity to nonstrong terminals; they do not already
  connect all terminals.

The closest mechanism is the exact path-edge XOR ledger. The corrected
near miss is to change to another Hamiltonian generator and pretend the
original labelled edge has been repaired. The less-used sidecar is the
action of a permutation supported on a proper interval of a circulant.
The live board is: marked paths, cyclic symmetry, proper-interval
permutations, block boundary arcs, endpoint score charges, and SCC basins.

Let B be the circulant tournament on Z/rZ with skew step set S:

    u -> v iff v-u mod r belongs to S,
    S contains exactly one of d,-d for every nonzero d.

Thus r is odd. Assume 1 belongs to S, so the marked path is
Q=(0,1,...,r-1), with its closing arc r-1->0. Take

    T = B^op -> Z -> B,

where the left endcap is read in converse-reversed coordinates and M is
any Hamiltonian path of Z, possibly empty. Write A=(left path,M), and
P=(A,Q). Put N=2r+|Z| for the TOTAL order. Endcap order r and total order
N must not be confused. P is blue exactly when the middle presentation
M is blue.

R_P reverses just the path arcs. It is the converse of nonpath complement,
so represents the same edge in the converse-merged line multigraph.
Every displayed word is evaluated in the tournament after the preceding
reversals. Our claim is literal equality R_P3 R_P2 R_P1(T)=R_P(T), not
an equality after a free isomorphism or change of Hamiltonian path.

## 2. No proper-interval rotation is a circulant automorphism

**Lemma (PROVED).** In a circulant tournament, a permutation rotating a
proper nonempty consecutive interval nontrivially and fixing its
complement is never an automorphism.

By translation, let the interval be I={0,...,L-1}, with 2<=L<r, and let
the nonzero rotation amount be ell modulo L. Put g=gcd(L,ell). Write
s(d)=[d in S], for 1<=d<r.

For every 1<=d<L, some ordered pair in I at difference d straddles the
rotation's wrap cut. Its image has difference d-L. If the permutation
were an automorphism, then

    s(d)=s(r-L+d)=1-s(L-d).                                (1)

If L were even, d=L/2 would contradict (1). Hence L is odd, and so is g.

Each outside vertex v=L,...,r-1 is fixed. Its relations to I force the
length-L word

    s(v-L+1),...,s(v)

to be invariant under cyclic shift ell. Consequently every such window
is constant on congruence classes modulo g. These windows imply global
g-periodicity on s(1),...,s(r-1): any two positions d,d+g in this range
belong to a common length-L window, because g<L. Equivalently consecutive
windows overlap on at least g positions and agree on every residue class.

Since g is odd, choose 1<=d<=g with 2d congruent to r modulo g. Both
d and r-d lie in {1,...,r-1}. Global periodicity gives s(d)=s(r-d),
contradicting skew symmetry. Thus the purported automorphism cannot exist.

This proof uses only the circulant orientation, not primeness of r,
triviality of stabilizers, or a classification of circulant automorphism
groups. It also handles noncoprime interval rotations: the retained g
coordinate is essential in the argument.

## 3. A guarded block extraction with an exact target

First allow B to be an arbitrary endcap with a chosen Hamiltonian word

    (a,C,K,D),

where C,K,D are NONEMPTY ordered blocks. Let c be first(C), x=last(C),
y=first(D), and z=last(K). Require the two specific arcs

    z -> a,       x -> y.                                 (2)

These replace the earlier single-vertex closing/chord guards; they are
not inferred merely from strong connectivity. Define

    P  = (A,a,C,K,D),
    P1 = (A,K,a,C,D),
    P2 = (a,rev K,rev A,rev D,rev C),
    P3 = (C,K,D,A,a).                                     (3)

**Block extraction lemma.** The successive words in (3) are directed
paths and have path-edge XOR E(P). The last two moves are black. The
first is black whenever the reordered right word (K,a,C,D) is not
obtained from (a,C,K,D) by an automorphism of B.

Directedness: P1 uses (2), plus the original internal block paths and
the cross-arcs A->B. In P2 the arc a->z reverses the last K-to-a edge;
rev(K) and rev(A) use the previously reversed internal paths; the last
A vertex in rev(A) has its unchanged arc to last(D); rev(D) then joins
rev(C) through the reversed x->y edge. In P3 the original directions
inside all four blocks have been restored. Its D-to-A edge is the one
just reversed by P2, and all other joins are untouched original path
edges. All vertex labels and all nonpath internal block arcs are retained.

For the parity ledger, the internal edges of A,C,K,D occur three times.
The unintended edges last(A)-first(K), z-a, first(A)-last(D), x-y each
occur twice. The remaining four joins are precisely the original joins
last(A)-a, a-c, x-first(K), and z-y. Hence

    E(P1) triangle E(P2) triangle E(P3) = E(P),             (4)

and the final matrix equals R_P(T).

The endpoint sums for both P2 and P3 are exactly

    d_B(a)+d_B(c) <= 2r-3 < N-1.                           (5)

P1 ends in D, so a and c are unchanged before P2. P2 lowers a by one and
raises c by one, cancelling in the P3 sum. A blue positional reversal
would require endpoint sum N-1, so the last two moves are black.

The left endcap path is unchanged in P1. A positional reversal pairing
it with the new right path would require the new right presentation to
have the same adjacency matrix as the original one, hence an automorphism
of B. This argument is local to the two endcaps and is valid even when M
is not blue.

### Applying the extraction to a marked circulant

Find indices satisfying

    t in S, k not in S, 2<=t<=k<=r-2.                     (6)

Put ell=t-1 and u=k-ell+1, and choose

    a=0, C=(1,...,u-1), K=(u,...,k), D=(k+1,...,r-1).

All three blocks are nonempty. The first guard in (2) is k->0, because
k is not in S. The second is u-1->k+1, whose step is ell+1=t in S.
The first path rotates the proper interval {0,...,k} by moving K to its
front. Section 2 excludes an automorphism. Thus (3) is three black moves.

For the usual cyclic interval tournament R_(2m+1), S={1,...,m}, take
t=2,k=m+1. The block K is a singleton; this repairs every standard marked
R_(2m+1) presentation, m>=2. It is not the earlier full-cycle rotation.

The single-vertex version is not enough for all circulants. For Paley11,
S={1,3,4,5,9}, no nonzero d satisfies both d in S and 2d in S. Nevertheless
t=3,k=6 gives K=(5,6) and a valid three-black detour. No change of the
marked Hamiltonian cycle is needed.

## 4. The unique missing step-set pattern and its repair

If no pair (6) exists, the membership word on indices 2,...,r-2 has no
present step preceding an absent step. It is therefore a block of absent
steps followed by a block of present steps. Skew symmetry, with 1 present
and r-1 absent, gives exactly (r-3)/2 of each. Thus r=2m+1 and

    S = {1} union {m+1,...,2m-1}.                          (7)

This is an exact classification of failure of the extraction grammar,
not a disconnectedness claim.

For m>=3 put

    F=(2,...,m-1), G=(m,...,2m-1), z=2m.

Both F and G are nonempty. The original path and replacement are

    P  = (A,0,1,F,G,z),
    P1 = (A,0,1,F,z,G),
    P2 = (rev G,z,rev F,0,rev A,1),
    P3 = (1,A,0,F,G,z).                                   (8)

The three special original arcs required for validity are

    m-1 -> z,       z -> m,       2 -> 0.

Their steps are m+1,m+1,2m-1, all in (7). The last arc supplies the
join from rev(F) to 0 in P2. After that move it supplies 0->2 in P3.
The internal A,F,G paths alternate direction exactly as displayed;
the cross-boundary A-to-right arcs provide the remaining joins.

The edges (m-1)-z, z-m, 0-2, first(A)-1 cancel in pairs; internal block
path edges occur three times. The surviving edges are exactly E(P).
Thus (8) again has the literal target R_P(T).

P1 rotates the proper interval {m,...,2m}, moving z to its front.
It is black by Section 2 and the unchanged left endcap. Every right
vertex originally has score m. P1 raises last(G) by one, so P2's endpoint
sum is (m+1)+m=r. P2 raises vertex 1 by one, so P3 has the same endpoint
sum r. Since N>=2r, both are strictly below N-1. All three moves are black.

The boundary m=2 is inherited rather than extrapolated across empty F.
For the five-vertex step set {1,3}, the third regular-pentagon row of
[moves.md, Section 3](../scc_basins_20261020d/moves.md) gives the path-reversal
words

    (A,0,1,4,2,3),
    (1,2,0,rev A,3,4),
    (3,A,0,2,4,1).

These are the earlier nonpath-complement words with the middle word
reversed; they retain the same literal reversal target. The checker
replays this inherited boundary explicitly.

## 5. The all-order consequence and its precise boundary

**Theorem.** For every skew S on odd r>=5 with 1 in S, every tournament Z
and every Hamiltonian path M of Z, the marked source B^op->Z->B above has
three successive black path reversals ending at exactly R_P(T).

Proof: either (6) holds, when Section 3 applies, or S is (7), when Section
4 applies. The proof is elementary and constructive, with no Hamiltonian-
connectivity theorem or census of large tournament classes as a dependency.

More intrinsically, suppose cyclic rotation of the chosen right-endcap
Hamiltonian path is an automorphism of B. It is a transitive r-cycle, so
numbering vertices along that path presents B as a circulant with 1 in S.
An even r-cycle cannot be a tournament automorphism, so r is odd. For
r>=5 this theorem supplies the bypass of that exact marked edge.

At r=3, B=C3 and empty Z gives C3->C3, whose nine Hamiltonian paths are
all blue. Thus a theorem asserting this black bypass for all odd r>=3
and arbitrary middle would be false.

For a Cayley presentation on Z/rZ and a UNIT step d, a marked step-d
Hamiltonian path can be numbered in its actual order; its normalized step
set is {j:jd in S}. The final route is transported back to the same labels.
Changing d after choosing the target instead changes the marked edge and
does not constitute a bypass of the original edge. Similarly this theorem
does not assert that every Hamiltonian path of a circulant has automorphic
cyclic rotation.

Under black-component contraction, every blue incidence covered here has
its two endpoints in the same component. This removes the full automorphic-
rotation family as a possible bridge obstruction. It neither identifies
the other terminal basins nor proves connectivity of the entire black graph.

## 6. Exact replay universe

Run from the repository root:

    python 05-knowledge/results/rotation_owners_spectrum_20261020e/black_check.py
    python -O 05-knowledge/results/rotation_owners_spectrum_20261020e/black_check.py

The checker is self-contained, uses no removable assertions, and checks
each directed path, all three exact positional blue predicates, an
independent path-edge XOR, endpoint charges, and the final labelled matrix.
Its [receipt](black_check.out) records:

* All 2,046 marked circulant step sets containing 1 at orders 5,7,...,23.
  At each order exactly one step set needs the threshold/base branch.
* All 2,034 proper-interval rotations starting at 0 for every skew step set
  at orders 3,5,...,11, including step sets not containing 1: none is an
  automorphism. Other interval starts are conjugate by translation.
* All 1,101 fixed-path middle presentations through order six, without a
  blue filter, with carousel7, Paley11, and threshold15 endcaps: 3,303
  actual three-black replays.
* The variable-length threshold formula at every odd order 7 through 101,
  with a three-vertex middle; the Paley11 single-vertex obstruction; and
  the pure-blue r=3 boundary.

These are finite controls of independently proved variable-size identities.
They are not a tournament isomorphism census, exhaustive arbitrary-path
enumeration, or numerical evidence substituted for the all-order proof.
