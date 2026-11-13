# The inverse problem for involutions: recover the two-sheet wiring

The [inherited pairing enumerator](../joint_metagraph_dyadic_20261113/PAIRING_WIRING.md)
asks for a permutation of ordinary classes compatible with **merged** weighted
wiring. Here we retain also the actual weighted dart matrix on ordinary
classes before merging converse. This is stronger, and still weaker than
asking for a specified lift through marked tilings.

Let R be converse. Write i+,i− for the two representatives of a chiral
fibre i, and keep each singleton as its own representative. Between two
double fibres (including i=j) the ordinary dart matrix, by R symmetry, has
block [[a,b],[b,a]]. Its even coordinate a+b is observed by the merger;
its odd coordinate d_ij=a−b is the missing signed wiring. Links from a
singleton into a double fibre have equal weights on both sheets. These
facts follow just from equivariance of the dart matrix. Parallel lines,
including opposite sheets in a common fibre, are not discarded.

**Inverse theorem.** Let g be an involution of the merged weighted diagram
preserving fibre sizes. For a lift to preserve the ordinary dart matrix,
first require |d_ij|=|d_g(i)g(j)| for every i,j and equality of the signed
diagonal d_ii=d_g(i)g(i). For off-diagonal nonzero entries make a constraint
edge ij with binary label w_ij=0 when their signs agree and 1 otherwise.
Add an edge i--g(i) labelled 0 for every moved chiral fibre. Parallel
constraints are retained. Then:
1. An involutive ordinary-wiring lift exists iff every closed walk in the
   constraint graph has label sum zero mod 2.
2. If feasible, each connected component independently supplies two flip
   assignments. With one chosen solution, let a_C and b_C count the
   g-fixed chiral fibres in C assigned no flip and flip respectively.
   The fixed-count contribution of this g is
   x^(number of g-fixed singletons) times the product over C of
   (x^(2a_C)+x^(2b_C)). A component with no fixed fibres contributes 2.
3. In particular at g=1 every component of the nonzero off-diagonal odd
   support must flip as a unit; if it has size j, its factor is 1+x^(2j).

**Proof.** A proposed lift has a flip bit e_i per double fibre. The signed
matrix condition is d_ij=(-1)^(e_i+e_j)d_g(i)g(j), because flipping one
sheet exchanges a and b. The first requirements and edge equations are
exactly equivalent to this condition. Singleton and mixed entries are
already determined by merged wiring and fibre sizes. Squaring the lift
requires e_i=e_g(i), precisely the added zero-label edges. Such binary
equations have a solution exactly when path integration from an arbitrary
root is consistent. A solution is unique given the root bit, hence two per
component. Fixed ordinary vertices come from fixed singleton fibres or
fixed double fibres with flip bit zero. This proves the enumerator.
For simultaneous colors repeat the matrix test for every color, and put all
constraints into the same graph. The sign of a diagonal is an invariant,
not an edge that can be corrected by a sheet change.

This is a representation-theoretic inverse problem: the merger retains the
trivial character block and forgets the sign-character block. Even with
both blocks known only up to gauge, choosing a lift requires solving an
actual F_2 coboundary equation. In inverse mechanics terms, the observable
is not the completed field; loop compatibility decides whether local sign
choices can be integrated. The formula gives a genuine obstruction without
assuming any geometric or K-theoretic equivalence.

## Exact first comparison: orders 5 and 6

The literal replay starts from Hamilton-spine presentations, canonicalizes
ordinary classes by permutations, and builds the uncolored complement dart
matrix. Weighted vertex refinement by the multisets of
(neighbor color, dart multiplicity), starting from degree, distinguishes
all merged vertices in orders 4,5,6. Thus g must be the identity at these
orders, even before colors are required. At order 5 the two chiral fibres
have no off-diagonal ordinary coupling, so the stronger uncolored test
still gives x^8(1+x^2)^2. At order 6 there are 12 singleton and 22 double
fibres. The nonzero odd support on the latter is connected. Hence:

| order | merged-wiring-compatible | ordinary-wiring-compatible |
|---|---|---|
| 5 | x^8(1+x^2)^2 | x^8(1+x^2)^2 |
| 6 | x^12(1+x^2)^22 | x^12(1+x^44) |

The only ordinary-wiring lifts at order 6 are identity and converse,
which also preserve the colored and literal reflection data. No analogous
all-order rigidity claim is made. The point of the comparison is not
that 8 is periodic. It is that the independence of two choices at order
5 is real at that scale, while the odd wiring synchronizes all 22 at the
next scale. That is a concrete reason to retain the marked odd channel
before comparing growth with stable module equivalences.
