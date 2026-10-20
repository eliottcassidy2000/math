# Eulerian Fourier modes and bridge-free cut aggregates

2026-10-20. **PROVED elementary / FINITE-EXACT checks.** This addendum concerns only the weighted cut-partition aggregate of the path-complement incidence graph. The prior trijections/tree-gauge and order-five distance-four packets remain unchanged.

The inherited carrier is Theorem 3, Section 8 of [the trijections packet](../tournament_graph_trijections_20261020/TRIJECTIONS_AND_TREE_GAUGES.md). Let E=E(K_n), X=F2^E, B be the cut space, and Q=X/B. For an ordered Hamiltonian vertex word p with underlying undirected path P, put g_p=[E\P]. The aggregate matrix acts by

    (C f)(q) = sum_(ordered words p) f(q+g_p).          (1)

For n>=3 its positive entries all have a common weight w_n, it has no diagonal entries, and its row sum is n!. The support graph has adjacency A=C/w_n. The known values at n=3,...,8 are w_n=6,4,8,4,4,2; for every n>=8 the weight is2.

The closest inherited spectral mechanism is ordinary character diagonalization of a finite binary translation operator. The repo's [THM-4433, signed Hamilton second minimum](../../../01-canon/theorems/THM-4433-signed-hamilton-second-minimum-and-cross-scale-stability.md) concerns signed **cycles**. The transform below concerns ordered **paths**, so its constants cannot be copied from that cycle theorem. No priority claim for Fourier diagonalization or inclusion–exclusion is intended.

## 1. The complete Fourier transform is indexed by Eulerian graphs

Call an edge set Z Eulerian when every vertex has even degree in Z, including isolated vertices. The annihilator of the cut space is exactly the binary cycle space:

    B^perp = {Eulerian edge sets Z}.                    (2)

Indeed the scalar product with the vertex cut delta(v) is deg_Z(v) mod2; these cuts span B. Thus

    chi_Z(q)=(-1)^|Z intersect x|,      q=[x],           (3)

is well-defined exactly for Eulerian Z. These characters form an orthogonal basis of functions on Q. For n>=2 there are 2^binom(n-1,2) of them; at n=1 there is only the empty character.

**Theorem 1 (signed Hamilton-path transform).** The complete weighted spectrum of C is indexed by Eulerian Z and is

    lambda_Z = (-1)^|Z| sum_(ordered words p)
                            (-1)^|Z intersect P|.      (4)

The unweighted support eigenvalue is lambda_Z/w_n when n>=3.

**Proof.** Insert chi_Z into(1): each summand is chi_Z(q)chi_Z(g_p). Since g_p is represented by E\P,

    chi_Z(g_p)=(-1)^(|Z|-|Z intersect P|).

Factoring out the q-dependent character gives(4). Orthogonality and the correct number of characters make this the entire spectrum, not just a list of some eigenvalues. QED.

The Eulerian condition belongs to the **dual edge set Z**, not to the tournament itself. Equation(4) is not a formula for the tournament's raw H. Nor does it assert that signed path sums are invariant under switching the negative-edge signing Z. Path endpoints make that a different operation.

## 2. A finite linear-forest formula, not a formal spectral restatement

For any edge set Z, not necessarily Eulerian, define the signed path sum

    S_n(Z)=sum_(ordered words p) (-1)^|Z intersect P|.

A linear forest is a graph whose nontrivial connected components are paths. Write r(F)=|E(F)| and c_+(F) for its number of nontrivial components; isolated vertices contribute no orientation factor.

**Theorem 2 (explicit forest expansion).**

    S_n(Z) = sum_(F subset Z, F a linear forest)
                 (-2)^r(F) 2^c_+(F) (n-r(F))!.         (5)

**Proof.** Expand the product over Z:

    (-1)^|Z intersect P|
       = product_(e in Z) (1-2 [e in P])
       = sum_(F subset Z) (-2)^|F| [F subset P].

A subset of a Hamiltonian path cannot contain a cycle or a vertex of degree at least three, so only linear forests remain. Contract every nontrivial path component of F to a block. There are n-r(F) blocks in total. Arrange them in (n-r(F))! orders and orient each of the c_+(F) nontrivial blocks in two ways. Every ordered Hamiltonian word containing F is obtained exactly once. This proves(5), including the empty forest contribution n!. QED.

For n>=2 there is a precise optional connection to signed cycles: adjoin a new vertex with positive signs on all its incident edges. Deleting that distinguished vertex from an unoriented Hamiltonian cycle gives an undirected Hamiltonian path. Therefore S_n(Z) equals twice the signed **unoriented** Hamilton-cycle sum of this positive-apex graph on n+1 vertices. The extra vertex, its positive signs, the factor2, and the change of order are necessary; same-order cycle formulas are not interchangeable with(5).

## 3. Triangle modes and the exceptional zero at order nine

Take Z to be one triangle, with all other vertices isolated. Its linear-forest subsets consist of the empty set, three single edges, and three two-edge paths. The whole triangle is forbidden. Hence

    S_n(C3)=n!-12(n-1)!+24(n-2)!
           =(n-2)!(n-4)(n-9),

and its odd edge count gives

    lambda_C3 = -(n-2)!(n-4)(n-9),       n>=3.          (6)

This proves the exact sign and vanishing pattern:

- n=3: lambda=-6.
- n=4 and n=9: lambda=0.
- n=5,6,7,8: lambda=24,144,720,2880, respectively.
- n>=10: lambda is negative.

Every distinct labelled triangle is a distinct Eulerian character, and these characters are orthogonal. Therefore at either n=4 or n=9 the adjacency nullity is at least binom(n,3). At n=4 this gives4 zero modes and is exact by the inherited K8-minus-a-matching description.

### Complement characters distinguish full and component multiplicities

When n is odd, E itself is Eulerian. Since every generator has the same edge parity k=binom(n-1,2),

    lambda_(E xor Z)=(-1)^k lambda_Z.                  (7)

At n=9, k=28 is even. The84 triangle characters and their 84 distinct complementary characters are all zero modes, giving the stronger statement

    nullity(C on the full quotient) >=168,
    nullity(C on either parity component) >=84.        (8)

For the component statement, the annihilator of the even displacement subspace consists of0 and E. Thus restrictions identify precisely the pairs Z,E xor Z. Distinct triangle pairs remain distinct, so there are84 independent zero characters on each component. These are lower bounds, not a complete order-nine nullity census.

At n=5 the same mechanism gives 20 full-quotient modes of weighted eigenvalue24, pairing into 10 modes on each 32-vertex component. Dividing by the edge weight8 gives the eigenvalue3 of multiplicity 10 in [the distance-four six-bit graph](../tournament_graph_distance_four_20261020/DISTANCE_FOUR_SIX_CUBE.md). This is an exact cross-check of the new transform against the prior structural graph identification, with no full/component multiplicity ambiguity.

## 4. Exact peripheral spectrum and second moment

The main packet proves that differences g_P+g_Q span the image of the even-edge-weight subspace. For even n that image is all Q, because a one-vertex cut has odd size. For odd n it has codimension one, because all cuts are even.

An eigenvalue in(4) has absolute value n! exactly when its character is constant on every generator. Such a character annihilates all generator differences. Consequently:

- For even n>=4, lambda_empty=n!, and every other eigenvalue has absolute value strictly less than n!.
- For odd n>=3, the only possible peripheral characters are empty and E. They have values n! and (-1)^k n!, respectively; all others have strictly smaller absolute value.

Equivalently, for n>=3 the top eigenvalue n! has multiplicity 2 precisely when n=1 mod4, and otherwise multiplicity 1. The bottom eigenvalue -n! occurs precisely when n=3 mod4, then with multiplicity 1. This recovers the component and bipartiteness boundary directly in the complete spectral coordinates.

For n>=3, looplessness gives

    sum_(Eulerian Z) lambda_Z = 0.                     (9)

Character orthogonality, or the diagonal of C^2, also gives

    sum_(Eulerian Z) lambda_Z^2
       = |Q| sum_g C[0,g]^2
       = |Q| n! w_n.                                 (10)

From n>=8, distinct undirected paths give distinct generators and w_n=2. Thus the exact character-average squared eigenvalue is

    (1/|Q|) sum_Z lambda_Z^2 = 2 n!.                   (11)

The factorial normalization and incidence weights are part of these identities. Replacing C by its unweighted support adjacency divides the first power by w_n and the second moment by w_n^2.

## 5. Every aggregate support edge lies on a commuting square

**Theorem 3 (bridges and the sharp small boundaries).** For every n>=4, every edge in the aggregate support graph belongs to a four-cycle, so the support has no bridges, including inside either parity component when n=1 mod4.

**Proof.** A nonzero generator g defines an arbitrary support edge q--(q+g). Choose a distinct nonzero generator h. Then

    q, q+g, q+g+h, q+h, q                             (12)

is a four-cycle. Its four vertices are distinct because g and h are distinct nonzero vectors over F2, and its successive edge differences are g,h,g,h. The cycle need not be induced; possible diagonals do not matter. For n>=4 the inherited generator-span dimension is at least two, so a distinct h always exists. An edge on a cycle is not a bridge. QED.

The small cases must not be silently included:

- n=3: the support is K2, whose unique support edge **is** a bridge. Its aggregate weight is6. If this weight is instead represented by six parallel incidence edges, no individual parallel edge is a bridge; deleting the support adjacency means deleting all six at once.
- n=1: Q has one vertex and the incidence matrix is[1].
- n=2: Q again has one vertex and the incidence matrix is[2].

The last two entries count identity-return incidences. After discarding loops, both support graphs have one vertex and no edges, hence no bridges. The diagonal convention here is literal incidence count, not the convention of counting an undirected loop twice in a degree.

### Why these squares do not lift automatically

The matrix C aggregates over every tournament representative of a switching class. Consecutive quotient edges may be realized from different representatives of their common intermediate class. Equation(12) guarantees a quotient square, not a sequence of four legal physical path complements from a specified tournament.

There is an explicit hostile at every n>=3. A transitive tournament has exactly one directed Hamiltonian word. By unique path witnesses, its degree in the **original labelled path-complement graph** is therefore1. Its complement move is not an identity when n>=3, so the incident edge is an actual bridge. At n>=4 its image in the aggregate nevertheless lies on the square(12).

Thus the original labelled graph can have bridges while the cut aggregate has none. This neither proves nor refutes a bridge assertion about a separately defined ordinary-isomorphism, converse-merged, black-edge, or regular-return metagraph. Those have different node/edge sets and require their own maps.

## 6. Reproducibility and scope

Run

    python3 check_fourier_squares.py > verification.json
    python3 -O check_fourier_squares.py > verification.optimized.json
    cmp verification.json verification.optimized.json

Both modes produce identical receipts with 2,744 explicit live checks. The checker imports no earlier producer code. It builds the aggregate kernels and full Walsh spectra at orders 1–6; compares every character against a literal signed Hamilton-word sum; checks the linear-forest expansion for every Eulerian signing through order 5; verifies the two spectral moments, peripheral multiplicities, all generator-square certificates at orders 4–6, and the transitive-leaf hostile. A separate signed Hamilton-path subset dynamic program verifies(6) at orders 3–11, including the order-nine zero. The168 order-nine masks are checked directly for distinctness and even degrees; no 2^28-state quotient is enumerated.

The bounded spectrum receipts at n=4,5 reproduce the prior cross-polytope and distance-four models. The n=6 receipt is only a control, not a new all-order spectral classification. No additional spectral family, mixing theorem, or global original-metagraph conclusion is claimed.
