# From binary Moore matrices to elliptic curves, through three-sheet tournaments

2026-11-26. **PROVED** graph-to-tournament transform and Moore extremal test,
using the **CITED** Monsky formula for the arithmetic interpretation.
**FINITE-EXACT** graph controls and prime witness. No rational rank,
density statement, or existence of a degree-57 Moore graph is asserted.

## 1. A lift that retains an arbitrary binary adjacency matrix

Let A be the symmetric, zero-diagonal binary adjacency matrix of any
simple graph on n ordered vertices. Make a tournament T(A) on
(i,x), with i a graph vertex and x in Z/3. Define the arcs by:

- If i<j and ij is a nonedge, direct all nine cross arcs from i to j.
- If i<j and ij is an edge, direct (i,x) to (j,x) and the other six
  cross arcs from j to i.
- Inside fibre i direct x to x+1 if n-1-i is odd (zero-based i), and
  x to x-1 otherwise.

Let L be the binary tournament matrix with the diagonal set to make
each row sum zero. There is an exact operator relation:

    dim_F2 ker(L²+L+I) = 2 dim_F2 ker(A²+A).                 (1)

**Proof.** Extend scalars to F4=F2(ζ), with ζ²+ζ+1=0, and split by the
three simultaneous-shift modes: functions constant on sheets, and
functions c_i ζ^x or c_i ζ^(2x). The cross blocks are I or the
all-ones 3-by-3 matrix J3, or their complements. On constant sheets
the quotient off-diagonal is the transitive tournament, since the
forward count is either one or three and the reverse count zero or
two. The corrected restriction of L is triangular over F2, so the
polynomial L²+L+I has diagonal entries one and is invertible there.

On c_i ζ^x the all-ones cross block vanishes and I remains I.
Consequently the cross operator is A in both directions. The diagonal
correction equals 1+q_i, where q_i=n-1-i mod 2 is the cross degree parity.
The internal arc contributes ζ if q_i=1 and ζ² if q_i=0. In either case
the sum with the correction is ζ. Thus the operator on this mode is
A+ζI and its polynomial is A²+A. The conjugate mode gives the same
kernel. Extending scalars does not change matrix rank, proving (1).

On a **marked** shift mode L itself retains A after subtracting ζI:
this is not merely a coincidence of its nullity. The scalar nullity
does not recover a graph in general, nor have we claimed that
forgetting the fibre marking recovers the graph from the isomorphism
class of the tournament.

### The arithmetic realization

The [October Selmer note](../bsd_selmer_growth_20261013a.md) and
[one-apex three-state calculation](../chanpak_affine_sieve_metagraph_20261122/MONSKY.md)
already use quadratic-reciprocity tournaments. Every finite
tournament is realized on distinct primes p_i ≡ 3 mod 8, with
i→j meaning (p_j/p_i)=-1. Choose p_j successively in prescribed
nonzero square or nonsquare residue classes modulo earlier p_i,
and in class 3 modulo 8, using CRT and Dirichlet. Reciprocity sets
the reverse arrow. For N the product of the primes, put

    E_N: y² = x³ - N² x.

[Monsky's formula, as explicated by Watkins et al.](https://www.dpmms.cam.ac.uk/~taf1000/papers/rankcongr.pdf)
gives dim Sel₂(E_N) - dim E_N(Q)[2] = dim ker(L²+L+I). Indeed on these
primes its two diagonal matrices for 2 and -2 are I and zero, giving
block rows [L+I,I] and [I,L]. Eliminating the second vector gives
the polynomial. The pure dimension is therefore exactly the
right side of (1), independent of which primes realize the marked
tournament. The curve itself and its rational rank are not determined
by unlabelled T. A large pure Selmer dimension is not a large rank
certificate; it may live in the Tate–Shafarevich direction. This
is the distinction in the [rank-moment follow-up](RANK-MOMENT-BARRIERS-AND-LONG-SEAMS.md).

## 2. Exact Moore criterion in a restricted arithmetic family

For a k-regular diameter-two Moore graph of order n=k²+1, the standard
route-count identity is A²+A=(k-1)I+J over the integers.
For odd k it reduces to J over F2, so (1) gives pure dimension 2(n-1).

Conversely, start from **any k-regular graph on k²+1 vertices with k odd**.
The diagonal of A²+A over F2 is identically one, hence its rank is at
least one. Equality in the associated Selmer bound ν≤2(n-1) is
equivalent to rank one, which for a symmetric binary matrix with all
diagonal entries one means J: a rank-one symmetric binary matrix
is v v^T and here v is the all-ones vector. Therefore every distinct
pair admits an odd (hence positive) number of length-one or length-two
routes, with the latter counted through a common neighbor. Each vertex
has precisely k² total such routes to other vertices. There are k²
targets, so each is reached exactly once. This proves the Moore condition.

An additional low-defect case has a direct interpretation. If A²+A has
binary rank two, its odd-route relation partitions the vertices into
two nonempty even sets: it is all ones within each set (including the
diagonal), and zero across. To see this, pick a zero off-diagonal entry.
Its two vertices form a principal identity matrix, giving coordinates
in the rank-two form. Every diagonal is one, so there are only the
two coordinate types (1,0) and (0,1). The cell sizes are even because
(A²+A)1=0 for a regular graph. The graph K4 disjoint union K(3,3)
is a 3-regular order-10 example, with pure dimension 16 rather than 18.

| Input graph | Graph order | Lifted tournament order | Pure 2-Selmer dimension |
|---|---:|---:|---:|
| Pentagon | 5 | 15 | 2 |
| Petersen | 10 | 30 | 18 |
| Hoffman–Singleton | 50 | 150 | 98 |
| Degree-57 Moore graph, if it exists | 3250 | 9750 | 6498 |

The even-valency pentagon has A²+A=I+J over F2; its kernel is the
constant line. The [Hoffman–Singleton classification](https://www.mathpuzzle.com/MAA/29-HoffmanSingleton/Venn%20Diagrams/ibmrd0405H.pdf)
supplies the historical context. The converse here concerns curves
obtained by this encoding of regular graphs at the Moore order, not
arbitrary high-Selmer curves.

For the pentagon an explicit prime realization is

    3,19,83,131,523,1051,1307,12539,13451,128987,
    155723,205427,280883,612331,1708411.

The replay verifies primality, all residue signs, and dimension 2.
It also gives the tournament Hamiltonian-path count 45,385,395.

## 3. Why 4, 8 and 9 have different jobs

- **4=2² is the coefficient field of the lift.** A three-entry
  binary circulant block has eight possible words. They correspond
  bijectively to their value at 1 in F2 and at ζ in F4, by the
  factorization X³-1=(X-1)(X²+X+1). One binary quotient coordinate
  and one four-valued coordinate coexist without loss.
- **9=3² gives a graph, not a Paley tournament.** In F9, -1 is a
  square, so square differences are symmetric. The square connection
  set is the union of the nonzero vectors in two distinct F3-lines
  (pair a nonzero square with its negative). In their coordinate
  basis the Paley graph is the three-by-three rook graph. It has
  degree four; adjacent pairs have one common neighbor, nonadjacent
  pairs two. Thus A²+A=0 over F2, so applying (1) gives pure
  dimension 18 on 27 tournament sites. It is not a Moore graph.
- **The Catalan pair is 8 and 9, not 4 and 9.** [Mihăilescu's theorem](https://doi.org/10.1515/crll.2004.048)
  says that these are the only consecutive positive nontrivial
  powers. On the *different* elliptic curve y²=x³+1 this coincidence
  is the point P=(2,3), where x³=2·4=8. Doubling has slope 2,
  gives 2P=(0,1), and tripling gives 3P=(-1,0). So the adjacent
  powers yield a point of order six, not a large-rank example.
  This arithmetic link is exact; it does not imply the matrix
  lift or prove any new adjacent-power theorem. The repo's
  [Pillai/dyadic cone work](../../../01-canon/theorems/THM-5109-dyadic-cones-pillai-chain-law-and-critical-fixed-port-sheets.md)
  asks about larger bounded positive gaps between powers of 2
  and 3 in a particular Hamiltonian gadget. The Catalan gap-one
  condition is much narrower.

There is one literal 4-and-9 overlap with the tournament counts.
For these three-sheet lifts, simultaneous sheet rotation acts freely
on ordered Hamiltonian paths, so 3 divides H. This family cannot
supply an exact +2 pair. The edgeless graph on two vertices lifts to
the ordered join of two directed triangles: H=9, ν=4. The earlier
[one-apex examples](../chanpak_affine_sieve_metagraph_20261122/MONSKY.md)
have H=9 with ν=0 and 2; append a terminal vertex to each to keep
both statistics unchanged. Thus on **six vertices with H=9** the
exact pure Selmer spectrum is {0,2,4}. There is no larger value:
ν is even by the two conjugate roots, and the constant vector is
never in the kernel because L1=0. This gives a sharp illustration
that equal Hamiltonian totals do not encode the arithmetic reading.

Replay: `python moore_monsky_check.py`.
