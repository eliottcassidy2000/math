# Strong components split the tournament Selmer obstruction exactly

2026-10-13. **PROVED binary matrix identities; CITED Monsky arithmetic
interface; FINITE-EXACT complete small-order controls.** General BSD is
not proved here. The separately grounded E17 and E627 terminals are in
[the Tunnell packet](bsd_tunnell_terminals_20261013c.md).

## Inheritance, hostile and board

The closest mechanism is the row-zero tournament matrix in
[finite Selmer growth](bsd_selmer_growth_20261013a.md), where all prime
factors are 3 modulo 8 and the full pure 2-Selmer dimension is
nullity(A^2+A+I). The canonical hostile is the cyclic triple for
627=3*11*19: it has dimension two, but this is not a rank-two conclusion.
The corrected near miss is substituting a different curve after adding a
prime. The little-used sidecar is ordered strong-component factorization,
already used for path counts in
[THM-1936, signed Redei join multiplicativity](../../01-canon/theorems/THM-1936-signed-redei-join-multiplicative.md).
We use its elementary join proof, and repair its incidental strong-four
unsigned path count from three to five, preserving the signed result.

Board: ordered joins, quadratic characters, F4 obstruction, retained prime
labels, local conditions, and finite versus divisible Selmer growth.

## 1. The entire polynomial matrix splits

For a tournament T, let A_ij be its binary orientation indicator for i!=j,
and choose the diagonal to make every row sum zero over F2. Put
Phi(X)=X^2+X+1. For an ordered join T=T1->T2, let B,C be the row-zero
matrices of the two induced tournaments, let b=|T2|, and let J be the
all-one rectangular matrix. Then

             [ B+(b mod2)I   J ]
    A(T) =   [                 ].                       (1)
             [      0       C ]

**PROVED, every pair of finite tournaments.**

    Phi(A(T)) = diag(Phi(B),Phi(C)).                     (2)

The row-zero condition gives BJ=0. In a b-vertex tournament, each column
of C sums to b-1 modulo two: its row sum and reciprocity give this directly.
Hence JC=(b-1)J. The upper-right block of Phi(A(T)) is therefore

    (B+bI)J+JC+J=(b+b-1+1)J=0.

The diagonal blocks are unchanged because Phi(X+1)=Phi(X) over F2.
This is an identity of matrices, stronger than equality of their nullities.

Every tournament has a unique ordered chain of strong components C_i.
After grouping its vertex labels by that chain, iteration gives

    Phi(A(T)) = direct_sum_i Phi(A(C_i)),
    ker Phi(A(T)) = direct_sum_i ker Phi(A(C_i)).         (3)

No cancellation between different ordered components is possible. The
all-one cross blocks of A remain relevant to other observables; (2) is a
specific cancellation for this polynomial and this orientation law.

In particular, the shift invariance of X^2+X+1 is not numerology about
three: it is exactly the coefficient identity killing both the parity
shift on diagonal blocks and the rank-one cross block. The irreducible
quadratic also supplies the inherited F4 structure on each kernel.

## 2. Arithmetic meaning and its failure boundary

Take distinct primes p_i=3 mod8 and E_n:y^2=x^3-n^2x, n=product p_i.
Set A_ij=[p_j/p_i], where the bracket is the nonsquare bit. Reciprocity
gives an intrinsic tournament, and the cited Monsky matrix specializes to

          [ A+I  I ]
    M =   [        ],      ker Phi(A) ~= ker M,
          [  I   A ]

via x |-> (x,(A+I)x). Thus (3) exactly decomposes the **full pure
2-Selmer dimension** in this family. Prime labels, residue classes, and
the cited identification with full local conditions remain required.
The binary matrix by itself is not an arbitrary Selmer-group certificate.
The source and theorem are documented in the inherited packet, citing
[Watkins et al., section 4.1](https://dms.umontreal.ca/~andrew/PDF/RankQuadTwist.pdf).

For an ordered chain of r directed triangles, dimension is exactly 2r.
Every such finite orientation pattern occurs for distinct primes all
3 mod8 by CRT and Dirichlet's theorem, choosing the prescribed square or
nonsquare class at each prior prime. This constructs arbitrarily large
first-descent obstructions; it does not prove arbitrarily large rank.

Conversely, every tournament whose strong components are singletons or
strong four-vertex tournaments has dimension zero. The four-vertex claim
can be obtained from the inherited even switching/deletion theorem:
the unique strong four-vertex class switches, after a marked sink
normalization and deletion, to a transitive triple. The explicit matrix
and full enumeration are also retained here. This is an infinite
low-corank family after the same arithmetic realization. Under the
owner's supplied low-corank BSD theorem, every realized curve has full
BSD. No density of those prime patterns is asserted.

**Concrete failure of monotone obstruction.** E627 has dimension two.
Adjoining the prime 43 gives n=26961 and dimension zero, with all four
primes still 3 mod8. The fourth vertex mixes the old component, rather
than forming an ordered new component. This explains why (3) does not
apply to arbitrary prime extension. It certifies the new curve only.
Grounding the original E627 requires an unchanged-source argument; the
Tunnell terminal packet now supplies one by central-value nonvanishing.

The strong four-vertex representative has H=5 Hamiltonian paths and
signed path sum -1, yet Phi(A) is invertible. Thus a directed triangle
inside a graph does not automatically survive as an arithmetic kernel.
The precise operation is strong-component decomposition, not listing
cyclic subgraphs or confusing signed with unsigned path counts.

## 3. A counting recursion retaining the obstruction

Let F_n(t) count all labelled n-vertex tournaments, weighted by
t^(dim ker Phi(A)), and let S_n(t) count the strongly connected ones.
Set F_0=1. Choosing the unique source strong component gives

    F_n(t)=sum_(k=1)^n binom(n,k) S_k(t) F_(n-k)(t).     (4)

Each cross edge is forced by the ordered join, and (3) adds exponents.
Equivalently, their exponential generating functions satisfy F=1/(1-S).
This is a lossless component recursion when the labelled components and
their order are retained. Scalar polynomial counts do not recover an
individual prime-labelled source.

| n | Strong count at dimension 0 | dimension 2 | dimension 4 |
|---|---:|---:|---:|
| 1 | 1 | 0 | 0 |
| 2 | 0 | 0 | 0 |
| 3 | 0 | 2 | 0 |
| 4 | 24 | 0 | 0 |
| 5 | 264 | 240 | 40 |
| 6 | 13920 | 7680 | 720 |

The all-tournament row at n=5 is 624+360t^2+40t^4, totaling1024.
The total number on orders one through five is1099. These are the same
Boolean orientation counts occurring in the owner's tournament/even-graph
questions, now refined by an actual arithmetic invariant. No claim that
Hamiltonian count, chromatic number, and Selmer dimension are equivalent
is made.

## Reproduction and next target

Run `python 04-computation/experiments/bsd_component_recursion_20261013c.py`.
It exhausts all labelled tournaments through six vertices, computes
components by an independent reachability closure, checks the entire
matrix identity before nullities, and tests (4). A separate loop tests
all joins with sizes1..4 and1..3. Actual prime controls authenticate both
the Phi matrix and the full Monsky matrix. The literal five Hamiltonian
paths authenticate the inherited arithmetic correction.

Normal and optimized runs agree at **102441 checks**, with output saved
in the paired `.json`. There are no rank or favourable-component filters.
The useful next arithmetic target is a genuine higher-Selmer receipt on
a source where the central-value terminal also vanishes. A growing
first-level kernel alone cannot decide that problem.

**Pre-freeze API repair.** Independent review caught a source-order bug:
an in-place `list.sort` key iterated over the same list, which Python
temporarily exposes as empty. Kernel decomposition still passed because
it is insensitive to block order. The reader now returns a separately
sorted list, and every census item checks all cross edges point forward.
The two-vertex reverse-labelled tournament is the minimal hostile; the
actual E1947 component order is `[0],[2],[1]`. No theorem changed.
