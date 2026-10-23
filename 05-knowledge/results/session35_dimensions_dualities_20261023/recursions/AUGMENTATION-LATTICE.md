# Third-round theorem: augmentation Pfaffian and integral Heisenberg data

Status: PROVED identities; explicit finite-exact hostiles. This is an arithmetic refinement of the rational Heisenberg construction, not a Hamiltonian-count reconstruction.

## Inheritance

`01-canon/theorems/THM-174-pfaffian-determinant-identity.md` already proves that det(I+2A) is a perfect square, with odd-order value (sum signed principal Pfaffians)^2. `THM-213-paley-pfaffian.md` supplies the Paley principal-minor spectral computation. The augmentation-lattice presentation below identifies that existing square root as a concrete integral alternating volume.

Do not confuse this square root with signed Rédei count R(T)=sum over Hamiltonian paths of permutation sign. `THM-1936-signed-redei-join-multiplicative.md` defines R; `THM-1966-signed-redei-count-independent-invariant-n7.md` proves that |R| can vary inside an adjacency-cospectral class. Our Pfaffian absolute value is adjacency-spectral and therefore cannot equal |R| universally. Paley7 provides a clean numerical hostile: 49 versus147.

## 1. The augmentation form is always nondegenerate

Let n=2m+1≥3, M=A−A^T, and B be the n×(n−1) matrix whose column i is e_i−e_{n−1} for i=0,...,n−2. Thus B is an integral basis of the augmentation lattice L={v in Z^n:sum v_i=0}. Set K=B^T M B. Explicitly

K_ij=M_ij−M_{i,n−1}+M_{j,n−1}.

For i≠j each entry is odd, and diagonal entries are zero. Modulo2, K=I+J in even dimension n−1. Since J²=(n−1)J=0 mod2, (I+J)²=I. Consequently detK is odd and nonzero.

The resulting Lie lattice L⊕Zz with [u,v]=(u^T K v)z is two-step nilpotent with one-dimensional center over Q. All nondegenerate alternating forms of fixed dimension are equivalent over Q, so these rational Lie algebras are the same Heisenberg algebra. The rational isomorphism forgets the integral alternating elementary divisors.

Under integral congruence, K has blocks d_i[[0,1],[-1,0]], where positive d_1 divides ... divides d_m. All d_i are odd. The ordinary Smith normal form has each d_i twice. Their product is |Pf(K)|; they retain the integral commutator pairing, not merely its rational rank.

For n≥5, d_1=1. Indeed among at least four old vertices at least two lie on the same side of the last vertex. For such i,j, M_{i,n−1}=M_{j,n−1}, so K_ij=M_ij=±1. At n3 the cyclic tournament instead has d_1=3.

## 2. Exact Pfaffian and determinant identity, with signs fixed

Use zero-based deleted index i, natural order on every remaining principal minor, and the ordered columns of B above. Exterior multilinearity / the Pfaffian Cauchy–Binet formula gives

Pf(K)=sum_{i=0}^{n−1}(-1)^i Pf(M with row/column i removed).

The sign is fixed by det of B after deleting row i, which equals (-1)^i for odd n. Write w_i=(-1)^i Pf(M deleting i); the skew cofactor identity is adj(M)=ww^T. The generalized determinant lemma then gives

detK = Pf(K)^2 = (sum_i w_i)^2
     = 1^T adj(M)1 = det(M+J) = det(I+2A).

Thus the proposed identity is exact with no scale factor n or n². The augmentation basis is unimodular in L, not orthonormal; switching to an orthonormal real basis would introduce a covolume factor and must not be confused with this formula.

In the earlier marked pair p_M(z), q_M(z)=p_M(z)F_M(z), odd n gives

detK=q_M(0).

Indeed adj(−M)=adj(M) since n−1 is even. Therefore the product of alternating elementary divisors is already determined by the marked spectral pair, equivalently by adjacency spectrum. Smith profile may retain more arithmetic than that product; any strict-refinement assertion needs a separate witness.

## 3. Small exact examples and a minimal odd-order obstruction

With cyclic orientations i→j when j−i mod5 is1 or2:
cyclic5: Pf(K)=−5, alternating Smith(1,5), H=15, signed R=−5.

With Paley7 orientations j−i mod7 in{1,2,4}:
Paley7: Pf(K)=−49, alternating Smith(1,7,7), H=189, signed R=−147.
Thus detK=25 and2401 respectively. The signs depend on the stated labelled basis; the absolute values do not.

The first Smith-profile/H collision among odd orders occurs already at n5. Exact outneighbor bitmasks are
T: [18,4,9,19,6],
U: [30,4,8,18,6].

Both have Pf(K)=−1 and alternating Smith(1,1), but H(T)=11 versus H(U)=5. Their signed Rédei counts are −3 and−1. Their augmentation matrices are

K_T = [[0,−1,−3,−1],[1,0,1,1],[3,−1,0,3],[1,−1,−3,0]],
K_U = [[0,−1,−1,1],[1,0,1,1],[1,−1,0,3],[−1,−1,−3,0]].

Independent integer Smith normal form calculation returns diag(1,1,1,1) in both cases. Thus even the full integral Heisenberg alternating type cannot recover Hamiltonian count. At n3 the two tournament isomorphism types are transitive (d=1,H=1) and cyclic (d=3,H=3); n1 is trivial. This proves minimality among odd orders, not a statement about a differently defined even-order lift.

Full fixed-path n5 scan (64 presentations) gives:
(1,1): H in{1,5,11};
(1,3): H in{3,9};
(1,5): H in{9,15};
(1,7): H=13;
(1,9): H=15.
Fixed-path enumeration meets every isomorphism class because every tournament admits a Hamiltonian path. Presentation multiplicities are not class counts.

## 4. Verification

`augmentation.py` computes Pfaffians recursively with memoization, verifies the signed cofactor sum, computes alternating divisors through gcds of principal Pfaffian minors, and calculates H and signed R by subset/end-vertex DP. The two n5 witness matrices also receive independent SymPy Smith normal forms. `augmentation.out` records all numbers above. `augmentation7.py` performs a bounded fixed-path order-seven profile census without Hamiltonian enumeration.

Scope: this gives a concrete arithmetic sidecar and a proved obstruction to its completeness. It does not equate Pfaffian volume with signed path count, nor infer a Heisenberg/Riemann zero connection from shared integers.

### Order-seven census outcome (FINITE-EXACT)

All32768 fixed-path presentations were checked. Alternating profiles are (1,1,d) with
 d in{1,3,5,7,9,11,13,15,17,19,21,23,25,27,31,33,35},
and (1,7,7). The latter appears in9 fixed-path presentations. No two different Smith profiles share the same |Pfaffian| at this order. Thus through odd order7 this census gives NO strict refinement of the Pfaffian product by the Smith profile. Do not advertise additional information beyond adjacency spectrum from these cases. The general integral classification is finer in principle, but a tournament-family realization of that strictness remains unestablished here.
