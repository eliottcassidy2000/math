# Exact score-flow decomposition and regular completion

Status: PROVED algebraic identities; accompanying floating-point checks are
controls, not replacements for the proofs. No uniform Hamiltonian asymptotic
is asserted. Session30, 2026-10-22.

## 1. Types and the complete-graph flow decomposition

Let S be the skew sign matrix of an n-vertex tournament, d=S1, M=S/n,
e=1/sqrt(n), P=I-ee^T and C=PMP. Then

    M = C + (d1^T - 1d^T)/n^2.

This follows from 1^T S1=0 and expanding PMP. The two terms are orthogonal
in Frobenius inner product, C1=0, and

    ||M-C||_F^2 = 2||d||^2/n^3.

This is an exact divergence-free / gradient decomposition of an antisymmetric
flow on the complete graph. It is not a bijection to a regular tournament:
nC generally has fractional, rather than sign, entries. The restoration sidecar
is d. Throwing away d loses endpoint-sensitive Hamiltonian information.

For almost-regular even n, d_i is +1 or -1, with n/2 of each sign, so
||d||^2=n and the squared Frobenius change is exactly2/n^2.

## 2. A uniform resolvent gap, without regularity

For every tournament ||M||_F^2=(n-1)/n. The nonzero skew eigenvalues occur
in pairs ±ix. Hence

    ||M||_op^2 <= (n-1)/(2n) < 1/2.

Compression gives ||C||_op<=||M||_op. Thus I+M^2 and I+C^2 are positive
definite, their inverses have norm at most2n/(n+1), and I-M^2,I-C^2 have
inverse norm at most1. These estimates are uniform over tournament classes.
They control each fixed resolvent contraction, NOT the infinite sum of cores.

Writing G=-M^2(I+M^2)^(-1), H=M^2(I-M^2)^(-1), O=M(I-M^2)^(-1),
Cauchy-Schwarz applied to the rows/columns of M gives

    |G_ij| <= 2(n-1)/(n(n+1)),
    |H_ij| <= (n-1)/n^2,
    |O_ij| <= 1/n + sqrt((n-1)/(2n)) (n-1)/n^2.

For O use O=M+M(I-M^2)^(-1)M^2. All are O(1/n) entrywise.
For almost-regular δ=d/n, (I+M^2)^(-1)δ=(I+G)δ is O(1/n)
entrywise as well. This supplies a uniform fixed-kernel leaf bound.

## 3. Exact determinant effect of the score-flow component

Define rho(B)=det(I-B)/det(I-iB) for a real skew B with ||B||<1.
The denominator is positive real: each conjugate skew pair contributes1-x^2.
This is the same positive ratio as det(nI-S)/|det(inI-S)|.

Set a=Me=d/(n sqrt(n)). In the orthogonal splitting span(e)⊕e-perp,

    M = [[0,-a^T],[a,C|e-perp]].

Schur complementation yields

    det(I-tM)=det(I-tC)(1+t^2 a^T(I-tC)^(-1)a).

Skew odd powers vanish in a quadratic form. Therefore, putting

    u=a^T(I-C^2)^(-1)a,
    v=a^T(I+C^2)^(-1)a,

we obtain the exact identity

    rho(M)/rho(C) = (1+u)/(1-v).

For almost-regular n, ||a||^2=1/n^2, so

    1 <= rho(M)/rho(C)
      <= (1+1/n^2)/(1-2/[n(n+1)]) = 1+O(n^-2),

uniformly. By contrast the path-tree correction to logΦ is already O(n^-1).
Thus a small spectral change can coexist with a first-order endpoint defect.
This is a structural reason not to substitute spectral normalization for the
full forest expansion. It does not alone prove an H asymptotic.

## 4. Unique regular completion

An almost-regular even tournament has a unique labelled completion by one new
vertex * to a regular tournament: orient S_i*=-d_i, S_*i=d_i. Then

    S+ = [[S,-d],[d^T,0]],        S+1=0.

Each old row needs exactly this final entry, proving uniqueness. Conversely,
deleting a vertex from a regular odd tournament gives an almost-regular even
one, and completion restores it. Consequently almost-regular isomorphism
classes correspond to rooted regular isomorphism classes. Equivalently their
number is the sum, over regular classes of order n+1, of the number of vertex
orbits. Unrooted regular isomorphism classes alone are not enough.

## 5. Deleted doubly regular tournaments: an exact hostile/control family

Suppose the completion is doubly regular of order q=n+1, i.e.
(S+)^2=J_q-qI_q. Block multiplication gives

    S^2=J_n-(n+1)I_n+dd^T,       Sd=-1,       S1=d.

Thus span(1,d) has squared eigenvalue -1, and its orthogonal complement has
squared eigenvalue -(n+1). The normalized skew squared spectrum has -1/n^2
with multiplicity2 and -(n+1)/n^2 with multiplicity n-2. Therefore

    rho = [(n^2+1)/(n^2-1)]
          *[(n^2+n+1)/(n^2-n-1)]^((n-2)/2).

The path-tree penalty t=nδ^T(I+M^2)^(-1)δ is exactly

    t=n^2/(n^2-1).

The centered C has two zero eigenvalues, and otherwise the second squared
value above. These statements hold for every deleted doubly regular tournament,
not only Paley examples. They fix a parameter-free flat-spectrum control.
For general almost-regular tournaments the same norm estimates prove

    1 <= t <= 2n/(n+1) < 2.

## Verification

Run python3 check_centering.py. Independent direct-matrix/determinant comparisons
cover all labelled tournaments n=2..5, all2640 labelled almost-regular n=6,
and deleted Paley tournaments q=3,7,11,19,23,31,43,47,59,67,71,79,83.
centering-results.json records PASS. Algebraic proofs above provide generality.

Open: uniform loop-expansion remainder; arbitrary second-order nonregular
prediction; transport of these rooted/spectral invariants to full metagraph
wiring. No equality between score deficit and blue-orbit deficit is asserted.
