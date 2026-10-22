# Exact deletion identity and nonzero lollipop limits

Use the notation of nonregular_forest_extension.md. Let T be almost-regular of even order n and U its unique regular completion. Set δ=M1, where M=S_T/n. Normalize the larger matrix at the same n:

K=S_U/n = [[M,−δ],[δᵀ,0]].

For any skew matrix A, write G(A)=−A²(I+A²)^−1 and d(A)=diag(A²(I−A²)^−1−G(A)). Put a=δᵀ(I+M²)^−1δ. Then

G(M)1=−(1−a) G(K)_(old,*).                                           (1)

Proof. For s=±i, the off-diagonal block of (I−sK)^−1 is

−s(I−sM)^−1δ/[1+s²δᵀ(I−sM)^−1δ].

The scalar quadratic form is a at both signs because odd powers of a skew matrix vanish in a quadratic form. Average the two signs. The numerator becomes M(I+M²)^−1δ=−G(M)1, and the common denominator is 1−a. This proves (1).

For n sufficiently large, all resolvents used below have uniformly bounded operator norm. In fact for M the tournament Frobenius norm and paired skew eigenvalues give ||M||op²≤(n−1)/(2n)<1/2. The regular completion K differs only in its harmless same-n normalization and still has an asymptotic gap from one. Moreover δ_i=±1/n and ||δ||₂=1/√n. Since every row of M has norm O(n^−1/2), each entry of (I−sM)^−1δ is O(1/n), for s=±1,±i: separate its identity contribution and bound the remaining row-resolvent product by C/n.

The Schur correction to the old-old resolvent block is rank one and entrywise O(n^−2). Therefore

d(M)_i=d(K)_i+O(n^−2), uniformly in i.

Also a=O(1/n), d(K)_i=O(1/n), G(K)_(i,*)=O(1/n), and g_i=(G(M)1)_i=O(1/n). From (1), including the omitted single diagonal term at cost O(n^−2),

L(T)=−2d(M)ᵀg=2(G(K)d(K))_*+O(n^−2).                               (2)

Now take a fixed regular base tournament B on m vertices, and any sequence of vertex-transitive regular tournaments C_t on odd t vertices. Rotational C_t suffices. Form U=B[C_t], delete a vertex in base block b, and put n=mt−1. Let P_t=J_t/t. Then

K=(tB/n)⊗P_t + I_m⊗(C_t/n).

These two summands annihilate each other, since B1=0 and C_t1=0. For every analytic matrix function f with f(0)=0,

f(K)=f(tB/n)⊗P_t+I_m⊗f(C_t/n).

Vertex transitivity makes the diagonal vector d(C_t/n) constant. Both G(tB/n) and G(C_t/n) annihilate the all-ones vector. Consequently, exactly,

(G(K)d(K))_*=(1/t)[G(tB/n)d(tB/n)]_b.

Insert this into (2), multiply by n, and let t tend to infinity. Since t/n→1/m and n/t→m,

lim nL(T)=2m[G(B/m)d(B/m)]_b.

No Hamiltonian-count asymptotic is used. Thus any base with a nonzero displayed vector proves that the lollipop is not uniformly a second-order term in the almost-regular class.

The explicit regular seven-vertex B in lollipop_obstruction.json gives, by exact rational matrix arithmetic,

2mG(B/m)d(B/m) = (−12607/2980692,−12607/2980692,12607/496782,
                 −12607/2980692,−12607/2980692,−12607/2980692,−12607/2980692).

Its sum is zero, as required by G(B/m)1=0. The exceptional coordinate is the base vertex deleted in that family; different labelings permute this vector.
