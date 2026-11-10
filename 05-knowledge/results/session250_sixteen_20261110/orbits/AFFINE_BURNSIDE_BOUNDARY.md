# Affine Burnside retains the missing switching sidecar

Status: PROVED-IN-NOTE, with FINITE-EXACT examples at orders 2 through 5. No canon ID, independent audit, priority claim, or Chan–Pak closure.

## Inheritance and intended gain

The inherited mechanisms are ordinary Burnside, the free cut-switching action, and the endpoint Fourier bound proved in [the switching boundary note](../geometry/SWITCHING_BOUNDARY_RECURRENCE.md). The prior [sequence synthesis](../../sequence_gf_20261110/SYNTHESIS.md) gives an odd-cycle Burnside compiler but explicitly lacks Hamiltonian-path data. The present combination identifies the precise affine stabilizer sidecar that converts switch packets back to **ordinary unlabelled tournaments**, then intersects its fixed-coset constraints with the endpoint degree bound. This is an exact alternative carrier, not a faster general enumeration algorithm.

Concept board: free switch torsors / affine stabilizers / odd permutation cycles / endpoint Fourier support / unlabelled H moments. Anchor: recover A000568 and its H-refinement; niche: remove fixed-coset Fourier terms before computing moments; wildcard: descend the low-order walk recurrence to the unlabelled quotient. Hostile control: even-order signed symmetries exist but have no fixed switching state. Another: packet-average H does not equal the unlabelled average. The exact provenance of the requested 1-5-4-5-1 partition is unresolved; no identification with this construction is used.

## 1. The affine sidecar and exact class census

Fix n>=2. Write K=F_2^n/<1>, N=|K|=2^(n-1), and T^x for cut-switching a tournament by x. K acts freely: a cut fixes every edge only if all its vertex bits agree. Let G=K semidirect S_n. For one representative T of each G-orbit, let

    P_T={pi in S_n : pi T=T^b for some b in K}.

The bit class b=b_pi is unique, and

    b_(pi sigma)=b_pi+pi b_sigma.

Thus P_T is a group and acts affinely on K by alpha_pi(x)=b_pi+pi x. The map x -> T^x identifies P_T-orbits with ordinary isomorphism classes contained in the G-orbit of T: pi(T^x)=T^(b_pi+pi x). The projection of the G-stabilizer onto S_n is injective, with image P_T. Record both P_T and its cocycle b, not merely its order.

Let f(x)=H(T^x), where H counts directed Hamiltonian paths. It is affine-P_T-invariant. Define

    A_C(t)=sum_(ordinary classes U in C) t^H(U).

Weighted Burnside gives the exact identity

    A_C(t)=1/|P_T| sum_(pi in P_T) sum_(x : (I+pi)x=b_pi) t^f(x),
    A_n(t)=sum_(G-orbits C) A_C(t),
    A_n(1)=A000568(n).

Indeed, each affine orbit contributes |P_T| copies of its common weight to the double sum, by orbit-stabilizer. No division by a guessed switching orbit size is involved.

### Odd-order collapse

An even-order permutation pi in P_T has **no affine fixed point**. Otherwise it would be an automorphism of T^x. A tournament automorphism has odd order: an even-order element has a nontrivial involutory power; an involution swaps a pair of vertices and would reverse their edge, impossible.

An odd-order pi always has affine fixed points. Its affine action has the same odd order (the signed stabilizer projects injectively), so if its order is m, the sum in the F_2 affine space of the m points alpha_pi^j(0), 0<=j<m, is a fixed point. Write pi's c cycles as C_1,...,C_c; all have odd length. A homogeneous fixed class has a representative with pi x=x: the alternative pi x=x+1 contradicts pi^m=I for odd m. Such vectors are constant on each cycle, modulo the global constant. Therefore the fixed locus has 2^(c-1) points.

Consequently the **unweighted class contribution loses the cocycle completely**:

    a_C=A_C(1)=1/|P_T| sum_(pi in P_T, ord(pi) odd) 2^(c(pi)-1).

The cocycle remains necessary for the H-weighted refinement. This is a cycle-index formula over each packet's actual permutation group, not over all S_n. It provides a local route to A000568 from switching classes and their permutation projections.

For a general affine permutation, a fixed locus is empty unless b belongs to im(I+pi) on K; when nonempty it is x0+ker(I+pi). Its dimension is c-1 if some permutation cycle is odd, and c if all cycles are even. The latter extra dimension comes from solutions of pi x=x+1, possible exactly when every cycle is even. Tournament-compatible even-order lifts never realize either nonempty case.

## 2. Endpoint degree meets fixed-cycle parity

The characters of K are indexed by even subsets S of [n]: chi_S(x)=(-1)^(sum_(i in S) x_i). Write

    f(x)=sum_(S even) h_S chi_S(x),
    h_S=N^(-1) sum_x f(x) chi_S(x).

The boundary theorem gives h_empty=n!/N and h_S=0 for |S|>d=2 floor(n/3). For odd-order pi, its fixed locus is x0+V with V the cycle-constant vectors modulo 1. Averaging a character over this coset is zero unless it is trivial on V. Precisely,

    V^perp={S : |S intersect C_i| is even for every i}.

Therefore the exact first-moment contribution of pi is

    sum_(x fixed by alpha_pi) H(T^x)
      =2^(c-1) sum_(|S|<=d; all cycle intersections even) h_S chi_S(x0).

The answer is independent of the chosen fixed point. This is the useful intersection of two different restrictions: low endpoint degree from path cancellation and separate even parity on each permutation cycle from affine Burnside. A fixed vertex can never belong to a surviving S; its coefficient need not be evaluated for this contribution. Identity pi retains only the constant coefficient, so its contribution is always n!.

For higher moments put

    h_S^[r]=sum_(S_1 xor ... xor S_r=S) product_(i=1..r) h_(S_i).

These are the Fourier coefficients of f^r. The same fixed-cycle formula holds with h^[r] and cutoff

    d_r=2 min(floor(n/2), r floor(n/3)).

Thus fixed moment order r needs only cycle-even endpoint tensors through degree d_r, not the full exponential weight t^f as an unstructured function. The possible surviving supports for a permutation with cycle lengths ell_i are counted by coefficients through d_r of

    product_i ((1+u)^ell_i+(1-u)^ell_i)/2.

This is a pruning bound, not a claim of independent or nonzero coefficients. For example, a single 3-cycle and otherwise fixed points allows only its three internal pairs and the empty set, at every n where d>=2.

To recover the complete H-enumerator one can retain all moments or directly use weighted Burnside. A fixed number of moments does not determine support uniformly as n grows. The bound also becomes vacuous in degree once r floor(n/3)>=floor(n/2).

## 3. A recurrence that survives the ordinary-isomorphism quotient

Let A be the star-switch adjacency operator: (Ag)(x)=sum_v g(x+e_v). It commutes with the affine P_T action because translations commute with star switches and permutations permute the n generators. Consequently A descends to the finite set of ordinary classes inside C, with an integer quotient matrix Q_C. Its entry U,V counts vertex-star switches from a representative of U landing in V. This is independent of representative, includes loops and multiplicities, and every row sums to n. It need not be symmetric: the orbit-size sidecar satisfies detailed balance |O_U|Q_UV=|O_V|Q_VU.

For r>=1 define

    B_C,r,k=sum_(P_T-orbits O) (A^k f^r)(x_O).

This is an ordinary-class-weighted sum of switch-walk endpoint H^r, not the number of distinct classes reachable in k steps. In particular B_C,r,0=sum_U H(U)^r. A character indexed by |S|=2j has A-eigenvalue n-4j. With J_r=min(floor(n/2),r floor(n/3)),

    product_(j=0..J_r)(E-(n-4j)) B_C,r,k=0,
    sum_(k>=0) B_C,r,k z^k
      =sum_(j=0..J_r) beta_C,r,j/(1-(n-4j)z),

where E shifts k and beta_C,r,j is the sum over affine orbits of the size-2j spectral projection of f^r. The components are themselves affine invariant, since A commutes with the action and its different weight eigenvalues are distinct. A zero eigenvalue gives a constant-in-z summand, with the usual k=0 convention.

The same recurrence holds after summing C. For r=1 the constant-mode coefficient is exactly (n!/N) A000568(n), whereas the k=0 value is A_n'(1). Hence this bridge retains both the census and its first H moment in a common walk generating function. It does **not** produce an n-recurrence for either sequence, and the quotient operation has not removed the need to find P_T and b_pi.

At n=5, every packet satisfies

    B_C,1,k=(15/2)a_C 5^k + (sum_U H(U)-(15/2)a_C).

Thus the global unlabelled switch-walk H GF at order five is

    90/(1-5z)+2/(1-z).

The first H moment is 92; the class count is 12. This is stronger than assigning the packet mean 15/2 to each ordinary class, which would incorrectly give 90.

## 4. Exact small examples and failure boundaries

Run `python3 05-knowledge/results/session250_sixteen_20261110/orbits/check_affine.py` from the repository root. It refreshes `output.json`. The independent direct route canonizes all labelled tournaments under S_n and computes H by all vertex permutations; the sidecar route first partitions by K semidirect S_n and computes affine fixed sets. The exhaustive universes have 2, 8, 64, 1024 labelled tournaments at n=2,3,4,5. It checks weighted Burnside, odd/even fixed-set rules, Fourier parity filtering, the boundary cutoff, factorial packet sum and the quotient walk recurrence.

Results:

- n=2: one packet, |P|=2; its nonidentity even-order element has no fixed point. A_2(t)=t.
- n=3: one packet, |P|=3; A_3(t)=t+t^3. Four switch states have H values 1,1,1,3, so their mean 3/2 differs from the ordinary-class mean 2.
- n=4: two packets, |P|=4 and 12, contributing t+t^5 and 2t^3. Both contain even-order fixed-point-free symmetries.
- n=5: two packets, |P|=5 and 3, contributing respectively t+t^9+t^11+t^15 and 3t^3+2t^5+t^9+t^13+t^15. Their sum is t+3t^3+2t^5+2t^9+t^11+t^13+2t^15.

The n=5 example separates ordinary classes with equal H and keeps the two H=15 classes distinct. Scalar switch spectra alone would not recover these multiplicities. The finite computation does not identify the unexplained 1-5-4-5-1 partition and does not prove any all-order realizability claim.

## 5. Concrete next question

The map is labelled tournament -> switching torsor with endpoint Fourier coefficients -> affine quotient. Exact H and ordinary isomorphism multiplicities survive only when P_T and its cocycle are retained. Internal path-forest addresses and labelled tournament identity are discarded by the Fourier aggregation. The missing practical step is a vertex-extension update for **both** the endpoint coefficients and signed stabilizer sidecar. A scalar recurrence on switching graphs cannot supply it. A useful next experiment is to extend the two order-five packets by one vertex and seek a small marked refinement that predicts their child cocycles, rather than fitting a recurrence to A000568 values.
