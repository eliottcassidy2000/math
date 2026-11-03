# Modular centers of the apex wedge lattice

Status: elementary deductions from G1–G4 in `apex_wedge_constructor_20261103/GEOMETRY.md`, with independent small finite checks in `center_verify.py`. These are properties of the response lattice, not of an individual tournament's automorphism group. The source uses J=D/2; all modular statements reduce this *integer* J, not D/2 after division in a modular ring.

Let m>=1, b_j=2^max(0,m-2-j), and S=S_m from G1. Put C=Z^m/S Z^m. For m=2d+2>=4, G3 gives C ≅ (Z/d)^2 and G4 gives a unique rational solution a=S^(-1)b, with a=g_bits/d and g_bits containing -1. Thus [b] has exact order d in C: the order is the least t such that t S^(-1)b is integral. Let λ([x],[y])=x^T S^(-1)y mod Z. This is a well-defined perfect alternating pairing on C (changing sign gives an equivalent convention). By |C|=d² and ord([b])=d, the cyclic subgroup <[b]> is a canonically specified Lagrangian, including the trivial case d=1. It has exactly d complementary Lagrangians: a primitive vector of (Z/d)^2 extends to a symplectic basis (e,f); complements are exactly <f+t e>, t mod d. These are integral response-lattice statements; they do not assert that an integral basis change preserves Boolean words.

**Full modular theorem.** Let N>=2 and t∈Z. If m=2d+2>=4, the congruence Sx≡t b (mod N) has a solution exactly when gcd(d,N) divides t. When it is soluble, it has exactly gcd(d,N)^2 solutions modulo N. In particular Sx≡b is uniquely soluble precisely when gcd(d,N)=1, and is otherwise insoluble.

Proof. The obstruction is the image of t[b] in C/NC. Since [b] has exact order d in (Z/d)^2, it is a primitive vector: for each prime dividing d its coordinates are not both divisible by that prime. Hence its image has order gcd(d,N), proving solvability. The Smith invariants are 1,...,1,d,d, so the modular kernel has size gcd(d,N)^2, and a nonempty solution set is its affine coset.

For odd m, B_m is unimodular by G2 and S has a primitive one-dimensional radical by G3. Let k generate it. Expanding the determinant of B_m after extending k to an integral basis yields (b^T k)^2 times an integer equal to 1, so b^T k=±1. Thus Sx≡t b (mod N) can hold only when t≡0 mod N; conversely then x=0 works. G3 shows N solutions in that case. At m=2, S=0 and b=(1,1), so exactly N² solutions exist when N|t and none otherwise. At m=0 the source has no bit coordinates and the question is vacuous.

**Binary distinction.** An actual center of J modulo 2 exists exactly when m≥4 and m≡0 mod 4; it is unique. For m≡2 mod 4 (including m=2) and odd m it does not exist. For m=6 the rational center has denominator 2 but no binary center; at m=8 the rational center has denominator 3 and a unique binary center. If gcd(d,N)=1, use the unique center a_N and for every actual Boolean pair w,u evaluate

    J(w,u) ≡ (w-a_N)^T S (u-a_N) (mod N).

This congruence is a valid response formula, while for general N the center need not be Boolean and translation need not preserve the Boolean cube. Modulo 2, every translation does preserve the full Boolean cube. To determine J modulo N from raw Hamiltonian counts, compute D modulo 2N before exact division by 2.

**Verification.** `center_verify.py` exhaustively enumerates all x and t for m=2..8 with declared moduli 2..5 as listed in `center_verify.out`, and checks all exact counts, including m=6,N=2,4 and m=8,N=2,3. These are finite controls; the proof above establishes all sizes and moduli.

**Exact source-presentation parity counts.** For all 2^(2m) unrestricted Boolean word pairs, the number with D ≡ 2 (mod 4), i.e. J odd, is 2^(2m−1) when no binary center exists; for m divisible by 4 it is 2^(2m−1)−2^(m−1). For a fixed w, the coefficient of u is b+Sw modulo 2. If nonzero, half of u give J odd; if zero, w is a center and the constant is w^T S w=0.

For the strong chart, let m=k+2 and take w=(1,x,0), u=(1,y,0) with k free bits. If k>=1 the number of word pairs with D ≡ 2 (mod 4) is (2^k−t_k)2^(k−1), where t_k=0 for odd k, t_k=1 for k divisible by 4, and t_k=4 for k≡2 mod 4. At k=0 the count is 0. Indeed J=c^T(y−x)+x^T S_k y. Modulo 2, c=e_1+e_(k−1) in zero-based coordinates when k>=2, and c=e_0 for k=1. S_k is the adjacency matrix of two paths on the even and odd positions. If k is odd, exactly one path has odd length and c pairs nontrivially with its radical, so the center equation has no solution. For even k, both paths have length k/2: when this length is even there is one solution; when odd, c pairs trivially with each radical and there are four. The same affine-function count as above finishes. These count marked Boolean presentations, not tournament isomorphism classes; the strong guard is as in G3.

`boolean_jump_verify.py` independently enumerates the unrestricted word counts through m=8; `strong_jump_verify.py` directly enumerates the guarded pairs through k=7.

Independent review checked the modular theorem and both parity formulas, including k=1,2 and the path radical argument. All three verifier scripts use explicit conditional raises, and their normal and `python -O` outputs agree; finite controls do not replace the all-size proofs.
