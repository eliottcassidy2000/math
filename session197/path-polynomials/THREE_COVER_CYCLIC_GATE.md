# A positive three-cover cyclic gate, a complete restricted inverse, and a self-converse 7p witness

Session 197. **PROVED exact physical kernel and target-bounded inverse within the stated construction family.** A16-vertex strong self-converse graph has H=37050461=7*5292923, with prime quotient. No infinitude or full-spectrum coverage is asserted.

This is an explicit specialization of the inherited colored-fragment algebra, not a new general composition theorem. The pinned S195 proof is https://github.com/eliottcassidy2000/math/blob/ee07a5aac9cdfd8d3834fa834df4fb357fcaf769/session195/pathcover-algebra/PROOF.md . Its Section4 already gives the one-input C3[A,Tr2,Tr2] ordered-cover row11,17,4. The bilinear kernel below recovers that row and exposes an additional paired construction and inverse bound.

## 1. Exact positive kernel

For a nonempty tournament T let c_j(T) count unordered covers by exactly j nonempty directed paths. Singleton paths are allowed. Put v(T)=(c_1,c_2,c_3)^T, with absent coordinates zero and c_1=H(T).

For any two nonempty tournaments T,U form the actual cyclic quotient

    G(T,U)=C3[T,U,Tr2],     T->U->Tr2->T.

Then

    H(G(T,U)) = v(T)^T K v(U),
    K = [5  6  0
         6 28 24
         0 24 72].                                      (1)

### Literal maximal-run proof

Cut a global Hamilton path at every interblock edge. Its maximal block runs are distinct actual fragments of local unordered path covers. Adjacent run colors must follow the directed cycle T,U,V,T, where V=Tr2; adjacent equal colors would not be maximal. The number k of V-runs is1 or2. Thus each of T and U has at most3 runs.

For given run counts i,j,k, a legal color word is uniquely determined by its starting color and length. Let N_(i,j,k) be the number of the three starting colors yielding those counts. Assign the actual local cover fragments to their named run positions in i! j! k! ways. Both unordered cover counts of Tr2, for k=1 and k=2, equal1. Therefore

    K_ij = i! j! sum_(k=1,2) k! N_(i,j,k).               (2)

For clarity, before multiplication by i!j!, the two N matrices are

    N_1 = [3 1 0       N_2 = [1 1 0
           1 1 0                1 3 1
           0 0 0],              0 1 1].

Equation(2) gives exactly (1). Cutting and assigning fragments are inverse because every joining edge crosses a module boundary. In particular no two adjacent fragments of the same module are allowed and no cycle of fragments is admitted. The factorials assign distinct actual paths; identical endpoint types do not identify fragments.

Taking U=Tr2 gives v(U)=(1,1,0), so the coefficients on unordered covers of T are11,34,24, equivalent to S195's coefficients11,17,4 on ordered covers. This checks both inheritance and normalization.

## 2. Bin-evaluation form and an uncancellable binary mode

Write h=H(T), b=B_T(2), c=B_T(3). The exact labelled-bin relations are

    c_1=h, c_2=(b-2h)/2, c_3=(c-3b+3h)/6.

Thus (1) is equivalently a symmetric bilinear form on (h,b,c), with matrix

    M = [15 -11  2
        -11  13 -4
          2  -4  2].                                    (3)

For fixed nonempty U, put h_U=H(U), b_U=B_U(2), c_U=B_U(3). Then

    H(G(T,U))=A_U B_T(3)+D_U B_T(2)+C_U H(T),
    A_U=2c_U-4b_U+2h_U,
    D_U=-4c_U+13b_U-11h_U,
    C_U=2c_U-11b_U+15h_U.                                (4)

In particular for T=Tr_m, m>=1, this is A_U*3^m+D_U*2^m+C_U. The coefficient D_U is always odd: b_U is even because swapping the two labels acts freely on nonempty two-bin configurations, and h_U is odd by the inherited odd-cycle identity. Therefore D_U is never zero. This entire family with one exterior transitive two-vertex block cannot produce a pure ternary clone ray. The positive form (1) is the counting proof; negative entries in (3) are only a change of coordinates.

## 3. The paired self-converse gate

Take U=T^op. Reversing every local path preserves each c_j, so the Hamilton count becomes

    F(h,u,v)=5h^2+12hu+28u^2+48uv+72v^2,                 (5)

where (h,u,v)=v(T). The graph has2|T|+2 vertices. Swapping the T and T^op blocks by corresponding vertex labels and reversing the two-vertex block is a fixed-point-free anti-involution. The cyclic quotient is strong and all blocks are nonempty, so the entire graph is strong regardless of whether T itself is strong.

As h is odd,

    F(h,u,v) = 5 (mod 8),                                (6)

because12hu+28u^2=4u(3h+7u) and u(h+u) is even; the remaining terms vanish modulo8.

There is also a useful prime-scope condition. The parity identity for bin profiles gives B_T(2)=B_T(-1)=(-1)^n h modulo3, where n=|T|. Consequently u=0 modulo3 if n is odd, and u=h modulo3 if n is even. Formula(5) then gives

- even n: F=0 modulo3;
- odd n: F=-h^2 modulo3.

For n>=2 the minimum bound proved below is45 or larger. Hence a paired prime output must have odd n and3 not dividing h, and is5 modulo24. A paired7p output likewise must have odd n,3 not dividing h, and p=11 modulo24. These are necessary conditions, not representation claims.

## 4. A sharp order bound and genuine restricted inverse

Every tournament induced on a nonempty vertex set has a Hamilton path, by the elementary insertion proof of Redei's theorem. Therefore every set partition into j nonempty parts contributes at least one unordered j-path cover:

    c_j(T)>=S(n,j),                                     (7)

where S(n,j) is a Stirling number of the second kind. The transitive tournament attains equality simultaneously for every j.

The kernel K has nonnegative entries. For the paired gate, (7) therefore gives the sharp order-specific bound

    H(G(T,T^op)) >= Phi_n,
    Phi_n = F(1,2^(n-1)-1,(3^n-3*2^n+3)/6)
          = 2*9^n-8*6^n+13*4^n+4*3^n-22*2^n+15.       (8)

The first values for n=1,...,8 are

    5,45,509,6069,69485,744405,7545389,73529109.

They are strictly increasing for n>=1: S(n,2) increases strictly, the other relevant Stirling numbers are nondecreasing, and (5) has a positive u coefficient and nonnegative coefficients. Equality in (8) holds precisely for transitive T. Indeed a nontransitive tournament contains a directed triangle and has H(T)>=3 by the odd-cycle identity, so the positive5h^2 term makes the inequality strict.

For a prescribed M, enumerate exactly those orders n with Phi_n<=M, then enumerate every labelled tournament of those orders, compute its first three cover counts, and test (5). A hit constructs the actual paired graph; exhaustion proves absence from this paired language at all orders. This is a complete finite inverse algorithm for this named gate, not a test of the unrestricted Hamilton spectrum. No efficiency beyond finiteness is claimed.

The unpaired family also has a complete target-bound: because v(U)>=(1,0,0) and (7) holds,

    H(G(T,U)) >= 5H(T)+6c_2(T) >= 3*2^|T|-1,

and symmetrically for U. Thus both input orders are bounded by the largest n with3*2^n-1<=M. Enumerating both finite tournament lists gives a complete inverse for the entire two-input gate as well.

The paired inverse through M=2000 has exactly the four outputs5,45,509,621. This follows from Phi_4>2000 and exhaustive inputs of orders1 through3. It illustrates that a complete restricted inverse need not cover all admissible global targets.

## 5. An explicit16-vertex7p graph

Use the seven-vertex tournament

    T=[78,124,120,112,97,65,0].

It is the transitive-sink padding of C3[K1,Tr3,Tr2], the same physical graph used in the separately reviewed non-TN counterexample. Direct path-cover enumeration gives

    (h,u,v)=(23,273,611).

The paired construction therefore has

    H=37050461=7*5292923.

The quotient5292923 is prime: trial division by every integer2 through2300 gives a nonzero remainder, and2300^2<=5292923<2301^2. The verifier records the complete remainder digest and the bound; it also directly counts the actual16-vertex graph, so the displayed Hamilton count is not only obtained from (5). Its anti-involution, every arc and strongness are verified on that graph. The exact rows and all certification data are written to the witness file.

This is one explicit new paired-gate witness, not a claim of least order for its count or of infinitely many7p outputs. The order bound Phi_8>37050461 proves only that any paired representation of this same target has an input of order at most7.

## 6. Reproduction and dependencies

Run `python3 verify_three_cover_cyclic_gate.py` and its `-O` form. The checker independently enumerates the finite run-color kernel, counts actual cyclic graphs for small arbitrary input pairs, verifies the coefficient transform and binary parity, checks the paired sharp order bounds, runs the complete paired inverse through2000, and directly counts/certifies the16-vertex witness. `--target M` runs the exact paired inverse; it may require a very large finite enumeration.

The counting proof in Section1 is self-contained and specializes the pinned S195 path-cover algebra linked above. Computational helpers are the frozen `verify.py`, SHA256 `10c2a35ce3f506aefb031645e81cffc1a2fb4d73bf9bbb38a163ccbf8e35dfcb`. The odd-cycle/parity convention is recorded in the frozen `PROOF.md`, SHA256 `b88c86ac0afa64d3a1bf399c56fa7f8aad046fb93424c70a9ace8d164708183b`. The non-TN packet identifies the same input graph but is not a dependency of this proof or count.
