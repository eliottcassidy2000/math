# An unbounded class of block-uniform modulo-eight blindness

**PROVED by a finite exact polynomial certificate:** Let T=C3[A1,A2,A3] and U=C3^op[A1,A2,A3] on the same physical child sets, retaining every child interior. If each child has even order at least four, then every two-port query that is uniform on each child satisfies

    F_T(a,c)=F_U(a,c) (mod 8).

The children may have different orders, different profiles, and arbitrary nontransitive interiors. There are 64 such block-uniform queries. This is not an exact collision and not equality of full named responses: nonuniform masks remain unrestricted. For equal-profile children, the independent exact positive-square theorem already provides a nonzero integer gap, which here is divisible by eight.

## 1. Universal finite profile reduction

The companion `mod8_transfer/PROOF.md` proves that for an even-order child of order m>=4,

    q_0=0,  q_1=H,  q_2=2H-4,  q_3=4,  q_k=0 for k>=4   (mod 8),

with H odd. For q_3, use 6S(m,3)=3^m-3·2^m+3=4 modulo eight. Thus write H=1+2s, with unrestricted integer s. The entire relevant profile reduces to

    (q_0,q_1,q_2,q_3)=(0,1+2s,-2+4s,4) (mod 8).

The hypothesis m>=4 is essential: an order-two child has profile (0,1,2), and the three-TT2 quotient contrast is not blind modulo eight.

## 2. Exact physical-word functional

Put R=1-XYZ and

    N=[X,XY,XYZ; XYZ,Y,YZ; XZ,XYZ,Z],       G=N/R.

For a,c in {0,1}^3, write o=1-a and p=1-c. Define

    L_a=R+1^t N a,      R_o=R+o^t N 1,      M_(o,c)=o^t N c,
    W_(a,c)(N)=L_a(R+M_(o,c))R_p+L_c M_(p,a)R_o.

Then W/R^3 counts all words using both ports exactly once: the first term puts x before y, the second y before x, with the actual x->y direction. Prefixes and suffixes may be empty, and the middle segment is counted in full. This is the same maximal-child-run bijection as in `covered_response/PROOF.md`, now with arbitrary physical block masks and possibly unequal child profiles.

Set D_(a,c)=W_(a,c)(N)-W_(a,c)(N^t), and apply the functional

    L(X^rY^sZ^t)=q_r(A1)q_s(A2)q_t(A3).

The literal response contrast is L(D/R^3). For calculation modulo eight, indices above three vanish, so only finitely many monomials and terms of 1/R^3 are needed. The checker uses

    1/R^3=sum_(k>=0) binom(k+2,2) X^kY^kZ^k,

with k=0,1,2,3 sufficient after the profile reduction.

## 3. Exact coefficient certificate

Insert the three reduced profiles, with independent integer parameters s1,s2,s3. For each of the 64 physical mask pairs, the resulting contrast is a multilinear polynomial in these parameters: each child contributes one profile entry per monomial, and every entry is affine in its own parameter.

The retained certificate lists all eight integer coefficients of each polynomial, indexed by subsets of {1,2,3}. Every coefficient is divisible by eight. The checker obtains them by exact Möbius inversion from the eight evaluations s_i in {0,1}, then checks divisibility coefficient by coefficient. This proves the identity for every integer parameter value, not only these evaluations. The arithmetic is standard-library integer polynomial multiplication and finite coefficient extraction; no numerical approximation or external solver is used.

Together with Sections 1–2, this proves the claimed blindness at all child orders. It is a finite algebraic certificate for an unbounded theorem, not an extrapolation from a small tournament census.

## 4. Ambient and full-response scope

If the cyclic core is a module of a larger tournament, and the two ports are lifted with that module's exterior incidences, the port core is a common module of even order. The even-order transfer theorem therefore preserves the zero residue for every such canonical lifted block-uniform query in every ambient quotient.

This does not assert blindness for arbitrary exterior port attachments that break the augmented module, or for masks nonuniform within a child. It does not assert an exact integer collision. For example C3[TT4,TT4,TT4] and its quotient converse have

    F_T(15,240)=16075817,  F_U(15,240)=16483385,

an exact gap -407568, divisible by eight but nonzero. The complete named problem is not replaced by this restricted query class; cubic and quartic terms are included in the physical-word functional, not dropped.

## 5. Reproduction

Run `python verify.py > receipt.json` and `python -O verify.py > receipt.optimized.json`. The sibling `covered_response` and `response` packets must be available. All 512 symbolic coefficients are retained, along with independent literal queries on three twelve-name examples with transitive and strongly connected children. No claim is inferred from those finite literal controls beyond checking the implementation.
