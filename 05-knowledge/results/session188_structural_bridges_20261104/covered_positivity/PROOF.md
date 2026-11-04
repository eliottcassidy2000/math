# A positive cubic functional for every finite nonnegative profile

**PROVED.** The equal-profile cyclic-quotient contrast functional is strictly positive without any log-concavity hypothesis. The proof is a finite sum-of-squares identity over the integers. Its tournament application uses the separately derived exact response formula; no general response-reconstruction theorem is asserted.

## 1. Statement

Let q_0=0, let q_1,...,q_m be nonnegative real numbers, and extend q_j=0 for every j>m. Define

    P(q)=sum_(k>=0) (k+1)(q_(k+2)+2q_(k+1))
                          (q_(k+1)^2-q_k q_(k+2)).       (1)

All sums below are finite after zero extension. Then the following exact identity holds, even without nonnegativity:

    3P(q)=q_1^3
       +sum_(j>=1) j(q_(j-1)-q_(j+1))^2
                       (q_(j-1)+3q_j+2q_(j+1))
       +sum_(j>=1) (q_j-q_(j+1))^2
                       ((j+2)q_j+(2j+1)q_(j+1)).         (2)

The first sum in (2) includes j=m+1. Its last potentially nonzero term is (m+1)q_m^3. Omitting that boundary term would make the identity false.

Every term on the right is nonnegative when q is nonnegative. If q is nonzero, take m to be its last positive index. Then the first sum's j=m+1 term is strictly positive, so P(q)>0. In particular q_1>0 alone gives P(q)>=q_1^3/3>0. Keeping both upper boundary contributions also gives

    P(q)>=(q_1^3+(2m+3)q_m^3)/3                     (3)

when m is the last nonzero index.

## 2. Derivation, with all shifts explicit

Put j=k+1 and abbreviate a=q_(j-1), b=q_j, c=q_(j+1). Formula (1) is

    P=sum_(j>=1) j(2b^3+b^2c-ac^2-2abc).

Use the exact elementary identity

    -2abc=b(a-c)^2-a^2b-bc^2.

The shifted adjacent terms satisfy

    sum_(j>=1) j(q_j^2 q_(j+1)-q_(j-1)^2 q_j)
       =-sum_(j>=1) q_j^2 q_(j+1),

because q_0=0. Define

    A=sum_(j>=1) j q_(j-1) q_(j+1)^2,
    B=sum_(j>=1) j q_j q_(j+1)^2,
    C=sum_(j>=1) q_j^2 q_(j+1),
    D=sum_(j>=1) j q_j(q_(j-1)-q_(j+1))^2.

Thus

    P=2 sum_(j>=1) j q_j^3-A-B-C+D.                 (4)

The cubic comparison sums are, exactly,

    sum_(j>=1) j(q_(j-1)^3+2q_(j+1)^3)/3
       =sum_(j>=1) (j-1/3)q_j^3,
    sum_(j>=1) j(q_j^3+2q_(j+1)^3)/3
       =sum_(j>=1) (j-2/3)q_j^3,
    sum_(j>=1) (2q_j^3+q_(j+1)^3)/3
       =sum_(j>=1) q_j^3-q_1^3/3.

In the first equality, the contribution from q_m in the shifted q_(j-1) term occurs at j=m+1. These equations also hold for m=1.

Finally retain, rather than merely bound, each cubic Young deficit:

    u^3+2v^3-3uv^2=(u-v)^2(u+2v),
    2u^3+v^3-3u^2v=(u-v)^2(2u+v).

Substitute these identities for A,B,C in (4). The cube sums cancel except for q_1^3/3. Combining D with the A deficit and combining the B and C deficits gives precisely (2). This derivation is polynomial algebra; nonnegativity is used only after the identity is established.

## 3. Consequence for the equal-profile C3 comparison

The companion cyclic-quotient calculation states that, for three nonempty child tournaments with one common ordered nonempty path-cover profile q, the actual query with a equal to all of block 0 and c equal to all of block 1 has

    F_(C3[A0,A1,A2])(a,c)
      -F_(C3^op[A0,A1,A2])(a,c)=-P(q).              (5)

Assuming that exact formula, (2) proves a strict negative query difference for every common tournament profile. No log-concavity, transitivity, equality of labelled child orientations, or isomorphism between children is needed; equality of their q profiles suffices. The children retain their interior orientations and physical labels while only the quotient triangle is conversed.

This does not by itself cover unequal arbitrary profiles or lift an exact nonzero contrast through an arbitrary ambient module context. The formula (5) and its full named-query interpretation are independently audited in the companion packet.

## 4. Exact controls

The standalone checker expands both sides of (2) as cubic polynomials in formally independent variables q_1,...,q_m, with integer coefficients. It compares coefficient dictionaries for every m=1,...,20, explicitly retaining the j=m+1 boundary. It also checks all 97,648 nonzero sequences of lengths one through seven over {0,1,2,3,4}, including internal zeros, for strict positivity and the bound (3), after removing trailing zeros. It checks the zero sequence separately and signed controls to ensure the identity itself does not secretly rely on nonnegativity.

These finite checks support the indexing. The proof of the unbounded identity is Section 2. All checks use explicit exceptions and remain active under Python -O.
