# Adjacent-cut height control and an infinite diagonal repair

Session 201. **PROVED, independently audited in the stated scope:** every full C1 target n has an exact physical pair inverse with both cover heights at most floor(sqrt(3(n-1))). The same sharper bound holds for balanced targets and for the synchronized joint image, without increasing a supplied J budget. The coefficient 3 is sharp for the adjacent-cut inequality and for re-cutting a specified word, not claimed optimal among all representations of a target.

**PROVED:** all even-index diagonal Mersenne targets have explicit balanced witnesses and at least four distinct numerical source profiles. **OPEN:** odd-index diagonal repair, universal nonempty fibers, and the global {7,21} spectrum conjecture.

## 1. Inheritance and connection contract

The inherited source orbit S starts at t=(1,0,2,1) and uses chronological moves

    O(D,E,U,d)=(D,0,2U,d),
    Z(D,E,U,d)=(D-E+U,d,2U-d,d),
    J(D,E,U,d)=(D-E+U,D-E,D+U,D).

All orbit members are actual rooted start-one tournament responses. A final J closes a source into a strong tournament with exactly one Hamilton cycle and F=D-E+U Hamilton paths. The balanced closing OJ has count G=D+2U. Only J changes the root.

The physical adjoint uses b=(3,1,3,1), z=(1,0,2,0) and

    Q=[[0,0,0,1],[0,1,0,-1],[0,0,1,-1],[1,-1,-1,1]],
    K(v,w)=v^T Q w, H(s)=2U-d.

Reversing a chronological word and exchanging O with Z while fixing J is dagger. Then K(c v,w)=K(v,sigma(c) w), K(t,s)=F(s), K(b,s)=G(s), and K(z,s)=H(s). All these identities hold on formal coordinates, but physical existence uses actual words. The old cutoff was floor(2 sqrt(n-1)).

Anchor: improve the exact physical inverse rather than assert occupancy. Niche: repair a specified infinite part of the diagonal Carmichael hostile. Wildcard, stopped: fixed-letter pumping from the smallest odd diagonal does not automatically produce a uniform balanced family.

Concept board: physical positive cone; two adjacent cuts; synchronized auxiliary boundary; minimum J costs; diagonal Mersenne source fibers. The source-to-target map is word reversal under K; it preserves the target and J budget. The joint version must transport b and z with the same word. Forgetting that word or the auxiliary vector still loses H. The new ingredient is a height inequality involving BOTH sides of the crossed letter, not an independent occupancy argument.

## 2. A one-letter inequality on the physical positive cone

Use nonnegative coordinates

    y=(a,f,e,p)=(D-E-d,U-D+E-d,E,d-E),
    (D,E,U,d)=(a+2e+p,e,a+f+2e+2p,e+p).

Write P for this raw-coordinate substitution matrix. Set

    u=(1,1,2,2), delta=(0,0,1,1),
    A_c(y)=U(c P y), B_c(y)=U(sigma(c) P y).

The following are row vectors in y coordinates:

    A_O=(2,2,4,4), B_O=(2,2,3,3),
    A_Z=(2,2,3,3), B_Z=(2,2,4,4),
    A_J=B_J=(2,1,4,3).

For each letter put

    R_c=3 P^T c^T Q P - A_c^T B_c - 3 delta^T delta.

Direct integer multiplication gives

    R_O=[[2,2,3,3],[2,2,0,0],[4,1,3,3],[4,1,0,0]],
    R_Z=R_O^T,
    R_J=[[2,1,4,3],[1,2,2,3],[4,2,2,3],[3,3,3,3]].

All entries are nonnegative. Consequently, for any two actual sources v,w and any c in {O,Z,J},

    U(c v) U(sigma(c) w) <= 3(K(c v,w)-d_v d_w).       (1)

This is an all-order bilinear inequality. It follows by multiplying R_c by two nonnegative coordinate vectors. It also holds on the entire stated positive cone; no cone point is declared to be a physical source.

## 3. The sharper cut theorem

**Theorem 1.** Every full C1 count n has a presentation K(v,w)=n with actual v,w in S and

    U_v,U_w <= C(n):=floor(sqrt(3(n-1))).              (2)

Every balanced target n other than the separate singleton value 3 has such a presentation with v from t and w from the physical orbit of b. Every joint pair (m,h)=(G(s),H(s)) has such a presentation

    m=K(v,w), h=K(v,q),
    v=t.alpha, (w,q)=(b,z).beta,
    U_v,U_w <= C(m).                                 (3)

The auxiliary q remains bounded coordinatewise by w, as proved in the frozen synchronized inverse. Minimum total J budgets and literal reconstruction are retained.

**Proof.** Start with a physical word gamma, and move its suffix to the other boundary under dagger. Let v_i and w_i denote the two physical sources at cut i. The left cover heights increase strictly and the right cover heights decrease strictly. The scalar K(v_i,w_i) stays n; in the joint case moving z by exactly the same suffix also preserves h.

Suppose a crossing occurs between consecutive cuts. Write their common letter c and their profiles as

    previous cut: (v, sigma(c) w),
    next cut:     (c v, w).

Before the first crossing, U_v<U(sigma(c) w); at the next cut U(c v)>=U_w. Thus the maximum of the two heights at the previous cut is B=U(sigma(c) w), and at the next cut it is A=U(c v). Choose the cut with smaller maximum. Equation (1) and d_v d_w>=1 give

    min(A,B)^2 <= AB <= 3(n-d_v d_w) <= 3(n-1).

The chosen cut proves (2), or (3) with the synchronized sidecar.

For full C1 the initial left height is 2 and the terminal right height is 2, so a crossing always occurs, allowing an initial equality. If the initial cut already has left height at least the right height, both are 2 and n>=3 gives (2). For balanced or joint targets the terminal right height is 3. Any nonempty source word has terminal left height at least 3, so it crosses; the empty word has n=5 and heights (2,3), satisfying 3^2<=3(5-1). This exhausts the boundary cases.

The chosen cut only repartitions the original word. It neither inserts J letters nor changes any physical recipe. Thus the two J costs sum to the original J cost. The synchronized q bound is inherited for every shared suffix, not inferred from the scalar m. QED.

## 4. Exact membership criteria with the smaller inventories

Let S[C] be the physical source inventory from t through U<=C, B[C] the physical inventory from b, and R[C] the synchronized inventory of pairs (b,z).beta through U_w<=C. All recipes are retained.

Theorem 1 and the inherited reconstruction give the exact equivalences:

* Full C1: n belongs iff some v,w in S[C(n)] have K(v,w)=n.
* Balanced: n=3, or some v in S[C(n)], w in B[C(n)] have K(v,w)=n.
* Joint: (m,h) belongs iff m>=5 and some v in S[C(m)], (w,q) in R[C(m)] satisfy both bilinear equations in (3).

The minimum J-cost versions are also exact. For the joint inverse the cost belongs to the full pair (w,q), not merely to w. The inherited synchronized collision ((27,9,27,9),(15,5,15,5)) versus ((27,9,27,9),(23,7,23,7)) still forbids merging only on w.

For the inherited nonbalanced descent write L_k=2(2^k-1) and n=m+L_k h with a joint parent. Since 3<=h<m<2h, one has n>2^k m. Its physical inverse therefore needs only

    C_k(n)=floor(sqrt(3(floor(n/2^k)-1)))              (4)

on both unweighted physical inventories, instead of the previous factor 2 outside the square root. The condition 2^k<=(n+1)/6 still bounds k. Accepted pairs compile the literal source alpha beta-dagger O Z^k and then the final J.

The improvement reduces the required physical height by a factor approaching sqrt(3)/2. It is an algorithmic and structural improvement, not a proof that any prescribed absent target acquires a parent. The canonical 443/617 no-k=1 hostiles remain intact.

## 5. Sharpness of the cut argument, with its precise limitation

**Theorem 2.** No coefficient c<3 can replace 3 in (1) uniformly on actual sources. Also no coefficient c<3 bounds the best pair height obtainable solely by re-cutting every specified physical word, even with only one J.

For k>=0 put d=2^(k+1)+1 and take

    v=t.Z.O^k=(3,0,3*2^k,1),
    w=t.Z^(k+1).J=(2d-1,d-1,2d,d).

For c=O one computes

    n=K(O v,w)=3d^2,
    A=U(O v)=3(d-1), B=U(Z w)=3d,
    AB/(n-d_v d_w)=9(d-1)/(3d-1) -> 3.

The whole source word

    gamma_k=Z O^(k+1) J O^(k+1)

has F=3d^2. At its crossing the two candidate maxima are 3d and 3(d-1). All earlier cuts have right height at least 3d; all later cuts have left height at least 3d. Hence the exact smallest maximum over all its cuts is 3(d-1), and

    [3(d-1)]^2/(n-1) -> 3.

For k=0, the displayed expression for v must retain E=1 (v=t.Z); the expressions for Ov, w, n and all height conclusions are unchanged. The verifier checks this boundary explicitly rather than using the E=0 expression at k=0.

Removing gamma_k's final O gives a balanced precursor with G=3d^2. The same cuts prior to that final letter give exactly the same physical pairs with the triangle boundary b, so the cut-only sharpness applies to balanced and synchronized inverses too.

This does NOT prove that 3 is the optimum cutoff coefficient after changing the representation of the target. Such a different source word may have a smaller pair inverse.

## 6. An explicit repair of half the diagonal Mersenne hostile

**Theorem 3.** For every even a>=2, the target

    n_a=(2^(a+1)-1)^2+2

is balanced with a zero-J precursor. It has at least four distinct numerical source profiles in its F fiber.

Put p=(a-2)/2, q=a+2 and

    x=t.O^p Z^q O^p.

Direct iteration gives, for every p>=0 and q>=1,

    G(t.O^p Z^q O^p)=3+2^q(2^(2p+2)-1).          (5)

Indeed with A=2^(p+1), B=2^q, the middle profile after O^p Z^q is

    D=2+(B-1)(A-1), E=d=1, U=1+B(A-1).

The final O^p multiplies U by 2^p, while leaving D unchanged; when p=0 it is the identity. Thus G=3+B(A^2-1). Our p,q make this
3+2^(a+2)(2^a-1)=(2^(a+1)-1)^2+2. The balanced graph is x.OJ. Its order is 2a+4, one more vertex than the original diagonal closing t.O^a Z^a J.

Here are four explicit F recipes:

    A: O^a Z^a,
    B: O^p Z^q O^(p+1),
    C: O^p Z^q O^p J,
    D: Z^(p+1) O^q Z^p.

A is the original diagonal. B and C are the inherited O/J twins of x; D is the dagger of B. All have F=n_a. C has deletion D_x=2+(2^q-1)(2^(p+1)-1)>1, while A,B,D have deletion 1, so C is distinct. B ends O and has E=0; A ends Z and has E=1, so A and B differ.

It remains to distinguish D. Let r=2^p>=1, so a=2p+2, q=2p+4 and 2^q=16r^2. Direct iteration yields

    U_A=32r^4-4r^2+1,
    U_B=64r^4-32r^3+2r,
    U_D=32r^4+16r^3-r+1.

Then U_D-U_A=16r^3+4r^2-r>0. Also

    U_B-U_D=32r^4-48r^3+3r-1.

At r=1 this is -14; for r>=2 it is positive because
16r^3(2r-3)+3r-1>0. Therefore A,B,C,D are pairwise distinct. This counts numerical source profiles, not graph isomorphism classes.

At a=2 the four profiles are (23,1,29,1), (17,0,34,1), (33,16,34,17), (3,0,48,1), all with F=51. At a=4 the balanced witness has H=963. The odd-index diagonal, including a=1 giving 11 and a=5 giving 3971, is not settled by (5).

## 7. Verification and remaining obligation

The verifier checks the exact one-letter matrices by raw-coordinate basis multiplication, their entrywise positivity, and all common-word cuts through a declared length. It compares the smaller full, balanced, joint and J-budget inverses against a complete forward source inventory. It independently counts Hamilton paths and cycles in reconstructed small graphs by subset dynamic programming. Normal and optimized Python runs must agree.

The same script checks the sharpness family algebraically on many k, all four diagonal profiles on many a, and literal C1 balanced graphs at the first even diagonal indices. These finite checks support the implementation; the all-order claims rest on the displayed positive matrices and proofs. Universal odd coverage remains open.
