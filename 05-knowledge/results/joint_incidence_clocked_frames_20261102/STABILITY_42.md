# Joint triangle incidence proves sharp weighted stability from order 42

**Status: PROVED, same-agent proof audit plus exact rational certificates.
No independent-agent or formal audit is claimed.** This extends the incoming
[order-48 theorem](../stability_748_guarded_encodings_20261102/STABILITY_48.md)
by a new triangle surplus, not by adding two overlapping lower bounds.
The accepted order-48 proof supplies all orders at least 48; the argument
below supplies the six remaining orders 42 through 47.

**Concurrent provenance:** while this proof's exact checks were running,
mainline `4e3e096444` published an independent
[capacity-potential proof of the same range](../tensor_primorial_ports_20261102/STABILITY_42.md).
Both routes start from order48. No priority over that published extension is
claimed. The triangle surplus (2), its overlap-safe aggregation (4), and the
full-quarter capacity corollary below are the distinct contributions here.

For a weighted tournament A (zero diagonal, A_uv+A_vu=1, entries in [0,1]),
write H for its weighted Hamilton-path sum,

    Delta=H-det(I+A)-abs(det A),  V=sum_e q_e,
    q_e=p_e(1-p_e),  R=floor(n/2),
    kappa_n=(4/R)[1-(3/4)^R-1_(n even)(1/4)^R].

**Theorem.** For every n>=42, Delta>=kappa_n V, with optimal coefficient.
For V>0 equality consists exactly of maximal ordinal fair pairs (and one
singleton if n is odd); at V=0 it is transitivity. The accepted n=3..9
results remain separate, so the unresolved interval after this result is
**10..41**. No Hamilton-path spectrum or A000568 extension is asserted.

## 1. Inheritance and the budget that must not be spent twice

Use the SAME independent binary rounding throughout, and put B=E Delta(T).
The accepted inputs from the order-48 proof and its cited predecessors are

    B>=tau>=W,  W>=2D^2/n,  W>=2V^2/n-V/4,
    Delta>=B+(1-exp(-4V))/4,  Delta>=B+kappa_n M.

Here tau is the expected number of cyclic triangles, W=sum q_e q_f over
unordered incident edge pairs, M is a greedy matching's weight and
D=V-M. Assign each discarded edge to the first selected edge incident to
it. A selected weight w and assigned mass S have a disjoint wedge budget

    G(w,S)=max(wS,wS/2+S^2/4).

The tensor matching theorem in
[Session139](../session139_incidence_moduli_and_stability_20261102/tensor_matching/PROOF.md)
also gives, padding a matching to R slots by zeros and applying AM-GM,

    Delta >= B+1-(1-M/R)^R-4^(-R).                    (1)

The last subtraction is only needed for an actual perfect matching;
using it universally weakens the bound safely. A retained nonperfect
r-edge matching of mass M_r gives Delta>=B+a_r M_r, where
a_r=(4/r)[1-(3/4)^r] decreases in r.

The prior [joint-profile reserve](../stability_priority_moduli_20261102/PAIR_RESERVE.md)
proved tau_uvx>=alpha(q_a+q_b+(a-b)^2). It could not simply be added to W.
The next lemma explicitly leaves the entire triangle's W budget intact.

## 2. A triangle surplus above all three wedges

Put p=A_uv, a=A_vx, b=A_ux, alpha=min(p,1-p), w=p(1-p), and

    z(p)=min(alpha,7/20),   rho(p)=z(p)-w.

**Lemma J.** On the entire real cube [0,1]^3,

    p a(1-b)+(1-p)(1-a)b >= z(p)[q(a)+q(b)]+q(a)q(b)
                           = W_triangle+rho(p)[q(a)+q(b)].       (2)

In particular rho>=w^2>=0. For alpha<=7/20 the exact surplus used here
is rho=alpha^2. This is a lower bound, not a claim that z is optimal.

**Proof.** Set c=7/20. Interchanging a,b replaces p by 1-p. The difference
in (2) is affine in p on each of [0,c], [c,1-c], [1-c,1]. Thus it suffices
to check p=0 and p=c. At p=0 the difference is
(1-a)b[1-a(1-b)]>=0. At p=c put r=a-b, s=a+b-1 and A=s^2. Sixteen times
the difference is

    1/5+(4/5)A+(44/5)r^2-(A-r^2)^2-(12/5)r.          (3)

For fixed r, A ranges over [0,(1-|r|)^2]; (3) is concave in A. At the
upper endpoint a or b is 0 or 1, and the inequality follows from
z<=min(p,1-p). At A=0 the worse sign has r>=0, and r^4<=r^2 gives

    (3) >= 1/5+(39/5)r^2-(12/5)r
         = (39/5)(r-2/13)^2+1/65 >0.

This proves (2). Finally alpha-(w+w^2)=alpha^3(2-alpha)>=0 and
w+w^2<=5/16<7/20, proving rho>=w^2. QED.

Every triangle contains at most one edge of a matching. Subtract its
three wedge terms FIRST, and apply (2) only when it has a selected edge.
All other triangles retain their nonnegative inherited tau-W surplus.
Summing yields the joint-incidence inequality

    B >= W+sum_(selected i) rho(p_i)
                        sum_(x outside i)[q(ix at first end)+q(ix at second end)]
      >= sum_i [G(w_i,S_i)+rho(p_i) S_i].              (4)

An edge connecting two selected pairs can occur in two DIFFERENT triangles;
neither triangle is counted twice. Assigned discarded edges are only a
subset of these nonnegative terms. This is the missing overlap audit.

Hostile: the tempting uniform replacement tau>=4W/3 is FALSE. The triangle
(p,a,b)=(7/8,7/8,1) has tau=1/64 and W=49/4096, ratio 64/49<4/3.
Thus the distinguished matching edge, rather than a uniform factor, matters.

## 3. A sharper local dual, paid by the triangle surplus

Fix s>=1/6 with s<7/20 and put t=s(1-s). For 0<=w<=t let
r=(1-sqrt(1-4w))/2, so 0<=r<=s and rho=r^2. Then

    G(w,S)+rho S
       =max(rS,(r+r^2)S/2+S^2/4)
       >= sS-(s/4)(4s-4w),                 S>=0.       (5)

**Proof.** Set u=r/s in [0,1], v=S/s. The normalized difference is

    max(uv,u(1+s u)v/2+v^2/4)-v+1-u+s u^2.

It increases with s for fixed u,v, so take s=1/6. The two branches meet
at v=2u(1-u/6). On the linear branch the minimum is at this boundary.
The quadratic branch has its free minimum at v=2-u-u^2/6. If that lies
in the quadratic branch, 2-3u+u^2/6>=0, hence u<3/4, and its value is

    u^2(12-12u-u^2)/144 >=0.

Otherwise the boundary value is

    1-3u+(5/2)u^2-u^3/3
       >=(47-13sqrt(13))/12 >0,

as follows by differentiating on [0,1] and using 47^2>169*13.
This proves (5). For a heavy w>t, z(p)>s since s<7/20; consequently
G(w,S)+rho S>=z(p)S>=sS, strictly if S>0. QED.

## 4. Six exact finite endpoints

Assume 42<=n<=47, V>0, and Delta<=kappa V. Then B<=kappa D.
Write g(V)=(1-exp(-4V))/4. The decreasing ratio g(V)/V excludes
V<=13/10: the positive degree-16 Taylor sum proves

    exp(26/5)>1250/7,
    g(13/10)>1243/5000,  (1243/5000)/(13/10)>4/21>=kappa.

Thus kappa V>B+c_0 for c_0=1243/5000. If kappa V>B+c, the discarded
degree bound gives, by completing a square,

    M>m_c:=c/kappa-n kappa/8.                          (6)

Using (1), advance c to any smaller rational number than

    F_n(c)=1-(1-m_c/R)^R-4^(-R).                       (7)

For each of the SIX explicit values of n, the exact verifier checks the
same chain

    1243/5000 -> 1/4 -> 51/200 -> 27/100 -> 8/25
               -> 12/25 -> 3/4 -> 9/10 -> 49/50.      (8)

Every m_c lies in (0,R/4), every advance is strict, and hence
kappa V>B+c with c=49/50. The argument does not extrapolate a floating
numerical iteration or a finite tournament census to all weights.

Take s=4/R and t=s(1-s), so 1/6<s<=4/21<7/20. Suppose L>=1 selected
edges are light (w<=t), with mass U; let r_h be the heavy count.
Since U<=tL and L+r_h<=R,

    r_h(1-4t)>4m_c-4tR.                               (9)

The verifier checks (4m_c-4tR)/(1-4t)>9 in all six cases, so r_h>=10.
Set m=R-L. Then 10<=m<=R-1 and the retained heavy matching is nonperfect.
Applying a_(r_h)>=a_m while retaining B gives

    cL-4U < 4D-RB+4V(3/4)^m.                         (10)

For clarity, multiply B+a_m(V-D-U)<=kappa V<=4V/R by m; then use
4V/R>B+c. This is exactly the retained-baseline rearrangement, not a
new independent payment of B.

The degree bound and B<=kappa V imply

    4V<=Q_R:=R+33/2+8/R.

For each R=21,22,23 and each 1<=L<=R-10, exact rational arithmetic checks

    Q_R(3/4)^(R-L)/L < c-4s.                         (11)

Equivalently one can check the two endpoints, since (3/4)^(-L)/L is
convex on L>0. Equations (10),(11) imply

    4s L-4U <4D-RB.                                  (12)

Apply (5) to light groups and the heavy bound after (5) to all others.
Using (4) and lambda=s/4=1/R gives

    B >= sD-lambda(4sL-4U)
      > (s-4lambda)D+lambda RB = B,

a contradiction. Thus L=0. Every selected edge has z>s>=kappa; if D>0,
(4) gives B>kappa D, another contradiction. Therefore D=B=0.
The inherited matching-support and H=1 equality theorems give exactly
the stated fair-pair equality locus and show the coefficient is attained.
This completes orders 42..47; the incoming theorem supplies n>=48.

## 5. What this does and does not compress

**Corollary J-capacity (PROVED).** For every weighted tournament, every
descending greedy matching as above, and every 0<=k<=1/4,

    B >= kD-k^2(R-4M).                                (13)

This extends the concurrent wedge-only capacity proof's parameter interval
k<=(sqrt(2)-1)/2. To prove it, consider one group. If alpha>=k, then
z>=k and (4) pays kS. If alpha=r<k, then rho=r^2 and

    G+rho S >= (r+r^2)S/2+S^2/4.

Writing b=(r+r^2)/2, one has 0<=k-b<=k(1-2r), because k<=1/4.
Completing the square gives

    G+rho S-kS >= -(k-b)^2 >=-k^2(1-2r)^2
                             =-k^2(1-4w).

Sum over the selected groups and pad to R slots with zero groups. In
particular B>=V/4-R/16. Optimizing (13) gives B>=D^2/[4(R-4M)] when
R-4M>0 and D<=(R-4M)/2, and B>=D/4-(R-4M)/16 otherwise. At zero unused
capacity the latter formula gives B>=D/4. These are additional scalar
consequences, not an all-order solution: the low-matching-mass bootstrap
obstruction at40,41 remains.

The new usable coordinate is the attachment of TWO external incidences to
the SAME selected pair. After paying the ordinary wedge sum, enough of that
coordinate survives as rho(p)S to improve the local dual. This is a sufficient
certificate for this inequality, not a complete attachment classifier.

Five-vertex stability remains the accepted local base case; its full joint
feasible region has not been compiled here. The present three-vertex surplus
is compatible with it and strengthens the all-size matching method.
At n=40,41 the specific bootstrap (8) fails; that is a failure of this
certificate, not a weighted-tournament counterexample. Orders 10..41,
general strong-atom minimal carriers, and arbitrary complement/re-marking
closure remain OPEN. In particular this result does not erase the
twelve-vertex 81/27 joint-response hostile.
