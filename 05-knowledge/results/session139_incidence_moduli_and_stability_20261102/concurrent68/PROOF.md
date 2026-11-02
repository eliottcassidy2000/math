# Combining concurrent degree cancellation with the tensor reserve

Status: PROVED, independently audited for every n>=68.
The earlier n>=140 proof and its verified checkpoint remain unchanged.

## Provenance and scope

The degree-cancellation lemma and retained negative-RB light budget below
are recovered from the concurrent packet
`../../stability_priority_moduli_20261102/STABILITY_N500.md` and
`../../stability_priority_moduli_20261102/STABILITY_N276.md`, read at commit
c16226249f1f27a335cf7c84afdeaa90879b9580. They are credited concurrent
mechanisms, not claimed discoveries of this note. We independently reprove
the parts used here. The new tensor product reserve is Session139's
`../tensor_matching/PROOF.md`, and the geometric-tail method is retained
from `../stability/PROOF_N526.md`. This combination is not an order census.

Let n>=68 and R=floor(n/2). Define the fair-pair candidate by

    kappa_n=(4/R)[1-(3/4)^R-1_(n even)(1/4)^R].

Write kappa=kappa_n.

For a weighted tournament A, put Delta=H-det(I+A)-abs(det A),
q_e=A_ij(1-A_ij), V=sum_e q_e, and W=sum_(incident unordered e,f)q_e q_f.
Let B=E Delta(T) for independent binary rounding T; the inherited
proper-triangle bound gives B>=tau>=W. We prove Delta>=kappa V, with the existing positive-V
fair-pair equality classification. The case V=0 is the inherited binary
bound, with deterministic transitive equality. Orders10..67 are not claimed.

## 1. The sharper degree budget and matching mass

Take the descending-weight greedy matching, total weight M, discarded mass
D=V-M. Give each discarded edge f the selected edge that first removed it.
If d_v are weighted degrees in the discarded-edge graph, discarded/discarded
wedges total (1/2)sum_v d_v²-sum_f q_f². The assigned selected/discarded
wedge pays at least q_f², because its selected weight is at least q_f.
These wedges are distinct, and do not overlap discarded/discarded wedges.
Thus

    W >= (1/2)sum_v d_v² >= 2D²/n.                    (1)

Assume Delta<=kappa V and V>0; otherwise no work is needed. Retaining the
binary baseline in the matching bound gives B<=kappa D. Whenever a reserve
c satisfies Delta>B+c, square completion gives

    M=V-D > c/kappa+B/kappa-D
      >= c/kappa-n kappa/8
      >= cR/4-1-1/(2R),                              (2)

using kappa<=4/R and n<=2R+1. Also

    V<=n/8+n kappa/2,    4V<=Q_R:=R+33/2+8/R.        (3)

## 2. Bootstrap the reserve using the new tensor theorem

Since kappa<=2/17<g(1), where g(V)=(1-exp(-4V))/4, the inherited nonlinear
bound Delta>=B+g(V) rules out V<=1. As g(1)>6/25, start with c=6/25.

The tensor reserve applied after rounding outside the greedy matching says

    Delta>=B+1-exp(-M)-4^(-R).                       (4)

The last subtraction is needed only when the selected matching is perfect,
in which case its size equals R. Otherwise it is zero, so (4) is a valid
uniform weakening. For R>=34, 4^(-R)<1/10000.

The lower bound in (2) increases with R. At R=34:

* c=6/25 gives M>1;
* c=63/100 gives M>4;
* c=49/50 gives M>7.

These bounds boost the valid reserve successively to63/100,49/50,99/100.
Exact elementary exponential controls suffice:

    exp(1)>65/24,
    exp(4)>sum_(j=0)^7 4^j/j!=16319/315>51,
    exp(7)>sum_(j=0)^4 7^j/j!>110.

They respectively give
1-24/65-1/10000>63/100,
1-1/51-1/10000>49/50,
and1-1/110-1/10000>99/100.
Consequently every candidate equality or counterexample satisfies

    kappa V-B > c:=99/100.                           (5)

This reserve is conditional on Delta<=kappa V; it is not a universal
positive constant gap for arbitrary weighted tournaments.

## 3. Retain the complete light-mass budget

Fix

    eta=2/17,   delta=1/21,   a=1979/4200.

Let L be the number of selected edges of weight at most eta. Suppose L>=1,
and retain the other r_h selected edges with total mass M_h. By (2),

    M_h>=M-eta L
        >R(c/4-eta)-1-1/(2R).

At R=34, four times the right side is greater than13. It increases with R,
so r_h>=4M_h>13 and r_h>=14. Put m=R-L>=r_h>=14. Since at least one edge
was removed, this retained matching is nonperfect. Its sharp coefficient
is at least a_m=(4/m)(1-(3/4)^m), by monotonicity. Let epsilon=(3/4)^m.
Keeping B, not dropping it, gives

    (kappa V-B)(R-L)>=4(1-epsilon)(V-D-eta L).

Since kappa R<=4, rearrangement yields

    L[kappa V-B-4eta(1-epsilon)]
       <=4epsilon V+4D-RB.                           (6)

## 4. Uniform geometric-tail payment

For1<=L<=R-14, convexity of log((3/4)^(R-L)/L) puts its maximum at an
endpoint. Therefore

    4epsilon V/L <=max{Q_R(3/4)^(R-1),
                       Q_R(3/4)^14/(R-14)} <delta.    (7)

At R=34 both inequalities are exact rational checks. Both endpoint bounds
decrease for all larger R: the first has successive ratio less than
(3/4)(1+1/R)<1, since Q_(R+1)<Q_R+1 and Q_R>R. For the second,

    Q_R/(R-14)=1+(61/2+8/R)/(R-14),

whose positive excess decreases. Thus no finite extrapolation is used.

The constants were chosen so that

    c-4eta-delta-a = a-4eta =43/71400>0.

Combining (5)--(7) gives the strict bound

    aL<4D-RB.                                        (8)

In particular D>0.

## 5. A shared budget contradicts L>=1

Split discarded mass into D_l from light groups and D_h from the others.
Disjoint group wedges give

    B>=W>=eta D_h+D_l²/(4L)>=eta D-L eta².            (9)

The last comparison completes a square; its two budgets are disjoint.
Combining (8),(9) gives

    D eta(a-4eta)<B(a-R eta²).                        (10)

The left side is positive. If a-R eta²<=0, this is impossible. Otherwise
use B<=kappa D<=4D/R and R eta>=34*(2/17)=4:

    B(a-R eta²)<=D(4a/R-4eta²)
                   <=D(eta a-4eta²),

again contradicting (10). Hence L=0.

All selected weights now exceed eta>=kappa. If D>0, assigned wedges give
W>=sum w_i S_i>eta D>=kappa D>=B>=W, impossible. Thus D=B=0 and
Delta>=B+kappa M=kappa V. Every positive-V equality configuration has
matching uncertainty and all supported binary completions transitive;
the accepted H=1 scalar theorem forces maximal fair pairs, with a singleton
in odd order. Conversely that construction attains kappa.

The cutoff68 is a conservative proof cutoff, not a claimed phase change
or optimal threshold. The argument uses the concurrent stronger graph
inequality together with the independently accepted tensor product theorem.
