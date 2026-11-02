# Eventual exact weighted stability for n >= 1500

Status: PROVED, independently audited and accepted. This quantitative
refinement uses the accepted [asymptotic proof](../asymptotic_stability/PROOF.md),
[matching-support theorem](../stability/MULTISTAR_STABILITY.md),
[triangle-variance bounds](../stability/TRIANGLE_VARIANCE.md), and
[rounding framework](../rounding_graph/PROOF.md). Equality uses the accepted
[H=1 scalar minimization](../stability/H1_STABILITY.md). No step depends on
the earlier n>=10,000,000 eventual-exactness conclusion; that historical
proof remains unchanged. This file is outside earlier frozen manifests
and makes no Git changes.

## Statement

For every integer n >= 1500, K_n = kappa_n, where, writing R=floor(n/2),

    kappa_(2R)   = (4/R)(1-(3/4)^R-(1/4)^R),
    kappa_(2R+1) = (4/R)(1-(3/4)^R).

The equality classification is inherited from the accepted matching-support
and H=1 scalar minimization: equality with V>0 consists exactly of ordinal
fair pairs, with one singleton in odd order.

## Accepted inputs

For q_e=p_e(1-p_e), put V=sum q_e, B=E Delta(T), and let W be the weighted
sum of unordered incident edge pairs. We use

* B >= W, W >= 2V^2/n-V/4;
* Delta >= B+(1-exp(-4V))/4;
* for every nonempty matching of m>=1 edges and weight M,
  Delta >= B+a_m M when n>2m,
  where a_m=(4/m)(1-(3/4)^m), and Delta >= B+kappa_n M in all cases;
* the descending-weight greedy matching has discarded mass D, with
  D <= sqrt(2nW); for each selected edge of weight w its assigned discarded
  edge group has mass S and wedge contribution at least both wS and S^2/4.
  The contributions of different groups are disjoint.

All these inputs are already recorded in the accepted asymptotic proof or
reviewed matching-support framework. Only positive-weight edges are selected.

## 1. Candidate-beating configurations

Fix V>0 and suppose Delta <= kappa V, with kappa=kappa_n. We have

    kappa R <= 4,           kappa <= 8/(n-1) < 1/187.

Write c0=(1-exp(-4))/4. The elementary bound

    exp(4) > sum_(j=0)^7 4^j/j! = 16319/315 > 50

shows c0>49/200. Thus kappa<c0. If V<=1, the nonlinear input gives
Delta/V >= c0, a contradiction. Therefore V>1, and

    kappa V >= c0 > 49/200,
    V > 49(n-1)/1600.                                      (1)

The wedge lower bound and W<=B<=Delta<=kappa V give

    V <= n/8+n kappa/2 <= n/8+4n/(n-1) <= n/7.              (2)

The last inequality follows from n-1>=224. Retaining the binary baseline
for the greedy matching of total weight M=V-D gives B<=kappa D. Hence

    B <= 2n kappa^2,          D <= 2n kappa <= 16n/(n-1).  (3)

For B=0 the same statements follow directly from W=0 and D<=sqrt(2nW).

## 2. Delete all light selected edges

Set eta=1/44. Let L be the number of selected greedy edges of weight at
most eta. Assume L>=1, and delete these edges from the matching. Let r_h
and M_h be the size and weight of the retained matching. Since the original
matching has at most R edges,

    r_h <= R-L,
    M_h >= V-D-eta L
        > 49(n-1)/1600-16n/(n-1)-n/88 >= n/120.             (4)

Here the last bound holds uniformly for n>=1500: after subtracting n/120,
the expression is

    f(n)=(577/52800)n-49/1600-16n/(n-1).

Its derivative is 577/52800+16/(n-1)^2>0, and
f(1500)=9253539/26382400>0. Thus r_h>=4M_h>n/30, so r_h>0.
Because L>=1, n>2r_h, so the coefficient is a_(r_h).

Put epsilon=(3/4)^(n/30). Then (3/4)^(r_h)<=epsilon. For n>=1500,

    epsilon <= 1/(1000n).                                 (5)

Indeed, n(3/4)^(n/30) decreases there: its logarithmic derivative is
1/n-log(4/3)/30 < 1/n-1/120<0, using log(4/3)>1/4. At n=1500, (5) is
the exact integer inequality 1,500,000*3^50<4^50.

Apply the retained-matching inequality and multiply by R-L>=r_h>0:

    kappa V(R-L) >= 4(1-epsilon)(V-D-eta L).

Rearrangement, followed by kappa R<=4, gives

    L[kappa V-4eta(1-epsilon)]
      <= [kappa R-4+4epsilon]V+4(1-epsilon)D
      <= 4epsilon V+4D.                                   (6)

Unlike use of n/2, use of the integer R has no parity-error term. By
(2) and (5), delta:=4epsilon V<=1/1750. Since L>=1, delta<=delta L.
The left bracket in (6) is strictly greater than 49/200-1/11. Therefore

    L[49/200-1/11-1/1750] < 4D.

The square bracket exceeds 3/20, with exact excess 271/77000. Consequently

    L < (80/3)D.                                          (7)

In particular D>0.

## 3. The light-edge contradiction

Write D_h and D_l for discarded mass assigned to heavy and light selected
edges. On heavy groups, W_group>=wS>=eta S; on light groups,
W_group>=S^2/4. Since the respective group collections are disjoint,

    D_h <= W/eta <= 44B,
    D_l <= 2sqrt(LW) <= 2sqrt(LB).

Thus D <= 44B+2sqrt(LB). If B=0 this already contradicts D>0.
Using B<=kappa D, (7), and D>0 yields

    1 < 44kappa+2sqrt((80/3)kappa)
      <= 44/187+2sqrt(80/561)
       = 4/17+sqrt(320/561) < 1.                          (8)

The last strict comparison is rational after squaring positive numbers:

    (13/17)^2-320/561 = 137/9537 > 0.

This contradiction proves L=0.

## 4. Closure and equality

Every greedy selected edge now has w>eta>kappa. If D>0, at least one
assigned group has S>0, so

    B >= W >= sum_groups wS > kappa D,

contradicting B<=kappa D. Hence D=0, and then B=0. The uncertainty graph
is a matching, and Delta>=B+kappa(V-D)=kappa V. This proves the universal
inequality, since all configurations with a larger ratio required no work.
The fair ordinal construction gives equality and proves K_n=kappa_n.

For any equality configuration with V>0, the same argument gives D=B=0.
Every supported binary completion then has zero strong gap and is transitive
(at n>=4); the accepted H=1 scalar minimization gives precisely maximal
fair-pair ordinal configurations, with the required singleton in odd order.

## Quantitative audit

The companion `verify_n1500.py` checks each nontrivial rational endpoint or
squared inequality exactly using Python's standard-library `Fraction` and
integer arithmetic. Uniform extension is supplied by the monotonic arguments
above, not by finite numerical sampling. This proof does not claim the
threshold 1500 is optimal and makes no assertion at unproved smaller orders.
