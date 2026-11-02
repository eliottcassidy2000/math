# Sharp weighted tournament stability for every n >= 500

**Intermediate proof, superseded in range by [n>=276](STABILITY_N276.md).**
The open-range statements below describe this intermediate argument, not the
current frontier. The graph lemma and its proof remain dependencies of the refinement.

Status: **PROVED** (derivation below; exact arithmetic and independent graph
identities checked by `verify.py`). This extends the accepted n>=1500 theorem,
not an enumeration at order 500. The interval 10<=n<=499 remains OPEN.
No claim of an optimal threshold or of an independent human/agent audit is made.

## Setting and inherited inputs

For a real weighted tournament A, A_ij+A_ji=1, A_ii=0, 0<=A_ij<=1, put

    H = sum over vertex permutations of the product of successive arc weights,
    u = det(I+A),  v = (-1)^n det(A),  Delta = H-u-|v|,
    q_e = A_ij(1-A_ij),  V = sum_e q_e,  R = floor(n/2),
    kappa_(2R)   = (4/R)(1-(3/4)^R-(1/4)^R),
    kappa_(2R+1) = (4/R)(1-(3/4)^R).

Then Delta>=kappa_n V for all n>=500. Equality with V>0 consists exactly
of ordinal maximal fair pairs, with one singleton in odd order. Consequently
K_n=inf_(V>0) Delta/V=kappa_n.

Dependencies, all in
`../session138_stability_and_actual_contexts_20261102/`:

* `asymptotic_stability/PROOF.md` and `threshold_refinement/PROOF_N1500.md`:
  independent binary rounding T has B=E Delta(T)>=W, where W is the sum
  of q_e q_f over unordered incident distinct edge pairs;
  W>=2V^2/n-V/4, and Delta>=B+(1-exp(-4V))/4.
* `stability/MULTISTAR_STABILITY.md`: a matching of m edges and weight M
  gives Delta>=B+kappa_n M. If n>2m, its coefficient improves to
  a_m=(4/m)(1-(3/4)^m).
* The same asymptotic proof: descending-weight greedy matching assigns every
  discarded edge to the selected edge that removed it. A group of discarded
  mass S attached to an edge of weight w contributes at least wS and S^2/4
  to W. Groups' counted wedges are disjoint.
* `stability/H1_STABILITY.md` and the binary gap theorem in Session137:
  zero binary gap at n>=4 means transitive; the H=1 scalar problem has the
  equality classification stated above.

The argument keeps the binary baseline B in the retained-matching inequality.
It does not strengthen the false termwise Pfaffian bound identified by
`stability/PFAFFIAN_UNIT_COEFFICIENT_HOSTILE.md`.

## 1. A stronger greedy loss bound, sharp as a graph inequality

This lemma holds for every nonnegative weighted simple graph on n vertices,
without a tournament or q<=1/4 assumption. Let E_d be the discarded edges,
D their total weight, and d_v their weighted degrees. Separate W into
discarded/discarded and selected/discarded wedges. Selected edges form a
matching. Each discarded edge f was removed by a selected edge e of weight
q_e>=q_f. Its own selected/discarded wedge pays at least q_f^2; these wedges
are distinct for distinct f. Therefore

    W >= [ (1/2) sum_v d_v^2 - sum_(f in E_d) q_f^2 ]
         + sum_(f in E_d) q_f^2
      = (1/2) sum_v d_v^2 >= 2D^2/n.                     (1)

This improves the previously used D^2/(2n) by a factor four. The constant 2
is the best uniform constant: in the unit-weight K_(r,r), greedy matching
discards D=r(r-1), W=r^2(r-1), and nW/D^2=2r/(r-1)->2.

There is also a useful *joint*, rather than separately spent, wedge budget.
Fix eta>0 and let L selected edges have weight <=eta. Denote the discarded
mass assigned to these groups by D_l, and the rest by D_h. Disjointness gives

    B >= W >= eta D_h + D_l^2/(4L) >= eta D - L eta^2,  (2)

when L>0; the last step is (D_l-2L eta)^2>=0. The same last bound holds
when L=0, directly. No wedges have been counted twice.

## 2. Uniform numerical bounds on a candidate equality or counterexample

Assume V>0 and Delta<=kappa V. Set

    c=24999/100000, eta=1/62, a=1/15.

Throughout n>=500, kappa<=8/(n-1). Since exp(12)>25000 (already its
Taylor sum through degree 10 exceeds 25000), (1-exp(-12))/4>c.
The function (1-exp(-4V))/(4V) decreases for V>0. If V<=3, then
Delta/V>c/3>8/499>=kappa, a contradiction. Thus V>3 and

    kappa V>c,  V>c(n-1)/8.                              (3)

As in the inherited argument, the lower bound for W gives

    V<=n/8+n kappa/2<=n/7.                               (4)

Let M=V-D be the weight of the greedy matching. Crucially,

    Delta>=B+kappa(V-D),  hence B<=kappa D.

Together with (1) this yields the sharper bounds

    D<=n kappa/2<=4n/(n-1),
    B<=n kappa^2/2<=32n/(n-1)^2.                         (5)

They also hold when D=0. These are variance-mass bounds, not bounds on
the number of discarded edges.

## 3. Retain heavy matching edges, including the baseline

Suppose L>=1. Remove the L selected edges of weight <=eta. The remaining
matching has size r_h<=R-L and mass

    M_h >= V-D-eta L
        > c(n-1)/8 - 4n/(n-1) - n/124
        >= n/64 - 5/16.                                 (6)

The last inequality holds at n=500 with slack
654055969/12375200000>0; its difference has derivative
c/8-1/124-1/64+4/(n-1)^2>0. Since q_e<=1/4,
r_h>=4M_h>n/16-5/4>0. Put

    epsilon=(3/4)^(n/16-5/4).

Then (3/4)^r_h<epsilon and, uniformly for n>=500,

    epsilon<=1/(11n),  4epsilon V<=4/77.                 (7)

At n=500 the first inequality is exactly 5500*3^30<4^30.
The function n(3/4)^(n/16-5/4) decreases on this range since
log(4/3)>1/4 and 1/n-1/64<0.

Now n>2r_h, so retain B in the matching inequality:

    (kappa V-B)(R-L) >= 4(1-epsilon)(V-D-eta L).

Here R-L>=r_h>0, the right-hand mass is positive by (6), and
kappa V-B>=0. Rearranging, using kappa R<=4, gives

    L[kappa V-B-4eta(1-epsilon)]
       <= 4epsilon V + 4D - RB.

Absorb 4epsilon V<=4/77 into L(4/77), since L>=1. By (3), (5),

    kappa V-B-4eta(1-epsilon)-4/77
      > c-32n/(n-1)^2-4eta-4/77 >= a.

The last difference is increasing, and at n=500 its excess over a is
464029788839/178309616100000>0. Consequently

    aL < 4D-RB.                                         (8)

In particular D>0. Combining (2) and (8) gives

    D eta(a-4eta) < B(a-R eta^2).                        (9)

But a-4eta=1/465>0. If a-R eta^2<=0, (9) is impossible.
Otherwise B<=kappa D<=4D/R makes its right side at most

    D(4a/R-4eta^2) < D(eta a-4eta^2),

because R>=250>248=4/eta. This again contradicts (9).
Thus L=0. This step is where keeping the baseline pays for the whole
discarded mass; independently bounding the heavy and light budgets loses it.

## 4. Closure and equality

Every selected weight now exceeds eta>kappa. If D>0, the assigned wedges
give B>=W>=sum wS>eta D>kappa D, a contradiction. Hence D=B=0.
The uncertainty graph is a matching. The inherited matching-support theorem
gives Delta>=kappa V; every positive-probability binary completion is
transitive. Its Hamilton-path expectation is therefore H=1. The accepted
H=1 scalar classification yields precisely the asserted equality family.
That family attains kappa and proves the infimum statement.

The proof is uniform in all real edge weights and all n>=500. Its finite
checks authenticate rational constants and graph identities; they are not
being extrapolated to an infinite assertion. Orders 10 through 499, stronger
cycle-sensitive inequalities there, and the best eventual threshold remain
open. Even a sharper graph-only constant in (1) cannot help: its coefficient
2 is already optimal as n varies.
