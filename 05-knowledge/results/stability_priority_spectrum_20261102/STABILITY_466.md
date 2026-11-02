# Exact weighted tournament stability from order 466

Status: **PROVED**, using the accepted Session138 inputs below; same-agent
proof audit and exact rational certificate. No independent-agent audit or
formal proof is claimed. The threshold is sufficient, not claimed optimal.

Later refinement: [STABILITY_320.md](STABILITY_320.md) combines this retained
baseline with concurrent Session139's local dual and reaches n>=320. The
466 proof and its exact receipt remain unchanged as the first checkpoint.

## Statement and inherited mechanism

Let A be a real weighted tournament: A_ii=0, 0<=A_ij<=1, and
A_ij+A_ji=1. Let H be its weighted directed Hamiltonian-path sum, put
Delta=H-det(I+A)-|det A|, q_e=A_ij A_ji, and V=sum_e q_e.
For R=floor(n/2), set

    kappa_n=(4/R)[1-(3/4)^R-1_{n even}(1/4)^R].

**Theorem.** For every integer n>=466, Delta>=kappa_n V. Equality with
V>0 occurs precisely at ordinal sums of R fair two-vertex blocks and, in
odd order, one singleton. Thus K_n=inf_{V>0} Delta/V=kappa_n.
At V=0, equality occurs precisely at transitive binary tournaments.

This extends the [accepted n>=1500 proof](../session138_stability_and_actual_contexts_20261102/threshold_refinement/PROOF_N1500.md).
Together with its separate small-order results, the unresolved orders are
10..465 at this checkpoint. No assertion at those orders follows from this
note alone; the later320 refinement leaves10..319.

The closest mechanism is binary rounding with a matching retained. The
hostile family is the ordinal sum of fair pairs, which rules out a
dimension-free positive stability constant. The corrected near miss is
discarding the binary baseline before estimating light groups; the
least-used useful coordinate is the allocation of wedge mass between the
heavy and light groups. This proof keeps that allocation and the baseline.

## Accepted inputs

Let T be independent binary rounding of A, B=E Delta(T), and W the sum of
q_e q_f over unordered incident edge pairs. The accepted
[asymptotic argument](../session138_stability_and_actual_contexts_20261102/asymptotic_stability/PROOF.md),
[matching theorem](../session138_stability_and_actual_contexts_20261102/stability/MULTISTAR_STABILITY.md),
[triangle bounds](../session138_stability_and_actual_contexts_20261102/stability/TRIANGLE_VARIANCE.md), and
[rounding framework](../session138_stability_and_actual_contexts_20261102/rounding_graph/PROOF.md)
give, for n>=4:

1. B>=W>=2V^2/n-V/4, and Delta>=B+f(V), where f(V)=(1-exp(-4V))/4.
2. Retaining any nonempty matching of h edges and mass M gives
   Delta>=B+kappa_n M. If n>2h, also Delta>=B+a_h M, where
   a_h=(4/h)[1-(3/4)^h].
3. A descending-weight greedy matching selects weights w_i and assigns
   discarded mass S_i to each selected edge. Its disjoint group wedge
   contributions W_i satisfy W_i>=w_i S_i and W_i>=S_i^2/4.
   With D=sum S_i, W>=D^2/(2n).

Only positive-weight edges are selected. Group wedge contributions are
disjoint subsets of the actual W; they are not two independent budgets.

## 1. Keep the baseline in the putative extremizer

Assume V>0 and Delta<=kappa V, writing kappa=kappa_n. Let

    c=999/4000, eta=1/58, delta=11/100,
    Astar=c-4 eta-delta=8211/116000.

We have R>=233, kappa<=4/R, and 4/233<c/2. Since
exp(8)>sum_{j=0}^7 8^j/j!=47259/35>1000, f(2)>c.
The ratio f(V)/V is decreasing for V>0. If V<=2, then
Delta/V>=f(2)/2>c/2>kappa, a contradiction. Hence

    V>2,                  kappa V>B+c.                    (1)

The greedy matching has mass V-D, so input 2 gives

    B<=kappa D.                                           (2)

The wedge lower bound gives

    V<=n/8+n kappa/2<=n/8+5,                              (3)

using n<=2R+1 and n kappa/2<=4+2/R<=5.

## 2. A light selected edge forces a large retained matching

Suppose L>=1 selected edges have weight at most eta. Delete these L
edges, retaining h<=R-L edges of mass M_h>=V-D-eta L. From (1),
W>=D^2/(2n), and completion of a square,

    V-D-eta R
      >(B+c)/kappa-D-eta R
      >=c/kappa + D^2/(2n kappa)-D-eta R
      >=c/kappa-n kappa/2-eta R.

Every retained weight is at most 1/4. Consequently

    h>=4 M_h>(c-4 eta)R-16-8/R=:g(R).                    (4)

The endpoint g(233)-26=2390619/27028000 is positive, and
g(R+1)-g(R)>c-4 eta>1/6, since c-4 eta-1/6=4913/348000.
Put j=floor((R-233)/6). Then

    h>=27+j, R<=238+6j, n<=477+12j.

In particular h>0 and n>2h. Set epsilon=(3/4)^(27+j).
By (3),

    4 epsilon V<=(517/2+6j)(3/4)^(27+j)<delta.             (5)

The last expression decreases with integer j>=0: its next/previous ratio
is (3/4)(517/2+6j+6)/(517/2+6j)<1. At j=0 the strict endpoint is exactly
25850*3^27<11*4^27. Thus (5) holds at every order, without numerical
extrapolation.

## 3. Deletion keeps its negative baseline term

Input 2 applied to the retained matching yields

    (kappa V-B)(R-L)>=4(1-epsilon)(V-D-eta L).

Here a_h>=4(1-epsilon)/(R-L), and the last mass bound is positive by (4).
Expanding and using kappa R<=4 gives

    L[kappa V-B-4 eta(1-epsilon)]
      <=[kappa R-4+4 epsilon]V+4(1-epsilon)D-BR
      <=4 epsilon V+4D-BR.

By (1), (5), and L>=1, this implies

    L Astar<4D-BR.                                        (6)

In particular D>0. The negative term -BR is essential to the improvement.

## 4. Spend each group contribution once

Write W_h,W_l for the sums of the disjoint heavy and light group wedge
contributions; W_h+W_l<=W<=B. The input 3 bounds and Cauchy--Schwarz give

    D_h<=W_h/eta,
    D_l<=2 sqrt(L W_l)<=W_l/eta+eta L.

Thus D<=B/eta+eta L. Combining this with (6),

    D(1-4 eta/Astar)<B(1/eta-eta R/Astar).                 (7)

The left side is positive: Astar-4 eta=211/116000>0.
If the right coefficient is nonpositive there is already a contradiction.
Otherwise substitute (2), divide by D>0, and rearrange:

    1<kappa/eta + eta(4-kappa R)/Astar.                    (8)

But kappa/eta<=232/233, while

    0<4-kappa R
       =4[(3/4)^R+1_{n even}(1/4)^R]
       <=8(3/4)^233<1/1000,
    eta/Astar<1.

The last two comparisons use 8000*3^233<4^233 and
Astar-eta=6211/116000. Hence the right side of (8) is less than
232/233+1/1000<1, the final contradiction. Therefore L=0.

## 5. Closure and equality

All selected weights now satisfy w_i>eta>kappa. If D>0, some S_i>0 and
B>=W>=sum_i w_i S_i>kappa D, contradicting (2). Therefore D=B=0.
The uncertainty graph is a matching and the universal inequality follows.

Every supported binary completion has Delta(T)=0 and is transitive at
n>=4. The accepted [H=1 classification and scalar minimization](../session138_stability_and_actual_contexts_20261102/stability/H1_STABILITY.md)
then gives precisely the maximal fair-pair equality configurations stated
above. Those configurations attain kappa, completing the sharpness proof.

## Verification and limitation

`python -B verify_stability.py` checks every rational endpoint and both
integer exponential comparisons, the deletion identity independently by
polynomial arithmetic, and exact fair-pair/triangle controls. Its finite
checks are supplements: the inequalities above establish the infinite
range. The threshold 466 is not an obstruction at 465; the chosen constants
simply stop certifying there. The unresolved range requires new work.
