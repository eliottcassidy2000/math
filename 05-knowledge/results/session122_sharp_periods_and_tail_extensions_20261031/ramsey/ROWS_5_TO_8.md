# Complete capacity rows five through eight

2026-10-31. **PROVED with exact finite certificates and an analytic tail cutoff.**
The class D and function F are exactly those in `DYADIC_PERIOD.md`: no empty
C5 slots, and no claim about unrestricted Ramsey numbers.

## Exact formulas

For every positive integer B, write B=8t+r with 0<=r<8. Then

    F(5,B)= 62t + (0,5,12,19,29,35,43,51)[r],
    F(6,B)= 80t + (0,6,15,23,35,43,55,65)[r],
    F(7,B)= 98t + (0,7,17,27,41,51,65,79)[r],
    F(8,B)=125t + (0,8,20,32,50,62,80,98)[r].

Each formula holds from B=1, so its least eventual-translation threshold is
1. Each has least period 8, also following the all-A theorem in
`DYADIC_PERIOD.md`. Complementation supplies the corresponding columns.
These extend the inherited exact first four rows; they are not fitted tails.

## Finite-cutoff theorem behind the certificates

Fix A and a proposed H(B)=sL floor(B/L)+h_(B modL), where s=w(A), L is
positive even, h_0=0, sL integral, and D=max_(0<=r<L)(sr-h_r)>0.
Assume all smaller clique rows are already known exactly. The following
finite checks certify H as the exact row for every B.

First list all coordinatewise maximal positive integer clique patterns
p=(a_0,...,a_4) with a_i+a_(i+1)<=A. Coordinatewise maximality is equivalent
to every vertex having at least one tight incident edge. Any other positive
pattern can be increased to a maximal one, and the child capacities are
monotone in clique budget. Thus checking these patterns loses no root.

Reorder positions as (0,2,4,1,3), so independence constraints are now
b_i+b_(i+1)<=B. Put c_i=w(a_i) in this order. The maximum slope of the
continuous port polytope is

    m=max(0, max_i c_i, max_i(c_i+c_(i+2)), (sum_i c_i)/2).

This is the complete vertex list, not merely a chosen dual bound.

### Nonsharp C5 patterns

If the half-vector slope is below s, the strict-superadditivity lemma makes
all other vertex slopes strictly below s too. Put delta=s-m>0.
A counterexample n>H(B) would satisfy

    delta B <= sB-n < sB-H(B) <= D.

Consequently B<=ceil(D/delta)-1. Exhaustively check the exact root objective
up to this explicit cutoff. This argument does not require the half-vector
dual to be positive for nonsharp patterns; some row8 patterns have negative
half-vector dual coordinates and are correctly handled by the full vertex max.

### Sharp C5 patterns

Suppose sum c_i=2s. Solve the cycle dual

    y_i+y_(i-1)=c_i,
    y_i=(c_i-c_(i-1)+c_(i-2)-c_(i-3)+c_(i-4))/2.

The certificate verifies every y_i>0. For port slacks
s_i=B-b_i-b_(i+1)>=0, the exact linear identity is

    sum_i c_i b_i = sB - sum_i y_i s_i.

A counterexample therefore has sum_i y_i s_i<D. Writing y_min=min y_i,
its total integer slack obeys

    sum_i s_i <= S := ceil(D/y_min)-1.                       (1)

The invertible odd-cycle system gives 2b_i=B-t_i, where t_i is an alternating
sum of the five slacks. In particular |t_i|<=S. The certificate checks that
every child row has an exact translation identity at EVERY positive argument
(with threshold 1), with period dividing L/2. Therefore increasing
B by L and every b_i by L/2 preserves all slacks, changes the objective by
(L/2)sum c_i=sL, and changes H by the same amount.

Given a violating port vector, reduce B in its residue class modulo L to the
smallest positive B' satisfying B'>=2+max_i t_i. It exists, satisfies B'<=B,
and has

    1 <= B' <= S+L+1.                                      (2)

Because L is even, all reconstructed b'_i=(B'-t_i)/2 remain integers and
positive. Exact child translation preserves the violation in either
direction along positive arguments. Thus it suffices to check each sharp
pattern through (2). This is a proof of the infinite tail, not an inference
from observing a finite plateau.

### Joins, unions, and lower witnesses

At a join splitting A into x+(A-x), strict superadditivity gives positive
gap delta=s-w(x)-w(A-x). The same counterexample inequality bounds B by
ceil(D/delta)-1, so exact join checks through that cutoff suffice. Singleton
roots are trivially bounded by H(B); the supplied seed values are positive
and H increases by sL between residues.

For disjoint unions, check the L^2 residue inequalities

    h_r+h_t <= sL floor((r+t)/L)+h_((r+t) modL).

They imply full superadditivity of H. Structural induction on the graph
construction now proves the upper bound F(A,B)<=H(B).

For the lower bound, the certificate supplies realizable root constructions
for B=1,...,L, with recursively known smaller-row children or smaller-B
union children. Their orders are exactly H(B). Disjoint union with copies
of the seed at B=L extends these witnesses to every B. No empty child is
introduced and no scalar LP point is substituted for an actual graph.

## Executable certificate and independent controls

`check_rows.py` uses only the Python standard library. It enumerates every
maximal positive clique pattern for A=5,6,7,8, checks all slope gaps and
sharp duals exactly, computes the analytic cutoffs, and verifies every
finite root upper bound through them. The root optimizer is a max-plus path
DP with fixed first port and an explicitly checked closing edge, retaining
all positive independence ports. It also checks all 64 union residues per
row and all 32 seed construction statistics.

The output lists every pattern, dual or gap, cutoff, finite upper result,
and seed construction. Full (not maximal-only) clique-pattern checks for
2<=A<=6 and 2<=B<=6 supply a separate pruning control. Positive strictness
checks cover 16,384 pairs; row6's integral slope with a nonzero residue
penalty is a denominator-only hostile control. The all-A dyadic proof is
analytic and is not inferred from this finite range.

Reproduce:

    python3 -B check_rows.py > rows.normal.json
    python3 -B -O check_rows.py > rows.optimized.json
    cmp rows.normal.json rows.optimized.json

The exact first four rows used as the induction base are G3 in
`ramsey_neighbourhood_lifts_20261030/GRAPH_CAPACITY.md`. The known sharp
slope and polytope vertex list are its G2. There is no claim here to have
implemented the general Presburger engine, nor to have computed all later
nonzero residue laws.
