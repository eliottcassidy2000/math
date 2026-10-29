# A native exterior payment coordinate and an exact orientation selector

**Status: PROVED scoped elementary constructions and inherited obstructions;
FINITE-EXACT controls. Universal Collatz payment, six-colour forcing and
the planar Jacobian conjecture remain OPEN.**

## 1. Four coordinates, six signed edges, one actual Collatz debt

The earlier tournament counter has `Lambda^2 C` integrally conjugate to
`C-I`; its cubic spectrum and modulo-seven law are inherited from
[WEDGE](../atomic_descent_wedge_20261029/WEDGE.md). They are statements
about that counter, not a Collatz spectral estimate. The six-state
legality machine is not the same exterior representation; the actual
native two-run antisymmetric sector supplies a different genuine wedge.

Here is a direct exterior representation of the **native Collatz map**.
Fix the original positive odd source n0. After J odd steps and D divisions
write

    N=3^J*n0+C,   Q=2^D,   current integer n=N/Q.

In Z^4 put

    u=(n0,1,0,0),   v=(0,0,N,Q),   p=u wedge v.

Index the six coordinates by the edges 01,02,03,12,13,23 of K4. They are

    p=(0,n0*N,n0*Q,N,Q,0),
    p01*p23-p02*p13+p03*p12=0.                              (1)

The exact native valuation a acts by

    G_a=diag(I_2, [[3,1],[0,2^a]]),
    p'=(Lambda^2 G_a)p.                                    (2)

Every coordinate remains integral. The marked state includes J and D;
Q itself retains D, but no claim is made that p alone recovers J.
The source n0 is recoverable as p03/p13, so this construction has not
replaced the input by an unlabeled endpoint.

The payment functional is the difference of two opposite K4 edge slots:

    L=p03-p12=n0*Q-N=Q(n0-n).                              (3)

Consequently L>0 means exactly that the trajectory has returned below
the original source. Its exact update is

    L'=3L+Q*((2^a-3)n0-1).                                 (4)

For a=1 the increment relative to 3L is `-Q(n0+1)`; for a=2 it is
`Q(n0-1)`. The sign is not forced by the Pluecker relation (1), by the
dimension six, or by an average over phases. Resetting n0 to the current
endpoint replaces the actual account; projectivizing away Q without
retaining D loses the integer division clock. Neither is a payment proof.

**An exact normalization is nevertheless available.** Since N=Q*n at
every native integer source, all six coordinates have common divisor Q,
and p13=Q shows it is exactly the gcd. Thus

    pbar=p/Q=(0,n0*n,n0,n,1,0),
    pbar'=(Lambda^2 G_a)pbar / 2^a                            (5)

is integral on the exact valuation guard. Indeed its only nonconstant
numerators are `3n+1` and `n0*(3n+1)`, both divisible by 2^a; the
remaining nonzero coordinates are multiplied by 2^a. Store J and D
as separate integer counters, updating J by1 and D by a. This gives an
exact source- and clock-preserving projective gauge, with payment
`Lbar=n0-n` and

    Lbar'=[3Lbar+(2^a-3)n0-1]/2^a.

The construction compresses a redundant scale, not the mathematical
debt. It does not apply arbitrary rational matrices to all of Z^6:
integrality depends on the actual native guard. All divisions, normalized
coordinates and integer clocks are checked in the stated native universe.

This was sharpened after reading incoming Session105's
[integer projective recurrence](../session105_metagraph_recursion_20261029/enumeration/INTEGER_PROJECTIVE.md).
There, a combinatorial integer gauge pays negative binary shifts even
after histories merge; here the native valuation pays the division in
(5). The shared proved mechanism is exact divisibility under a specified
operation and retained sidecar, not an identification of the two systems.

**Normalization does not commute with an unmarked weighted mixture.**
The two genuine native one-step trajectories

    45 ->17 with a=3,     59 ->89 with a=1

have the same odd-step clock J=1 and the same original Elias weight1/512
(source indices23 and30). Nevertheless their weighted unnormalized
payment sum is

    [8*(45-17)+2*(59-89)]/512=41/128>0,

while the weighted original-integer comparison is

    [(45-17)+(59-89)]/512=-1/256<0.                         (5a)

Thus an aggregate positive scaled score need not even have the sign of
the aggregate actual debt. Pointwise division preserves each source's
sign, but varying division clocks reweight a mixture. Keep each Q (or D)
through normalization and apply the original source weight afterward;
one common scalar cannot repair an already merged unmarked sum. These
are native one-step paths, not claimed completed auxiliary queue returns.
The hostile is especially relevant to ROOT cancellation: the original
weights, source labels and common odd clock are all authenticated.

Thus the K4 connection is an exact map with a preserved predicate:

* Source: `(n0,J,D,N,Q)` with exact native valuations.
* Target: the marked decomposable bivector (1), acted on by (2).
* Preserved: the current integer, original source and sign of repayment.
* Missing in an unmarked quotient: the odd-step clock and possibly Q.
* Cheapest test: literal native execution versus (2)--(4), including
  both positive and negative values of L. The checker performs twelve
  steps from each of the two hundred odd sources 1,...,399.

This is not a theorem identifying Collatz with tournament enumeration.
It supplies an exact carrier on which a proposed cancellation argument
can be required to preserve the payment functional.

## 2. A rational cross-edge chooses the geometric branches

Fix four distinct boundary ports

    A=(-3/5,4/5),   B=(-3/5,-4/5),
    C=(8/5,4/5),    D=(8/5,-4/5).

The common unit-distance neighbours of A,B are
`X0=(0,0), X1=(-6/5,0)`. Those of C,D are
`Y0=(1,0), Y1=(11/5,0)`. Of the four pairs Xi,Yj, exactly X0,Y0 have
distance one: their possible distances are 1,11/5,11/5,17/5.
Thus the cross-edge XY selects a unique pair of contact orientations.
All coordinates are rational, hence also in Q(zeta55), but these are
newly specified ports, not a relocation of a previously pinned boundary.

At the selected realization the complete unit-distance graph has exactly

    AX, BX, CY, DY, XY.

It is a six-vertex tree. No accidental unit edge is omitted.

Use half the derivatives of squared distances as the constraint rows.
With boundary velocities fixed at zero, the pinned matrix is

    [ 3/5 -4/5    0    0 ]
    [ 3/5  4/5    0    0 ]
    [   0    0 -3/5 -4/5 ]
    [   0    0 -3/5  4/5 ]
    [  -1    0    1    0 ].

It has rank four and unique left stress, up to scale,

    (5,5,5,5,6).                                            (6)

Therefore a requested infinitesimal change b in the five constraints
has a solution iff

    5(b1+b2+b3+b4)+6b5=0.                                   (7)

This is a first-order compatibility condition, not a nonlinear lifting
theorem. One extra scalar equation can distinguish four discrete root
pairs even while its differential row is dependent. Discrete orientation
information and differential rank are different quantities.

**A reusable rational family.** Replace3/5 by a and4/5 by b, where
`a,b>0` and `a^2+b^2=1`; use ports `(-a,+/-b)` and `(1+a,+/-b)`.
The candidate roots are `0,-2a` on the horizontal axis and `1,1+2a`
on that axis. Again only the pair0,1 has unit separation. The pinned
stress becomes `(1,1,1,1,2a)`. Every rational parameter

    0<t<1,   a=(1-t^2)/(1+t^2),   b=2t/(1+t^2)

gives the same faithful six-vertex tree. The only possible extra
port-pair unit edge would require b=1/2, impossible for rational a,b
on this circle. Thus this supplies infinitely many exact rational
orientation selectors, all with the same six-colour non-obstruction.

The nonselected root-pair distances are `1+2a,1+2a,1+4a`, so the
smallest selection gap is Delta=2a. The determinant squared of a
two-contact block is

    4a^2 b^2=Delta^2*(1-Delta^2/4).

As a approaches0, the discrete separation and the differential
conditioning both deteriorate; as a approaches1 the ports coalesce
even though the selection gap stays large. This exact tradeoff is
useful when composing selectors: a unique branch alone supplies no
uniform inverse bound. The checker covers all21 reduced rational
parameters in(0,1) with denominator at most8.

## 3. What the same gadget does and does not force in colour space

For every prescribed six-colouring of A,B,C,D, each interior vertex has
at least four available colours. If its lists are L_X,L_Y, the number
of extensions is

    |L_X|*|L_Y|-|L_X intersect L_Y| >=12.

Equality occurs when both available lists are the same four colours.
All 6^4 boundary assignments and all 6^2 interior pairs are checked.
So the geometry selector imposes **no six-colour boundary obstruction**.
It does obstruct the three-colour boundary `(1,2,1,2)`: both X and Y
would be forced to colour 0, violating XY. The unpinned graph itself
is two-colourable.

This separates three predicates: selecting coordinate roots, extending
a specified boundary colouring, and the chromatic number of the whole
graph. The inherited
[six-colour interface profile](../six_colour_interface_profile_20261015f.md)
already makes the second predicate the correct one for forcing.

There is a useful exact specialization for larger contact systems.
Suppose h interior vertices each have two fixed port contacts with
linearly independent unit normals, and there are e authenticated
cross-edges between interiors. The 2h port-contact rows form an invertible
block-diagonal matrix B. Appending cross-edge rows C leaves column rank
2h, so the pinned stress space has dimension exactly e. Each row c of C
has a stress with coordinates `(-c B^(-1), unit_c)`.

If port distance is d, the determinant squared of its two unit-normal
rows is `d^2*(4-d^2)/4`. Thus the inverse carries a genuine conditioning
cost near d=0 or d=2; exact tangency violates the hypothesis.

For six colours the residual lists have size at least four. A minimal
nonextendible induced interior subproblem must have every vertex of
cross-edge degree at least four; otherwise delete a vertex, extend the
rest, and then colour it. It needs at least 2h cross-edges on that minimal
subproblem, hence at least 2h independent pinned stress coordinates in
this two-contact model. A three-degenerate cross-edge graph always
extends every such boundary assignment.

This is a **specialization**, not an improvement, of the stronger
stress/list bounds already in
[six-colour games, Sections 2--5 and 7](../strong_blue_games_20261020f/six_colour.md).
The current gadget has h=2,e=1: it is far below the necessary stress
budget for a six-colour obstruction. Repeated geometric branch selection
alone cannot be counted as accumulating colour force.

## 4. The planar Jacobian comparison retains the return operator

For an actual polynomial Keller pair `Jac(f,g)=1`, the inherited
[inverse-clock theorem](../keller_inverse_clock_20261015f.md) gives

    D=f_x partial_y-f_y partial_x,
    E=g_y partial_x-g_x partial_y,
    M=K[x,y]/D(K[x,y])        (an additive quotient),
    E^(-1)[r]=[R_f r],  R_f r=f_x integral_0^x r(s,y) ds.

On a fixed compact polydisc the displayed representatives R_f^m r decay
factorially, while `[E^m R_f^m r]=[r]`. The recovery operator must be
charged before this becomes a statement about the original class.
The elementary regular-singular control

    nabla=d/dz+1/(2z),
    nabla^(-m)1=z^m/(3/2)_m,
    nabla^m z^m=(3/2)_m

shows exactly reciprocal decay and recovery. It is a control example,
not a Keller counterexample.

The [genus-two transfer-word classification](../genus_two_transfer_words_20261011.md)
and [quartic/sextic source anchor exclusion](../shape_ports_20261025/anchors_quartic_sextic.md)
remain stronger inherited restrictions than any claim that rationalizing
a comparison curve alone produces a polynomial Keller pair. The latter
retains the actual scalar Jacobian, not just a rational parameter.

The common structure used here is narrowly specified: after a change of
representation, retain the functional or operator that returns to the
original obligation. In Collatz it is L and the clock; in the geometric
selector it is the pinned stress equation and boundary list relation;
in the Keller module it is E^m and the original additive class. There is
no asserted equivalence among the three conjectures, and no common
numeric dimension is used as a substitute for these maps.
