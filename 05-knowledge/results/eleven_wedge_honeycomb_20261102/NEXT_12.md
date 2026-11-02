# Next even order: a different payment problem

Status: **OPEN**, with **FINITE-EXACT hostile probes**, not order12 stability.

The next candidate sharp coefficient is561/1024. Six ordinal fair pairs
have H=1, u=(3/4)^6, |det A|=(1/4)^6 and V=3/2, giving

    [1-(3/4)^6-(1/4)^6]/(3/2)=561/1024.

The perfect-matching term must be retained. Reusing the odd-order formula
would incorrectly give3367/6144. The worse signed full-cover matching group
is

    P12=V-M22+M222-M2222+M22222-2M222222.

A cheap next graph obligation, with no assumed wedge reserve, is

    P12+P21 >= (561/1024)V.                           (O12)

At q endpoints, multiplying the difference by4096 gives

    J12=463m-256m2+64m3-16m4+4m5-2m6+64p21.

This is separately affine, so a proof for ALL twelve-vertex simple graphs
would prove(O12) on the variance box. But(O12) alone would still leave the
order12 directed-cycle ledger and its payments to prove. The order11
ledger does not automatically lift across an even full-cover boundary.

`next12_probe.py` checks all12,346 eight-vertex cores, augmented with
k=0,1,2 disjoint fair edges and enough isolated vertices to have12 vertices.
There are37,038 **cases**, not necessarily distinct isomorphism classes.
The minimum J12 is respectively0,319,0. It also checks all77 disjoint-clique
supports of order12, with minimum0 only for the empty graph and6K2.
All values are exact integers; normal and optimized replays agree.

These tests include the previous P3+4K2 hostile padded by an isolated vertex.
They do not cover arbitrary twelve-vertex supports. The first decisive
next work is either a complete graph certificate for(O12), or a new hostile
showing precisely which directed-cycle reserve it needs. A graph failure
would refute this payment route, not the weighted-tournament conjecture.

For a one-vertex augmentation of an eleven-vertex base G, the same incidence
expansion gives nonnegative quadratic terms64Qij, with
Qij=m-d_i-d_j+Aij. Thus the negative-linear-coordinate pruning still applies,
but some quadratic coefficients can vanish and the base universe is much
larger. Keep zero coordinates and equality analysis honest; the old56
strict-positivity margin is not available here.
