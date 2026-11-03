# Incoming native branches, cubic norm shells, and live-port precedence

Read through `031b85ab30`. **PROVED scoped transport below; incoming results
attributed, not rediscovered.**

## 1. A literal two-wing version of the incoming tribonacci bridge

Incoming [NATIVE_BRANCHES](../six_roots_hidden_bits_20261103/NATIVE_BRANCHES.md)
constructs all three65 branches with one shared output port and proves the
all-size native map `E_n -> TT2 -> S_(n-2)`, with S the March interval
tournament. Our two-wing source has the following literal companion identity,
for every m,n>=1:

    native(T_(m,1,1,n), natural path)
           = S_(m+1) -> TT2 -> S_(n+1).              (1)

The left factor is block0 plus singleton block1; the middle consists of
blocks2,3; the right is singleton block4 plus block5. Inside either exterior
factor, precisely the consecutive arcs survive forward. Between factors,
each arc was either a backward quotient arc that flips or one of the two
protected boundary arcs1->2 and3->4. This checks every pair and proves (1).

The native intermediate therefore has H=`h_(m+1) h_(n+1)`. At m=3,n=3 this
explains the25-path intermediate of the485 escape. All old vertices and path
incidences survive. When a parameter is1, its exterior S2 is not strong;
the ordinal-join identity and product formula still hold. `incoming.py`
also checks all16 cases1<=m,n<=4 literally.

Both native constructions expose actual interval factors in their original
frames. They do not preserve dyadic growth or prove a common reset controller.

## 2. The escaped scalar streams occupy different norm shells

Incoming [Session167 norm theorem](../session167_recurring_structures_20261103/norm/NOTE.md)
uses `C(s,u,d)=(u,s+u+d,s)` and a cubic integral norm N invariant under C.
Its positive norm shells have finite descent-root catalogues; its primitive
norm-two shell is one chain. Our H recurrence has the exact observation map

    (H_n,H_(n+1),H_(n+2))^T
       = [[1,1,0],[1,2,1],[3,3,2]] (s,u,d)^T.        (2)

The matrix has determinant2. Its inverse requires

    d=(H_(n+2)-3H_n)/2,
    u=H_(n+1)-H_n-d, s=H_n-u.

The actual odd H values satisfy this parity guard. At n=2:

| Family | State | Descent root | Invariant norm |
|---|---|---|---:|
|A|(97,194,0)|(97,0,0)|2*97^3=1825346|
|B|(85,170,62)|(0,23,39)|85838|

The B descent passes through(62,85,23),(23,62,0). These roots and norms
classify the **scalar streams** within that arithmetic dynamical system.
The H3007 collision does not identify their recurrence orbits. The matrix
intertwining, norm identity, and exact descents are in `incoming.py`.

Equation(2) does NOT prove that these numbers are the guard profile of a
chosen vertex in A_n or B_n. It transports three observed H values, not
tournament incidence. Even matching the canonical norm shell supplies an
arithmetic recognition test, not permission for a class-valued replacement.

## 3. Retained vertices and live ports remain different

Incoming [joint-port geometry](../joint_port_geometry_20261103/JOINT_EFFECTS.md)
corrects an earlier C3/C5 legalizer claim: the old second vertex survives,
but its live port is replaced. NR uses a different policy. Our four new
branch certificates retain old root2 through their displayed complement/birth.
Subsequent multiplier continuations follow the inherited explicit
port-replacement traces; root2 is not claimed to remain their live guard.

The incoming [incidence synthesis](../joint_port_geometry_20261103/INCOMING.md)
also distinguishes H0 reflection fibers from H1 edge-voltage classes. Our
three-sheet clock has nontrivial transport around its specified base cycle,
the H1 position. It is not a fiber count from measuring independent vertex
reflections. Tournament realization remains separate in each construction.
