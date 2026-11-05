# Even-path capacity expansion and the circulation quotient

**PROVED derivation from OCF, with a separate finite verifier; independent audit pending.**
This continuation does not modify the frozen first packet.

For ordered distinct i,j in U define

    P_ij = sum_{directed simple paths p:i to j,
                 |V(p)| positive and even} H(U-V(p)).        (1)

Path order is retained and paths contain no repeated vertex. The shortest
possible path has two vertices. Let c_ij be the first packet's capacity.
Then

    c_ij/2 = P_ij+P_ji.                                     (2)

## Proof and topological interpretation

Remove the new root x from an odd cycle through x in the marked OCF formula
(THM-070). It leaves a directed even-vertex path with first endpoint in S
and last endpoint outside S. Consequently

    Delta_U(S)/2 = sum_{i in S,j outside S}P_ij.             (3)

Both sides are quadratic multilinear functions of the cut indicators.
The coefficient of x_i x_j in (3), compared with the first packet's
Delta_U(S)=sum_{i->j}c_ij x_i(1-x_j), is exactly (2).

There is also a direct cycle-space explanation for the oriented gauge.
On each old arc i->j put q_ij=P_ji. The reverse endpoint path j to i,
followed by the old arc i->j, is a directed EVEN cycle. Conversely every
such cycle C contributes H(U-V(C)) to q on each of its arcs. Hence

    q = sum_{even directed cycles C} H(U-V(C)) 1_C.          (4)

Cycles are counted modulo cyclic rotation, as usual. Thus q is a
nonnegative circulation, with zero divergence at every vertex. Move the
reverse P_ji weight onto the forward old arc i->j, replacing the pair
(P_ij,P_ji) by (P_ij+P_ji,0). The change in every outgoing cut sum is
q(delta^+(S))-q(delta^-(S))=0. Equation (3) therefore becomes the positive
old-tournament cut network of the first packet.

This is an exact quotient by the circulation space. It preserves every
cut value and destroys the path/cycle decomposition, which must be kept
for more than one birth or for structural growth arguments. No equality
of the complete path tensor is asserted.

## Complete positivity of the capacity presence cube

Expand each H(U-V(p)) in (1) by OCF. Each contribution consists of:

- a directed even-vertex path with endpoints {i,j};
- a collection of disjoint directed odd cycles avoiding that path;
- weight 2^(number of cycles), with an overall factor 2 for c_ij.

If U[A] is an induced subtournament containing i,j, exactly the structures
whose whole vertex support lies in A survive. Therefore the capacity
function on the presence cube A has nonnegative Boolean Möbius atoms and
ALL its mixed presence differences are nonnegative. In particular

    c_ij(U[A]) <= c_ij(U[B]) whenever {i,j} subset A subset B.  (5)

The length-two path alone proves c_ij>=2H(U-{i,j}); this independently
explains the block-insertion lower bound in the first packet.

Under a fixed induced enlargement, old-coordinate zonotope generators
therefore grow. This does not make the next full zonotope fixed: new
vertices create new generators and new possible cut states. For a cut S
of a larger parent W containing an induced U, one immediate consequence is

    Delta_W(S)>=Delta_U(S intersection V(U)).               (6)

When the restricted cut is trivial this gives only zero. One must not use
(6) to assume every low-energy cut keeps a nontrivial old restriction.

## Fewer flow certificates for strong parents

For a strong parent U and any fixed vertex r, the minimum nontrivial cut
is the minimum among the 2(n-1) terminal-pair problems (r,v) and (v,r),
v!=r. Every nontrivial cut separates r from some v, in one of those two
orientations. Thus this smaller family already covers all cuts for purposes
of the minimum. The first packet's n(n-1) checks remain valid and provide
extra redundancy. For weak parents retain the source-SCC by sink-SCC
constraint of that packet; the unrestricted minimum would otherwise be zero.

## Verification scope

`verify_paths.py` independently enumerates all even-vertex simple paths and
even cycles on every labelled parent of orders 2 through 5. It checks (2),
(4), flow conservation, the block lower bound and the 2(n-1) minimum cover.
It does not import a cycle bank or a precomputed path polynomial.
