# Changing feedback inventory by four black moves

Status: **PROVED**, all orders. Global black-core connectivity remains
**OPEN**; the theorem covers the entire active three-atom SCC language.

## Inheritance and the missing operation

The closest mechanism is the fixed-inventory normalization in
[`scc_basins_20261020d/arithmetic_interior.md`](../scc_basins_20261020d/arithmetic_interior.md),
Section 6. It connects active words in `{1,C3,R4}` with fixed vertex
number and feedback-arc number. The missing coordinate was that inventory:
R4 dissolution preserves it at both endpoints. The hostile is a pure-blue
palindromic rigid word; such a vertex cannot lie on a black path. The
least-used sidecar is the score of a *nonendpoint* pair exchanged by a
putative path-reversal anti-automorphism.

The inherited seven-vertex certificate in
[`black_descent_20261020c/scc_probe.json`](../black_descent_20261020c/scc_probe.json)
suggested a three-move local operation. The even four-move operation below
restores an arbitrary exterior literally and changes the inventory.

Our graph is the converse-merged Hamiltonian-path-complement line
multigraph. Reversing all arcs of a directed Hamiltonian path P, denoted
R_P, gives the same merged line. It is blue exactly when positional
reversal of P is an anti-automorphism. In particular every pair of
opposite positions on a blue path has outdegrees summing to n-1.

## A protected triangle creation

Let A be any nonempty tournament having a source vertex u, and choose
any Hamiltonian path p of A. Necessarily p starts at u; write v for its
last vertex and m=|A|. The case m=1, u=v, is included. Let x,y,z be three
successive singleton SCCs, and let a,b,c induce the directed triangle
a→b→c→a. Start at

    T = A → x → y → z → C3(a,b,c).

Use these four path reversals in order:

    P1 = (p, x,y,z, b,c,a),
    P2 = (a,b, z,y,x, reverse(p), c),
    P3 = (x,b,a,c, p, y,z),
    P4 = (c,a,b, x,z,y, reverse(p)).

**Theorem.** All four are directed black Hamiltonian paths at the time
they are used, and the literal final adjacency matrix is

    A → C3(x,y,z) → C3(a,c,b).

Every internal arc of A and every cross arc is restored. The first new
triangle has x→y→z→x. The old guard triangle has all three arcs reversed,
so it remains a triangle of the same isomorphism type.

### Directedness

P1 follows the source path and a cyclic rotation of the guard. After P1,
the links a→b and b→z are respectively an untouched triangle edge and a
reversed boundary edge; z,y,x,reverse(p) follow reversed P1 edges, and
u→c is an untouched cross arc. These give P2.

After P2, x→b is untouched, b→a is the reversed first edge of P2, a→c
is the reversed last edge of P1, and c→u is the reversed last edge of
P2. The internal edges of p have been reversed twice, v→y is untouched,
and y→z has been reversed twice. Thus P3 is directed.

After P3, c→a, a→b, and b→x are the reversed guard links, x→z is
untouched, z→y is reversed, and y→v is reversed. The internal path p has
now been reversed three times, giving the final reverse(p). Thus P4 is
directed. All statements remain valid when u=v.

### Exact edge ledger

Every internal edge of p is reversed four times. The boundary edges
vx, uc, and vy, and the extra edge xb, occur twice. The edge xy occurs
twice and yz four times; zb occurs twice. What remains modulo two is
exactly

    {xz, ab, bc, ca}.

This creates the first triangle, reverses the guard, and changes nothing
else. In particular this proof does not collapse A to a vertex or assume
anything about its nonpath edges.

### Blackness without a size exception

Path reversal lowers the first vertex's outdegree by one, raises the
last vertex's by one, and changes no other outdegree. Here n=m+6.

* P1 has endpoint-degree sum m+6, since u is a source of T and a has
  degree one. It differs from n-1=m+5.
* At P2 the endpoint degrees of a,c are 2,1, with sum 3<n-1.
* At P3 the second vertex b has degree 1, while its opposite-position
  partner y has degree 4. Their sum is 5<n-1. This also handles m=3,
  where the *endpoint* sum happens to equal n-1.
* At P4 the endpoint degrees of c,u are 2,m+4, with sum m+6≠n-1.

Each displayed imbalance excludes the required anti-automorphism.
The nonempty-context assumption is substantive: m=0 removes the strict
inequality used at P3, and erasing the context from the paths invalidates
their boundary ledger. The proposed target C3→C3 is pure blue, so no
all-black route to that target could exist.

## All feedback inventories in the golden language connect

**Corollary (PROVED).** For each n≥4, every black-active ordered SCC word
over `{1,C3,R4}` lies in one component of the ambient black metagraph.
In particular every such word can reach the merged type

    1^(n-3) → C3.

This is an ambient connectivity statement. Intermediate vertices need
not have SCCs in the three-atom alphabet, and the induced subgraph on
nonstrong terminal words has no edges.

The inherited fixed-inventory theorem normalizes a word with feedback
number f to

    N_(n,f)=1^(n-3f) → C3^f,
    1≤f≤floor((n-1)/3).

For any adjacent admissible values f=t and f=t+1, set

    r=n-3(t+1)≥1,     A=1^r → C3^(t-1).

A has a source. The new theorem joins

    A → 1 → 1 → 1 → C3
    to A → C3 → C3 = N_(n,t+1).

The first endpoint has feedback number t and belongs to the preceding
fixed-inventory component. Both endpoints start with a singleton and
end with a triangle; hence neither is self-converse, and both are active.
This joins every adjacent inventory, proving the corollary. For n=4,5,6
there is just one nonempty inventory, already covered by the inherited
theorem. For n≤3 there are no active words in this language.

The new macro changes the feedback-arc number by exactly one, and the
directed-triangle number by exactly one, at its endpoints. Feedback count
and triangle count need not be monotone or fixed along its intermediate
states. In the Fibonacci matching model this operation creates one
matching edge, while the guard triangle retains one edge.

## Exact controls and remaining obstruction

`inventory.py` checks the literal four-path identity and every blackness
inequality on all standard-Hamiltonian-path presentations A of orders
1–7 with a source. It also independently constructs the finite abstract
word graph through order 21 using the inherited active replacement and
swap rules and the new inventory edges. These finite controls accompany
the all-order proof; they are not its basis.

`inventory_probe.py` records the small discovery searches at orders 7–9.
The separately implemented `outer_inventory_audit.py` tests the identity
on all 1,100 fixed-path A presentations of orders 1–6, without imposing
a source, and verifies the distinct role of that blackness hypothesis.

Still open: normalize arbitrary larger strong SCC atoms, connect every
terminal basin, or find a genuine mixed-to-mixed blue bridge. This result
rules out feedback inventory as a separator *within the three-atom
language*, not as a global invariant that has been fully classified.
