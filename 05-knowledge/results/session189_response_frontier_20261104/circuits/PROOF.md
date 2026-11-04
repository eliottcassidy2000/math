# Constant-region ports, cycle budgets, and shielded circuit cores

**Status: proved in the explicit shared-DAG model below.** This packet extends
Session 188's constant-region and sequential fanout identities. Its new
accounting retains the connected-component partition of the outside gate
ports. It proves a stronger *consumed topology* budget, and then proves why
that budget is not extra gate deletion. A universal shielding construction
preserves an arbitrary permitted core, including unbounded grid cores,
through arbitrarily many tight two-gate restriction steps.

No unrestricted circuit lower bound, minimum-circuit theorem, elimination
algorithm certification, or literature-priority claim is made.

## 1. Inheritance, model, and the exact sidecar

Pinned sources, both at commit
`8f14b96ab020f7d5ceeb88e7bda413bb6f8bc9b7` in
`eliottcassidy2000/math`:

- `05-knowledge/results/session188_structural_bridges_20261104/circuits/PROOF.md`:
  the certified deletion identity `|C|+|B|=2r+f+E`.
- `05-knowledge/results/session188_structural_bridges_20261104/circuits/SEQUENTIAL_BUDGET.md`:
  the potential `Phi=sum_g max(d(g)-1,0)`, its monotonicity, and the
  alias-ladder hostile to repeated boundary fanout charges.

Use a finite acyclic single-output circuit with at most two pins per gate.
Pins read primary inputs, allowed boundary values, or earlier gates. Count
one gate edge per distinct parent/child pair, not per pin occurrence. Every
old gate lies on a directed path to the output. A legal restriction and any
semantic constant/literal claims require separate certification. Constants
and signed literals at the output are free for the deletion statements.

The allowed simplifications introduce no gate: substitute certified
constants; replace a gate by a signed literal of one of its current input-pin
sources; redirect its users and absorb signs into their truth tables; and
remove dead gates. Input restrictions may substitute constants or literals
of other *primary inputs*, not uncharged new internal computations. Arbitrary
resynthesis and aliases to nonparents are outside the minor assertion unless
a separate valid minor certificate is supplied.

Let C be a nonempty set of certified constant gates excluding the output,
closed under gates whose every input pin reads C. Let B be the outside gates
with a parent in C. Thus every B gate has exactly one distinct C parent.
Write r and f for the numbers of vertices of C with respectively zero and
one C-parent, and

    E = sum_(g in C) (d(g)-1).

All these degrees refer to the current live graph before this step.

Let G be the undirected simple *gate-only* graph of the circuit. Primary
inputs and output references are not vertices. Let G-C mean vertex deletion,
retaining isolated outside gates. Define

    a = number of directed gate edges from outside C into C,
    kappa = number of connected components of G-C,
    Q = E + a - kappa + 1.

The output remains outside, so kappa>=1. The new sidecar is the partition of
outside gate ports by connected components of G-C, together with each
physical incidence. Scalar source/fanout counts alone discard this data.

## 2. Exact topological loss and its equality case

For an undirected finite graph H, including disconnected or empty graphs,
write beta(H)=|E(H)|-|V(H)|+components(H), the binary cycle-space dimension.
The empty graph has beta=0. This dimension formula and the spanning-forest
construction are classical; see Kavitha, Liebchen, Mehlhorn, Michail,
Rizzi, Ueckerdt, and Zweig, *Cycle Bases in Graphs: Characterization,
Algorithms, Complexity, and Applications*, Theorem 2.3:
https://people.mpi-inf.mpg.de/~mehlhorn/ftp/SurveyCycleBases.pdf .
We give all circuit-specific arguments here.

**Theorem 1 (exact constant-region port loss).**

    beta(G) - beta(G-C) = Q = E+a-kappa+1 >= E >= 0.       (1)

Moreover Q=E if and only if every component of G-C other than the output's
component has exactly one edge directed from that component into C, and
the output's component has no such edge.

**Proof.** G is connected because every gate reaches the output. Put c=|C|,
and let e_CC count edges with both endpoints in C. The total number of
edges leaving C is c+E. These comprise e_CC and the edges from C to outside.
Thus the number of edges incident with C is (c+E)+a. Deleting C removes c
vertices and those c+E+a edges, while changing the component count from 1
to kappa. Subtracting the two cycle-rank formulas gives (1)'s equality.

Each outside component K not containing the output has at least one edge
from K into C: start at any gate of K and follow a directed path to the
output; the first step leaving K must enter C. These kappa-1 edges are
distinct. Hence a>=kappa-1. Equality occurs exactly in the stated pattern.
Finally every C gate is live and is not the output, so d(g)>=1 and E>=0.
QED.

**Port quotient and what must not be collapsed.** Contract each component
of G-C to one vertex, retain the C vertices, and retain every edge incident
with C separately. Delete the internal edges of the contracted components.
The resulting connected *multigraph* P_C has cycle rank Q. In fact

    beta(G) = sum_K beta(K) + beta(P_C).

This is also a cycle-space statement: the cycle spaces of the outside
components form a direct-sum subspace of the cycle space of G, and the
quotient maps isomorphically onto the cycle space of P_C. One can prove
surjectivity by joining the incident ports within a spanning tree of each
outside component; its kernel consists exactly of cycles internal to those
components. Parallel edges in P_C must be retained. Identifying them changes
the claimed invariant. This is a connectivity interface, not a table of
Boolean input/output behavior.

## 3. Minor certificates and a stronger sequential no-reuse budget

**Lemma 2 (ordinary simplification is a gate-graph minor).** Each allowed
simplification makes the new undirected gate graph a minor of the old one.

**Proof.** Constant and dead-gate deletion are vertex/edge deletion. If g is
rerouted to an input pin p which is itself a gate, first delete g's edges to
its other input-pin sources, then contract the edge pg, and remove duplicate
edges. The resulting adjacencies are exactly the rerouted distinct child
adjacencies, with p retaining its previous users. Acyclicity precludes a
reverse directed adjacency creating an extra gate edge. If p is an external
source, delete g; all its new external wires are absent from the gate-only
graph. Signs change gate truth tables, not this topology. Further dead
cleanup is vertex deletion. Compositions are again minor operations. QED.

This proof supplies branch sets: a surviving gate's set contains itself and
the alias gates contracted into it; discarded vertices belong to no set.
The sets are disjoint and connected in the original graph, and every final
gate edge has an original edge between the corresponding sets.

On a live one-output DAG with at least one gate, only the output has zero
gate children. Therefore

    Phi(D) = sum_g max(d(g)-1,0) = |E(G)|-|V(G)|+1 = beta(G).  (2)

The identity also holds with value zero for an empty final gate graph.
Cycle rank cannot increase under a graph minor. For deletion this follows
by cycle-space inclusion; contracting a nonloop edge preserves cycle rank
before duplicate/loop deletion, and deletion cannot increase it.

**Theorem 3 (port-loss no-reuse).** At step j take a current live DAG,
a region C_j as above, and actually remove C_j, simplify its unary boundary,
and clean up using the declared operations. The result is the next live DAG.
Then

    Phi(D_j)-Phi(D_(j+1)) >= Q_j = E_j+a_j-kappa_j+1,
    sum_j Q_j <= Phi(D_0)-Phi(D_final) <= Phi(D_0).          (3)

**Proof.** Delete C_j first. The later boundary reroutes and cleanup make
the final graph a minor of G_j-C_j by Lemma 2. Its cycle rank is at most
beta(G_j-C_j). Apply (1) and (2), then telescope. If additional constants
outside C_j appear, their deletion is also a permitted minor operation.
QED.

The new charge a_j-kappa_j+1 can be strictly positive when E_j=0. It measures
reconvergence lost at the *incoming* constant-region ports. It is not a new
lower bound on how often a circuit must exhibit such a region.

## 4. Exact decomposition of the deletion certificate

Define the binary gate-input capacity deficiency

    lambda = sum_(v in C) (2 - number of distinct gate parents of v).

This is nonnegative. It counts non-gate input capacity, including primary or
boundary inputs and any repeated-pin or unused-pin deficiency. In particular,
it is not the number of independent primary inputs or restriction parameters.

**Proposition 4.**

    |C|+|B| = Q + lambda + kappa - 1.                     (4)

Thus 0<=Q<=|C|+|B|, and Q equals the whole certificate exactly when
lambda=0 and kappa=1.

**Proof.** The C-induced edges number 2c-2r-f. Hence
lambda=2c-(2c-2r-f)-a=2r+f-a. Substitute this and Q into the inherited
identity |C|+|B|=2r+f+E. Both correction terms are nonnegative. QED.

This is a decomposition of the existing certificate, not a source of extra
deletions. In the signed-output full B2 model all gates C union B are
removable as in Session 188. Strict B2 without free terminal NOT can lose
one net gate saving; the graph topology statements are unchanged, but adding
that new NOT is outside their no-new-gate operation list.

**Sharp four-gate hostile.** On independent inputs a,b,c,y take

    P = b AND c,
    A = P AND a,
    B = A XOR y,
    O = P AND B.

Fix a=0. The full constant region is C={A}; its boundary is {B}. P and O
remain nonconstant genuine binary gates; the reduced circuit is
P=b AND c, O=P AND y. Exactly A and B disappear under constant/unary/dead
cleanup. Here r=1,f=E=0,a=1,kappa=1,lambda=1,Q=1. Initial Phi=1 and final
Phi=0. The certificate is still two gates, not 2+Q=3. In P_C the two physical
incidences P-A and A-B become parallel edges between A and the same outside
component; treating that quotient as a simple graph would incorrectly give
Q=0. This one example tests both the port multiplicity and the non-additivity
of topological consumption.

## 5. Universal core shielding

The Session 188 ladder preserves one diamond. The following construction
preserves any declared core, without turning it into a formula.

Let H be a live single-output B2 DAG, let p be any nonoutput gate with d>=1
distinct children, and let k>=1. Introduce independent control inputs
(a_i,b_i) and add gates

    A_i = a_i AND b_i,
    G_0 = p,
    G_i = G_(i-1) XOR A_i,     i=1,...,k.

Redirect *every original use* of p in H to G_k, leaving all other old gate
and input identities unchanged. Place the inserted gates immediately after
p in a topological ordering. Call the new circuit H[k,p].

**Theorem 5 (arbitrarily long zero-consumption shielding).**

1. H[k,p] has |H|+2k live gates and Phi(H[k,p])=Phi(H). It contains H as a
   gate-graph minor; no expansion of shared gates into copies is involved.
2. Restrict a_k=0,a_(k-1)=0,...,a_1=0 in that order. At each step the local
   rewrite A_i=0, G_i=G_(i-1) deletes two gates and restores the same form with
   one shorter ladder. Its region C_i={A_i} has r_i=1, f_i=E_i=0,
   a_i^(ports)=0, kappa_i=1, and Q_i=0. The final circuit is exactly H.
3. If every gate of H is nonconstant and is not a signed literal of one of
   its input-pin sources on the original primary-input domain, then these
   are exactly the gates removed at each step by exhaustive certified
   constant/direct-pin-literal/dead cleanup. The output remains nonconstant.
4. The deleted boundary gate G_i has d children at every step, so the sum
   of its contemporaneous fanout excesses is k(d-1), while the preserved
   gate-graph minor and Phi are unchanged throughout.

**Proof.** All inserted gates reach p's old children and then the output.
p's child degree changes from d to one; the excess d-1 moves to G_k. Every
other inserted gate has one child, and all old child degrees except p are
unchanged. This proves the Phi identity. Deleting A_i and contracting each
G_i into its preceding gate gives H as a minor.

At reverse-order restriction i, higher stages have already vanished. A_i
is a leaf of the gate-only undirected graph because both of its pins are
new primary inputs. Thus G-C_i is connected, a_i^(ports)=0 and E_i=0.
The two local simplifications restore exactly the prior form, proving 2
and 4.

For 3, any old gate becoming constant or a signed input-pin literal on the
current restricted domain would remain so after further specializing all
remaining a_j to zero. But that specialization recovers H and the original
pin sources, contradicting the hypothesis. Every remaining A_j is a genuine
AND of two fresh free bits. Every remaining G_j is the XOR of a nonconstant
preceding value and a fresh independently variable A_j, so it is neither
constant nor a signed literal of either pin. All gates remain structurally
live. Thus no additional permitted simplification exists. QED.

The conclusion is about a specific legal restriction schedule and a specified
simplification model. It does not exclude better chosen restrictions,
arbitrary resynthesis, or stronger semantic rewriting.

### Unbounded reconvergence and grid cores, not just one diamond

For m>=2 construct distinct independent boundary gates

    L_i = x_i XOR x'_i,      i=1,...,m,
    T_j = y_j XOR y'_j,      j=1,...,m,
    V_(i,0)=L_i, V_(0,j)=T_j,
    V_(i,j)=V_(i-1,j) XOR V_(i,j-1),  1<=i,j<=m,

with output V_(m,m). The core has m^2+2m gates, 2m^2 gate edges,
and Phi=(m-1)^2. Its interior gate graph is the m-by-m square grid.
Every gate is live. Each V_(i,j), as a binary linear form, contains L_i and
T_j with coefficient one, has no L index greater than i or T index greater
than j, and its two parents differ because one has L_i and the other lacks
it, while one has T_j and the other lacks it. Hence the parents are distinct
nonzero linear forms and the gate is neither constant nor a signed literal
of either. The boundary gates also satisfy the irreducibility condition.

Shield p=V_(1,1), whose degree is two. For arbitrary k this gives k exact
two-gate steps, Q_j=0 throughout, initial and final Phi=(m-1)^2, and an
unchanged m-by-m grid minor. Thus an arbitrarily large existing cycle budget
or grid minor, without an additional hypothesis linking the chosen seeds to
that core, does not force any consumed-surplus event. These XOR grid circuits
are not claimed to be minimum circuits for their output functions.

## 6. Connection contract and limitations

The useful cross-domain connection is graph-cycle-space/port connectivity:

    shared circuit + constant region
        -> outside-component partition + retained physical port edges
        -> exact cycle-space quotient of dimension Q.

It preserves exactly the topology consumed by deleting C. It forgets the
Boolean truth functions, input-domain dimension, semantic certificates,
minimality, and the cost of computing a restriction. As in named tournament
response interfaces, shared identities cannot be replaced by unrelated
marginals: here parallel physical ports must survive the quotient. No
transfer of a tournament nonvanishing theorem to Boolean circuits is asserted.

The four-gate hostile is decisive against adding Q to a gate-count guarantee.
The shielded grid family is decisive against forcing a Q-event from the
mere presence of a large surviving topology. Neither example resolves the
Session 118 GP/P2 bottleneck, nor provides an asymptotic circuit lower bound.

## 7. Reproduction

Run `python3 verify.py` and `python3 -O verify.py`. The checker uses only the
Python standard library and explicit failures, not optimization-removable
assertions. It enumerates stated small live circuit topologies and all their
closed nonempty constant regions, checks the port identity, capacity identity,
equality case, sequential bound, and explicit branch-set certificates. It
also checks complete small truth-table universes and named shielding/grid
controls. `RECEIPT.json` records exact tested universe sizes and examples.
These finite checks corroborate the displayed all-size proofs; they do not
replace them or certify all circuit functions at larger sizes.
