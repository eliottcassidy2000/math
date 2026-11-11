# Two-sided endpoint rigidity and canonical neighbours of branching cuts

Date: 2026-11-11. **PROVED:** the endpoint rigidity lemma and canonical
closed-tree neighbours below. **Subsequent closure:**
[exchange/STARS.md](../exchange/STARS.md) resolves the intermediate
fixed-pair boundary, proves every closed tree is a star, and obtains the
all-order quantitative domination inequality. Section 3 below preserves
the earlier reduction and identifies the coordinate used by that proof.
The full Chan–Pak spectrum conjecture remains open.

## Independent audit of the exchange carrier

The defect identity, component matching construction, and transitive
cut-leaf lemma in [exchange/PROOF.md](../exchange/PROOF.md) were checked
independently. The identity partitions each literal newborn path exactly
once. For a cyclic component, use a cycle plus a spanning forest rooted
on that cycle; for an open component, a stub plus a rooted spanning tree.
These give distinct path assignments as stated. The cut-leaf proof correctly
uses its domination condition to rule out a last strong component with
only one B vertex. Empty sides are compatible with the endpoint counts.
This preliminary audit does not itself establish the component bound;
the subsequent all-order proof is in exchange/STARS.md.

## 1. A two-sided endpoint lemma

Let U be a nonempty tournament whose vertices are partitioned into nonempty
sets A,B. Suppose

    every a∈A has an outgoing neighbour in B,
    every b∈B has an incoming neighbour in A.             (1)

Then the number of Hamilton paths starting in A and ending in B equals
one **if and only if U is transitive**. In particular a nontransitive U
under (1) has at least two such paths.

### Strong case: explicit second path

Take a directed Hamilton cycle, using the same Hamilton-cycle consequence
already employed in the inherited exchange proof. Each B→A transition
on the cycle gives an A-starting/B-ending path by deleting that edge.
Two or more transitions prove the assertion immediately.

Otherwise write the cycle in its unique-block form

    a_1,...,a_r,b_1,...,b_s,a_1.

There is already the path a_1,...,a_r,b_1,...,b_s. Since b_s→a_1,
condition (1) forces r,s≥2: a_1 must beat a B other than b_s, and b_s
must receive from an A other than a_1.

Let j be the first index with a_1→b_j.

* If j>1, use

      a_2,...,a_r,b_1,...,b_(j−1),a_1,b_j,...,b_s.

  The two new boundary arcs are b_(j−1)→a_1 and a_1→b_j.
* If j=1 and a_r→b_s, use

      a_2,...,a_r,b_s,a_1,b_1,...,b_(s−1).

* If j=1 and b_s→a_r, choose the largest i with a_i→b_s.
  Condition (1) gives such an i, and i<r. Hence b_s→a_(i+1).
  Use

      a_1,...,a_i,b_s,a_(i+1),...,a_r,b_1,...,b_(s−1).

Every displayed path visits all vertices, starts in A, ends in B, and
differs from the original path. All other adjacencies come from the
Hamilton cycle. This proves the strong case constructively.

### Several strong components

Order the strong components from source to sink. The first meets A:
a B vertex there needs an incoming A neighbour in the same component.
The last meets B by the dual argument. Paths of the components concatenate,
so the endpoint count is the A-starting count of the first component,
times the unrestricted counts of intermediate components, times the
B-ending count of the last component.

A nontrivial first component contains at least two A vertices. If it had
only one, every B vertex in that component would receive from this A,
contradicting strong connectivity. Breaking a Hamilton cycle at each of
two A vertices gives two A-starting paths. The last component is dual;
any nontrivial intermediate component gives at least three paths. Thus
uniqueness forces every component to be a singleton, hence U transitive.
Conversely a transitive tournament has one Hamilton path, and (1) forces
its first vertex into A and its last into B. QED.

The cheap finite control enumerated every labelled tournament through
order five and every nontrivial A,B partition. Among the nontransitive
cases satisfying (1), 136 order-four and 7,680 order-five cases were
checked; none has a unique constrained path. These finite counts are
controls, not the proof.

## 2. Canonical neighbours in a closed tree

Use the exact incidence graph from exchange/PROOF.md and the branching
cut lemma from [classification/EXCHANGE_IDEAS.md](../classification/EXCHANGE_IDEAS.md):
in a closed tree, left endpoint multiplicity forces the whole right side
to lie in A, and right multiplicity forces the whole left side into B.

Let L be a **proper upper cut**, B⊊L, in a closed tree. Put R=V\L⊆A.
The same branching lemma implies U[R] is transitive: otherwise its
A-starting multiplicity exceeds one and forces L⊆B, a contradiction.
The empty R has its one empty path, as usual.

Define

    B0={b∈B : no a∈L∩A points to b},
    W=L\B0.                                          (2)

Every B0 vertex beats every A in L. Domination at L implies every A
in W points to some B in W, and the definition of B0 says every B in W
receives from some A in W. Both parts of W are nonempty: L contains A,
and each such A has a B out-neighbour by domination.

Every Hamilton path p of U[B0] can be concatenated with every A-starting,
B-ending Hamilton path q of U[W]: the boundary arc goes B0→A. This is a
B-ending path of U[L]. Closedness at L forces its ending vertex to point
to the first vertex of the unique path r of U[R], because otherwise
`p q z r` would be a deletion-invalid stub at L. Thus `p q r` is an old
Hamilton path with legal newborn slots at both B0 and L.

The cut B0 is dominating: its left side is all B, and every B outside
B0 receives from A in W on the right. If there were two choices of the
pair (p,q), they would form a four-cycle between the distinct cut vertices
B0 and L. Therefore a closed tree forces

    U[B0] transitive,    U[W] transitive.              (3)

For W this follows from the new endpoint lemma; for B0 it is the standard
unique-Hamilton-path characterization. In particular B0 and L have a
canonical unique connecting path. This remains valid if B0 is empty.
The special central cut L=B is excluded: then the two proposed cut
vertices would coincide and the four-cycle argument would be invalid.

The lower-cut statement is dual. A proper lower cut L⊊B has a canonical
upper neighbour obtained by deleting the A vertices which have no
outgoing neighbour in the remaining B part; the corresponding two
induced pieces are transitive in a closed tree.

## 3. The intermediate monotone reduction, now completed by STARS.md

Starting from a proper upper L, first form B0 by (2), then define its
dual upper neighbour

    A0={a∈A : a has no outgoing neighbour in B\B0},
    L′=V\A0.                                         (4)

Every a∈L∩A has an outgoing neighbour in B\B0 by domination at L.
Consequently

    A0⊆V\L,  hence L′⊇L.                             (5)

Thus canonical alternation expands upper cuts; it does not make an
arbitrary walk through the incidence graph. Strict expansions can occur
only finitely often. This reduction by itself leaves the fixed-pair case
L′=L; the subsequent star proof replaces canonical alternation alone
with a maximal-backward-arc splice that forces a strictly smaller lower
branching cut and then a strictly larger upper branching cut.

At a fixed pair, writing V=B0⊔W⊔R, the preceding results make each block
transitive. B0 beats every A in W, and every B in W beats every vertex
of R. Backward arcs that can still create multiplicity inside the two
canonical cut sides have the forms

    W∩B → B0,        R → W∩A.

Those are same-colour arcs. Cross-block arcs R→B0 may also remain; they
already occur in path-stars and do not by themselves create multiplicity
inside either canonical cut side. The same-colour arcs are not excluded
by domination or by scalar endpoint counts. They were the precise missing
construction at this intermediate stage. STARS.md now uses such an arc
to produce a new path and an explicit directed triangle certifying the
next branching cut, which completes the finite strict-progress argument.

Source→target: labelled domination fibres→canonical cut neighbours.
Preserved: actual simultaneous old-path witnesses, deletion validity,
and closed-tree cycle exclusion. Lost by endpoint cardinalities alone:
the same-colour backward arcs at the fixed pair. This is the precise
remaining local coordinate after the positive reduction.
