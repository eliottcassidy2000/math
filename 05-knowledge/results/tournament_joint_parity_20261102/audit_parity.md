# Independent audit of the ordinary-class involution

**Verdict: PASS.** This audit checks the mathematical construction in
[parity.md](parity.md), its implementation in [parity_probe.py](parity_probe.py),
and its ordinary-class/H data against the independent
[joint_audit.json](joint_audit.json). No substantive correction to the parity
proof was required. The all-order parity theorem is inherited; the finite
checks below are not used to promote its scope.

## Mathematical checks

1. **Correct carrier and odd fibers.** The fixed forward-spine cube has
   `2^binom(n-1,2)` presentations. It surjects onto ordinary classes because
   every tournament has a Hamiltonian order. The fiber is exactly the set
   of Hamiltonian orders modulo the free automorphism action, so its size is
   `H(T)/|Aut(T)|`. Both numerator and denominator are odd, and the quotient
   is an integer. The empty tournament is treated separately with one empty
   order. This does not claim the spine cube is closed under arbitrary
   `S_n` relabellings.
2. **The two involutions have the declared meanings.** `sigma` permutes
   chord coordinates by reversal and represents converse downstairs; it
   does not negate bits. `tau` negates every nonspine bit while retaining
   the forward spine. It is free precisely when there is at least one
   chord, namely `n>=3`. Coordinate permutation commutes with addition of
   the all-ones vector, so `tau` preserves `Fix(sigma)`.
3. **Self-converse fibers are odd and nonempty.** On a self-converse class,
   `sigma` acts on its odd presentation fiber, hence has an odd number of
   fixed presentations. A fixed presentation necessarily lies over a
   self-converse class. This proves both directions needed for the map
   `Fix(sigma) -> {self-converse ordinary classes}`.
4. **The cancellation graph is valid with parallel edges.** The local
   matching leaves exactly one representative per odd fiber. Every vertex
   has one `tau` edge and at most one local edge; therefore components are
   paths or cycles, counting two coincident matching edges as parallel
   edges. A representative has degree one and cannot belong to a cycle.
   Its component has two distinct degree-one endpoints belonging to
   distinct fibers. Reversing the component swaps the same endpoints, so
   the induced map is a free involution. No parity-of-cardinality shortcut
   is substituted for this algorithm.
5. **Extension to all ordinary classes.** Self-converse classes are paired
   by these blue paths. Remaining classes are paired with their distinct
   converses. The two domains are disjoint and exhaustive; the resulting
   involution preserves order and the self-converse predicate. Reversing
   Hamiltonian orders proves H-preservation on converse pairs only.
6. **Gauge is explicit.** Integer mask ordering chooses the one unpaired
   representative and the local matching within each fiber. The resulting
   output depends on this convention but not on the names used to present
   the input tournament: its full normalized fiber is the same. It is a
   specified, choice-dependent class involution, not a claim that cube
   complement itself descends or that a local intrinsic class reflection
   has been found.
7. **The H-obstruction is all-order.** A directed triangle has three
   Hamiltonian orders; successive insertion extends each to the whole
   tournament while preserving the old relative order. Thus any
   nontransitive tournament has at least three orders. Transitive
   tournaments, including the empty one under the stated convention, have
   exactly one. The unique `H=1` class prevents a free H-preserving class
   involution at every order. At order three the two self-converse classes
   already force `H=1 <-> 3`.

The independent arithmetic inequalities in parity section 2 were also
checked algebraically. The two equal odd parts require the separate
valuation `n-2` (only possible for even `n>=6` in that case), exactly as
written. The normalized minimal terms at even orders are the coprime
unordered odd pairs, whose count is `phi(n)/2`. The resulting equality
classification retains the exceptional `n=4` and excludes `n=2`.

## Independent finite checks

[audit_parity.py](audit_parity.py) imports no producer code. It reads the
published receipts and reconstructs each blue matching as a complete
undirected multigraph. Component breadth-first search independently gives
the same endpoint class pairs. Adjacency-matrix canonicalization over every
vertex permutation checks that every local matching remains in one ordinary
class, that distinct blue fibers are distinct classes, and that each is
self-converse. Literal vertex permutations verify each displayed H label.

At orders `3,4,5,6`, there are respectively `2,4,16,64` blue masks and
`1,1,4,6` induced class pairs. At order six there are also **two closed
components without any class representative**; ignoring them in endpoint
traversal is correct, not a hidden termination assumption. Ordinary and
self-converse joint H enumerators agree exactly with `joint_audit.json` at
every order `0..6`.

The abstract cancellation lemma is separately replayed for **all set
partitions into odd blocks** at presentation-set sizes `2,4,6,8,10`, with
fixed free matching `(0,1),(2,3),...`. The numbers of tested partitions are
`1,5,37,457,8169`. This includes matching-edge coincidences and disconnected
cycles. All representatives are paired correctly. Relabelling a free
matching to this standard form explains why the universe is a useful
hostile control; the written graph proof supplies the all-order claim.

Reproduction:

```bash
python3 -I -B 05-knowledge/results/tournament_joint_parity_20261102/audit_parity.py
python3 -I -B -O 05-knowledge/results/tournament_joint_parity_20261102/audit_parity.py
```

Both modes use explicit checks, and write the same
[audit_parity.json](audit_parity.json). A fresh independent invocation of the
producer also matched its saved receipt exactly. The SVG is not used as
proof; its mask trails and H labels are checked as mathematical data.

## Cross-lane wording repair

The first toothpick tangent repair had described an unqualified spine-cube
quotient too broadly. The parent review caught this, and T1544 now states the
precise surviving fact: the fixed-spine cube surjects onto ordinary classes
with odd fibers `H/Aut`, while the full `S_n` action lives on the ambient
labelled tournament carrier. The current parity construction is linked
there. This correction changes no toothpick growth proof or incidence
identity.
