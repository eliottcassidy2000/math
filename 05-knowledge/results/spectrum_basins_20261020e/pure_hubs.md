# Inherited pure-blue hubs and the remaining connectivity obstruction

Status: **PROVED** integration consequences, with the combined-connectivity
input explicitly **CITED**. The pure-blue hub theorem is inherited, not a
new result of this packet. The elementary rigid-strong argument below is
a new dependency reconciliation of two already proved mechanisms.

Throughout, F is the actual converse-merged blue-and-black line multigraph,
and K its literal black subgraph. A black-active vertex has positive black
degree; a pure-blue vertex has none. Dividing the parallel black twin
multiplicities by two preserves black connected components, but changes
bridge questions. The paired-black graph and the black-component quotient
are different objects.

## 1. Recovery of the exact pure-blue common hub

The result proposed during integration was already proved in
[defect_transitions_20261020/blue.md, Theorem 2](../defect_transitions_20261020/blue.md):
the three distinct mixed neighbors of every trivalent pure-blue node have
a common black neighbor. That source includes the empty-outer boundary
and literal controls through total order 29. We retain its provenance.

Here is the mechanism, including the singleton outer case. The proved
pure-blue classification writes every trivalent source as

    U=A -> C5 -> A^op,

where A is a rigid word in 1 and C3. In C5 with increments {1,2} modulo 5,
the three path orbits have representatives

    (0,1,2,3,4), (0,2,3,4,1), (0,2,4,1,3).             (1)

All start at vertex 0. If A is nonempty, fixing the outer paths gives
three full source paths with the same first, second, and last vertices
`a,b,z`. This includes |A|=1: a is the initial singleton and b is the
common first central vertex. For |A|>=2 the first two vertices lie in A.

For any nonstrong self-converse U of order n>=4, the inherited rotation
lemma says that the complement partner of a source path P has a black
path to the exact labelled hub

    V=U triangle {ab,az}.                              (2)

Indeed the partner retains P and reverses the original endpoint arc a->z,
so its cyclic rotation is directed. The two K moves differ on exactly ab
and az. The source scores satisfy d(a)+d(b)>n; hence the rotated path's
endpoint sum at V is d(a)+d(b)-1>n-1, proving blackness. For |A|=1,
the direct values are d(a)=n-1 and d(b)=3, so there is no missing boundary.

At A empty, U=C5 and the nonstrong-source lemma does not apply. The
inherited three explicit five-vertex routes instead meet, after the
recorded relabellings, at the common hub with outgoing bitmasks

    (26,20,25,18,0).

This hub is R4 followed by a sink. Its Hamiltonian-path count is five.
The original source partners' intrinsic central counts are 11,15,13,
distinguishing the three neighbors even after converse merging.

Thus every pure-blue vertex at order n>=4 has all its mixed neighbors in
one black component: this is immediate for a leaf, and follows from (2)
or the C5 boundary for a trivalent vertex. The n>=4 restriction matters:
at n=3 the two pure-blue vertices are adjacent to each other.

## 2. An exact reduction of active black connectivity

**Theorem (PROVED relative to CITED Reid connectivity).** Fix n>=4. The
following are equivalent:

1. All black-active vertices belong to one black connected component.
2. Every blue line between two mixed vertices has a black detour between
   its endpoints.
3. Every blue line between two strong mixed vertices has such a detour.

Here a detour means a black walk between the same two merged endpoints;
it need not reproduce an arbitrarily chosen labelled target matrix.

**Proof.** (1) implies (2), and (2) implies (3). Assume (3). The new
[outer theorem](outer.md) supplies the missing detours whenever at least
one mixed endpoint is nonstrong, so every mixed-to-mixed blue line stays
inside a black component.

Contract every active black component in F, leaving the pure-blue vertices
as individual vertices. Every mixed-to-mixed blue line is now a loop.
Section 1 and the proved pure-blue branch classification say that each
pure-blue vertex is adjacent to exactly one contracted active component
(possibly by three parallel lines). There are no pure-to-pure lines at
n>=4. Consequently this quotient cannot join two distinct active black
components.

Reid's theorem makes F connected, hence so is its contraction. Therefore
there is only one active black component. QED.

The CITED input is Kenneth B. Reid's connectivity under directed
Hamiltonian-path reversal, up to isomorphism. The exact translation and
primary-source access boundary are recorded in
[BRIDGES_THROUGH_NINE_AND_LOCAL_CYCLE_LAWS.md, Section 6](../paired_black_bridges_order_nine_20261020/BRIDGES_THROUGH_NINE_AND_LOCAL_CYCLE_LAWS.md).
Reversing a path and taking converse is K, so converse merging permits the
translation. The indexed author announcement was read there; the full
journal proof was not independently audited in that packet or here.

This equivalence is stronger than the assertion that no mixed blue line
is a bridge. For example, an abstract quotient with three active black
components joined by a blue triangle has no mixed blue bridge, yet has
three black components. Pure-blue leaves or triple attachments can be
added without changing that distinction. This is a graph-theoretic
countermodel, not a claimed realizable tournament metagraph.

## 3. An elementary replacement for the rigid-strong citation

The incoming
[rotation_owners_spectrum_20261020e/black.md, Sections 2–4](../rotation_owners_spectrum_20261020e/black.md)
proves two particularly useful statements:

* a nontrivial rotation of a proper consecutive interval is never an
  automorphism of a circulant tournament;
* every marked circulant of odd order r>=5 with step 1 has a directed
  Hamiltonian word obtained from the marked cyclic word by such a proper
  interval rotation.

Both have elementary proofs there. The first uses the interval rotation's
gcd period and the skew connection set. The second uses the guarded block
extraction, or the unique threshold step set when extraction fails. The
r=5 threshold has the explicit word `(0,1,4,2,3)` for step set {1,3}.
The right-hand word itself is the needed nonautomorphic Hamiltonian path;
the surrounding metagraph detour is not required for this consequence.

**Corollary (PROVED).** The rigid strong alphabet is exactly {1,C3},
without a Benhocine–Wojda or pure-blue-alphabet dependency.

**Proof.** Suppose a nontrivial strong S has one Hamiltonian-path orbit.
Choose an HC by Camion's theorem. Mapping one cut word to its next cut
must be an automorphism; the coordinatewise map is exactly the full cyclic
successor. Tournament automorphisms have odd order, so r is odd, and S
is a circulant in these coordinates. If r>=5 the incoming construction
supplies a nonautomorphic Hamiltonian word, contradicting one orbit.
The only remaining nontrivial order is three, whose strong tournament C3
is rigid. Include the singleton. QED.

This gives an alternative dependency path for Lemmas 1–2 of `outer.md`.
The frozen outer proof remains valid with its original cited rigidity
input; the present reconciliation removes that stronger literature input
from the outer-order-at-least-four theorem. Camion remains the classical
strong-tournament HC input. The earlier singleton-padding proof has its
own separately recorded cited dependencies.

## 4. Other incoming theorems: compatible, but different targets

The [double-descent note](../bridge_hubs_and_double_descent_20261020/BOUNDARY_HUB_AND_DOUBLE_DESCENT.md)
proves that every strong nonselfconverse vertex has two distinct descending
**paired-black** edges. Its one-cycle orbit argument, prime-orbit block
construction, two-prime-orbit exception, and bridge-cut consequence survive
an independent read here. This strengthens descent and constrains a
paired-black bridge branch without nonstrong terminals to a self-converse
single-descent bottleneck. It does not connect all basins or imply anything
about individual literal black bridges, which already have parallel twins.

The [Hamiltonian relocation theorem](../hamilton_path_relocation_connectivity_20261020/THEOREM.txt)
is also sound on its stated objects: legal vertex relocation connects all
labelled Hamiltonian words of a fixed tournament. Its insertion-edge
lifting and front-promotion proofs correctly handle endpoint gaps and
the strongly connected initial block. After quotienting, it establishes
the connectivity of ordinary/converse-merged **fibres**, making the local
relocation/reflection forest certificate necessary as well as sufficient
under the existing complete-wiring hypotheses. Those auxiliary edges are
not K moves. It is not the Reid input in Section 2 and does not itself
prove metagraph black connectivity.

## 5. Fresh finite replay

`pure_hubs.py` independently reconstructs all 39 rigid outer words of
weights 1–8, checks all 117 common-hub routes, and replays the three empty
outer C5 certificates. It verifies every actual path, blackness, the exact
target matrix, the intrinsic central score shell, and central path counts.
It also constructs nonautomorphic Hamiltonian words for all 126 marked
circulant step sets containing 1 at orders 5,7,...,15. Each order has
exactly one threshold step set. These controls audit the cited variable-
size identities; they do not replace their proofs.

    python 05-knowledge/results/spectrum_basins_20261020e/pure_hubs.py
    python -O 05-knowledge/results/spectrum_basins_20261020e/pure_hubs.py

Both runs passed. The exact universe and receipt are in `pure_hubs.json`.
