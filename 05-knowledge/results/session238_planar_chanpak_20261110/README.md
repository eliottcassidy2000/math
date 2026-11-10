# Session238: planar groups, actual symmetries, and positive tournament constructions

**PROVED, independently audited partial results. Chan–Pak completeness remains OPEN.**
No new unrestricted coverage cutoff is claimed. This session recovers the
user's actual planar-group work and connects it to Sessions236–237 through
three explicit maps, rather than identifying similarly named objects.

## 1. Cyclic symmetry has an exact orbit-profile cost

[THM-5094](../../../01-canon/theorems/THM-5094-cyclic-orbit-profile-cost-and-sharp-sevenfold-path-floor.md)
uses the vertex orbits of an actual cyclic subgroup A=C_m. If their sizes are
d_i and mu_d is the minimum count among d-vertex cyclic Cayley tournaments,

    H(T)/m >= product_i mu_(d_i)/m.

This is sharp for every possible orbit profile. Equality holds exactly for
ordinal sums of minimizing orbit blocks. The proof uses positivity of OCF
and classifies the no-mixed-triangle boundary.

In particular every specified C7 action has H(T)/7>=25. Equality is one
carousel Q7 with transitive padding; s nontrivial orbits cost at least175^s/7.
The eight-circulant input was already in Session190 and is explicitly credited.
A separate corollary, using the inherited THM-4812(c) classification, says a
cyclic A with fewer than23 path orbits must have order dividing15. It does not
re-audit that inherited large classification or classify all tournaments.

[Proof](cyclic_cost/PROOF.md) · [Independent audit](cyclic_cost_audit/AUDIT.md).

## 2. Graph symmetry lifts, but extra cycles do not disappear

[THM-5095](../../../01-canon/theorems/THM-5095-equivariant-odd-group-cycle-carriers-and-exact-completion-count.md)
strengthens Session237's graph-to-cycle-carrier construction equivariantly.
Every free odd-group action A on a graph G lifts to a faithful action on a
tournament whose selected odd cycles induce G corona K1. For m vertices,
e edges and o odd-degree vertices, the construction has N=e+5m+o vertices
and exactly

    2^((binom(N,2)-2e-6m-o)/|A|)

labelled invariant completions. These are neither isomorphism classes nor
claims that A is the full automorphism group.

The inherited positive completion tax remains:

    H(T) >= 3^m I(G,2/3)+2m.

Thus preserving cyclic planar symmetry is possible, but does not give a
count-preserving reduction. Noncyclic odd actions lift too: a planar selected
carrier does not imply a faithful planar geometric group. Even-order faithful
actions are impossible for the familiar involution reason.

[Proof](graph_groups/PROOF.md) · [Independent audit](graph_groups_audit/AUDIT.md).

## 3. Planar cancellation exposes precisely why inverses are dangerous

[THM-5096](../../../01-canon/theorems/THM-5096-planar-commutator-forbidden-orbits-and-positive-tail-word-decoder.md)
recovers the October19 planar Hamiltonian cancellation idea on a literal affine
plane of physical port states. For base x>0 and fixed deletion count d, the
sink-ending physical moves are A(h,R)=(h,2R) and C(h,R)=(h+xR/2,4R-2d).
Adjoining inverses permits a polynomial-Hamiltonian commutator translating h
by -3xd/8 while fixing R.

Inverse powers formally reach forbidden binary counts7,21 and ternary10,31
from genuine seeds. Sink-padded examples even satisfy the inherited necessary
physical cone h<=xR/2, positivity and the base congruence. They are not physical
because group-orbit membership forgets legal construction history.

In contrast the positive A,C semigroup is free. Each word has an exact binary
normal form R'=2^n R-2dm, with no adjacent one-bits in m. At the seed R=2d,
final R alone uniquely decodes the positive word. This supplies a concrete
realizability sidecar for a small alphabet, not a decoder for all tournaments.

[Proof](planar_actions/PROOF.md) · [Independent audit](planar_actions_audit/AUDIT.md).

## What was actually merged from earlier work

The [sourced inheritance report](provenance/INHERITANCE.md) separates:

- Literal S190/Maschke/cyclic group work: planar point groups, SO(3), planar
  Cayley groups, signed actions and graph planarity are different predicates.
- Actual subgroup path-orbit counts and stabilizer-corrected prefixes
  (THM-4987/4990), distinct from a numerical gcd with no group action.
- October19 planar source-flow groups: finite-word rigidity can change under
  generator enlargement or closure; cancellation is operation-relative.
- S236's exact port cones and fixed-alphabet all-word rigidity, and S237's
  arithmetic memory, reroot escape, endpoint slack and positive cycle carriers.

The old eight-circulant census, OCF, gap exclusions, positive carrier tax,
fixed-port formulas, and low-q classification remain credited dependencies.
The new claims are their precise action/profile, equivariance/completion, and
commutator/positive-word connections. No external literature priority is claimed.

## What this changes about the open problem

The three maps preserve different things. Actual cyclic actions give exact
lower bounds and equality cases. Equivariant graph completion preserves a
selected carrier and a group but adds positive counting terms. Formal inverse
closure produces invalid targets, whereas positive words retain a usable clock.

A future constructor must therefore control target magnitude AND its physical
positive history, or supply an independently proved legal operation that changes
roots or generators. Primitive content, full residues, a cone, or abstract group
membership alone cannot establish integer occupation. The full all-odd-except7,21
statement remains OPEN, with no new unrestricted coverage claim in this packet.

## Verification and close-out

[Reproduction](REPRODUCE.md) · [Consolidated replay](REPLAY.json) · [SHA256 manifest](SHA256SUMS).
Six producer/independent programs pass in normal and optimized Python; each
has a deliberate failure control tested in both modes. The cyclic audit uses
literal DFS, the planar audit independently reconstructs odd-cycle packings,
and the graph audit covers multi-orbit actions and all16 smallest completions.
The inherited large vertex-transitive classification and global exclusion
proofs are dependencies, not newly rerun computations.

Verified nonforce checkpoints include
`6cd9b65c319d4ee4f5ba6d4199e2e2480f716392` and
`cd6f29853943f10843322021a45e2d00d3b0bfd8`.
Publication checks compare the exact local Git tree with the remotely created
one and fetch the resulting main ref. Unrelated checkout changes are excluded.

**Session end:** research, independent audits and verification are complete.
These are partial theorems; Chan–Pak completeness remains OPEN.
