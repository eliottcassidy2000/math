---
id: THM-5095
title: "Free odd-group graph actions lift to tournament cycle carriers with an exact completion count"
status: "PROVED, independently audited. Positive completion tax persists; Chan-Pak completeness OPEN."
session: Session238, 2026-11-10
source: 05-knowledge/results/session238_planar_chanpak_20261110/graph_groups/PROOF.md
audit: 05-knowledge/results/session238_planar_chanpak_20261110/graph_groups_audit/AUDIT.md
---

# THM-5095: an equivariant carrier is not a count-preserving reduction

[Full proof](../../05-knowledge/results/session238_planar_chanpak_20261110/graph_groups/PROOF.md) · [Independent audit](../../05-knowledge/results/session238_planar_chanpak_20261110/graph_groups_audit/AUDIT.md)

Let a finite odd-order group A act freely on the vertices of a nonempty
finite simple graph G. Write q=|A|, m=|V(G)|, e=|E(G)|, and let o count
odd-degree vertices. The THM-5090 corona construction has an A-equivariant
version on N=e+5m+o<=e+6m tournament vertices. Its selected long odd cycles
and private triangles induce G corona K1 as their conflict graph.

There are exactly

    2^((binom(N,2)-2e-6m-o)/q)

labelled A-invariant completions of the prescribed arcs and fixed action.
This is not an isomorphism-class count and A need not be the full group.
The proof transports incidence orders from vertex-orbit representatives;
free odd actions have free unordered-pair orbits, each with two orientations.
Every completion has faithful free A-action and satisfies

    I(Omega(T),x) >=_coeff (1+x)^m I(G,x/(1+x)) + m*x,
    H(T) >= 3^m I(G,2/3)+2m.

The positive tax is inherited from THM-5090, now retained equivariantly.
Consequently q divides H(T); THM-4987 separately gives H(T)/q not in{7,21}
for this actual subgroup. The latter exclusion is an inherited dependency,
not a fresh proof from the carrier construction.

Even-order groups cannot act faithfully on tournaments. Odd cyclic groups
therefore give a literal cyclic planar-point-group specialization. But a
planar selected graph can also carry a free C3xC3 action, so its planarity
does not imply planar geometric symmetry for the tournament. Freeness is
sufficient here; arbitrary nonfree graph actions are not classified.

Independent controls cover73 graph/action cases, multi-orbit actions, all16
smallest completions and separate path/packing counts. They supplement the
all-order proof. Extra cycles and their positive interactions remain the
obstruction to importing arbitrary graph counts into Chan-Pak.
