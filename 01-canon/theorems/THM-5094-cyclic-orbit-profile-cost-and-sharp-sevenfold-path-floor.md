---
id: THM-5094
title: "Exact cyclic orbit-profile cost, a sharp sevenfold path floor, and the low-orbit planar symmetry boundary"
status: "PROVED from inherited OCF, independently audited; the low-q corollary additionally uses inherited THM-4812(c). Chan-Pak completeness OPEN."
session: Session238, 2026-11-10
source: 05-knowledge/results/session238_planar_chanpak_20261110/cyclic_cost/PROOF.md
audit: 05-knowledge/results/session238_planar_chanpak_20261110/cyclic_cost_audit/AUDIT.md
---

# THM-5094: the whole cyclic orbit profile has an exact price

[Full proof](../../05-knowledge/results/session238_planar_chanpak_20261110/cyclic_cost/PROOF.md) · [Independent audit](../../05-knowledge/results/session238_planar_chanpak_20261110/cyclic_cost_audit/AUDIT.md)

For odd d, let mu_d be the minimum directed Hamilton-path count among Cayley
tournaments on Z/d, with mu_1=1. If A=C_m acts faithfully on a tournament T
with vertex-orbit sizes d_1,...,d_t, then

    H(T)/m >= product_i mu_(d_i)/m.

Equality holds exactly when the orbit blocks form an ordinal sum and each
attains its cyclic minimum. Every odd orbit profile with lcm(d_i)=m attains
this bound. The proof retains internal odd-cycle packings in OCF; equality
forbids mixed triangles, forcing uniform incidence and a transitive quotient.

For an actual C7 action with s nontrivial orbits,

    H(T)/7 >= 175^s/7 >=25.

Equality in the first bound is an ordinal sum of s carousel Q7 blocks and
fixed singleton vertices. Equality25 has one Q7. Session190 already computed
the eight seven-circulants: six have H175 and two Paley orientations H189.
That finite input and the Q7/P7 comparison are inherited, not new here.
The same profile theorem gives mu_3=3 and mu_5=15. A p11 reproduction is
recorded separately, not used for the sevenfold theorem.

Using the inherited computer-assisted THM-4812(c) classification, a cyclic
A with H(T)/|A|<23 has order dividing15. The possible groups are exactly
1,C3,C5,C15. This is a group-type conclusion, not a classification of all T
or all values. The proof excludes a P7 full-group orbit because its cyclic
projection has order at most7, hence contributes at least27 to the orbit
product. Remaining full-group orbit groups have exponent dividing15.

These statements give genuine cyclic/planar-symmetry costs. They do not
permit changing the denominator from7 to a larger |A| in the sharp floor.
P7 has H/21=9 but H/7=27. The actual action and its orbit profile matter.
No unrestricted coverage extension or Chan-Pak completeness follows.
