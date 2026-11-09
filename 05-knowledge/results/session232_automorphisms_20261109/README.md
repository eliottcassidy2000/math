# Session 232: where the 7/21 gaps survive automorphisms

**Status: completed mathematical research; all four result packets independently
PASS. Chan–Pak completeness remains OPEN.** The final publication checks are
recorded in REPLAY.json and SHA256SUMS.

The same two numbers behave differently depending on the object being counted.
This session makes that distinction precise rather than inferring a graph
correspondence from matching numbers.

## 1. The global gaps survive every subgroup

[THM-4987](../../../01-canon/theorems/THM-4987-subgroup-path-orbit-gaps-and-small-symmetry-cost.md):
for every tournament T and every A<=Aut(T), the Hamilton-path orbit count
H(T)/|A| is neither 7 nor 21. This extends the inherited full-group theorem
[THM-4812](../../../01-canon/theorems/THM-4812-presentation-fibres-7-and-21-are-absent-at-every-order-orbit-product-lemma-vertex-transitive-lemma-closed-window.md),
not a new proof of its large finite certificates.

There is a sharper symmetry classification. If q=H/|A|<23 and
b=3^(v_5(|A|)+2v_7(|A|)), then q=b or q>=3b. Equality completely determines
the tournament's ordered strong components and forces A to be the full group.
In particular:

- Seven-symmetry below 23 forces q=9 and exactly one Paley P7 component,
  with any number of singleton and C3 components.
- Twenty-five-symmetry below 23 forces q=9 and exactly two R5 components.
- Five-symmetry below 9 forces q=3 and exactly one R5 component.

Here “p-symmetry” means p divides the acting group's order. These conclusions
are deliberately scoped. The cyclic seven-vertex carousel has automorphism
order 7 and quotient count 25, so the unrestricted reading is false.

[Proof](subgroups/PROOF.md) · [Sharp strata](subgroups/SHARP_STRATA.md) ·
[Independent audit](subgroup_transfer_audit/AUDIT.md).

## 2. The gaps disappear at an intrinsic endpoint boundary

[THM-4988](../../../01-canon/theorems/THM-4988-endpoint-721-boundaries-and-faithful-automorphism-graphs.md)
shows that the fibre over the orbit of a final vertex v has size

    e_T(v)/|Aut(T)_v|,

where e_T(v) counts paths ending at v. It is the vertex stabilizer, not the
whole group, that belongs in this denominator.

These fibres attain 7 and 21 at the **smallest possible orders 6 and 7**.
Their complete path counts are 19 and 55. Prepending k>=1 directed triangles
creates nontrivial groups C3^k. The full path-action graphs have 19 or 55
Hamming-graph components; the designated endpoint subgraphs have exactly
7 or 21. At k=1 these are literally unions of triangles.

The same theorem gives an explicit uncoloured simple graph encoding with
Aut(Gamma(T)) isomorphic to Aut(T). Its automorphisms are even, and its
specified traversal traces faithfully encode the tournament paths. Ordinary
spanning graph paths do **not** carry the statistic: the encoding has many
leaves and no Hamiltonian graph path. Automorphism orders 7 and 21 are
present, even though complete tournament path-orbit counts 7 and 21 are absent.

[Proof and exact witnesses](boundary_graphs/PROOF.md) ·
[Independent audit](graph_audit/REPORT.md).

## 3. A coarse quotient can recreate both forbidden numbers

[THM-4989](../../../01-canon/theorems/THM-4989-colored-cycle-packing-coarse-orbits-and-stabilizer-mass.md)
uses the OCF model: two-coloured collections of disjoint odd directed cycles
have the same total cardinality H as Hamilton paths, but a different group action.

- R5 has H=15 and group order 5. There are **7 ordinary packing orbits**,
  but their stabilizer-weighted mass is 3, equal to the Hamilton-path orbit count.
- R5=>C3 has H=45 and group order 15. It has **21 packing orbits**, again
  with weighted mass 3.
- For R5=>C3^k, the path-orbit count and packing mass stay 3 while the coarse
  packing count is 7*3^k. The path restriction maps are bijections, while
  packing restriction maps have three-point fibres.

Thus even equal weighted mass does not identify the underlying groupoids.
The packing tower's images stabilize immediately, but its inverse limit is
uncountable and the levels never become constant. The stabilizer and actual
transition map are indispensable data.

[Proof](packing_mass/PROOF.md) · [Independent audit](packing_audit/AUDIT.md).

## 4. An exact symmetry-aware counting instrument

[THM-4990](../../../01-canon/theorems/THM-4990-stabilizer-corrected-hamilton-prefix-orbit-transfer.md)
gives an integral recurrence on orbits of prefix states (used set, final vertex).
If N(b,a) counts legal extensions and k_b is the pointwise stabilizer order
of the used set, the correct orbit transfer is N(b,a)*k_a/k_b.

A source followed by C3 is the smallest hostile: three extensions produce
one orbit, not three. Pointwise stabilizers vanish once only two vertices
remain unused, and that terminal bound is sharp. The recurrence works for
any nonempty finite digraph; the terminal bound uses the tournament property.

The cut-stabilizer order divides all seven Session231 cut-state coordinates.
This explains why a shared marked core in the specified singleton-padded star
closures must be rigid to support a unit count change, but it is weaker than the inherited content filter. The twelve-bit
core is independently checked rigid; its placement obstruction remains.

[Proof](orbit_transfer/PROOF.md) · [Independent audit](subgroup_transfer_audit/AUDIT.md).

## The papers: precise inspiration, not an imported solution

Both supplied v1 HTML papers were accessible. Relevant theorem/proof passages
were read: [Shao](https://arxiv.org/html/2610.31289v1) and
[Zhang](https://arxiv.org/html/2610.28149v1). Their restricted summaries and exact
pins are in [SOURCE_NOTES](orbit_transfer/SOURCE_NOTES.md).

The finite results above retain stabilizer/transport data before quotienting
and distinguish kinds of stabilization. They do not assert a functor to either
paper's geometric objects, identify a genuine Nygaard stack with a coarse
orbit set, or import p-adic cohomology theorems into tournament theory.

## What remains open

These are omission, boundary, symmetry and counting theorems. They do not
construct every other odd Hamilton-path count. Nor does a small boundary fibre
automatically give a tournament with that total count.

[Session231](../session231_overnight_completion_20261109/README.md)'s unbounded
consecutive windows still need arbitrary target placement. Its complete
fixed-core context language is unchanged by this session. A useful next
construction would exploit the endpoint coordinates while controlling the
whole interface and baseline, rather than forgetting them after counting.

## Verification, inheritance and publication

- Producer subgroup checks: all 532 tournament classes through order 7 and
  all 640 literal class/subgroup pairs, plus sharp larger controls.
- Independent subgroup/transfer checks: all 1,099 labelled tournaments
  through order 5, 1,301 subgroup pairs, 74,521 endpoint-state orbits and
  328 corrected transfer edges; hostile stabilizer choices and failure sentinels.
- Independent packing checks: permutation-cycle covers, explicit path words,
  full automorphism groups and actual restriction maps, including order 11.
- Independent graph checks: 78 uncoloured graph controls, sharp endpoint
  witnesses, exact generator-action components and intrinsic decoding.
- Producer and independent programs agree between normal Python and -O.
  All load-bearing gates use explicit exceptions, not removable assertions.

[Reproduction](REPRODUCE.md) · [Synthesis and research lessons](SYNTHESIS.md).
Inherited THM-4812, OCF, classical free actions and prior R5/P7 constants are
credited in the proof packets. No external novelty is claimed. Neither the
inherited graph spanning-tree comparison nor the new graph encoding supplies
a count-preserving correspondence to ordinary graph Hamilton paths.

Verified nonforce checkpoints:
1. `1219194a7f3f6fe195e18ee0101e9b2de1f1b876`: inherited route and research plan.
2. `5cf220a72f286ac8f446167f5ef59e450d428153`: proof candidates and finite controls.
3. `98f781742bfc7d6884be157a84b2b2ddd1a2aafd`: first independent audit and honest
   empty theorem reservations.
4. The current final version promotes the reservations only after all audits
   passed and preserves all concurrent work.
