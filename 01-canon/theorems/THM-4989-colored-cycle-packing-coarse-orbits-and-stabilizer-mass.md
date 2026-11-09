---
id: THM-4989
title: "Coarse colored-cycle orbits realize 7 and21 while weighted mass retains the Hamilton-path gaps"
status: >
  PROVED, independently audited in Session232. FINITE-EXACT controls have
  their stated universes; inherited all-order theorem certificates remain
  explicit dependencies, not newly rerun claims. Chan-Pak completeness OPEN.
session: Session232-automorphisms-20261109
source: 05-knowledge/results/session232_automorphisms_20261109/packing_mass/PROOF.md
audit: 05-knowledge/results/session232_automorphisms_20261109/packing_audit/AUDIT.md
---

# THM-4989

[Full proof](../../05-knowledge/results/session232_automorphisms_20261109/packing_mass/PROOF.md) · [Independent audit](../../05-knowledge/results/session232_automorphisms_20261109/packing_audit/AUDIT.md)

Let P(T) be the ordered Hamilton paths and X(T) the two-coloured collections
of vertex-disjoint directed odd cycles, including the empty packing. Cycles
are identified by rotation, not reversal; colours stay fixed under vertex
relabeling. OCF gives |X|=|P|=H(T). For A<=Aut(T), orbit–stabilizer gives

    |P/A|=H/|A|=sum_[x in X/A] 1/|A_x|,
    |X/A|-H/|A|=sum_[x] (1-1/|A_x|).

The first sum is action-groupoid mass, not an ordinary number of packing
orbits. The second is nonnegative and strictly positive if A is nontrivial,
because the empty packing is fixed. Thus no nontrivial A-equivariant
bijection P->X exists, although OCF equates the cardinalities. By inherited
THM-4812 and the separately proved THM-4987, the weighted mass excludes7/21.

## Explicit coarse countercontrols

For R5=Cay(Z/5,{1,2}), the full group is C5. There are five directed triangles,
a free orbit under translation, and exactly two directed Hamilton cycles,
the all-step1 and all-step2 cycles, both fixed. With two colours, X has15
objects: empty,10 coloured triangles,4 coloured Hamilton cycles. Its coarse
orbit count is7, with stabilizer histogram {1:2,5:5}, while its mass is3.
The Hamilton-cycle completeness follows because a length5 step sum in{1,2}
must be5 or10; the group determination follows by fixing0 and its two ordered
neighbour pairs. LEM-003 already records the two rotation-fixed Hamilton cycles and path
quotient3; those facts are credited rather than rediscovered.

C3 has three fixed packing objects but one free Hamilton-path orbit. Ordered
joins multiply both actions because cycles cannot cross SCC boundaries,
Hamilton words concatenate, and automorphisms fix the SCC order. Hence
R5=>C3 has H45, full group order15, coarse packing count21 and mass3.
Neither example has7 or21 complete Hamilton-path orbits.

## Constant mass is not constant structure

For T_k=R5=>C3^k (ordered components), all k>=0 satisfy

    H(T_k)=15*3^k,  |Aut(T_k)|=5*3^k,
    |P(T_k)/Aut(T_k)|=mass(X(T_k)/Aut(T_k))=3,
    |X(T_k)/Aut(T_k)|=7*3^k.

Deleting the final C3 induces bijections on path-orbit sets, but maps with
exactly three-point fibres on coarse packing-orbit sets. The latter inverse
limit is a seven-element set times {0,1,2}^N and is uncountable. Every fixed-level
image stabilizes immediately, since the maps are surjective. Image stabilization
therefore does not imply eventual constancy.

The packing action groupoids are pairwise inequivalent: the empty object's
stabilizer at level k has order5*3^k, larger than every isotropy order at an
earlier level. Equal mass does not imply equivalence. This is elementary
finite-action mathematics, not an import of a p-adic inverse-limit theorem.

The independent audit uses permutation cycle covers and explicit Hamilton
words on all1099 labelled tournaments through order5, both order8 joins and
the order11 tower, plus actual restriction maps. The source notes identify
Shao/Zhang methodological inspiration without asserting a geometric functor.
No external novelty or Chan–Pak completeness is claimed.

Inherited routes in this directory: THM-002-ocf.md;
LEM-003-aut-acts-freely-on-ham-paths.md;
THM-4812-presentation-fibres-7-and-21-are-absent-at-every-order-orbit-product-lemma-vertex-transitive-lemma-closed-window.md;
THM-4892-tournament-sites-plus-construction.md (a different, colour-swap quotient).
