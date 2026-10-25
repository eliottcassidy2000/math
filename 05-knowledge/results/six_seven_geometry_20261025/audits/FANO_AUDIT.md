# Independent audit: Fano code, branch inertia, and Campedelli quotients

2026-10-25. Reviewer: root agent, independent of the Fano producer.
**PASS** for sections 1--6 of `fano_codes/RESULT.md` as read in this session,
subject to their explicit characteristic-zero, no-fourfold, and marking scope.
The universal existence statement may be promoted from CONDITIONAL to a
PROVED geometric deduction with an independently replayed FINITE-EXACT lemma.
There is no assertion of historical novelty or quotient isomorphism census.

## Proof review

1. The effective sign group is the quotient by the global sign. The simplex
   code consists of the zero word and seven weight-four words and meets the
   global-sign line trivially. Its image is the kernel of the residual syndrome:
   the full syndrome kernel is `S + <1>`. For an explicit calculation of the
   inclusion, `sum_(v!=0) v v^T=0` over F2: a diagonal entry counts four vectors,
   an off-diagonal entry counts two. Thus every evaluation word has syndrome0.
   Dimensions then prove equality. This verifies the actual exact sequence,
   rather than just the group orders.
2. A nontrivial code element has sign blocks of sizes four and three. A fixed
   projective point has support within one block. The no-fourfold assumption
   eliminates support of size at most three. The only allowed fixed points are
   nodes above a concurrent triple equal to a Fano line. Every such triple has
   eight nodes, and its stabilizer has exactly two elements. Different Fano
   lines cannot give the same vanishing triple. This proves the iff criterion.
3. Freeness persists on the minimal resolution because a fixed point there
   would map to a fixed point on the normal model. For a free action, the
   quotient is smooth and etale of degree eight; canonical intersection and
   Euler characteristic divide by eight. Nefness and bigness descend by the
   finite surjection, proving minimal general type, with `K^2=2`, `pg=q=0`.
   The twisted canonical character calculation independently proves pg0.
4. At a bad triple, use a reference root from the nonzero sign block. The
   three zero-root tangent coordinates all change sign. The A1 cone is the
   quadratic quotient of A2 by scalar -1; scalar i upstairs changes all three
   quadratic generators by -1. Thus the new quotient is `1/4(1,1)`, not A1.
   It has one (-4)-curve, discrepancy -1/2 and K-square drop1. There are two
   such quotient nodes per bad triple. The formula `K^2=2-2r` is for the
   **minimal resolution of singularities**, not its minimal model. The text
   correctly does not infer its Kodaira dimension or minimal model from K².
   The final incoming Session60 classification now supplies these through
   additional pencil and contraction arguments, cited separately in the text.
5. The parameter reversal really permutes the seven given linear forms.
   Under the adapted marking it is the old modified complement. It moves
   the Fano structure: an involution of GL(3,2) fixes a two-dimensional
   subspace, hence three nonzero vectors, whereas this reversal fixes one.
   Conjugating the code gives an isomorphism between marked quotient models;
   it does not preserve one Fano incidence or the monoid specialization order.

## Independent finite path

`check_fano_independent.py` imports no producer code. It enumerates compatible
triple systems on seven labels using direct set intersection, then constructs
the Fano structures as the systems with seven triples. Thus it does not use
the producer's 5040 vector-label permutations to generate the 30 structures.
It reconstructs each simplex subgroup from the seven complementary triples
and compares the entire set of groups to the producer receipt. It also checks
the canonical invariants and every relevant sign-support stabilizer directly.

Both normal and optimized runs give the same receipt:

- 5596 complete partial Steiner systems, including all nonrealizable hostiles;
- 30 Fano structures;
- at least eight disjoint Fano structures for every system;
- hourglass intersection census `r=0:16, r=1:12, r=3:2`;
- all 30 subgroup sets agree with the producer;
- no invariant canonical section for any selected simplex subgroup.

The search is exhaustive: any increasing list of triples with pairwise
intersection at most one is reached once by the recursion, and candidates
are excluded exactly on an already repeated pair. Each seven-triple system
uses all 21 pairs once and hence is a Fano plane. This is an explicit finite
proof by enumeration, not a sample or a classification of complex surfaces.

Reproduction from repository root:

    python3 -B 05-knowledge/results/six_seven_geometry_20261025/audits/check_fano_independent.py
    python3 -B -O 05-knowledge/results/six_seven_geometry_20261025/audits/check_fano_independent.py

The [receipt](fano_independent_receipt.json) pins both source hashes. The
producer separately enumerates all36 seven-subsets of the nine magic lines.
Their 12--20 free-subgroup counts were inspected as bounded consequence data;
the universal existence assertion already follows from the larger audited
5596-system universe.

## Attribution and remaining gate

Alexeev--Pardini [Definition3.1 and Lemma3.4](https://arxiv.org/html/0901.4431v3)
and Oudompheng [sections1.1--1.3](https://arxiv.org/html/1106.4846v1)
were opened as primary literature. The classical seven-line Campedelli
construction and its singularity types are recovered. The packet's explicit
hourglass census, gauge-kernel connection and marking-changing reflection
are scoped deductions. No arithmetic solution or fundamental group is inferred.

## Final shared-action addendum

The completed section6a was subsequently reviewed and is also **PASS**.
For center c and the three pair differences d_j, each d_j is nonzero and
their sum is c. A bad triple is d_j=c, equivalent to equality of the other
two. Three distinct nonzero vectors with nonzero sum in F2^3 form a basis.
This proves the free/regular equivalence. With one bad triple the rank is2;
with three it is1. The explicit marking `(7,6,5,4,2,3,1)` gives differences
`(3,6,2)` and no bad triple. Evaluating these differences gives the actual
translation on conic sign labels, hence the asserted Walsh-to-canonical
character map. This preserves a group character, not a divisor class or
differential form. The independent root audit checks all5040 markings and
the actual conic intersection configuration over F41. Its receipt and proof
review are [here](ROOT_DELETION_AUDIT.md).
