# Exact correction and scope note for the incoming checkpoint

**Audited source commit:** 1e4f40f2d5f9419ed97d1c516d0db0a83c100a38. All source paths below were read from that exact Git object, including computation files absent from the working-directory materialization. This packet does not rewrite the incoming theorem or hypothesis.

## 1. Numerical completion, with a failed preregistration

[HYP-9344, statement 2](https://github.com/eliottcassidy2000/math/blob/1e4f40f2d5f9419ed97d1c516d0db0a83c100a38/05-knowledge/hypotheses/HYP-9344-after-the-sylow-day-defect-growth-of-c15-tau-t-at-least-tau-e-affine-and-first-order-rules-for-prime-circulants-and-prime-power-circulants.md) explicitly preregistered defect 88 for C15 on 3+5+1⁵, against the alternative 64. It also explicitly allowed that a value outside {64,88} would refute both. The actual value is **84**, independently verified here. Therefore both proposed extensions from the initial sequence 1,4,16 fail at their first uncomputed case. The corrected sequence for k=0,…,5 is 0,0,1,4,16,84. No closed all-k formula is asserted.

The source count prediction R(5)=17824 was correct. The first failed inference was extrapolating the full-type defect from three earlier values of the colored-tournament sequence, not an error in the relative Burnside calculation. At k=5, full embedding multiplicities have additional structure: the even graph K_(3,5) disjoint union I5 has multiplicity four; K_(5,8) has multiplicity three; K_(5,5) disjoint union I3 and K_(5,5,3) each have two relatively-even embeddings. Treating distinct witnesses as their relative multiplicities would give the wrong answer. The complete balance is 84=60+2+282−260.

[The source note, section 3](https://github.com/eliottcassidy2000/math/blob/1e4f40f2d5f9419ed97d1c516d0db0a83c100a38/05-knowledge/results/odd_symmetry_beyond_sylow_and_affine_circulant_rules_20261019.md) and [THM-4728](https://github.com/eliottcassidy2000/math/blob/1e4f40f2d5f9419ed97d1c516d0db0a83c100a38/01-canon/theorems/THM-4728-hyp-9340-beyond-p-groups-defect-identity-pronormal-frattini-lemma-and-the-smallest-failure-c15-at-n-10.md) correctly identify that row as uncomputed. They also omit the numerical F21(7+1⁶) row, while correctly citing its equality mechanism. This successor supplies 3040 on each side.

## 2. Exact bookkeeping of the scope

The source classification counts 18 cumulative classes through eleven, ten at twelve, and twenty at thirteen. Their sum is **48**, not 44. The raw receipts contain 41+3=44 rows; adding the two borrowed degree-twelve rows gives a displayed combined table of 46. The two omitted degree-thirteen rows explain the remaining gap. Thus “44 computed rows” is valid provenance, but not a complete-through-thirteen numerical universe.

The source receipts share 26 lower-degree rows and four degree-thirteen rows with the three earlier audited packets, hence **30 total common rows**. The source note's “26 and two probes” wording does not include the two later double-triple overlaps. All thirty agree. The current packet independently reproduces all eighteen incoming degree-thirteen rows.

## 3. Pronormality wording: supported restricted claim, unsupported broader one

THM-4728 statement 2 proves the valid implication:

If H≤Aut(X) is pronormal there and X has an odd-sign automorphism, then some odd-sign automorphism normalizes H.

Therefore, for any **N(H)-even H-fixed representative** X of a witness type, H is not pronormal in Aut(X). The theorem's added phrase “for every H-fixed representative” is not established by that proof. A witness type is defined by existence of at least one relatively-even representative; relative evenness need not be invariant under arbitrary full isomorphism that changes the embedding of H. The supported repair is to insert “N(H)-even” into that consequence. This audit does not assert that the broader sentence is false; no counterexample to that stronger sentence is supplied here.

## 4. Double-triple explanation: result survives, cited mechanism is insufficient

THM-4728 statement 5 attributes the relation between synchronized C15 and C5×C3² on 5+3+3+1^f to the previous “same fixed cube, index-two normalizer” mechanism. But their fixed cubes do **not** coincide: at n=13 they have dimensions 16 and 14, and R=1024 versus256. An index-two normalizer relation cannot by itself supply the stated conclusion when the carrier changes.

The numerical relations happen to hold in the listed rows, and the two deficits of six are valid. Their proper structural justification is the independent predicate-restricted C5 descent in [DOUBLE_TRIPLE_DESCENT.md](../double_triple_symmetry_boundary_20261020/DOUBLE_TRIPLE_DESCENT.md) and [SHARP_N13_BOUNDARY.md](../double_triple_symmetry_boundary_20261020/SHARP_N13_BOUNDARY.md), not an unproved application of the same-cube lemma. This scope issue does not undermine the independently verified census.
