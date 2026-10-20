# Dependencies, scope and hostile boundaries

The complete D=z^4 result has four different evidence layers, which are not interchangeable.

1. PROVED DIRECTLY: the half-supported permutation obstruction, its four-cycle prefix clique, the regular-rotation-orbit connectivity lower bounds, and the m7 lowdegree refinement.
2. CITED: Moon's vertex-pancyclic theorem, including Hamiltonicity. Primary paper: J. W. Moon, "On subtournaments of a tournament", Canadian Mathematical Bulletin9 (1966), Theorem1,p.297. DOI https://doi.org/10.4153/CMB-1966-038-7 . The independent cutoff reviewer checked the statement in the primary PDF: https://www.cambridge.org/core/services/aop-cambridge-core/content/view/5D7CCBE66AB07BB3F90F771D1CD3CD08/S0008439500052875a.pdf/on-subtournaments-of-a-tournament.pdf . The k6 case supplies an induced strong six-set; the k=n case supplies a Hamilton cycle that can be cut at a prescribed starting vertex.
3. INHERITED PROVED/FINITE-EXACT: the cycle-cut quotient and congruence exclusions in the previous paired-black leaf packet, plus its complete exact center coefficients through9. These are not reclassified from samples in this packet. The unique strong S4 has H5, automorphism order1, blue count3 and Delta2.
4. NEW FINITE-EXACT: the complete local six-vertex configuration test and complete2,121,984 normalized residual presentation exclusion. Separate independent implementations rederived the presentation spaces and all rejection witnesses. Those counts are presentations, not class counts.

No use is made of a speculative all-blue classification to infer near-purity. No factorial asymptotic is needed to obtain the final cutoff, and no generic tournament automorphism-group size bound is needed. The previous bounded-defect finiteness theorem remains valid but is not a logical prerequisite of this sharper classification.

Important hostile boundaries:

- Regular cyclic C5 is pure-blue but nonrigid. It has even-order relative permutations between blue words when reflected supports overlap; hence the half-support assumption is indispensable.
- The original even q2 cycle may have zero blue cuts. The residual enumeration correctly imposes rotation only at n10 and14. It does not impose a reflection that was not proved.
- q1 is not silently assigned the stronger q>=2 connectivity bound. Its separate regular bound, congruence obstruction and exact small orders exclude it at every unknown order.
- The m7 degree4 restriction is derived from a minimum-cut argument, not guessed from tested degree sequences. Both degree inventories are equal by self-converseness; only then is it safe to enumerate low-outdegree orbit cases alone.
- The q>=2 connectivity bound is sharp in general even among self-converse cyclic presentations. Two regular odd cyclic blocks, with diagonal arcs from the first to the second and all other cross arcs reversed, have connectivity(m+1)/2. Any further improvement must use deficit information rather than symmetry alone.
- S4 remains the positive Delta2 control. Rigid wrappers produce unbounded-order nonstrong Delta2 examples, so strongness cannot be dropped from the theorem.
- A leaf is a vertex property, not a full bridge classification. At n4 the two leaves share one edge; no formula here automatically counts distinct bridges.

The historical session14 packets and shared math checkout were not modified by this research. The new packet is intended as an add-only theorem and certificate publication, with explicit supersession of the formerly unknown D coefficients.


Optional proof simplification, independently approved after the main census: HAND_LOCAL_CONFIGURATION_SUPPLEMENT.md replaces the local six-vertex enumeration dependency with six balanced extensions of S4, four explicit exit-cycle certificates, and elementary extensions of C5 and Q5. The original local scans remain validation rather than a needed proof dependency when this supplement is used. The complete residual presentation census is still required.
