# Independent audit: complete strong-C2 classification

**PASS**, no correction to the frozen snapshot. This report is bound to:

- C2_COMPLETE.md:4cc0bf11efc435d7a2746daf6d8f1a58c49c94e57e608911b24464facdd1057d.
- verify_c2_complete.py:c5c6347c6bcfbad0a2caabd098287c6f01a7f034e5b838dbdbf2ad44e9d4d736.
- c2_complete_certificate.json:7c247d8b15143f7e04dd296cc6cccaaecea7a05478602312f4cca362935f099f.

Source files were copied before replay and left unchanged. The completed C1 theorem, separating-vertex product, start-two theorem and full bridge-C2 theorem were independently reviewed previously; the immediately preceding bridge audit is ../ribbon-review-bridge-two/REPORT.md, SHA41da81b4ff3ed5bf3756b3fe282841e2e4ad2b94ee9c8bf626c2cac569334413.

## Approved all-order conclusion

The regular five-tournament is the only bridge-free strong tournament with exactly two Hamiltonian cycles. Consequently the complete unlabelled generating function is

    z^5 + z^5(1−z)(1−z^2)/(1−3z+z^2)^2.

Every other strong-C2 tournament belongs to the complete bridge grammar and is rigid. The labelled count is n! b_(n−1,2)+24[n=5]. This is an all-order theorem, not a census extrapolation. It does not prove the unrestricted Hamilton-path spectrum conjecture.

## Proof checks

1. Cycle insertion is genuinely injective. Strongness gives both incidence signs at the deleted vertex on every old Hamilton cycle, hence a cyclic in-to-out transition. Insert at a deterministic such edge. Removing the same named vertex recovers the old cyclic word modulo rotation; neither distinct old cycles nor distinct insertion positions on one cycle can collide.
2. A bridge-free C2 tournament is2-strong by the already reviewed separating-vertex product. In a minimal counterexample of order at least7, every deletion is strong and has one or two cycles. The completed C1 theorem or minimality plus the bridge-C2 equivalence gives an outdegree-one vertex in each deletion. The regular-five exception cannot occur because deletion order is at least6.
3. Such a low vertex has full outdegree exactly2 and points to the removed vertex, using the minimum outdegree2 of a2-strong tournament. Thus the set L of outdegree2 vertices beats every vertex. Its total outside arcs are2l−binom(l,2); the resulting bound n<=l(7−l)/2<=6 for1<=l<=5 is correct, with no distinct-endpoint assumption. This contradicts a minimal counterexample of order at least7.
4. The small-order proof is elementary. Minimum semidegree2 forces n>=5 and uniquely determines the regular order-five tournament. Its two cycles are the step-one and step-two cycles. At order6, a regular-five deletion supplies three inserted cycles: any2-set is consecutive on exactly one of those cycles, giving1+2 insertion gaps, and the3-set case follows by complement.
5. If no deletion is regular-five, the L-cover at order6 has l=3 or4. For l=4, minimum indegree demands three outside-directed L arcs although only two exist. For l=3, the cover makes L a triangle; the outside triple has total indegree6 and minimum indegree2, and each receives an L arc. This forces its internal indegrees to be1 and the cross L-arcs to be a perfect matching. The three displayed cycles are valid for each of the two remaining orientations. In particular the closing arcs, sometimes easy to overlook, all have the required directions.
6. The regular-five automorphism group has order5; a fixed vertex fixes its internally ordered out- and in-neighbor pairs. The generating function therefore adds one unlabelled class and24 labelled graphs at order5, with no change at other orders. The older residual polynomial is exactly z^5.

The order-six structural simplification was independently derived in this review lane and incorporated by the author before freezing. It removes the need for a census premise or Busch's quantitative bound from the mathematical classification. No literature-priority assertion is made.

## Computational controls

Both normal and Python -O replays reproduce the frozen certificate byte for byte. The source independently checks all33,867 labelled tournaments through order6, including the two-strong histograms C2:24 at order5 and C4:720,C5:1920 at order6, along with15,960 strong-deletion cases and504 physical degree-cover cases.

The separate bases.py was written independently before reading the frozen verifier. It checks all20 regular-five extension patterns and verifies exactly the1+2 insertion-gap split and three resulting cycles for each. It checks all arcs of the three explicit cycles in each matched-triangle orientation; complete literal enumeration gives five cycles in each orientation. See BASES.json.

The preceding independent bridge audit additionally checked139,968 individual cycle-insertion images on all strong single-vertex extensions through order6,55,440 actual arc deletions on every C2 graph there, and13,584 rooted start-two instances. These tests support the precise lemmas and bases; the minimal-counterexample and degree arguments prove the conclusion at arbitrary order.
