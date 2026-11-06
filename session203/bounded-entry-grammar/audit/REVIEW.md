# Independent review: PASS in frozen scope

Reviewed session203-bounded-entry-grammar/PROOF.md and verify.cpp; hashes are in SOURCE_SHA256SUMS. Candidate bytes were not changed. This review excludes the author's separate newer finite-interface extension.

The strong-tournament bound sum_{x in S} s_x >= floor(|S|²/4)+1 is proved correctly. The mixed-vertex insertion preserves each old start path injectively. For roots and cycle predecessors outside the new vertex's outneighbor set, insertion before the first outneighbor and terminal insertion give distinct extensions. Counting excluded cycle positions yields the marked increment. The two-vertex growth step is valid: absent a mixed vertex, strongness forces both uniform sides and an arc from the dominated side into the dominating side. Its new starts both equal the old H. The marked-set induction handles zero marks separately and correctly reaches the stated quadratic bound. Singleton strong tournaments satisfy the boundary; strong order two is impossible.

The ordered-SCC deletion formula, converse strongness conditions, bounded eligible set, finite nontrivial-tail catalogue, and strictly smaller budget except the literal J/O/Z case are correct. The construction remains physical and does not assert a finite arithmetic state space, occupancy, or F(8,137).

Additional direct corollary, not needed by the frozen proof: apply the quadratic bound to the actual entry set A, rather than all eligible vertices. Then |A|²<4k/q. In particular |A|<=10 when k<=26. This does not bound the whole eligible set by 10.

## Controls

1. Independently compiled and ran the author's C++ verifier. AUTHOR_RERUN.json agrees with the author's receipt: all 33,867 labeled graphs through order six, including 22,891 strong graphs; 99,648 mixed marked-size checks.
2. independent.py uses brute permutation enumeration, not endpoint subset DP, and checks every literal marked subset on every strong labeled tournament through order five. It separately reconstructs ordered SCCs for every root and verifies the exact start formula and eligible-set bound. It checks every mixed extension and marked subset for strong inputs through order four, and their two-vertex extensions.
3. INDEPENDENT.json: 571 strong inputs, 17,239 marked subsets, 5,472 mixed marked subsets, 26 paired extensions, and 2,822 SCC-root decompositions; all PASS.

Finite controls verify implementations and examples; the all-order proof is the combinatorial growth argument. No mathematical gap found.
