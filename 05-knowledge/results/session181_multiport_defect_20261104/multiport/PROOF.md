# A literal multi-in-port +4 return, with its exact defect boundary

**PROVED all-order multiport identities; FINITE-EXACT named strong +4 realization.** This enlarges the single-in-port attachment grammar and realizes the previously absent numerical child 512425. It does not preserve zero-readiness and does not establish reusable two-branch closure or the marked suffix target (6589,126459).

## Exact unbounded multiport grammar

Let A be any nonempty tournament, B any tournament, and W a nonempty set of actual named vertices of B. Form M(A;B,W) on A, B and a new named root r, orienting A→B, W→r, and r→A∪(B\W). If B is strong, M is strong: B reaches r via W, r reaches A, and A reaches every vertex of B. The suffix size, number of ports and prefix size are not bounded.

Write h=H(B), F_W=Σ_{w∈W}E_w(B,2), and

J_W=Σ_{u,v∈W}J_B(u,v),

where J_B(u,v) counts ordered spanning two-path covers with first path ending at physical u and second starting at physical v. Both paths in J are nonempty. Coincident endpoints have zero contribution. All paths are directed and vertex-disjoint. Then

H(M)=H(A)(h−J_W)+P_A(2)F_W.                    (1)
E_r(M,2)−H(M)=H(A)J_W.                         (2)

These equations retain the JOINT relation, not only separate marginal port responses. The generalized evaluator is (c,F)=(h−J_W,F_W); its intercept c need not be positive. For transitive A=TT_s, equation (1) becomes c+2^sF.

Proof of (1): separate Hamilton paths according to whether all A vertices lie before r, after r, or on both sides. When all A lie before r, the B path preceding r must end in W, and the B path following r must either be empty or start outside W. Its response is F_W−J_W. When all A lie after r, the response is h+F_W. If A is split into two nonempty slots around r, the response is F_W and the number of A choices is P_A(2)−2H(A). Adding gives (1). Empty B slots are already included in E, while J excludes them.

For (2), a two-slot cover ending at r in its first slot has either r alone in that slot, contributing H(A)h, or a B segment ending in W immediately before r. In the latter case A may be distributed between the two slots, contributing P_A(2)F_W. Thus E_r=H(A)h+P_A(2)F_W, and subtraction proves (2).

The associated all-clock affine criterion is still c'=3c+d, F'=3F. This is an evaluator criterion only. Literal realization and closure of an infinite source class remain necessary. In particular increasing the number of root in-ports is not a free zero-ready operation: equation (2) gives the exact defect introduced. Zero-readiness holds exactly when J_W=0.

The full squarefree named endpoint matrix of Session180 remains a sufficient state for arbitrary births, reversals, deletions and changed marks in this grammar. Nothing here compresses that matrix to (c,F), equates physically different ports, or substitutes F49/F343 observer coordinates or symmetry orders for adjacency.

## Actual 170807→512425 witness

Use the Session180 inherited 170807 source, physical root39, prefix50→51, and suffix names

[1,2,4,6,7,8,14,21,22,31,32,33,40,54].

Reverse suffix edges {14,21} and {21,54}. Replace root39's sole in-neighbor31 by the two physical in-neighbors4 and21, changing exactly {39,31}, {39,4}, {39,21}. No vertices are born or deleted; all seventeen old names remain. The new suffix is strong. Its exact responses are

h=18583, F_{4,21}=124383, J_{4,21}=3690,
c=14893.

Consequently

H=14893+4·124383=512425=3·170807+4.
E_39=516115, E_39−H=3690.

`candidate.json` records every final row, physical name and changed edge. This is an actual strong tournament, not merely an arithmetic solution. It is NOT a zero-ready return. Its five reversals now have separately certified native guards: ten actual Hamilton-path complements use five temporary helper vertices, each subsequently deleted, with no old vertex deletion. Every completed repair is strong; completed counts are170807,217551,1215017,1412269,1294011,512425. See [independent native certificate](../multiport_audit/native_audit.json). Joint endpoint and phase audits use these new rows, not inherited phases.

## All-size isolation of this concrete success

Apply the same suffix and attachment change to every common nonempty prefix A. The source evaluator is (2195,42153), whereas the target evaluator is (14893,124383). Therefore

H_out−3H_in=8308H(A)−2076P_A(2).               (3)

Every offset is divisible by4, excluding +2 throughout this family. For |A|≥4, the standard bound P_A(2)≥(|A|+1)H(A) makes (3) at most −2072H(A). The four prefix isomorphism classes on at most three vertices give offsets4156,4,−8300,12 for TT1,TT2,TT3,C3 respectively. Thus +4 occurs exactly at TT2. This is a new concrete +4 branch, not an induction step valid on its output.

## Complete finite suffix menu, unbounded transitive prefix clocks

`search.cpp` starts with the actual fourteen-name 170807 suffix. It exhausts 4187 edits with zero, one or two internal reversals, retaining3815 strong suffixes. For every such suffix it permits every nonempty physical in-port subset W, and every transitive prefix length s≥1. The target512425 has58 realizations. The chosen witness is the only two-port s=2 record; other records do not improve its zero-ready scope.

Completeness of clock truncation is proved rather than assumed: J_W≤F_W because J only restricts the second slot in E. Therefore H(M)≥h+(2^s−1)F_W. Any candidate at a fixed positive target must have h<target and F_W≤target−h. For each such subset, c+2^sF_W increases strictly with s, so the search stops at the first value above target. The loop limit32 is redundant here since F_W≥1 and target512425<2^19.

The stricter generalized all-clock target(c,F)=(6589,126459) has no hit in this entire menu. Its necessary h=c+J_W≤6589+126459=133048 ensures it is never removed by the target-based pruning. This is a bounded-suffix exclusion, not a general nonexistence theorem.

The search uses exact signed64-bit arithmetic. For fourteen vertices, all two-slot path cover counts are bounded by15·14!<2^63, and all tested clock products are safe. `verify.py` independently verifies (1) and (2) for3069 combinations: every labelled suffix through order4, every nonempty W, and each prefix TT1,TT2,C3. For suffix orders≤3 with TT1/TT2 it additionally counts full Hamilton paths by direct permutation. Checks remain active under Python -O.

## Preserved distinctions and open task

Physical names, all actual changed incidences, and the joint J correction are explicit. F49 and F343 actions remain separate observer codecs; they are not asserted unchanged. The inherited actual81/27 automorphism hostile pair remains a control against scalar-state replacement and is not an automorphism claim for this candidate. No five-presentation Q4 identification is made.

Reusable closure still requires an explicitly defined infinite class containing its outputs and actual two-branch transformations that return to it. The new multiport grammar gives exact unbounded state semantics and a literal +4 realization, but neither those facts nor the58 finite hits prove that closure.

## Second bounded multiport menu

The same search from the actual512423 fourteen-name suffix yields43 transitive-clock realizations of1537273=3·512423+4, among3900 strong suffixes in the4187-edit menu. None has prefix length s≥4, so no member of this exact menu retains all four existing prefix names through mere transitive extension. The two two-port records have s=3 and require deleting a recent prefix vertex. These records are mathematical suffix/attachment witnesses, not additional native-compiled programs. The preferred published actual+4 program remains the all-old-name-preserving170807→512425 candidate. `continuation_hits.jsonl` and `continuation_summary.json` retain the full finite evidence.

Reproduction: compile `search.cpp` with the existing C++17 compiler, run `./search 512425 2 < input.txt` and `./search 1537273 2 < continuation_input.txt`; stdout gives hit records and stderr gives each census. Run `python verify.py` and `python -O verify.py`. `build.py` reconstructs the selected named candidate from the sibling independent audit's frozen source.
