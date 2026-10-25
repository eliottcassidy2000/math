# Cycle-packing audit

Verdict: PASS, with scope as stated in the producer README.

Frozen input manifest: 5213de5e2193cd43e952d529277bb44eca5053fa4d86d9b9afcee944458bc116.

The HP blocker definition is arc deletion hitting every directed spanning path. Each directed Hamilton cycle supplies a spanning path after deletion of zero or one of its arcs. Therefore r arc-disjoint Hamilton cycles force every HP blocker to have at least 2r arcs. A Hamilton decomposition of an n-vertex tournament has r=(n-1)/2 cycles; the n-1 incident arcs of a vertex are a blocker for n>=3. Hence beta_HP=n-1. At equality precisely two arcs lie in each specified decomposition cycle, a necessary condition only.

For odd prime order, every nonzero additive step is a unit, so each difference class is one spanning cycle; the skew choice of differences partitions tournament arcs. Thus all prime-order circulant tournaments have beta_HP=q-1, including P7=6 and P11=10. The argument does not prove the decomposition hypothesis for arbitrary regular tournaments. For composite order 9, step 3 has three 3-cycles rather than one Hamilton cycle, so the direct difference-cycle extension fails without ruling out other decompositions.

Independent forward subset DP recomputed all 63 inherited P7 blockers and all 121 inherited P11 blockers, with distinctness, cardinality and cycle-profile checks. Exhaustion of the 9,261 two-per-cycle P7 candidates yielded exactly those 63 minima. Explicit P7/P11 profile-only nonblocker paths were independently checked edge by edge and by DP. Among the 55 ten-arc sets fixed by multiplier 3 on P11, exactly one blocks HPs. No full P11 minimum count or classification is asserted.

Normal and optimized independent runs are byte identical. Producer replay passes from arbitrary cwd with exact file-set/digest checking and file-add/content-mutation hostiles. All pinned copied source files and source-manifest digests agree with their producer locations. No producer or repository changes were made.

Original HYP-9168 credit is explicit: general tournament beta=Hall remains OPEN, regular and Z11/Z13 finite verifications were prior, and THM-4553's P7 census was prior. This supplement supplies a structural cycle-packing proof and a restricted all-prime theorem, not a literature-priority claim. The literal earlier source includes a different computational packing lower bound, so absence of prior structural argument is only the stated targeted repository-search finding, not independently established novelty.

Artifacts: audit.py, result.json, result-O.json, replay.txt. The independent checker reads locally bundled frozen producer and provenance bytes. The audit release is now standalone and closed; replay.py checks its exact file set, all digests, both normal/optimized runs, the producer replay, and file mutation/addition/deletion hostiles. Mathematical content is unchanged from the initial audit.
