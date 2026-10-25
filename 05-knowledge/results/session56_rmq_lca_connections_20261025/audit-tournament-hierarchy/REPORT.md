# Independent audit: tournament hierarchy packet

Verdict: PASS, with an implementation-versus-interface qualification. No mathematical correction required in the supplied REPORT.md. No claim of novelty, unrestricted enumeration, A000568 prediction, or native all-size wiring is warranted.

## Verification

Run:

    python3 /tmp/session56-audit-tournament-hierarchy/audit.py
    python3 -O /tmp/session56-audit-tournament-hierarchy/audit.py

Both passed and produced byte-identical normal.json and optimized.json. Every validation uses explicit exceptions, not Python assert. The independent implementation does not import the producer. It uses bit-mask modules, parent-by-minimal-superset hierarchy, a sparse-table Euler RMQ, and binary-lifting child-under-ancestor navigation. This verifies an actual indexed navigation implementation independently; the producer remains a linear-scan reference verifier. Computational complexity statements use the standard word-RAM model with stored quotient adjacency or transitive ranks supporting constant-time local access, not Python bit-operation performance guarantees.

All 1099 labelled tournaments at n=1..5 passed hierarchy/node/arity bounds, transitive-or-prime quotients, both directed versions of every leaf-pair LCA arc reconstruction, total quotient payload and width-dependent prime payload, and the signed-cut identity. Independent reachability found strongly connected labelled counts 1,0,2,24,544 and prime labelled counts 0,0,2,0,264. These verify the labelled SEQ and prime-substitution formulas through n=5, using exact rational coefficient arithmetic.

All three supplied prime-star examples independently have Hamilton-path counts 9,13,15 by permutation enumeration. Both cut witnesses have ordered scores (3,2,1,2,2), full endpoint-inclusive flux (0,2,2,0,0,0), and H=11,13. Their symmetric difference is exactly reversal of the directed cycle 0->4->1->3->0 in the first witness. Full fixed-Hamilton-order censuses of sizes 1,1,2,8,64 establish no H-separating flux collision through n=4 and separation at n=5. The independent n=5 audit exhausts all 64 inputs, exceeding the producer's stop-at-first-witness check.

The producer was copied to producer-replay/ and run normally and with -O. Both outputs match the supplied certificate.json byte-for-byte. Its assertions are stripped by -O, so producer replay alone is not an optimized-mode validation; audit.py supplies that guarantee.

## Proof review

Strong-module laminarity follows directly from non-overlap. The cited restriction theorem makes recursive local strong modules consistent with the global inclusion tree. Children partition each non-singleton node. Distinct child modules have uniform cross orientations, proving the exact LCA reconstruction in both directions.

The full quotient-pair payload injects into leaf pairs via representatives with that node as LCA. For prime arities d<=w, binom(d,2)<=w(d-1)/2, and summing d-1 over internal nodes gives n-1. The bound concerns prime adjacency payload only, not tree/port overhead; unrestricted tournaments need not have bounded w. Transitive orders require O(n) words. The singleton case is harmless with zero prime payload.

The SCC species identity uses ordered components. For prime substitution, a module meeting two blocks has a nontrivial touched-index module in the prime quotient; hence it meets all blocks. Omitting any vertex then forces a universal source or sink in the prime quotient, impossible at order >=3. Thus the blocks are precisely maximal proper strong modules and the inverse construction is unique. This proves S=X+P composed with T_nonempty. Unlabelled substitution must use the cycle index, not ordinary composition; the C3 term (t(x)^3+2t(x^3))/3 is correct. The arithmetic inverse through n=7 remains inherited-count consistency, not an independent census.

The signed-cut identity cancels internal arcs. Its successive differences recover the full ordered degree vector, so the information-loss claim and degree-preserving cycle mechanism are exact. Minimality is restricted to H separation under a fixed Hamilton order and the specified full flux signature.

## Complexity wording

The supplied REPORT.md says that sparse-table RMQ plus binary lifting gives O(log n) queries after O(n log n) additional preprocessing, while explicitly calling check.py a brute-force verifier. This is mathematically sound and does not explicitly claim an efficient producer implementation. Nevertheless, summaries should say "proved query interface and reference verifier," not "producer implements O(log n) queries." Producer RMQ uses min over an interval and port recovery scans children. Constant-time classical navigation is not implemented here and must not be promoted to a new result.

## Source pins

Every supplied manifest entry was hash-verified. Each inherited theorem copy was independently compared byte-for-byte with git show at commit 4f1760fd1eff346767c032524c932d5e317401d2 in the read-only clean repository; see git-source-pins.json. No producer file, repository file, or remote was modified.

Primary literature inspected on 2026-10-25:

- https://arxiv.org/pdf/1805.04607 : retrieved 33-page PDF, printed p.6/PDF index5. Theorem 16 states the transitive-or-prime quotient result; Theorem 19 characterizes trivial strong-module families; Definition 20 identifies the inclusion tree; Proposition 21 states restriction to a strong module. These support the exact cited structural package. They are classical results recalled in this paper, not new to this session.
- https://www.cs.colostate.edu/pubserv/pubs/McConnell-rmm-digraphDecomp.pdf : abstract states O(n+m) directed modular decomposition; opening quotient discussion explains uniform intermodule adjacency and reconstruction by quotient plus factors. An explicit tournament has m=binom(n,2), so O(n^2) input-scale decomposition is consistent. This pin does not by itself document optimal LCA/level-ancestor data structures.

No blockers remain for these bounded conclusions.
