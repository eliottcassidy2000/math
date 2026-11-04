# Named port replay and current field observers

2026-11-04. FINITE-EXACT independent endpoint replay using the inherited S180 algorithm and codecs on the new exact S181-source lineage. No new field construction or native legality claim.

`verify.py` reads the frozen named source and both literal child matrices. A separately compiled endpoint subset DP (`joint.cpp`, inherited from S180) recomputes H, E1, S1, deletion H, E2, S2 at every vertex, and every directed disjoint two-slot end/start joint. Thus17,18,18 physical ports and272,306,306 directed joints are retained. Counts agree with512425,1537277,1537279. Physical root39 has defects3690,0,43200 and precisely the declared inport sets. These computations are independent of discovery's suffix evaluator, but reuse a previously audited algorithm rather than claiming a new independent implementation.

All integers fit uint64: at order18 the generous two-slot bound19*18!=121645100408832000 is below2^64. No extrapolation beyond these checked orders is intended.

For each current physical port compute golden g=(H,(S2-E1)/2,(S2+E1)/2), checking integrality. The inherited codec from pinned4ad839b20dd280cdcfbadf86cb62b2be3c9615ad, `return6325_completion_20261103/transport.py`, is revalidated against all four frozen historical states. F343 uses X^3+X^2+1 over F7. Frobenius closure and inverse linear coordinates are checked at every port. No ancestor phase is copied to a child.

At physical39:
- Source: golden(512425,4499757,4500460), F49/charge(2,6;5), F343(4,3,6)->(3,6,5)->(1,5,3)->(4,3,6).
- +2: golden(1537277,26050888,26052597), F49/charge(1,4;3), F343(0,3,4)->(1,6,4)->(0,5,6)->(0,3,4).
- +4: golden(1537279,21197939,21199170), F49/charge(5,6;3), F343(2,0,6)->(3,0,3)->(0,0,5)->(2,0,6).

A current explicit alias: in the+4 state, old physical50 has golden(1537279,79029636,79030867), but exactly the same F49/charge and F343 codec as39. Their E2 fields are1595342 and1580479 respectively, so equal field observers do not identify physical endpoint response or actual defect. In the+2 state, physical50 is F49-Frobenius-fixed at(6,0;0), while its F343 codec still has degree3.

The original hostile control is recomputed from its actual matrices. Both H=462915; cut256 gives1061775 versus1395279. Exhaustive adjacency-preserving backtracking gives full automorphism orders81 and27, with361 and121 search nodes. These are graph automorphism groups of the hostile examples, not identifications with either finite-field Frobenius.

Normal and optimized Python replay produce identical mathematical receipts and stdout. Run:

    g++ -O2 -std=c++17 joint.cpp -o joint
    python verify.py
    python -O verify.py

`receipt.json` contains all values and hashes of the frozen input files. The native guards, permanent birth55, completed-state strongness and actual reversal schedules are checked in the separate native audit. Neither these codecs nor those finite certificates establish an all-clock or all-depth two-branch grammar. `PROVENANCE.md` distinguishes the older numerical children on a different zero-ready source from the present S181-source continuation.
