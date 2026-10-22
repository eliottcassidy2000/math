# Check of the sign-network engine (loop_engine5.py) against the explicit enumeration (decomp.c, a key-for-key port
# of loop_engine.decompositions_fast): for every core of excess 1, 2, 3 (3 + 15 + 111), the 5-matrix-basis weight
# tensor obtained by contracting the half-edge sign network, with closed runs removed by inclusion-exclusion over sets
# of closed trails, equals decomp.c's integer tensor exactly.  Output: sign_network_basis5_check.out.
# mac-mini-2026-10-22-loop4.
import sys, time
import numpy as np
sys.path.insert(0, '.')
import loop_engine4 as LE4
import loop_engine5 as LE5
ops = LE5.ops_basis5()
t0 = time.time(); bad = 0; tot = 0
for e in (1, 2, 3):
    nruns = []; nsets = []
    for i, (V, edges, aut) in enumerate(LE4.load_cores(e)):
        sc = LE5.SignCore(V, edges, aut)
        W = np.rint(sc.contract('basis5', ops)).astype(np.int64)
        ok = np.array_equal(W, LE4.weight_tensor(e, i))
        bad += not ok; tot += 1; nruns.append(len(sc.runs)); nsets.append(len(sc.sets))
        if not ok: print(f"excess {e} core {i}: MISMATCH")
    print(f"excess {e}: {len(nruns)} cores, closed runs per core {min(nruns)}..{max(nruns)}, compatible sets {min(nsets)}..{max(nsets)}")
print(f"{tot} cores, {bad} mismatches against decomp.c, {time.time() - t0:.1f}s")
