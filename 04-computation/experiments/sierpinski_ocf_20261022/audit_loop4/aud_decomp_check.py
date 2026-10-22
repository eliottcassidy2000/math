# AUDIT: decomp.c weight tensors (5-matrix basis) against the Python reference loop_engine_fast.CoreB(V, edges).Wb
# for all excess-1 and excess-2 cores and the 8 excess-3 cores with V <= 2 (the HYP's key-for-key claim is about the
# descriptor maps; this checks the basis tensors that the engine actually uses).  Also: the stored decomp_e3 tensors
# of those 8 cores equal a fresh decomp.c run.
import json, subprocess, sys, os, time
import numpy as np
HERE = '/Users/e/Documents/GitHub/math/.claude/worktrees/gifted-chatterjee-298788/04-computation/experiments/sierpinski_ocf_20261022'
sys.path.insert(0, HERE)
import loop_engine_fast as LEF
tmp = f'{HERE}/audit_loop4/_tmp_core.bin'
for e in (1, 2, 3):
    cs = json.load(open(f'{HERE}/cores_e{e}.json'))
    for i, c in enumerate(cs):
        V, edges = c['V'], [tuple(x) for x in c['edges']]
        if e == 3 and V > 2: continue
        inp = f"{V} {len(edges)}\n" + ''.join(f"{a} {b}\n" for a, b in edges)
        subprocess.run([f'{HERE}/decomp', tmp], input=inp, capture_output=True, text=True, check=True)
        Wc = np.fromfile(tmp, dtype=np.int64).reshape((5,) * len(edges))
        t0 = time.time(); Wp = LEF.CoreB(V, edges).Wb; tp = time.time() - t0
        same = np.array_equal(Wc, np.rint(Wp).astype(np.int64)) and np.allclose(Wc, Wp)
        stored = ''
        if e == 3:
            Ws = np.fromfile(f'{HERE}/decomp_e3/core{i:03d}.bin', dtype=np.int64).reshape((5,) * len(edges))
            stored = f", stored decomp_e3 tensor equal: {np.array_equal(Ws, Wc)}"
        print(f"excess {e} core {i:2d} (V={V}, E={len(edges)}): decomp.c tensor == CoreB.Wb: {same} (nonzero {np.count_nonzero(Wc)}, python {tp:.1f}s){stored}", flush=True)
os.remove(tmp)
