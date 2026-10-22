# Weight tensors (5-matrix basis of loop_engine_fast) of all excess-3 cores via decomp.c; writes decomp_e3/coreNNN.bin
# (raw int64, C order, shape (5,)*E; about 380 MB in all, regenerated in ~6 min, not committed).  mac-mini-2026-10-22-loop4.
import json, subprocess, time, sys
cs = json.load(open('cores_e3.json'))
t0 = time.time()
for i, c in enumerate(cs):
    V, edges = c['V'], c['edges']
    inp = f"{V} {len(edges)}\n" + ''.join(f"{a} {b}\n" for a, b in edges)
    r = subprocess.run(['./decomp', f'decomp_e3/core{i:03d}.bin'], input=inp, capture_output=True, text=True)
    if r.returncode: sys.exit(f"decomp failed on core {i}: {r.stderr}")
    print(f"core {i:3d} V={V} E={len(edges)} aut={c['aut']}: {r.stderr.strip()}", flush=True)
print(f"done in {time.time() - t0:.0f}s")
