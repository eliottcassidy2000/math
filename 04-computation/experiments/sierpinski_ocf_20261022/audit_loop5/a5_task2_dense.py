# Audit (HYP-9371 section 8), task 2: the sign-network engine (loop_engine5.SignCore.contract('basis5')) against the
# explicit enumeration decomp.c (dense 5-basis int64 tensor) on excess-4 cores, exact integer equality.
# decomp.c is rebuilt from the repo source into a scratch directory (cc -O2); nothing in the repo is modified.
# Usage: python3 a5_task2_dense.py DECOMP_BINARY TMPDIR PART NPARTS [vset]   (vset: e.g. "1,2,3,4,5" or "6")
# Cores are taken in index order, restricted to V in vset, every NPARTS-th starting at PART.  For V = 6 an evenly
# spaced sample (argument "6:K" = K cores).  One line per core is printed (flushed).
import sys, os, time, subprocess
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import loop_engine5 as LE5

decomp, tmp, part, nparts = sys.argv[1], sys.argv[2], int(sys.argv[3]), int(sys.argv[4])
vspec = sys.argv[5] if len(sys.argv) > 5 else "1,2,3,4,5"
cores = LE5.load_cores(4)
if vspec.startswith('idx:'):                                  # explicit core indices, e.g. "idx:741,797"
    pick = [int(x) for x in vspec[4:].split(',')]
elif ':' in vspec:
    V0, K = map(int, vspec.split(':'))
    idx = [i for i, c in enumerate(cores) if c[0] == V0]
    pick = [idx[int(round(j * (len(idx) - 1) / (K - 1)))] for j in range(K)]
else:
    vs = set(int(x) for x in vspec.split(','))
    pick = [i for i, c in enumerate(cores) if c[0] in vs]
pick = pick[part::nparts]
ops = LE5.ops_basis5()
nbad = 0
t00 = time.time()
for i in pick:
    V, edges, aut = cores[i]
    E = len(edges)
    fn = os.path.join(tmp, f'p{part}_core{i:04d}.bin')
    inp = f"{V} {E}\n" + ''.join(f"{a} {b}\n" for a, b in edges)
    t0 = time.time()
    r = subprocess.run([decomp, fn], input=inp, capture_output=True, text=True)
    t1 = time.time()
    if r.returncode != 0:
        print(f"core {i}: decomp failed: {r.stderr.strip()}", flush=True); nbad += 1; continue
    D = np.fromfile(fn, dtype=np.int64).reshape((5,) * E)
    os.remove(fn)
    sc = LE5.SignCore(V, edges, aut)
    W = sc.contract('basis5', ops)
    t2 = time.time()
    ok = (W.dtype.kind == 'i') and np.array_equal(W, D)
    nbad += not ok
    info = r.stderr.strip().split(':', 1)[-1].strip()
    print(f"core {i:4d} V={V} E={E} |Aut|={aut:6d}: closed runs {len(sc.runs):4d}, compatible sets {len(sc.sets):6d}, "
          f"nnz {np.count_nonzero(D):7d}, sum|W| {int(np.abs(D).sum()):12d}: {'EQUAL' if ok else 'MISMATCH'} "
          f"(decomp {t1 - t0:.1f}s, engine {t2 - t1:.1f}s; {info})", flush=True)
print(f"part {part}/{nparts}, V in {vspec}: {len(pick)} cores, {nbad} mismatches, {time.time() - t00:.0f}s", flush=True)
