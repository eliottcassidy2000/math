# Is the rotational tournament R_n = Cay(Z_n, {1..m}) a strict local maximum of H on the circulant flip graph?
# Computes H(R_n) and H of its m one-flip neighbours by inclusion-exclusion (ie_ham.c, exact in Z/2^128).
# Pre-registered partial test of HYP-9362(a) (n = 25, 27, 29, 31).  mac-mini-2026-10-21-fanopsp.
# usage: python3 rot_local_max.py N IE_HAM_BINARY OUTFILE [WORKERS]
import sys, subprocess, re
from concurrent.futures import ThreadPoolExecutor
n = int(sys.argv[1]); BIN = sys.argv[2]; OUT = sys.argv[3]; W = int(sys.argv[4]) if len(sys.argv) > 4 else 2
m = (n - 1) // 2
R = list(range(1, m + 1))
sets = {0: R}
for k in range(1, m + 1): sets[k] = sorted([s for s in R if s != k] + [n - k])
done = {}
try:
    for line in open(OUT):
        mm = re.match(r'flip=(\d+) S=\S+ H=(\d+)', line)
        if mm: done[int(mm.group(1))] = int(mm.group(2))
except FileNotFoundError: pass
def run(k):
    r = subprocess.run([BIN, 'circ', str(n)] + [str(s) for s in sets[k]], capture_output=True, text=True, check=True)
    return k, int(re.search(r'H=(\d+)', r.stdout).group(1))
with open(OUT, 'a') as f, ThreadPoolExecutor(W) as ex:
    for k, h in ex.map(run, [k for k in sets if k not in done]):
        done[k] = h; f.write(f"flip={k} S={{{','.join(map(str, sets[k]))}}} H={h}\n"); f.flush()
h0 = done[0]
print(f"n={n} H(R_n)={h0}")
for k in range(1, m + 1): print(f"  flip class {{{k},{n - k}}}: H={done[k]} ratio={done[k] / h0:.6f}")
best = max(range(1, m + 1), key=lambda k: done[k])
print(f"R_{n} strict local max: {all(done[k] < h0 for k in range(1, m + 1))}; best neighbour flip k={best}, ratio {done[best] / h0:.6f}")
