#!/usr/bin/env python3
"""procgen_price_20260925 -- runner.  Reproduces 05-knowledge/results/procgen_price_20260925.out:
   python3 04-computation/experiments/procgen_price_20260925_run.py > 05-knowledge/results/procgen_price_20260925.out
Parts: (1) exact undecided densities (procgen_price_20260925_density.py); (2) the explicit construction G_L
(procgen_price_20260925_greedy.c) on the Collatz sheet, the 3n-1 sheet (SHEET control) and the 5x+1 family (DRIFT
control), with flip density / rho_L; (3) lower bound, 111 criterion, DRIFT undecided density
(procgen_price_20260925_lower.py); (4) independent Python re-implementation of G_L (cross-check of the C densities).
Runtime about 2 minutes; peak memory about 330 MB (N = 10^7: dense table 160 MB + hash table 150 MB).
"""
import os, subprocess, sys, tempfile, time, math

HERE = os.path.dirname(os.path.abspath(__file__))
TMP = tempfile.mkdtemp(prefix="procgen_price_")
BIN = os.path.join(TMP, "greedy")
subprocess.run(["cc", "-O2", "-o", BIN, os.path.join(HERE, "procgen_price_20260925_greedy.c")], check=True)
sys.path.insert(0, HERE)
from procgen_price_20260925_density import bad_counts   # noqa: E402

t0 = time.time()
def py(script):
    out = subprocess.run([sys.executable, os.path.join(HERE, script)], capture_output=True, text=True, check=True).stdout
    print(out, end="", flush=True)

def greedy(L, N, off, q, maxf=3):
    out = subprocess.run([BIN, str(L), str(N), str(off), str(q), str(maxf)], capture_output=True, text=True, check=True).stdout
    return out

py("procgen_price_20260925_density.py")
cnt = bad_counts(64)
LOG32 = math.log(2) / math.log(3)
h = -(LOG32 * math.log2(LOG32) + (1 - LOG32) * math.log2(1 - LOG32))

print("=" * 100)
print("PART 2. The explicit construction G_L (sequential certification with low rescues), FINITE-EXACT")
print("=" * 100)
print("  Every n <= N descends within L steps in G_L (independent re-run inside the C program); 'stuck' = n for which the")
print("  rescue repertoire found no certificate (then G_L would not be a member).  Density = flipped pairs among pairs <= N/2")
print("  (final: later rescues touch only pairs of larger points).")
for label, off in (("COLLATZ SHEET (offset 0, T = 3n+1)", 0), ("3n-1 SHEET (offset 1, U = 3n-1) -- SHEET control", 1)):
    print(f"\n  {label}, N = 10^7")
    print("   L    density(G_L)   rho_L      density/rho_L   stuck   max descent   rescues")
    for L in (4, 5, 6, 7, 8, 10, 12, 13, 16, 20, 24, 28, 32, 36, 40):
        out = greedy(L, 10 ** 7, off, 3)
        lines = out.splitlines()
        d = float(lines[0].split("density")[1].split()[0])
        st = lines[0].split("stuck")[1].split()[0]
        resc = lines[1].split("rescues:")[1].split("hash")[0].strip()
        mx = lines[2].split("max descent time")[1].split(";")[0].strip()
        rho = cnt[L] / 2 ** L
        print(f"  {L:3d}   {d:.6f}     {rho:.6f}   {d / rho:.3f}           {st:>3}     {mx:>3}          {resc}")
print("\n  Lemma checks (theorem of section 2 of the note: L >= 8), N = 10^7:")
for off in (0, 1):
    for L in (8, 12, 16, 24, 32, 40):
        out = greedy(L, 10 ** 7, off, 3)
        lem = [x.strip() for x in out.splitlines() if "LEMMA" in x or "partner descent" in x]
        print(f"  offset {off}, L = {L}:\n    " + "\n    ".join(lem))
print("\n  full output for L = 16, Collatz sheet:")
print("  " + greedy(16, 10 ** 7, 0, 3).replace("\n", "\n  "))

print("  DRIFT control: the 5x+1 pairing family (up v -> (5v+[v odd])/2, down v -> floor(v/2)), N = 3*10^5, DFS <= 8 flips")
print("   L    density(G_L)   beta_L (5x+1 undecided)   density/beta_L   stuck")
sys.path.insert(0, HERE)
from procgen_price_20260925_lower import drift_beta   # noqa: E402
for L in (6, 8, 10, 12, 16, 20, 24, 30):
    out = greedy(L, 300000, 0, 5, 8)
    lines = out.splitlines()
    d = float(lines[0].split("density")[1].split()[0])
    st = lines[0].split("stuck")[1].split()[0]
    b = drift_beta(L)
    print(f"  {L:3d}   {d:.6f}       {b:.5f}                   {d / b:.3f}            {st}")
print()
py("procgen_price_20260925_lower.py")

print("=" * 100)
print("PART 4. Independent cross-check: pure-Python re-implementation of G_L (offset 0), N = 2*10^5")
print("=" * 100)
EPS = {}
def pair(v): return (v + 1) >> 1
def step(v, b):
    i = (v + 1) >> 1
    return v + i if ((v & 1) ^ b) else v - i
def run(n, ext, L):
    v = n; used = []
    for _ in range(L):
        i = pair(v); b = ext.get(i, EPS.get(i, 0)); used.append((i, b)); v = step(v, b)
        if v < n: return True, used
    return False, used
def build(L, N):
    EPS.clear(); stuck = 0
    for n in range(3, N + 1):
        ok, used = run(n, {}, L)
        if ok:
            for i, b in used: EPS[i] = b
            continue
        opts = []
        if n & 1:
            Tn = (3 * n + 1) >> 1
            oA = {pair(n): 1} if pair(n) not in EPS else None
            oB = {pair(Tn): 1} if (Tn & 1 and pair(Tn) not in EPS) else None
            oF = None
            v = n
            for _ in range(L):
                b = EPS.get(pair(v), 0); w = step(v, b)
                if ((v & 1) ^ b) and (v & 1) and (w & 1) and w < 2 * n and (w & 3) == 3 and pair(w) not in EPS:
                    oF = {pair(w): 1}; break
                v = w
            opts = [oA, oF, oB] if n % 8 == 3 else [oF, oB, oA]
        else:
            r = step(n, EPS.get(pair(n), 0))
            opts = [{pair(r): 1} if (r & 1 and pair(r) not in EPS) else None]
        done = False
        for o in opts:
            if o is None: continue
            ok, used = run(n, o, L)
            if ok:
                for i, b in used: EPS[i] = b
                done = True; break
        if not done:
            stuck += 1
    H = N // 2
    fl = sum(1 for i in range(1, H + 1) if EPS.get(i, 0) == 1)
    bad = sum(1 for n in range(3, N + 1) if not run(n, {}, L)[0])
    return fl / H, stuck, bad
for L in (4, 8, 12, 16, 20):
    d, st, bad = build(L, 200000)
    c = greedy(L, 200000, 0, 3).splitlines()[0]
    dc = float(c.split("density")[1].split()[0])
    print(f"  L = {L:2d}: Python density {d:.6f} (stuck without DFS: {st}, re-check failures {bad});  C density {dc:.6f}")
print(f"\n[total runtime {time.time() - t0:.0f} s]")
