"""Exploratory Python version of fp_cover.c. Fixed-port tails over every core (T0, u): values h0 + offsets, letters S (sink) and I_t (bridge block TT_t + w,
w -> u), port kept at u.  State (h, F, Delta), hu fixed.
  S  : (h, 2F, 0)
  I_t: offset x*(F+Delta)/2 if t = 0 else x*2^(t-1)*F ;  F' = 2^(t+1) F - hu ; Delta' = -hu
"""
import sys
ports = sys.argv[1]
X = int(sys.argv[2]); LIM = int(sys.argv[3]); NMAX = int(sys.argv[4]); TMAX = int(sys.argv[5]) if len(sys.argv) > 5 else 2
cores = set()
for line in open(ports):
    code, u, h2, E2, O2, hu2, h3, E3, O3, hu3 = map(int, line.split())
    if X == 2:
        cores.add((h2, E2 + O2, E2 - O2, hu2))
    else:
        cores.add((h3, E3 + O3, E3 - O3, hu3))
print(len(cores), "distinct cores", file=sys.stderr)
vals = set()
for (h0, F0, D0, hu) in cores:
    if h0 > LIM:
        continue
    seen = set()
    stack = [(h0, F0, D0)]
    while stack:
        s = stack.pop()
        if s in seen:
            continue
        seen.add(s)
        h, F, D = s
        vals.add(h)
        # can any bridge still fit?
        if X * (F - hu) // 2 > LIM - h and X * F // 2 > LIM - h:
            continue
        stack.append((h, 2 * F, 0))
        for t in range(0, TMAX + 1):
            off = X * (F + D) // 2 if t == 0 else X * (2 ** (t - 1)) * F
            if h + off <= LIM:
                stack.append((h + off, 2 ** (t + 1) * F - hu, -hu))
mod = 1 if X == 2 else 1
miss = [v for v in range(1, LIM + 1) if v % X == 1 % X and (X != 2 or v % 2 == 1) and v not in vals]
print("cores", len(cores), "misses", len(miss), miss[:60], "...", miss[-10:])
