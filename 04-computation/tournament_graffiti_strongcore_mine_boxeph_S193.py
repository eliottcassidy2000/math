#!/usr/bin/env python3
"""tournament_graffiti_strongcore_mine_boxeph_S193.py -- boxeph-2026-07-21-S193

Corollary of THM-1855: for join-monotone inequalities the whole content lives in the STRONGLY
CONNECTED core. So mine there. This enumerates the strong iso classes n=3..7 (1,1,6,35,353) and:
  (1) auto-mines TIGHT linear inequalities  target <= a*src + b  and  target <= a*s1 + b*s2 + c
      over the strong core, keeping only those tight on many strong classes (candidate theorems),
  (2) reports the c3<=H margin distribution on the strong core (kind-pasteur's open candidate,
      reduced to strong by THM-1855),
  (3) flags each survivor's join-monotonicity (so we know if it lifts to ALL tournaments).
"""
import sys
from itertools import combinations
from fractions import Fraction as Fr
import tournament_graffiti_coupling_boxeph_S193 as G   # reuse enumeration + invariants (no main run)

NMAX = 7
reps = G.iso_classes(NMAX)
DATA = {n: [(A, G.invariants(A, n)) for A in reps[n]] for n in range(3, NMAX+1)}
STRONG = {n: [(A, d) for (A, d) in DATA[n] if d["strong"]] for n in range(3, NMAX+1)}
SALL = [d for n in range(3, NMAX+1) for (_, d) in STRONG[n]]
ALLC = [d for n in range(3, NMAX+1) for (_, d) in DATA[n]]
print("strong iso classes: " + ", ".join("n=%d:%d" % (n, len(STRONG[n])) for n in range(3, NMAX+1)),
      "  total strong=%d, total all=%d" % (len(SALL), len(ALLC)), flush=True)

KEYS = ["c3", "H", "tr", "dom", "kings", "scc", "smax", "smin", "srange", "n"]

def order_join(A1, n1, A2, n2):
    n = n1 + n2
    B = [[0]*n for _ in range(n)]
    for i in range(n1):
        for j in range(n1): B[i][j] = A1[i][j]
        for j in range(n2): B[i][n1+j] = 1
    for i in range(n2):
        for j in range(n2): B[n1+i][n1+j] = A2[i][j]
    return B, n

JDATA = {1: [([[0]], G.invariants([[0]], 1))], 2: [(reps[2][0], G.invariants(reps[2][0], 2))]}
for n in range(3, NMAX+1): JDATA[n] = DATA[n]

def joinmono(fn):
    for n1 in range(1, NMAX):
        for n2 in range(1, NMAX-n1+1):
            for (A1, d1) in JDATA[n1]:
                if not fn(d1): continue
                for (A2, d2) in JDATA[n2]:
                    if not fn(d2): continue
                    B, m = order_join(A1, n1, A2, n2)
                    if not fn(G.invariants(B, m)): return False
    return True

# ---------------- (1) mine tight PAIR inequalities on the strong core ----------------
print("\n" + "="*92)
print("NEW TIGHT PAIR INEQUALITIES on the STRONG core:  target <= a*src + b   (tight on >= 12 strong classes)")
print("="*92)
COEFS = [Fr(1, 2), Fr(1), Fr(2), Fr(3)]
seen = set()
rows = []
for t in KEYS:
    for x in KEYS:
        if t == x: continue
        for a in COEFS:
            need = max(Fr(d[t]) - a*Fr(d[x]) for d in SALL)
            b = need if need.denominator == 1 else Fr(int(need)+1)  # ceil
            b = int(b) if float(b) == int(b) else (int(float(b))+1 if float(b) > 0 else int(float(b)))
            b = -(-max(Fr(d[t]) - a*Fr(d[x]) for d in SALL)).__ceil__() if False else b
            # robust integer ceil of need:
            import math as _m
            b = _m.ceil(float(need) - 1e-9)
            if abs(b) > 6: continue
            tight = sum(1 for d in SALL if Fr(d[t]) == a*Fr(d[x]) + b)
            if tight >= 12:
                fn = (lambda a=a, b=b, t=t, x=x: (lambda d: Fr(d[t]) <= a*Fr(d[x]) + b))()
                rows.append((tight, t, a, x, b, fn))
rows.sort(reverse=True, key=lambda r: r[0])
shown = set()
for tight, t, a, x, b, fn in rows:
    if (t, x) in shown: continue
    shown.add((t, x))
    holds_all = all(fn(d) for d in ALLC)
    jm = joinmono(fn) if holds_all else False
    astr = ("%s*" % a) if a != 1 else ""
    bstr = (" + %d" % b) if b > 0 else ((" - %d" % (-b)) if b < 0 else "")
    tag = "ALL-n" if holds_all else "STRONG-only(fails on some reducible!)"
    lift = " [join-monotone: lifts to all tournaments]" if jm else (" [NOT join-monotone]" if holds_all else "")
    print("  %-6s <= %s%-6s%-7s  tight/%d strong=%d  %s%s"
          % (t, astr, x, bstr, len(SALL), tight, tag, lift))

# ---------------- (2) c3 <= H margin on the strong core ----------------
print("\n" + "="*92)
print("c3 <= H on the STRONG core (kind-pasteur's open candidate, reduced to strong by THM-1855)")
print("="*92)
margins = [(d["H"] - d["c3"], d["n"], d["c3"], d["H"]) for d in SALL]
minm = min(margins)
allpos = all(m[0] >= 0 for m in margins)
print("  holds on all %d strong classes: %s   min margin (H - c3) = %d at n=%d (c3=%d, H=%d)"
      % (len(SALL), allpos, minm[0], minm[1], minm[2], minm[3]))
# is it ever TIGHT (H==c3) on strong? and what is the max c3/H ratio?
tightc = [m for m in margins if m[0] == 0]
ratios = sorted(((Fr(m[2], m[3]), m[1], m[2], m[3]) for m in margins), reverse=True)[:5]
print("  tight (H==c3) strong classes: %d" % len(tightc))
print("  largest c3/H ratios (closest to the bound): " +
      ", ".join("n=%d c3/H=%d/%d=%.3f" % (r[1], r[2], r[3], float(r[0])) for r in ratios))

# ---------------- (3) mine a TRIPLE inequality that repairs srange tightly ----------------
print("\n" + "="*92)
print("Tighter srange repair search:  srange <= a*tr + b*c3 + c  minimal, tight on strong+reducible")
print("="*92)
best = None
for a in [Fr(1), Fr(2)]:
    for b in [Fr(0), Fr(1), Fr(2)]:
        need = max(Fr(d["srange"]) - a*Fr(d["tr"]) - b*Fr(d["c3"]) for d in ALLC)
        import math as _m
        c = _m.ceil(float(need) - 1e-9)
        if abs(c) > 4: continue
        fn = (lambda a=a, b=b, c=c: (lambda d: Fr(d["srange"]) <= a*Fr(d["tr"]) + b*Fr(d["c3"]) + c))()
        if all(fn(d) for d in ALLC):
            jm = joinmono(fn)
            tight = sum(1 for d in ALLC if Fr(d["srange"]) == a*Fr(d["tr"]) + b*Fr(d["c3"]) + c)
            score = float(a) + float(b)  # prefer smaller coefficients
            cand = (score, a, b, c, jm, tight)
            if best is None or cand < best: best = cand
            print("  srange <= %s*tr + %s*c3 + %d   HOLDS all n<=7, join-mono=%s, tight on %d classes"
                  % (a, b, c, jm, tight))
