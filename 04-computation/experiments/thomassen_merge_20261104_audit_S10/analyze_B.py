#!/usr/bin/env python3
"""Audit of THM-4870 (B) from the audit's own censuses (tourn_audit.c outputs).

Inputs (produced by ./ta.exe full 9 and ./ta.exe strongB 11 224):
  strong_classes.txt     : every strong class of order 1..9:  n H hc aut cyclic code
  strongB224_classes.txt : every strong class of order 1..10 with H <= 224:  n H aut cyclic code
"""
import itertools, sys
from collections import defaultdict, Counter

def pidx(n, i, j):
    return i * n - i * (i + 1) // 2 + (j - i - 1)

def decode(n, code):
    out = [set() for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            if code >> pidx(n, i, j) & 1:
                out[i].add(j)
            else:
                out[j].add(i)
    return out

full = []
for line in open('strong_classes.txt'):
    if line.startswith('#'):
        continue
    n, H, hc, aut, cyc, code = line.split()
    full.append((int(n), int(H), int(hc), int(aut), int(cyc), int(code, 16)))
low = []
for line in open('strongB224_classes.txt'):
    if line.startswith('#'):
        continue
    n, H, aut, cyc, code = line.split()
    low.append((int(n), int(H), int(aut), int(cyc), int(code, 16)))

print("=== consistency of the two censuses on orders <= 9, H <= 224 ===")
a = Counter((n, H, aut) for (n, H, hc, aut, cyc, code) in full if H <= 224)
b = Counter((n, H, aut) for (n, H, aut, cyc, code) in low if n <= 9)
print("multisets of (n, H, |Aut|) agree:", a == b, " total", sum(a.values()), sum(b.values()))
print("strong classes with H <= 224 by order (orders 1..10):", sorted(Counter(n for (n, *_r) in low).items()))
print("total strong classes with H <= 224 (incl. order 1):", len(low), " excl. order 1:", len(low) - 1)

# min H of strong classes by order (Busch 2006 pattern 3,5,9,15,25,45,75,125,225,...)
minH = {}
for (n, H, *_r) in full:
    minH[n] = min(minH.get(n, 10**9), H)
for (n, H, *_r) in low:
    minH[n] = min(minH.get(n, 10**9), H)
print("min H of strong classes by order:", sorted(minH.items()))

print("\n=== asymmetric spectrum below 225 (complete: strong order >= 11 has H >= 225, verified by the census) ===")
LIM = 224
asymH = sorted({H for (n, H, aut, cyc, code) in low if aut == 1})
allH = sorted({H for (n, H, aut, cyc, code) in low})
def closure(vals, lim):
    reach = {1}
    changed = True
    while changed:
        changed = False
        for x in list(reach):
            for v in vals:
                y = x * v
                if y <= lim and y not in reach:
                    reach.add(y); changed = True
    return reach
ra = closure([h for h in asymH if h > 1], LIM)
rf = closure([h for h in allH if h > 1], LIM)
odd = set(range(1, LIM + 1, 2))
print("holes of the FULL spectrum (all tournaments) in odd [1,224]:", sorted(odd - rf))
print("holes of the ASYMMETRIC spectrum in odd [1,224]:", sorted(odd - ra))
symonly = sorted(set(allH) - set(asymH))
print("strong H values <= 224 (orders <= 10) attained only by symmetric strong classes:", symonly)
print("|Aut| values among strong classes with H <= 224:", sorted({aut for (n, H, aut, c, code) in low}))
print("|Aut| | H violations (H<=224 census):", sum(1 for (n, H, aut, c, code) in low if H % aut))
print("strong classes with H = 7 or 21 (orders <= 10):", [(n, H) for (n, H, aut, c, code) in low if H in (7, 21)])

print("\n=== orders <= 8 and <= 9 (complete, all H): values needing symmetry ===")
for NM in (8, 9):
    sv = defaultdict(set)
    for (n, H, hc, aut, cyc, code) in full:
        if n <= NM:
            sv[H].add(aut == 1)
    need = sorted(H for H, s in sv.items() if True not in s)
    print("orders <= %d: strong H values attained only with nontrivial Aut: %s" % (NM, need))
    need1125 = [H for H in need if H <= 1125]
    print("   ... of which <= 1125:", need1125)
    print("   |Aut| values among strong classes of order <= %d (all H): %s; with H <= 1125: %s" % (
        NM, sorted({aut for (n, H, hc, aut, cyc, code) in full if n <= NM}),
        sorted({aut for (n, H, hc, aut, cyc, code) in full if n <= NM and H <= 1125})))
print("max H among classes of order 8 (all strong classes):", max(H for (n, H, *_r) in full if n == 8))

print("\n=== automorphism groups of small strong classes ===")
for n in range(1, 10):
    c = Counter(aut for (m, H, hc, aut, cyc, code) in full if m == n)
    print("order %d strong: |Aut| distribution %s" % (n, dict(sorted(c.items()))))
nonc = [(n, H, hc, aut, code) for (n, H, hc, aut, cyc, code) in full if not cyc]
print("strong classes with NON-CYCLIC Aut, orders <= 9:")
for n in range(1, 10):
    lst = [(H, hc, aut) for (m, H, hc, aut, code) in nonc if m == n]
    print("  order %d: %d classes %s" % (n, len(lst), sorted(lst)[:20]))

def is_selfconverse(n, out):
    # brute force over permutations (n <= 7)
    conv = [set(u for u in range(n) if v in out[u]) for v in range(n)]
    for p in itertools.permutations(range(n)):
        if all((p[w] in conv[p[v]]) == (w in out[v]) for v in range(n) for w in range(n) if v != w):
            return True
    return False

for (n, H, hc, aut, code) in nonc:
    if n == 7:
        out = decode(n, code)
        scores = sorted(len(o) for o in out)
        print("  order-7 non-cyclic: H=%d hc=%d |Aut|=%d scores=%s self-converse=%s" % (
            H, hc, aut, scores, is_selfconverse(n, out)))

# also non-strong order-6 classes with non-cyclic Aut are not in the strong file; the note's 'C3 => C3' example:
c3 = [{1}, {2}, {0}]
T = [set(s) for s in c3] + [set() for _ in range(3)]
for i in range(3):
    T[i] |= {3, 4, 5}
T[3] |= {4}; T[4] |= {5}; T[5] |= {3}
def hpaths(n, out):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for m in range(1, 1 << n):
        for v in range(n):
            c = dp.get((m, v), 0)
            if not c:
                continue
            for w in out[v]:
                if not m >> w & 1:
                    dp[(m | 1 << w, w)] = dp.get((m | 1 << w, w), 0) + c
    return sum(dp.get(((1 << n) - 1, v), 0) for v in range(n))
print("C3 => C3: H =", hpaths(6, T))
