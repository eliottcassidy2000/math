#!/usr/bin/env python3
"""truncation_ladder_thm868_macmini_S110.py -- mac-mini-2026-07-15-S110.
THM-868: THE GEOMETRIC LADDER OF THE SUPPORT FILTRATION, + the new-sequence factory.

The support filtration T_j (k-subsets of [r] with <= j runs; T_2 = polygonal, T_inf =
Pascal) has row-sum and diagonal-sum generating functions that are PARTIAL SUMS OF
GEOMETRIC SERIES in one substituted variable each:

  rows:      v = x^2/(1-x)^2      R_j(x) = 1/(1-x) + (1/(x(1-x))) * (v + ... + v^j)
  diagonals: u = x^3/((1-x)^2(1+x))   G_j(x) = 1/(1-x) + (1/(x(1-x))) * (u + ... + u^j)

and the ladders TELESCOPE EXACTLY:
  (1-x)^2 - x^2 = 1-2x            =>  R_inf = 1/(1-2x)         (all subsets)
  (1-x)^2(1+x) - x^3 = 1-x-x^2    =>  G_inf = 1/(1-x-x^2)      (Fibonacci)

with exact DEFICIT (tail) generating functions -- the "missing regions" laws:
  2^r  - R_j:  x^(2j+1) / ((1-x)^(2j+1) (1-2x))
  F(n+1)-G_j:  x^(3j+2) / ((1-x)^(2j+1) (1+x)^j (1-x-x^2))
whose leading exponents 2j+1 / 3j+2 are the first-miss law proved in S109
(Moser's missing 32nd region = j=2, r=5 = coefficient 1 of x^5).

Also here: layer GFs; the d-dimensional Moser law 1 + C(r+d-1,d) + C(r+d-1,d+2);
pyramidal/centered triangles; odd/even-support split (period-8 via (1+i)^N);
the level-multiplicity triangle L(n,k) of THM-866; majorization-gradedness test;
false-peak cross-check (kps cont.22); OEIS candidate strings.
All series arithmetic exact (integer polynomial ops, truncation order NORD).
"""
import sys
from math import comb, isqrt
from itertools import combinations
sys.stdout.reconfigure(line_buffering=True)

NORD = 41   # series truncation order

# ---------- exact truncated power series over Z (lists of ints, index = exponent)
def pmul(a, b):
    c = [0] * NORD
    for i, ai in enumerate(a):
        if ai:
            for j, bj in enumerate(b):
                if bj and i + j < NORD: c[i + j] += ai * bj
    return c
def padd(a, b): return [x + y for x, y in zip(a, b)]
def psub(a, b): return [x - y for x, y in zip(a, b)]
def pinv(a):                      # 1/a, needs a[0] = +-1
    assert a[0] in (1, -1)
    inv = [0] * NORD; inv[0] = a[0]
    for n in range(1, NORD):
        s = sum(a[k] * inv[n - k] for k in range(1, n + 1))
        inv[n] = -a[0] * s
    return inv
def pmono(c, e):
    p = [0] * NORD
    if e < NORD: p[e] = c
    return p
def ppow(a, k):
    r = pmono(1, 0)
    for _ in range(k): r = pmul(r, a)
    return r
ONE = pmono(1, 0); X = pmono(1, 1)
ONE_MINUS_X = psub(ONE, X); ONE_PLUS_X = padd(ONE, X)

# ---------- the triangles
def runs_entry(r, k, j):
    if k == 0: return 1
    m = r - k + 1
    return sum(comb(m, s) * comb(k - 1, s - 1) for s in range(1, j + 1))

def seq_rows(j, R=NORD):    return [sum(runs_entry(r, k, j) for k in range(r + 1)) for r in range(R)]
def seq_diags(j, N=NORD):   return [sum(runs_entry(n - k, k, j) for k in range(n // 2 + 1)) for n in range(N)]

fib = [1, 1]
for _ in range(NORD + 2): fib.append(fib[-1] + fib[-2])

print("=" * 78)
print("(1) THM-868 -- THE GEOMETRIC LADDER (all series exact to order %d)" % (NORD - 1,))
inv1x = pinv(ONE_MINUS_X)
den_rows  = pinv(psub(pmul(ONE_MINUS_X, ONE_MINUS_X), pmul(X, X)))               # 1/(1-2x)
den_diag  = pinv(psub(pmul(pmul(ONE_MINUS_X, ONE_MINUS_X), ONE_PLUS_X), ppow(X, 3)))  # 1/(1-x-x^2)
print("   telescope algebra: (1-x)^2 - x^2 = 1-2x:",
      psub(pmul(ONE_MINUS_X, ONE_MINUS_X), pmul(X, X))[:3] == [1, -2, 0],
      ";  (1-x)^2(1+x) - x^3 = 1-x-x^2:",
      psub(pmul(pmul(ONE_MINUS_X, ONE_MINUS_X), ONE_PLUS_X), ppow(X, 3))[:4] == [1, -1, -1, 0])
okR = okG = okDR = okDG = True
for j in range(1, 7):
    # R_j = 1/(1-x) + sum_{s<=j} x^(2s-1)/(1-x)^(2s+1)
    Rj = inv1x[:]
    for s in range(1, j + 1):
        Rj = padd(Rj, pmul(pmono(1, 2 * s - 1), ppow(inv1x, 2 * s + 1)))
    okR &= Rj == seq_rows(j)
    # G_j = 1/(1-x) + sum_{s<=j} x^(3s-1)/((1-x)^(2s+1)(1+x)^s)
    Gj = inv1x[:]
    inv1px = pinv(ONE_PLUS_X)
    for s in range(1, j + 1):
        Gj = padd(Gj, pmul(pmono(1, 3 * s - 1), pmul(ppow(inv1x, 2 * s + 1), ppow(inv1px, s))))
    okG &= Gj == seq_diags(j)
    # deficits
    defR = pmul(pmono(1, 2 * j + 1), pmul(ppow(inv1x, 2 * j + 1), den_rows))
    okDR &= defR == [2 ** r - v for r, v in enumerate(seq_rows(j))]
    defG = pmul(pmono(1, 3 * j + 2), pmul(ppow(inv1x, 2 * j + 1), pmul(ppow(inv1px, j), den_diag)))
    okDG &= defG == [fib[n] - v for n, v in enumerate(seq_diags(j))]
print("   R_j ladder == T_j row sums (j=1..6):", okR)
print("   G_j ladder == T_j diagonal sums (j=1..6):", okG)
print("   row-deficit GF x^(2j+1)/((1-x)^(2j+1)(1-2x)) exact (j=1..6):", okDR)
print("   diag-deficit GF x^(3j+2)/((1-x)^(2j+1)(1+x)^j(1-x-x^2)) exact (j=1..6):", okDG)
print("   telescope totals: R_6 -> 1/(1-2x) head:", seq_rows(20)[:13] == [2**r for r in range(13)],
      "; G_13 == Fibonacci to order 40:", seq_diags(13) == fib[:NORD])
print("   [S109 ERRATUM: the S109 'Fibonacci defect' line used F(n+2); true deficit below]")
print("   deficit sequences (j=2): rows 2^r - Moser:", [2 ** r - v for r, v in enumerate(seq_rows(2))][5:12])
print("                            diags F - G:     ", [fib[n] - v for n, v in enumerate(seq_diags(2))][8:16])
print("   pure layer-s diagonal sequences (new; layer GF x^(3s-1)/((1-x)^(2s+1)(1+x)^s)):")
def layer_diag(s, N=NORD):
    out = []
    for n in range(N):
        tot = 0
        for k in range(1, n // 2 + 1):
            m = n - 2 * k + 1
            tot += comb(m, s) * comb(k - 1, s - 1)
        out.append(tot)
    return out
for s in (1, 2, 3):
    print(f"      layer {s}: {layer_diag(s)[2*s+1:2*s+13]}")

print()
print("(2) THE d-DIMENSIONAL MOSER LAW: row sums of the d-dim polygonal triangle")
def Nfig(s, d, m): return (s - 2) * comb(m + d - 2, d) + comb(m + d - 2, d - 1)
okd = True
dseqs = {}
for d in range(2, 6):
    rows = []
    for r in range(30):
        tot = 1
        for k in range(1, r + 1):
            m = r - k + 1
            tot += Nfig(k + 1, d, m)
        rows.append(tot)
        okd &= tot == 1 + comb(r + d - 1, d) + comb(r + d - 1, d + 2)
    dseqs[d] = rows
    print(f"   d={d}: rows = 1 + C(r+d-1,d) + C(r+d-1,d+2): head {rows[:9]}")
print("   law exact d=2..5, r<30:", okd)
ddiags = {}
for d in range(3, 5):
    dg = []
    for n in range(24):
        tot = 0
        for k in range(0, n // 2 + 1):
            r = n - k
            if k == 0: tot += 1
            else: tot += Nfig(k + 1, d, r - k + 1)
        dg.append(tot)
    ddiags[d] = dg
    print(f"   d={d}: diagonal sums head {dg[:14]}   (Fibonacci-of-d analog)")

print()
print("(3) CENTERED-POLYGONAL TRIANGLE (col k: 1 + k*C(m,2)); row sums = 1 + r + C(r+2,4)")
cok = True
crows = []
for r in range(24):
    tot = 1 + sum(1 + k * comb(r - k + 1, 2) for k in range(1, r + 1)) if r >= 1 else 1
    crows.append(tot)
    cok &= tot == 1 + r + comb(r + 2, 4)
cdiags = []
for n in range(24):
    tot = 1 if n >= 0 else 0
    for k in range(1, n // 2 + 1):
        m = n - 2 * k + 1
        tot += 1 + k * comb(m, 2)
    cdiags.append(tot)
print("   row law exact:", cok, "; rows head:", crows[:10])
print("   diagonal sums head:", cdiags[:14])

print()
print("(4) ODD/EVEN-SUPPORT SPLIT of the FULL triangle (Pascal): period-8 structure")
# rows: sum over s odd of C(r+1, 2s) -- involves Re (1+i)^(r+1): period 8
def split_rows(par, R=26):
    out = []
    for r in range(R):
        tot = sum(comb(r + 1, 2 * s) for s in range(1, r + 2) if s % 2 == par) + (1 if par == 0 else 0)
        out.append(tot)
    return out
odd_rows = split_rows(1); even_rows = split_rows(0)
def gauss_check(r):
    z = complex(1, 1) ** (r + 1)
    return round((2 ** r - z.real) / 2)
print("   odd-support rows:", odd_rows[:12])
print("   = (2^r - Re(1+i)^(r+1))/2 [Gaussian-integer / period-8 law]:",
      all(odd_rows[r] == gauss_check(r) for r in range(26)))
odd_diag = []
for n in range(26):
    tot = 0
    for k in range(1, n // 2 + 1):
        m = n - 2 * k + 1
        tot += sum(comb(m, s) * comb(k - 1, s - 1) for s in range(1, min(m, k) + 1) if s % 2 == 1)
    odd_diag.append(tot)
print("   odd-support diagonals:", odd_diag[:16], "(Fibonacci-split analog)")

print()
print("(5) LEVEL-MULTIPLICITY TRIANGLE L(n,k) (THM-866 made quantitative)")
def landau_sequences(n):
    total = n * (n - 1) // 2
    out = []
    def rec(prefix, s):
        k = len(prefix)
        if k == n:
            if s == total: out.append(tuple(prefix))
            return
        lo = prefix[-1] if prefix else 0
        for v in range(lo, n):
            ps = s + v
            if ps < (k + 1) * k // 2: continue
            if ps + (n - k - 1) * (n - 1) < total: continue
            if ps > total: continue
            rec(prefix + [v], ps)
    rec([], 0)
    return out
from collections import Counter
Lrows = {}
for n in range(3, 12):
    seqs = landau_sequences(n)
    floor = 0 if n % 2 else n
    cnt = Counter((sum((2 * s - (n - 1)) ** 2 for s in q) - floor) // 8 for q in seqs)
    ks = sorted(cnt)
    row = [cnt[k] for k in range(max(ks) + 1)]
    Lrows[n] = row
    uni = all(row[i] <= row[i + 1] for i in range(row.index(max(row)))) and \
          all(row[i] >= row[i + 1] for i in range(row.index(max(row)), len(row) - 1))
    print(f"   n={n:2d}: L-row ({len(row)} levels, total {sum(row)}): {row if len(row)<25 else str(row[:14])+'...'}  unimodal:{uni}")
print("   distance-1-below-transitive multiplicities (n=4..11):",
      [Lrows[n][-2] if len(Lrows[n]) > 1 else None for n in range(4, 12)])
print("   max multiplicities (n=3..11):", [max(Lrows[n]) for n in range(3, 12)])

print()
print("(6) MAJORIZATION GRADEDNESS: are ALL covering relations |Delta-x| = 8?")
for n in range(4, 9):
    seqs = [tuple(sorted(q, reverse=True)) for q in landau_sequences(n)]
    xs = {q: sum((2 * s - (n - 1)) ** 2 for s in q) for q in seqs}
    def majleq(a, b):   # a <= b in dominance (b majorizes a)
        return all(sum(b[:i + 1]) >= sum(a[:i + 1]) for i in range(n))
    covers = []
    for a in seqs:
        ups = [b for b in seqs if b != a and majleq(a, b)]
        minimal_ups = [b for b in ups if not any(c != b and majleq(c, b) for c in ups)]
        covers += [(a, b) for b in minimal_ups]
    bad = [(a, b, xs[b] - xs[a]) for a, b in covers if xs[b] - xs[a] != 8]
    print(f"   n={n}: {len(covers)} covers, non-8 covers: {len(bad)}"
          + (f"  e.g. {bad[0]}" if bad else "  -- POSET IS GRADED BY x/8"))

print()
print("(7) FALSE-PEAK CROSS-CHECK (kps cont.22) at n=5,6,7: descent obstructions")
for n in (5, 6, 7):
    stuck_unit = stuck_all = 0
    total = 0
    pairs = list(combinations(range(n), 2))
    mm = len(pairs)
    floor = 0 if n % 2 else n
    stuck_scores = Counter()
    stuck_x = Counter()
    for code in range(1 << mm):
        adj = [[0] * n for _ in range(n)]
        for i, (u, v) in enumerate(pairs):
            if (code >> i) & 1: adj[u][v] = 1
            else: adj[v][u] = 1
        s = [sum(row) for row in adj]
        d = [2 * si - (n - 1) for si in s]
        x = sum(di * di for di in d)
        if x == floor: continue
        total += 1
        has_any = any((adj[u][v] and d[u] - d[v] > 2) or (adj[v][u] and d[v] - d[u] > 2)
                      for u, v in pairs)
        if not has_any:
            stuck_all += 1
            stuck_scores[tuple(sorted(s))] += 1
            stuck_x[x] += 1
        if n < 7:
            has_unit = any(adj[u][v] and d[u] - d[v] == 4 or adj[v][u] and d[v] - d[u] == 4
                           for u, v in pairs)
            if not has_unit: stuck_unit += 1
    print(f"   n={n}: above-floor {total}; NO unit-descent: {stuck_unit if n < 7 else '(skipped)'}; "
          f"NO descent at all (x-local-min above floor): {stuck_all}")
    if stuck_all:
        print(f"        stuck score census: {dict(stuck_scores)}   stuck x levels: {dict(stuck_x)}")

print()
print("OEIS CANDIDATE STRINGS:")
print("   row-deficit j=2 (2^r - Moser):", ",".join(str(2 ** r - v) for r, v in enumerate(seq_rows(2)))[:80])
cand = [fib[n] - v for n, v in enumerate(seq_diags(2))]
print("   diag-deficit j=2 (F - G):", ",".join(map(str, cand[8:24])))
print("   layer-2 diagonals:", ",".join(map(str, layer_diag(2)[5:20])))
print("   layer-3 diagonals:", ",".join(map(str, layer_diag(3)[8:24])))
print("   d=3 Moser rows:", ",".join(map(str, dseqs[3][:16])))
print("   d=4 Moser rows:", ",".join(map(str, dseqs[4][:16])))
print("   d=3 diagonals:", ",".join(map(str, ddiags[3][:16])))
print("   centered rows:", ",".join(map(str, crows[:16])))
print("   centered diagonals:", ",".join(map(str, cdiags[:16])))
print("   odd-support diagonals:", ",".join(map(str, odd_diag[3:18])))
print("   G_3:", ",".join(map(str, seq_diags(3)[:18])))
print("   G_4:", ",".join(map(str, seq_diags(4)[:18])))
print("   L distance-1 multiplicities:", ",".join(str(Lrows[n][-2]) for n in range(4, 12)))
print("   L max multiplicities:", ",".join(str(max(Lrows[n])) for n in range(3, 12)))
print("\nDONE")
