#!/usr/bin/env python3
"""
h_realizability_goldbach_macmini_0611s2.py  (mac-mini-2026-06-11-S2, THM-490 / HYP-2420)

THE GOLDBACH ANALOG FOR TOURNAMENTS.
Fermat polygonal: every natural number = sum of <= s s-gonal numbers (atoms = figurate numbers).
Goldbach/Lemoine: every even = p+q ; every odd>=7 = p+2q (atoms = primes; the SEGMENT/digon = 2 atoms).
Tournament analog: H(T) is ALWAYS ODD (Redei). Which odd numbers are realizable as H(T)?
And is there an "atom + composition" basis (the dominance-join as the additive/multiplicative operation)?

Computes, for n = 3..9 (gentourng over iso classes):
  - R_n = set of realizable H values; the realizable-odd-density and the GAPS (forbidden odds below max)
  - the JOIN law: H(T1 => T2) where every vertex of T1 beats every vertex of T2 -- is it multiplicative?
  - SC (self-complementary) classes' H-values: the "self-dual middle" slice
  - the Pfaffian ledger: for each class, Pf = sqrt(det(I+2A)) (integer), and Q = (H^2 - Pf^2)/8 (THM-442);
    is Q >= 0 integer always? the BSD-analog: H (analytic-like total) vs Pf (square-root/skeleton) vs Q (rank-like correction)
"""
import sys, subprocess, itertools
from math import comb, isqrt

sys.stdout.reconfigure(line_buffering=True)
GEN = "/opt/homebrew/bin/gentourng"

def gen_tournaments(n):
    """Yield adjacency matrices (n x n, 0/1, A[i][j]=1 means i->j) for all iso classes via gentourng."""
    out = subprocess.run([GEN, str(n)], capture_output=True, text=True)
    pairs = list(itertools.combinations(range(n), 2))
    for line in out.stdout.split():
        bits = line.strip()
        if len(bits) != len(pairs):
            continue
        A = [[0]*n for _ in range(n)]
        for (i, j), b in zip(pairs, bits):
            if b == '1': A[i][j] = 1
            else: A[j][i] = 1
        yield A

def H_count(n, A):
    out = [0]*n
    for i in range(n):
        r = 0
        for j in range(n):
            if A[i][j]: r |= 1 << j
        out[i] = r
    full = (1 << n) - 1
    # dp[mask][last]
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n): dp[1 << v][v] = 1
    for mask in range(1 << n):
        row = dp[mask]
        for last in range(n):
            c = row[last]
            if not c: continue
            nxt = out[last] & ~mask
            while nxt:
                b = nxt & (-nxt); j = b.bit_length()-1
                dp[mask | b][j] += c
                nxt ^= b
    return sum(dp[full][v] for v in range(n))

def det_int(M):
    """Exact integer determinant via fraction-free Bareiss."""
    from fractions import Fraction
    n = len(M); A = [[Fraction(x) for x in row] for row in M]
    det = Fraction(1)
    for i in range(n):
        if A[i][i] == 0:
            sw = next((k for k in range(i+1, n) if A[k][i] != 0), None)
            if sw is None: return 0
            A[i], A[sw] = A[sw], A[i]; det = -det
        det *= A[i][i]
        piv = A[i][i]
        for k in range(i+1, n):
            factor = A[k][i]/piv
            for j in range(i, n):
                A[k][j] -= factor*A[i][j]
    return int(round(float(det))) if det.denominator==1 else None

def pfaffian_data(n, A):
    """det(I+2A) where 2A is the 0/1 doubled; project convention THM-174: det(I+2A)=Pf(S)^2, S=A-A^T skew."""
    S = [[A[i][j]-A[j][i] for j in range(n)] for i in range(n)]   # skew +-1 off-diagonal
    I_plus = [[(1 if i==j else 0) + S[i][j] for j in range(n)] for i in range(n)]
    d = det_int(I_plus)
    return d

def join_tournament(A1, A2):
    """T1 => T2: block, every vertex of T1 beats every vertex of T2."""
    n1, n2 = len(A1), len(A2); n = n1+n2
    A = [[0]*n for _ in range(n)]
    for i in range(n1):
        for j in range(n1): A[i][j] = A1[i][j]
    for i in range(n2):
        for j in range(n2): A[n1+i][n1+j] = A2[i][j]
    for i in range(n1):
        for j in range(n2): A[i][n1+j] = 1   # T1 beats T2
    return A

print("="*72)
print("H-REALIZABILITY SPECTRUM  (the Goldbach analog: which ODD numbers are H(T)?)")
print("="*72)
Rcum = set()
per_n = {}
sample_class = {}  # H -> one adjacency for join experiments (use n<=5 reps)
for n in range(3, 9):
    Hs = []
    for A in gen_tournaments(n):
        h = H_count(n, A)
        Hs.append(h)
        if n <= 5 and h not in sample_class:
            sample_class[h] = [row[:] for row in A]
    Rn = sorted(set(Hs))
    per_n[n] = set(Rn)
    Rcum |= set(Rn)
    mx = max(Rn)
    odds = [v for v in range(1, mx+1, 2)]
    forbidden = [v for v in odds if v not in set(Rn)]
    allodd = all(h % 2 == 1 for h in Hs)
    print(f"\nn={n}: {len(Rn)} distinct H values, all odd: {allodd}, max H = {mx}")
    print(f"   R_{n} = {Rn}")
    print(f"   forbidden odds < max: {forbidden}")

print("\n" + "="*72)
print("CUMULATIVE realizability R = union R_n  (n<=8): which odds are EVER realized")
print("="*72)
mxc = max(Rcum)
forb_cum = [v for v in range(1, mxc+1, 2) if v not in Rcum]
print(f"max cumulative H = {mxc}")
print(f"odds NEVER realized at any n<=8, below {mxc}: {forb_cum}")
# which n first realizes each odd
print("first n realizing each small odd:")
firsts = {}
for n in range(3, 9):
    for h in sorted(per_n[n]):
        if h not in firsts: firsts[h] = n
for h in sorted(firsts)[:30]:
    print(f"   H={h:5d}: first at n={firsts[h]}")

print("\n" + "="*72)
print("THE JOIN LAW  (dominance composition as the 'addition' operator)")
print("="*72)
print("Testing H(T1 => T2) vs H(T1), H(T2) for small reps:")
reps = {}
for n in (3,4,5):
    for A in gen_tournaments(n):
        h = H_count(n, A)
        reps.setdefault(n, {}).setdefault(h, [row[:] for row in A])
tested = 0
for (n1, h1), A1 in [( (n,h),A) for n in (3,4) for h,A in reps[n].items()]:
    for (n2, h2), A2 in [( (n,h),A) for n in (3,) for h,A in reps[n].items()]:
        J = join_tournament(A1, A2)
        hj = H_count(len(J), J)
        rel = "=prod" if hj == h1*h2 else ("=sum" if hj==h1+h2 else "?")
        print(f"   H({h1})⇒H({h2}) [n={n1}+{n2}] = {hj}   (prod={h1*h2}, sum={h1+h2})  {rel}")
        tested += 1
        if tested >= 12: break
    if tested >= 12: break

print("\n" + "="*72)
print("PFAFFIAN LEDGER  (THM-442: H^2 - Pf^2 = 8Q ; the BSD-shaped split)")
print("="*72)
print("For each n: is det(I+S) a perfect square (Pf integer)? is Q=(H^2-Pf^2)/8 a nonneg integer?")
for n in range(3, 8):
    bad = 0; qneg = 0; cnt = 0; pf_set = set()
    examples = []
    for A in gen_tournaments(n):
        h = H_count(n, A)
        d = pfaffian_data(n, A)   # det(I+S)
        if n % 2 == 1:
            # odd n: S singular, det(I+S) -- project uses det(I+2A); recompute Pf via det(I+2A)
            twoA = [[(1 if i==j else 0) + 2*A[i][j] for j in range(n)] for i in range(n)]
            d = det_int(twoA)
        sq = isqrt(d) if d is not None and d >= 0 else -1
        is_sq = (sq*sq == d) if d is not None and d>=0 else False
        if not is_sq: bad += 1
        else:
            pf = sq
            q8 = h*h - pf*pf
            if q8 % 8 != 0 or q8 < 0: qneg += 1
            pf_set.add(pf)
            if len(examples) < 3: examples.append((h, pf, q8//8 if q8%8==0 else None))
        cnt += 1
    print(f"   n={n}: {cnt} classes, det not a perfect square: {bad}, Q not nonneg-int: {qneg}; "
          f"Pf values seen: {sorted(pf_set)[:12]}{'...' if len(pf_set)>12 else ''}")
    print(f"      examples (H, Pf, Q): {examples}")

print("\nDONE.")
