# The two peaks of the circulant H-landscape against the fourth spectral moment M4 = sum_{j != 0} |Im lambda_j|^4
# (equivalently c_4, c_5: on circulant tournaments Sum y_j^2 is fixed, c_5 = const - M4/2, c_4 = const + M4/4).
# Prior art: THM-162 (kind-pasteur S60: Paley minimises M4, corr(M4, c_5) = -1 at p = 11).  New here: which class
# MAXIMISES M4, and how H correlates with M4 across the complete class tables n <= 23.  mac-mini-2026-10-21-fanopsp.
import re, glob, cmath, math
from fractions import Fraction
from math import gcd
from collections import defaultdict
rows = defaultdict(dict)
for f in sorted(glob.glob('../tournament_pseudoprimes_20261021/circ_*.out')):
    for line in open(f):
        m_ = re.match(r'n=(\d+) S=\{([\d,]+)\} mult=(\d+) H=(\d+)', line)
        if m_: rows[int(m_.group(1))][tuple(sorted(int(x) for x in m_.group(2).split(',')))] = int(m_.group(4))
def canon(n, S):
    return min(tuple(sorted((a * s) % n for s in S)) for a in range(1, n) if gcd(a, n) == 1)
def M4(n, S):  # exact via integer cycle counts would need traces; float is fine for ranking (values well separated)
    tot = 0.0
    for j in range(1, n):
        y = sum(math.sin(2 * math.pi * j * s / n) for s in S)
        tot += y ** 4
    return tot
def ranks(v):
    order = sorted(range(len(v)), key=lambda i: v[i]); r = [0] * len(v)
    for k, i in enumerate(order): r[i] = k
    return r
def spearman(a, b):
    ra, rb = ranks(a), ranks(b); n = len(a); ma = sum(ra) / n; mb = sum(rb) / n
    num = sum((x - ma) * (y - mb) for x, y in zip(ra, rb))
    den = math.sqrt(sum((x - ma) ** 2 for x in ra) * sum((y - mb) ** 2 for y in rb))
    return num / den if den else float('nan')
for n in sorted(rows):
    m = (n - 1) // 2
    H = {canon(n, S): h for S, h in rows[n].items()}
    R = canon(n, range(1, m + 1))
    Ss = sorted(H); hv = [H[S] for S in Ss]; mv = [M4(n, S) for S in Ss]
    imax = max(range(len(Ss)), key=lambda i: mv[i]); imin = min(range(len(Ss)), key=lambda i: mv[i])
    hmax = max(range(len(Ss)), key=lambda i: hv[i])
    tagmax = 'R' if Ss[imax] == R else str(list(Ss[imax]))
    print(f"n={n:2d} classes={len(Ss):4d}  argmax M4 = {tagmax:>6}  argmin M4 = {list(Ss[imin])}  "
          f"argmax H = {'R' if Ss[hmax] == R else list(Ss[hmax])}  spearman(H, M4) = {spearman(hv, mv):+.3f}")

# Exact check (integers): additive energy E(S) = #{s1+s2 = s3+s4 mod n}, 4-cycle count c4 = #{s1+s2+s3+s4 = 0}/... per
# vertex, and 5-cycle count c5 = n*#{s1+..+s5 = 0 mod n}/5; uniqueness of the extremes and the affine relations.
print()
def zero_sum_count(n, S, k):
    cnt = [0] * n; cnt[0] = 1
    for _ in range(k):
        new = [0] * n
        for r in range(n):
            if cnt[r]:
                for s in S: new[(r + s) % n] += cnt[r]
        cnt = new
    return cnt[0]
for n in sorted(rows):
    if n < 7: continue
    m = (n - 1) // 2
    H = {canon(n, S): h for S, h in rows[n].items()}
    R = canon(n, range(1, m + 1))
    data = []
    for S in sorted(H):
        Sset = set(S)
        E = sum(1 for a in S for b in S for c in S if ((a + b - c) % n) in Sset)
        c4 = n * zero_sum_count(n, S, 4) // 4
        c5 = n * zero_sum_count(n, S, 5) // 5
        data.append((E, c4, c5, H[S], S))
    data.sort()
    Es = [d[0] for d in data]
    # affine relations across the class table: c4 - nE/4 and c5 + nE/2 constant?
    rel4 = {d[1] - Fraction(n * d[0], 4) for d in data}; rel5 = {d[2] + Fraction(n * d[0], 2) for d in data}
    top, bot = data[-1], data[0]
    print(f"n={n:2d}: E range [{Es[0]}, {Es[-1]}]; unique max {'R' if top[4] == R else list(top[4])} (next {Es[-2]}), "
          f"unique min {list(bot[4])} (next {Es[1]}); c4 - nE/4 constant: {len(rel4) == 1}; c5 + nE/2 constant: {len(rel5) == 1}; "
          f"H-max at E-{'max' if max(data, key=lambda d: d[3])[4] == top[4] else ('min' if max(data, key=lambda d: d[3])[4] == bot[4] else 'interior')}")
