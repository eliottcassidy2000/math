#!/usr/bin/env python3
"""audit_C item 4: Mersenne orphan statistics, recomputed from audit_C's own orbit data (c4_mersenne_12800_audit.txt,
produced by c4_mersenne_orbits.py; byte-identical to the session's mersenne_sigma_12800.txt).

Partner criterion (HYP-9242): odd K has a deletion partner iff some K' < K has the same odd-step count and Terras
difference sigma_T(K) - sigma_T(K') = K - K'.

Outputs: window table; exponent fits (Bernoulli MLE with profile-likelihood intervals; model comparison); the 'last three
doublings' local exponent with a Poisson interval; nearest-partner distribution; shells; length terciles; frac(K log2 3)
quarters; residues mod 6 of orphans.
"""
import math, bisect
from collections import Counter, defaultdict

data = {}
for line in open('c4_mersenne_12800_audit.txt'):
    K, o, t = map(int, line.split())
    data[K] = (o, t)
KMAX = max(data)

by_odd = defaultdict(list)
orphan = {}
partners = {}
for K in range(2, KMAX + 1):
    o, t = data[K]
    ps = [Kp for Kp in by_odd[o] if t - data[Kp][1] == K - Kp]
    if K % 2 == 1 and K >= 3:
        orphan[K] = (len(ps) == 0)
        partners[K] = ps
    by_odd[o].append(K)

print("== windows (odd K) ==")
wins = [(3, 50), (50, 100), (100, 200), (200, 400), (400, 800), (800, 1600), (1600, 3200), (3200, 6400), (6400, 12801)]
print(" window          #oddK  orphans  fraction  x sqrt(geo-mid)  x sqrt(arith-mid)   Poisson 95% (count)")
for a, b in wins:
    ks = [K for K in orphan if a <= K < b]
    no = sum(orphan[K] for K in ks)
    f = no / len(ks)
    g = math.sqrt(a * b)
    print(f" [{a:5d},{b:5d})   {len(ks):5d}   {no:5d}   {f:.4f}    {f*math.sqrt(g):.3f}           {f*math.sqrt((a+b)/2):.3f}"
          f"            [{max(0, no - 1.96*math.sqrt(no)):.1f}, {no + 1.96*math.sqrt(no):.1f}]")

# ---------------- exponent fits ----------------
import numpy as np

def fit_C(shape, Ks):
    """1-parameter Bernoulli model p = C*shape(K); maximise the log-likelihood over C (golden search on log C)"""
    s = np.array([shape(K) for K in Ks], dtype=float)
    y = np.array([orphan[K] for K in Ks], dtype=bool)
    def nll(lc):
        p = np.clip(np.exp(lc) * s, 1e-12, 1 - 1e-12)
        return -(np.log(p[y]).sum() + np.log1p(-p[~y]).sum())
    gr = (math.sqrt(5) - 1) / 2
    a, b = -12.0, 8.0
    c, d = b - gr * (b - a), a + gr * (b - a)
    for _ in range(90):
        if nll(c) < nll(d):
            b = d
        else:
            a = c
        c, d = b - gr * (b - a), a + gr * (b - a)
    lc = (a + b) / 2
    return math.exp(lc), -nll(lc)

for lo_K in (200, 400, 800):
    Ks = [K for K in orphan if lo_K <= K <= KMAX]
    print(f"\n== Bernoulli MLE on odd K in [{lo_K}, {KMAX}] (n = {len(Ks)}, orphans = {sum(orphan[K] for K in Ks)}) ==")
    prof = []
    for i in range(0, 161):
        al = 0.1 + i * 0.01
        C, ll = fit_C(lambda K, al=al: K ** (-al), Ks)
        prof.append((al, C, ll))
    best = max(prof, key=lambda x: x[2])
    inside = [al for al, C, ll in prof if ll >= best[2] - 1.92]
    print(f"  power law C K^-alpha: alpha_hat = {best[0]:.2f} (C = {best[1]:.2f}), 95% profile interval [{min(inside):.2f}, {max(inside):.2f}]")
    ll_half = [x for x in prof if abs(x[0] - 0.5) < 1e-9][0][2]
    print(f"  log-likelihood drop at alpha = 1/2: {best[2] - ll_half:.2f} (1.92 = 95% threshold)")
    models = {
        'C K^-1/2': lambda K: K ** -0.5,
        'C K^-0.6': lambda K: K ** -0.6,
        'C K^-0.7': lambda K: K ** -0.7,
        'C log K / K^0.6': lambda K: math.log(K) / K ** 0.6,
        'C log K / K^0.7': lambda K: math.log(K) / K ** 0.7,
        'C / (K^0.5 log K)': lambda K: 1 / (math.sqrt(K) * math.log(K)),
        'C (log K)^2 / K^0.75': lambda K: math.log(K) ** 2 / K ** 0.75,
    }
    for name, shape in models.items():
        C, ll = fit_C(shape, Ks)
        print(f"  {name:22s}: C = {C:8.4f}, log-lik = {ll:9.2f} (best power law {best[2]:9.2f}; difference {best[2]-ll:+.2f})")

# local exponent over the last three doublings, as in HYP-9242
def cnt(a, b):
    ks = [K for K in orphan if a <= K < b]
    return sum(orphan[K] for K in ks), len(ks)
n1, N1 = cnt(800, 1600)
n4, N4 = cnt(6400, 12801)
f1, f4 = n1 / N1, n4 / N4
alpha_loc = math.log(f1 / f4) / math.log(8)
# crude interval: +-1.96 sqrt(1/n1 + 1/n4) on the log ratio
se = math.sqrt(1 / n1 + 1 / n4)
print(f"\nlocal exponent [800,1600) -> [6400,12800]: {alpha_loc:.3f}, 95% interval [{(math.log(f1/f4)-1.96*se)/math.log(8):.2f}, {(math.log(f1/f4)+1.96*se)/math.log(8):.2f}]")

# ---------------- nearest partner ----------------
c = Counter(K - max(partners[K]) for K in orphan if K >= 1000 and not orphan[K])
nn = sum(c.values())
print(f"\nnearest partner, odd K >= 1000: non-orphans {nn}, D = 1 for {c[1]} ({c[1]/nn:.3f}); top: {sorted(c.items())[:8]}")

# ---------------- shells (odd K in [600, 6400]) ----------------
def has_D(K, Ds):
    o, t = data[K]
    return any(K - D >= 2 and data[K - D] == (o, t - D) for D in Ds)
rows = defaultdict(lambda: [0, 0, 0, 0])
for K in range(601, 6401, 2):
    m = (K - 1 & -(K - 1)).bit_length() - 1
    key = min(m, 6)
    r = rows[key]
    r[0] += 1
    r[1] += has_D(K, (1, 2))
    r[2] += has_D(K, (3, 4))
    r[3] += not orphan[K]
print("\nshells m = v2(K-1), odd K in [600, 6400]:  #K, D in {1,2}, D in {3,4}, any D")
for m in sorted(rows):
    cc, a, b, d = rows[m]
    print(f"  m={m}{'+' if m == 6 else ' '}: {cc:5d}  {a/cc:.3f}  {b/cc:.3f}  {d/cc:.3f}")

# ---------------- length terciles ----------------
def terciles(lo, hi):
    rows = [((data[K][1] - K) / K, orphan[K]) for K in orphan if lo <= K <= hi]
    xs = sorted(x for x, _ in rows)
    q1, q2 = xs[len(xs) // 3], xs[2 * len(xs) // 3]
    out = []
    for a, b in ((-1, q1), (q1, q2), (q2, 1e9)):
        sel = [o for x, o in rows if a <= x < b]
        out.append((sum(sel), len(sel)))
    return q1, q2, out
for lo, hi in ((1000, 6400), (1000, 12800)):
    q1, q2, out = terciles(lo, hi)
    print(f"\nlength terciles of (sigma_T - K)/K, odd K in [{lo}, {hi}] (cuts {q1:.3f}, {q2:.3f}):",
          ", ".join(f"{a}/{b} = {100*a/b:.2f}%" for a, b in out))

# ---------------- frac(K log2 3) quarters ----------------
L23 = math.log2(3)
for lo, hi in ((1000, 6400), (1000, 12800), (200, 12800)):
    q = defaultdict(lambda: [0, 0])
    for K in orphan:
        if lo <= K <= hi:
            f = (K * L23) % 1.0
            qq = int(f * 4)
            q[qq][0] += orphan[K]
            q[qq][1] += 1
    print(f"frac(K log2 3) quarters, odd K in [{lo}, {hi}]:",
          ", ".join(f"Q{i+1}: {q[i][0]}/{q[i][1]} = {100*q[i][0]/q[i][1]:.2f}%" for i in range(4)))

# ---------------- residues mod 6 (depth-3 rescue class K = 5 mod 6) ----------------
for lo, hi in ((1000, 6400),):
    orph = [K for K in orphan if lo <= K <= hi and orphan[K]]
    allk = [K for K in orphan if lo <= K <= hi]
    print(f"\nK mod 6 among orphans in [{lo},{hi}]:", dict(sorted(Counter(K % 6 for K in orph).items())),
          " all odd K:", dict(sorted(Counter(K % 6 for K in allk).items())))
    print("orphans in [1000, 6400]:", orph)
print("\nall orphans K >= 6400:", [K for K in orphan if K >= 6400 and orphan[K]])
