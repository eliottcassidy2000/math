"""a3_carriers.py (AUDIT) - item 4: the minus carrier d = 31 and the plus carriers d = 3, 21, 35.

For each carrier: (i) visit matrix by plain walk enumeration (a2_census.exe, no cap on visits) -> G, max visits;
(ii) direct Hamiltonian path counts of K[TT_a, TT_b] for 1 <= a, b <= 5 (hp.exe, subset DP on the explicit
tournament of order up to 17) against the closed formula; (iii) arithmetic of the grid claims."""
import subprocess, os, sys
from aud_common import *

out = []


def P(s):
    print(s)
    out.append(s)
    sys.stdout.flush()


def walk_matrix(m, code):
    r = subprocess.run([os.path.join(HERE, "a2_census.exe"), str(m)], input="%d 0\n" % code, capture_output=True, text=True)
    res = {}
    for l in r.stdout.splitlines():
        t = [int(z) for z in l.split()]
        res[(t[2], t[3])] = t[4:]
    return res


carriers = [
    ("minus d = 31", 8, 17082497, (0, 7), lambda a, b: 31 * 2 ** (a + b + 1) - 2 ** (a + 4) - 1, [[215, 124], [108, 62]]),
    ("plus d = 3", 8, 172278145, (0, 7), lambda a, b: 1 + 2 ** (a + 4) + 3 * 2 ** (a + b + 3), None),
    ("plus d = 21", 8, 71647617, (7, 0), lambda a, b: 1 + 2 ** (a + 2) + 21 * 2 ** (a + b + 1), None),
    ("plus d = 35", 9, 17332994050, (8, 1), lambda a, b: 1 + 2 ** (a + 5) + 35 * 2 ** (a + b), None),
]
for name, m, code, (x, y), formula, Gclaim in carriers:
    K = decode(code, m)
    wm = walk_matrix(m, code)
    maxi, maxj, n11, n12, n21, n22, det, two = wm[(x, y)]
    al = n11 - 2 * n12 - 2 * n21 + 4 * n22
    be = n21 - 2 * n22
    ga = n12 - 2 * n22
    de = n22
    P("%s: order %d, code %d, ports (%d, %d): scores %s, H(K) = %d" % (name, m, code, x, y, sorted(len(s) for s in K), ham_c([K])[0]))
    P("   walk enumeration: max visits (%d, %d), G = [[%d, %d], [%d, %d]], det %d, (al, be, ga, de) = (%d, %d, %d, %d)%s"
      % (maxi, maxj, n11, n12, n21, n22, det, al, be, ga, de, ("; claimed G %s: %s" % (Gclaim, Gclaim == [[n11, n12], [n21, n22]])) if Gclaim else ""))
    jobs = []
    for a in range(1, 6):
        for b in range(1, 6):
            jobs.append((a, b, substitute(K, [x, y], [transitive(a), transitive(b)])))
    Hs = ham_c([j[2] for j in jobs])
    bad = [(a, b, h, formula(a, b)) for (a, b, T), h in zip(jobs, Hs) if h != formula(a, b)]
    bad2 = [(a, b) for (a, b, T), h in zip(jobs, Hs) if h != al + be * 2 ** a + ga * 2 ** b + de * 2 ** (a + b)]
    P("   direct counts of H(K[TT_a, TT_b]), 1 <= a, b <= 5 (25 tournaments of order up to %d): mismatches with the closed formula %d, with the visit-matrix form %d"
      % (m + 8, len(bad), len(bad2)))
    if bad:
        P("      first mismatches: %s" % bad[:5])
    P("      row a = 1: %s" % [h for (a, b, T), h in zip(jobs, Hs) if a == 1])
    # non-transitive modules as well (the formula with H, Pi_2)
    import random
    rng = random.Random(5)
    jobs = []
    for t in range(12):
        A = random_tournament(rng.randint(1, 4), rng)
        B = random_tournament(rng.randint(1, 4), rng)
        pa = P_tuple(A, 2)
        pb = P_tuple(B, 2)
        hA, hB, piA, piB = pa[1], pb[1], pa[2] + 2 * pa[1], pb[2] + 2 * pb[1]
        jobs.append((substitute(K, [x, y], [A, B]), al * hA * hB + be * piA * hB + ga * hA * piB + de * piA * piB))
    Hs = ham_c([j[0] for j in jobs])
    P("   random (non-transitive) modules, form in (H, Pi_2): mismatches %d of %d" % (sum(1 for (T, pred), h in zip(jobs, Hs) if h != pred), len(jobs)))

# ---- arithmetic of the grids
P("grid arithmetic:")


def plus_addr(N):
    if N == 1:
        return None
    s = 0
    q = N - 1
    while q % 2 == 0:
        q //= 2
        s += 1
    if q == 1:
        return None
    q -= 1
    r = 0
    while q % 2 == 0:
        q //= 2
        r += 1
    return (q, r, s)


def minus_addr(N):
    s = 0
    q = N + 1
    while q % 2 == 0:
        q //= 2
        s += 1
    q += 1
    r = 0
    while q % 2 == 0:
        q //= 2
        r += 1
    return (q, r, s)


ok = True
for N in range(1, 200001, 2):
    d, r, s = minus_addr(N)
    ok &= (d % 2 == 1 and r >= 1 and s >= 1 and d * 2 ** (r + s) - 2 ** s - 1 == N)
    pa = plus_addr(N)
    if pa is None:
        ok &= (N == 1 or (N - 1) & (N - 2) == 0)
    else:
        d2, r2, s2 = pa
        ok &= (d2 % 2 == 1 and r2 >= 1 and s2 >= 1 and 1 + 2 ** s2 + d2 * 2 ** (r2 + s2) == N)
        if N % 4 == 1:
            ok &= (s == 1 and s2 >= 2)
        else:
            ok &= (s >= 2 and s2 == 1)
# uniqueness by exhaustive search of small triples
seen_p, seen_m = {}, {}
uniq = True
for d in range(1, 400, 2):
    for r in range(1, 12):
        for s in range(1, 12):
            Np = 1 + 2 ** s + d * 2 ** (r + s)
            Nm = d * 2 ** (r + s) - 2 ** s - 1
            if Np in seen_p:
                uniq = False
            if Nm in seen_m:
                uniq = False
            seen_p[Np] = (d, r, s)
            seen_m[Nm] = (d, r, s)
P("   every odd N <= 200001: minus address exists, plus address exists unless N = 2^s + 1 (or N = 1); level rule mod 4: %s; uniqueness of (d, r, s) (d < 400, r, s < 12): %s" % (ok, uniq))
# the minus grid of d = 31: N-(31, r, s) = 31 2^(r+s) - 2^s - 1 with s = a + 4, r = b - 3
cov = set()
for a in range(1, 40):
    for b in range(1, 40):
        cov.add(31 * 2 ** (a + b + 1) - 2 ** (a + 4) - 1)
grid_ok = all((31 * 2 ** (r + s) - 2 ** s - 1) in cov for s in range(5, 30) for r in range(1, 10))
not_cov = [(r, s) for s in range(1, 5) for r in range(1, 6) if (31 * 2 ** (r + s) - 2 ** s - 1) in cov]
P("   minus carrier: N-(31, r, s) covered for all s >= 5, r >= 1 (checked s < 30, r < 10): %s; levels s <= 4 covered: %s" % (grid_ok, not_cov))
lines = {}
for b in (1, 2, 3):
    vals = [31 * 2 ** (a + b + 1) - 2 ** (a + 4) - 1 for a in range(1, 6)]
    c, s0 = None, None
    v = vals[0] + 1
    s0 = 0
    while v % 2 == 0:
        v //= 2
        s0 += 1
    lines[b] = (v, s0, all((x + 1) == v * 2 ** (s0 + i) for i, x in enumerate(vals)))
P("   rows b = 1, 2, 3 are the lines c 2^s - 1 with (c, first s, is a line): %s  [claimed (27, 3), (29, 4), (15, 6)]" % lines)
P("   c = 31 2^r - 1 for r = 1..4: %s  [claimed 61, 123, 247, 495]" % [31 * 2 ** r - 1 for r in range(1, 5)])
P("   smallest members: N-(31,1,5) = %d, 27*8-1 = %d, 29*16-1 = %d, 15*64-1 = %d" % (31 * 64 - 32 - 1, 215, 463, 959))
open(os.path.join(HERE, "a3_carriers.out"), "w").write("\n".join(out) + "\n")
