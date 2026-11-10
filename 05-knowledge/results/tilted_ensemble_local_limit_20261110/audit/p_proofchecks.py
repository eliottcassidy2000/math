#!/usr/bin/env python3
"""AUDIT: machine checks of individual proof steps of the packet (written from scratch by the auditor).
 1. Theorem C (1): the insertion construction, executed literally on every labelled tournament of order 2..6, every
    pair {i,j} and every Hamiltonian path Q of R - {i,j} (including the empty Q at m = 2).
 2. Two-path-cover matrix M by brute force over permutations (m <= 5): mu even, mu >= 2, the stronger KNOWN bound
    mu >= 2 H(R - {i,j}) (Session 136 positive_cut_geometry (3); THM-5122 "Capacity"), and the reading
    "mu = paths of R through the arc + paths of the flipped tournament through the flipped arc".
 3. Lemma C4: both sides evaluated exactly/numerically for random R, U, V', theta; the theta = 1/4 and theta = a/q cases.
 4. Proposition 6.1 arithmetic for every odd n in [1025, 3*10^6].
 5. Proposition 5.1: monotonicity of mu_n, existence of the tilt in [1/2, p_n).
 6. Bridging tilts p_n and the Eulerian asymptotic (which normalisation is the right one).
 7. Rank distribution of alternating matrices over F_2: brute force (m <= 6) against the closed formula; table values.
"""
import itertools, math, random, cmath
from fractions import Fraction as Fr
from math import comb, factorial

def tournaments(m):
    prs = [(i, j) for i in range(m) for j in range(i + 1, m)]
    for code in range(1 << len(prs)):
        adj = [[0] * m for _ in range(m)]
        for k, (i, j) in enumerate(prs):
            if code >> k & 1: adj[i][j] = 1
            else: adj[j][i] = 1
        yield adj
def is_path(adj, P): return all(adj[P[t]][P[t + 1]] for t in range(len(P) - 1))
def ham_paths(adj, verts):
    return [P for P in itertools.permutations(verts) if is_path(adj, P)]
def H_of(adj, verts): return len(ham_paths(adj, verts))

# ------------------------------------------------------------------ 1. insertion argument of Theorem C (1)
print("== 1. Theorem C (1), insertion construction executed literally ==")
for m in range(2, 7):
    bad = 0; cases = 0; c_is_r = 0; c_is_0 = 0
    for adj in tournaments(m):
        for i in range(m):
            for j in range(i + 1, m):
                rest = [v for v in range(m) if v != i and v != j]
                Qs = ham_paths(adj, rest) if rest else [()]
                if not Qs: bad += 1          # Redei: there must be a Hamiltonian path of R - {i,j}
                for Q in Qs:
                    r = len(Q); cases += 1
                    c = next(c for c in range(r + 1) if c == r or adj[i][Q[c]] or adj[j][Q[c]])   # Q[c] is q_{c+1}
                    if not all(adj[Q[t]][i] and adj[Q[t]][j] for t in range(c)): bad += 1          # minimality claim
                    if c < r:
                        s = i if adj[i][Q[c]] else j; e = j if s == i else i
                    else:
                        e, s = i, j; c_is_r += 1
                    if c == 0: c_is_0 += 1
                    P1 = tuple(Q[:c]) + (e,); P2 = (s,) + tuple(Q[c:])
                    if not (is_path(adj, P1) and is_path(adj, P2) and sorted(P1 + P2) == list(range(m)) and len(P1) >= 1 and len(P2) >= 1): bad += 1
    print(f"   m={m}: {cases} (R, pair, Q) cases, c = r in {c_is_r}, c = 0 in {c_is_0}; invalid constructions: {bad}")

# ------------------------------------------------------------------ 2. M by brute force, parity, positivity, readings
print("== 2. two-path covers by brute force over permutations (m = 2..5) ==")
for m in range(2, 6):
    fail = 0; minratio = None; n_eq = 0; n_pairs = 0; minmu = 10 ** 9
    for adj in tournaments(m):
        M = [[0] * m for _ in range(m)]
        for P in itertools.permutations(range(m)):
            for t in range(1, m):
                if is_path(adj, P[:t]) and is_path(adj, P[t:]): M[P[t - 1]][P[t]] += 1
        Hs = ham_paths(adj, range(m))
        for i in range(m):
            for j in range(i + 1, m):
                mu = M[i][j] + M[j][i]; n_pairs += 1; minmu = min(minmu, mu)
                rest = [v for v in range(m) if v not in (i, j)]
                hrest = H_of(adj, rest) if rest else 1
                if mu % 2 or mu < 2 or mu < 2 * hrest: fail += 1
                if mu == 2 * hrest: n_eq += 1
                # reading: paths of R through the arc between i and j + paths of the flipped tournament through the flipped arc
                a, b = (i, j) if adj[i][j] else (j, i)
                thr = sum(1 for P in Hs if any(P[t] == a and P[t + 1] == b for t in range(m - 1)))
                adj2 = [row[:] for row in adj]; adj2[a][b] = 0; adj2[b][a] = 1
                thr2 = sum(1 for P in ham_paths(adj2, range(m)) if any(P[t] == b and P[t + 1] == a for t in range(m - 1)))
                if thr != M[a][b] or thr2 != M[b][a] or thr + thr2 != mu: fail += 1
                if (len(ham_paths(adj2, range(m))) - len(Hs)) != thr2 - thr: fail += 1        # flip increment = C - B
    print(f"   m={m}: {n_pairs} (R, pair) cases; failures {fail}; min mu {minmu}; equality mu = 2 H(R-{{i,j}}) in {n_eq} cases")

# ------------------------------------------------------------------ 3. Lemma C4
print("== 3. Lemma C4 (random R on 5..7 vertices, random disjoint U, V', several theta) ==")
def ham_dp(n, adj):
    out = [sum(1 << w for w in range(n) if adj[v][w]) for v in range(n)]
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n): dp[1 << v][v] = 1
    for S in range(1, 1 << n):
        for v in range(n):
            c = dp[S][v]
            if c:
                nx = out[v] & ~S
                while nx:
                    w = (nx & -nx).bit_length() - 1; nx &= nx - 1; dp[S | (1 << w)][w] += c
    return sum(dp[(1 << n) - 1])
def ext_values(m, adj):
    vals = []
    for X in range(1 << m):
        a2 = [row[:] + [0] for row in adj] + [[0] * (m + 1)]
        for u in range(m):
            if X >> u & 1: a2[m][u] = 1
            else: a2[u][m] = 1
        vals.append(ham_dp(m + 1, a2))
    return vals
def rank2(rows, ncols):
    rows = list(rows); r = 0
    for c in range(ncols):
        p = next((k for k in range(r, len(rows)) if rows[k] >> c & 1), None)
        if p is None: continue
        rows[r], rows[p] = rows[p], rows[r]
        for k in range(len(rows)):
            if k != r and rows[k] >> c & 1: rows[k] ^= rows[r]
        r += 1
    return r
rng = random.Random(4242); viol = 0; tests = 0; q14 = 0; qprime = 0; tight = 0; worst = 0.0
for trial in range(60):
    m = rng.choice([5, 6, 7]); adj = [[0] * m for _ in range(m)]
    for i in range(m):
        for j in range(i + 1, m):
            if rng.random() < 0.5: adj[i][j] = 1
            else: adj[j][i] = 1
    vals = ext_values(m, adj); H0 = vals[0]
    a = [vals[1 << i] - H0 for i in range(m)]
    mu = [[0 if i == j else a[i] + a[j] - (vals[(1 << i) | (1 << j)] - H0) for j in range(m)] for i in range(m)]
    perm = list(range(m)); rng.shuffle(perm); ku = rng.randint(1, m - 1); kv = rng.randint(1, m - ku)
    U = perm[:ku]; Vp = perm[ku:ku + kv]
    thetas = [Fr(1, 4), Fr(1, 3), Fr(2, 5), Fr(3, 7), Fr(1, 8), Fr(5, 12), rng.random(), rng.random()]
    for th in thetas:
        lhs = abs(sum(cmath.exp(2j * math.pi * float(th) * v) for v in vals) / len(vals)) ** 2
        rhs = 0.0; bound_q = 0.0
        for d in itertools.product((-1, 0, 1), repeat=len(Vp)):
            wgt = 1.0
            for x in d: wgt *= 0.5 if x == 0 else 0.25
            prod = 1.0; N = 0
            for i in U:
                s = sum(mu[i][j] * dj for j, dj in zip(Vp, d))
                prod *= abs(math.cos(math.pi * float(th) * s))
                if isinstance(th, Fr) and th.denominator % 2 == 1 and s % th.denominator != 0: N += 1
            rhs += wgt * prod
            if isinstance(th, Fr) and th.denominator % 2 == 1: bound_q += wgt * math.cos(math.pi / th.denominator) ** N
        tests += 1
        if lhs > rhs + 1e-9: viol += 1
        worst = max(worst, lhs - rhs)
        if abs(lhs - rhs) < 1e-9: tight += 1
        if th == Fr(1, 4):
            block = [sum((((mu[i][j] // 2) & 1) << k) for k, j in enumerate(Vp)) for i in U]
            if abs(rhs - 2.0 ** (-rank2(block, len(Vp)))) > 1e-9: q14 += 1
            full = [sum((((mu[i][j] // 2) & 1) << j) for j in range(m) if j != i) for i in range(m)]
            if lhs > 2.0 ** (-rank2(full, m)) + 1e-9: q14 += 1        # consistency with Theorem C (3)
        if isinstance(th, Fr) and th.denominator % 2 == 1 and rhs > bound_q + 1e-9: qprime += 1
print(f"   {tests} (R, U, V', theta) tests: violations of the inequality {viol} (largest LHS - RHS = {worst:.3e}); cases with equality {tight}")
print(f"   theta = 1/4: right side != 2^(-rank of the U x V' block of beta) or LHS > 2^(-rank beta): {q14} failures")
print(f"   theta = a/q (q = 3, 5, 7, q not dividing a): right side > E_d cos(pi/q)^N(d): {qprime} failures")

# ------------------------------------------------------------------ 4. Proposition 6.1 arithmetic
print("== 4. Proposition 6.1 arithmetic ==")
bad = 0; offs = set(); hmin = None
for n in range(1025, 3_000_001, 2):
    h = n // 3; h -= (h % 2 == 0)
    if not (h % 2 == 1 and 3 * h <= n and 3 * (h + 2) > n): bad += 1
    if not (h >= 341 and h >= 237 and n <= 3 * h + 5 and n < 5 * h and h < n and h not in (7, 21)): bad += 1
    offs.add(n - 3 * h); hmin = h if hmin is None else min(hmin, h)
print(f"   odd n in [1025, 3*10^6]: failures {bad}; n - 3h takes the values {sorted(offs)}; least h = {hmin} (at n = 1025: h = {1025 // 3 - (1 if (1025 // 3) % 2 == 0 else 0)})")
print(f"   base range: largest odd n below 1025 is 1023 = 3 * 341, so the base 'all odd numbers below 1024' is exactly what the step needs")

# ------------------------------------------------------------------ 5. Proposition 5.1
print("== 5. Proposition 5.1 ==")
mu_n = [None] + [Fr(factorial(n), 2 ** (n - 1)) for n in range(1, 40)]
print("   mu_1..mu_6 =", [str(x) for x in mu_n[1:7]], "; nondecreasing from n=1:", all(mu_n[n] <= mu_n[n + 1] for n in range(1, 38)),
      "; strictly increasing from n=2:", all(mu_n[n] < mu_n[n + 1] for n in range(2, 38)), "; mu_1 = mu_2:", mu_n[1] == mu_n[2])
def gamma_vec(n):
    A = [sum((-1) ** i * comb(n + 1, i) * (k + 1 - i) ** n for i in range(k + 2)) for k in range(n)]
    d = n - 1; c = A[:]; g = []
    for j in range(d // 2 + 1):
        g.append(c[j])
        for i in range(d - 2 * j + 1): c[j + i] -= g[j] * comb(d - 2 * j, i)
    assert not any(c)
    return g
ok5 = True
for n in range(3, 13):
    g = gamma_vec(n); mean = lambda w: sum(gj * w ** j for j, gj in enumerate(g))
    assert mean(Fr(1, 4)) == mu_n[n] and mean(Fr(0)) == 1
    hs = [h for h in range(1, int(mu_n[n]) + 1, 2) if mu_n[n - 1] < h <= mu_n[n]]
    for h in hs:                      # the tilt with mean h has w in (w_n, 1/4], i.e. p in [1/2, p_n)
        if not (mean(Fr(1, 4)) >= h > mu_n[n - 1]): ok5 = False
print("   for n = 3..12 and every odd h in (mu_(n-1), mu_n]: mu_(n-1) < h <= E_(1/2)[H_n], so by strict monotonicity the tilt lies in [1/2, p_n):", ok5)
print("   note: p_n is undefined at n = 2 (E_p[H_2] = 1 for all p) and p_3 = 1 (mu_2 = 1); the hypothesis needs n_0 >= 3")

# ------------------------------------------------------------------ 6. bridging tilt and the Eulerian asymptotic
print("== 6. bridging tilt p_n and the Eulerian asymptotic ==")
from decimal import Decimal, getcontext
getcontext().prec = 80
for n in (4, 6, 10, 16, 40, 80, 160):
    g = gamma_vec(n); target = Fr(factorial(n - 1), 2 ** (n - 2))
    lo, hi = Fr(0), Fr(1, 4)
    for _ in range(120):
        mid = (lo + hi) / 2
        if sum(gj * mid ** j for j, gj in enumerate(g)) > target: hi = mid
        else: lo = mid
    w = Decimal(lo.numerator) / Decimal(lo.denominator); p = Decimal("0.5") + (Decimal("0.25") - w).sqrt()
    pf = float(p); print(f"   n={n}: p_n = {pf:.5f}; (p_n - 1/2) sqrt(n/ln n) = {(pf - 0.5) * math.sqrt(n / math.log(n)):.4f}; sqrt(3 ln(n/2)/(4n)) = {math.sqrt(3 * math.log(n / 2) / (4 * n)):.4f}; p_n - 1/2 = {pf - 0.5:.4f}")
for n in (20, 60):
    A = [sum((-1) ** i * comb(n + 1, i) * (k + 1 - i) ** n for i in range(k + 2)) for k in range(n)]
    for t in (0.5, 0.8):
        exact = sum(A[k] * t ** k for k in range(n)) / factorial(n)
        asy = ((t - 1) / math.log(t)) ** (n + 1)
        print(f"   n={n}, t={t}: A_n(t)/n! = {exact:.6e}; ((t-1)/ln t)^(n+1) = {asy:.6e} (ratio {exact / asy:.6f}); with the factor 1/t: ratio {exact / (asy / t):.6f}")

# ------------------------------------------------------------------ 7. alternating matrices over F_2
print("== 7. rank distribution of alternating matrices over F_2 ==")
def formula(m):
    out = {}
    for r in range(m // 2 + 1):
        x = Fr(1)
        for i in range(1, r + 1): x *= Fr(2 ** (2 * i - 2), 2 ** (2 * i) - 1)
        for i in range(2 * r): x *= (2 ** (m - i) - 1)
        out[2 * r] = x
    return out
for m in (2, 3, 4, 5, 6):
    prs = [(i, j) for i in range(m) for j in range(i + 1, m)]; cnt = {}
    for code in range(1 << len(prs)):
        rows = [0] * m
        for k, (i, j) in enumerate(prs):
            if code >> k & 1: rows[i] |= 1 << j; rows[j] |= 1 << i
        r = rank2(rows, m); cnt[r] = cnt.get(r, 0) + 1
    f = formula(m)
    print(f"   m={m}: brute force {dict(sorted(cnt.items()))}; formula {({k: int(v) for k, v in f.items()})}; equal: {all(cnt.get(k, 0) == v for k, v in f.items())}")
census = {8: (6880, {2: 2, 4: 205, 6: 3795, 8: 2878}), 9: (191536, {2: 8, 4: 436, 6: 32535, 8: 158557})}
labelled = {8: (2 ** 28, {2: 53760, 4: 7728000, 6: 147784576, 8: 112869120}), 9: (2 ** 36, {2: 1612800, 4: 143781120, 6: 11631474304, 8: 56942608512})}
for m in (8, 9, 10):
    f = formula(m); tot = 2 ** (m * (m - 1) // 2); top = 2 * (m // 2)
    pr = {k: float(v / tot) for k, v in f.items()}
    line = f"   m={m}: random: top rank {pr[top]:.4f}, next {pr[top - 2]:.4f}, lower {sum(v for k, v in pr.items() if k < top - 2):.4f}; rank 2 probability {pr[2]:.3e}; E[2^(-rank/2)] = {sum(v * 2.0 ** (-k / 2) for k, v in pr.items()):.6f}"
    print(line)
    if m in census:
        N, c = census[m]; N2, c2 = labelled[m]
        print(f"         classes: top {c[top] / N:.4f}, next {c[top - 2] / N:.4f}, lower {sum(v for k, v in c.items() if k < top - 2) / N:.4f}; rank 2 fraction {c[2] / N:.3e} (ratio to random {c[2] / N / pr[2]:.1f}); "
              f"lower/random = {sum(v for k, v in c.items() if k < top - 2) / N / sum(v for k, v in pr.items() if k < top - 2):.2f}")
        print(f"         labelled (auditor, weights m!/|Aut|): top {c2[top] / N2:.4f}, next {c2[top - 2] / N2:.4f}, lower {sum(v for k, v in c2.items() if k < top - 2) / N2:.4f}; "
              f"largest |census - random| over the three bins: classes {max(abs(c[top] / N - pr[top]), abs(c[top - 2] / N - pr[top - 2]), abs(sum(v for k, v in c.items() if k < top - 2) / N - sum(v for k, v in pr.items() if k < top - 2))) * 100:.2f} pp, "
              f"labelled {max(abs(c2[top] / N2 - pr[top]), abs(c2[top - 2] / N2 - pr[top - 2]), abs(sum(v for k, v in c2.items() if k < top - 2) / N2 - sum(v for k, v in pr.items() if k < top - 2))) * 100:.2f} pp")
