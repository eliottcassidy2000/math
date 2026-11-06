# Audit S83, step 9: (a) the corollary constants under the PRESCRIBED split recipes of the note, versus the optimal split;
# (b) the analytic worst cases derived in the audit; (c) constants of the audit's two-block bound (U1).
import math
import numpy as np

def thmA_split(n, q, a, L, l):
    R = 2 * (2.0 ** L / q + 1) * math.sqrt((3 + 2.0 ** (3 - a) * q) * (3 + 2.0 ** (2 - l) * q))
    return 2.0 ** (-L / 2) * math.sqrt(1 + R)

def balanced(n, q, L):
    N = n - L; best = math.inf
    for a in range(max(2, (N + 1) // 2 - 1), min(N - 1, (N + 1) // 2 + 2) + 1):
        l = N - a
        if a >= 2 and l >= 1:
            best = min(best, thmA_split(n, q, a, L, l))
    return best

def U1(n, q):
    lq = 1 + math.log(q)
    return math.sqrt(min(6 / q + 7 * 2.0 ** (a - n) + 8 * lq * 2.0 ** (-a) + 8 * q * lq * 2.0 ** (-n) for a in range(2, n)))

# (a) recipes
infeasible = 0; feasible_max = 0.0; tested = 0
A1_floor = 0.0; A2_floor = (0.0, None); A2_floor1 = (0.0, None)
for n in range(6, 301):
    for x in np.arange(1.0, (n - 3) / 3 + 1e-9, 1 / 64):          # (A1) range
        q = 2.0 ** x; L = math.floor(x)
        need = math.ceil(x) + 2                                      # a, l >= log2 q + 2 (integers)
        tested += 1
        if L + 2 * need > n:
            infeasible += 1
        else:
            a = n - L - need; l = need
            feasible_max = max(feasible_max, thmA_split(n, q, a, L, l) * math.sqrt(q))
        A1_floor = max(A1_floor, balanced(n, q, L) * math.sqrt(q))
    xs = np.concatenate([np.arange((n - 6) / 3, n - 6 + 1e-9, 1 / 64), [(n - 6) / 3 + 1 - 1e-9, (n - 6) / 3 + 2 / 3 - 1e-9]])
    for x in xs:                                                     # (A2) range
        if x < 1 or x > n - 6:
            continue
        q = 2.0 ** x; norm = (q / 2.0 ** n) ** 0.25
        r0 = balanced(n, q, math.floor(x)) / norm
        r1 = balanced(n, q, math.floor(x) + 1) / norm
        if r0 > A2_floor[0]: A2_floor = (r0, (n, round(x, 4)))
        if r1 > A2_floor1[0]: A2_floor1 = (r1, (n, round(x, 4)))
print(f"(A1) note's recipe L = floor(log2 q), a, l >= log2 q + 2: infeasible in {infeasible} of {tested} grid points of the range;"
      f" where feasible, max bound*sqrt(q) = {feasible_max:.4f}")
print(f"(A1) L = floor(log2 q), balanced a: max bound*sqrt(q) = {A1_floor:.4f}   [audit analytic sup 5.369]")
print(f"(A2) note's recipe L = floor(log2 q), balanced a: max bound/(q/2^n)^(1/4) = {A2_floor[0]:.4f} at (n, log2 q) = {A2_floor[1]}   [claimed 12]")
print(f"(A2) L = floor(log2 q) + 1, balanced a: max = {A2_floor1[0]:.4f} at {A2_floor1[1]}   [audit analytic sup 11.30]")

# (b) analytic worst-case functions
c = 3 / math.sqrt(2)
phi1 = max(2 ** (f / 2) * math.sqrt(1 + 2 * (2 ** -f + 1) * math.sqrt(9 + 3 * c * 2 ** (1 - f / 2) + 2 ** (2 - f))) for f in np.linspace(0, 1, 100001))
psi_ceil = max(math.sqrt(2 ** (3 - g) * (1 + 2 * (2 ** g + 1) * math.sqrt(9 + 3 * c * 2 ** ((g - 1) / 2) + 2 ** (g - 1)))) for g in np.linspace(1e-9, 1, 100001))
psi_floor = max(math.sqrt(2 ** (3 + f) * (1 + 2 * (2 ** -f + 1) * math.sqrt(9 + 3 * c * 2 ** (-(1 + f) / 2) + 2 ** (-1 - f)))) for f in np.linspace(0, 1, 100001))
print(f"analytic: (A1) sup_f phi1(f) = {phi1:.4f}; (A2) with L = floor+1: {psi_ceil:.4f}; with L = floor: {psi_floor:.4f}")

# (c) U1 constants (vectorised; the a-dependent part 7 2^(a-n) + 8 (1 + ln q) 2^(-a) is convex in a, so the exact
# minimum over integers a in [2, n-1] is at floor/ceil of a* = log2(8 (1 + ln q) 2^n / 7) / 2, clamped)
def U1v(n, qs):
    qs = np.asarray(qs, dtype=np.float64); lq = 1 + np.log(qs)
    astar = 0.5 * np.log2(8 * lq * 2.0 ** n / 7)
    best = np.full(qs.shape, np.inf)
    for da in (-1, 0, 1, 2):
        a = np.clip(np.floor(astar) + da, 2, n - 1)
        best = np.minimum(best, 6 / qs + 7 * 2.0 ** (a - n) + 8 * lq * 2.0 ** (-a) + 8 * qs * lq * 2.0 ** (-n))
    return np.sqrt(best)
chk = max(abs(U1v(n, [q])[0] - U1(n, q)) for n, q in [(7, 3.0), (20, 1000.0), (40, 2.0 ** 17.3), (90, 2.0 ** 60.1), (150, 2.0 ** 140)])
print(f"U1 vectorised vs loop: max diff {chk:.2e}")
w1 = 0.0; w2 = 0.0; w3 = {}; w4 = 0.0
for n in range(7, 301):
    xs = np.arange(1.0, n - 1e-9, 1 / 64); qs = 2.0 ** xs; u = U1v(n, qs)
    m1 = xs <= (n - 3) / 3
    if m1.any(): w1 = max(w1, float((u[m1] * np.sqrt(qs[m1])).max()))
    m2 = (xs >= n / 2) & (xs <= n - 3)
    if m2.any(): w2 = max(w2, float((u[m2] / np.sqrt(qs[m2] * n / 2.0 ** n)).max()))
    w4 = max(w4, float((u * 2.0 ** (n / 4) / np.sqrt(n)).min()) if n >= 20 else 0.0)
    for eps in (0.05, 0.1, 0.2, 1 / 3, 0.5):
        m3 = (xs > eps * n) & (xs <= (1 - eps) * n)
        if m3.any():
            w3[eps] = max(w3.get(eps, 0.0), float((u[m3] * 2.0 ** (eps * n / 2) / np.sqrt(n)).max()))
print(f"U1: max over n in [7,300], 2 <= q <= 2^((n-3)/3) of U1*sqrt(q) = {w1:.4f}   [THM-4886 (A1): 6.5 claimed, 4.22 optimal]")
print(f"U1: max over n in [7,300], 2^(n/2) <= q <= 2^(n-3) of U1/sqrt(q n/2^n) = {w2:.4f}")
print(f"U1: best saving: max over n in [20,300] of (min_q U1) * 2^(n/4)/sqrt(n) = {w4:.4f}  (so min_q U1 <= {w4:.2f} sqrt(n) 2^(-n/4))")
print("U1 minor arcs: max over n in [7,300], 2^(eps n) < q <= 2^((1-eps) n) of U1 * 2^(eps n/2)/sqrt(n): " +
      ', '.join(f'eps={e:.3f}: {v:.4f}' for e, v in w3.items()) + "   [audit analytic: <= 4 for n >= 7]")
