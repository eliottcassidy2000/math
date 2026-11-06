# Audit S83, step 3: the constants of (A1), (A2), (A3), evaluated from the exact Theorem A bound minimised over ALL splits,
# for EVERY n (the session scanned only n = 9, 12, ..., 90 and q = 2^(k/4)), on a fine q grid that includes the range
# endpoints, and all integer q when n <= 30.
# For fixed (n, L) the factor (3 + 2^(3-a) q)(3 + 2^(2-l) q), l = n - L - a, is 9 + 3(alpha + beta) + alpha*beta with alpha*beta
# independent of a and alpha + beta convex in a, minimised at a* = (n - L + 1)/2; so the exact optimum over a is attained at
# floor/ceil of a* clamped to [2, n - L - 1]. We evaluate a* - 2 .. a* + 2 (clamped) and all L: this is the exact min over splits.
import math, sys
import numpy as np

def best_bound_vec(n, qs):
    qs = np.asarray(qs, dtype=np.float64)[:, None]
    best = np.full(qs.shape[0], np.inf)
    for L in range(1, n - 2):
        N = n - L
        amid = (N + 1) // 2
        for a in range(amid - 2, amid + 3):
            a = min(max(a, 2), N - 1)
            l = N - a
            if l < 1:
                continue
            R = 2 * (2.0 ** L / qs[:, 0] + 1) * np.sqrt((3 + 2.0 ** (3 - a) * qs[:, 0]) * (3 + 2.0 ** (2 - l) * qs[:, 0]))
            b = 2.0 ** (-L / 2) * np.sqrt(1 + R)
            best = np.minimum(best, b)
    return best

def best_bound_full(n, q):  # brute force over every split, for cross-checks
    best = math.inf
    for L in range(1, n - 2):
        for a in range(2, n - L):
            l = n - a - L
            R = 2 * (2.0 ** L / q + 1) * math.sqrt((3 + 2.0 ** (3 - a) * q) * (3 + 2.0 ** (2 - l) * q))
            best = min(best, 2.0 ** (-L / 2) * math.sqrt(1 + R))
    return best

def qgrid(n, lo, hi, step):
    if hi < lo:
        return np.array([])
    xs = np.arange(math.log2(lo), math.log2(hi), step)
    qs = np.concatenate([2.0 ** xs, [lo, hi]])
    if n <= 30:   # also every integer q <= 4096 in range (the bound is continuous in q; the log grid covers the rest)
        qs = np.concatenate([qs, np.arange(max(2, math.ceil(lo)), min(4096, math.floor(hi)) + 1, dtype=np.float64)])
    qs = np.unique(qs)
    return qs[(qs >= lo * (1 - 1e-12)) & (qs <= hi * (1 + 1e-12))]

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 200
# cross-check of the reduction against brute force
rng = np.random.default_rng(1)
mx = 0.0
for _ in range(300):
    n = int(rng.integers(6, 60)); q = float(2.0 ** rng.uniform(0.5, n - 1))
    mx = max(mx, abs(best_bound_vec(n, [q])[0] - best_bound_full(n, q)))
print(f'reduction check (300 random (n, q)): max |fast - brute force| = {mx:.2e}')

w1 = (0, None); w2 = (0, None); per_n_A2 = {}; per_n_A1 = {}
for n in range(6, NMAX + 1):
    qs = qgrid(n, 2.0, 2.0 ** ((n - 3) / 3), 1 / 128)
    if qs.size:
        r = best_bound_vec(n, qs) * np.sqrt(qs)
        i = int(np.argmax(r)); per_n_A1[n] = float(r[i])
        if r[i] > w1[0]:
            w1 = (float(r[i]), (n, float(qs[i])))
    qs = qgrid(n, max(2.0, 2.0 ** ((n - 6) / 3)), 2.0 ** (n - 6), 1 / 128)
    if qs.size:
        r = best_bound_vec(n, qs) / (qs / 2.0 ** n) ** 0.25
        i = int(np.argmax(r)); per_n_A2[n] = float(r[i])
        if r[i] > w2[0]:
            w2 = (float(r[i]), (n, float(qs[i]), float(math.log2(qs[i]))))
print(f'(A1) sup over n in [6,{NMAX}], 2 <= q <= 2^((n-3)/3) of bound*sqrt(q) = {w1[0]:.4f} at (n, q) = {w1[1]}   [claimed 6.5; session scan 4.213]')
print(f'(A2) sup over n in [6,{NMAX}], 2^((n-6)/3) <= q <= 2^(n-6) of bound/(q/2^n)^(1/4) = {w2[0]:.4f} at (n, q, log2 q) = {w2[1]}   [claimed 12; session scan 9.539]')
for r in (0, 1, 2):
    v1 = [v for k, v in per_n_A1.items() if k % 3 == r]
    v2 = [v for k, v in per_n_A2.items() if k % 3 == r]
    print(f'      n = {r} mod 3: max (A1) ratio {max(v1):.4f}; max (A2) ratio {max(v2):.4f}')
print('      (A2) per-n maxima, n = 6..40:', ' '.join(f'{k}:{per_n_A2[k]:.2f}' for k in range(6, 41) if k in per_n_A2))
# brute-force confirmation at the worst (A2) point
n0, q0 = w2[1][0], w2[1][1]
print(f'      brute-force bound at the worst (A2) point: {best_bound_full(n0, q0) / (q0 / 2.0 ** n0) ** 0.25:.4f}')

# (A3): q in (2^(eps n), 2^((1-eps) n)], claimed |S_n| <= 12 * 2^(n-1) 2^(-eps n/4); use min(bound, 1) (trivial bound allowed)
for eps in (0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 1 / 3):
    worst = (0, None); worst_nontriv = (0, None)
    for n in range(6, NMAX + 1):
        lo, hi = max(2.0, 2.0 ** (eps * n)), 2.0 ** ((1 - eps) * n)
        qs = qgrid(n, lo, hi, 1 / 64)
        qs = qs[qs > 2.0 ** (eps * n)]
        if not qs.size:
            continue
        b = best_bound_vec(n, qs)
        r = np.minimum(b, 1.0) * 2.0 ** (eps * n / 4)
        i = int(np.argmax(r))
        if r[i] > worst[0]:
            worst = (float(r[i]), (n, round(float(math.log2(qs[i])), 3), round(float(b[i]), 4)))
        if eps * n >= 6:
            r2 = b * 2.0 ** (eps * n / 4); j = int(np.argmax(r2))
            if r2[j] > worst_nontriv[0]:
                worst_nontriv = (float(r2[j]), (n, round(float(math.log2(qs[j])), 3)))
    print(f'(A3) eps = {eps:.3f}: sup_n<= {NMAX} of min(bound,1)*2^(eps n/4) = {worst[0]:.4f} at (n, log2 q, bound) = {worst[1]}; '
          f'over eps*n >= 6 of bound*2^(eps n/4) = {worst_nontriv[0]:.4f} at {worst_nontriv[1]}   [claimed 12; session 4.564]')
