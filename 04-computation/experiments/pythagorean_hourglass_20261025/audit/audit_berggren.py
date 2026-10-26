# Audit spot check of the results note section 2 (not in the core audit list): Berggren depth statistics by an
# independent method (the (m, n) generator form of Berggren's tree: children (2m-n, m), (2m+n, m), (m+2n, n); the
# parent of (M, N) is (N, 2N-M), (N, M-2N) or (M-2N, N) according as M < 2N, 2N < M < 3N, M > 3N; root (2, 1)),
# and Collatz total stopping times with the T-map (n -> n/2, (3n+1)/2).  Audit 2026-10-25.
import math, statistics
def depth(M, N):
    d = 0
    while (M, N) != (2, 1):
        if M < 2 * N: M, N = N, 2 * N - M
        elif M < 3 * N: M, N = N, M - 2 * N
        else: M, N = M - 2 * N, N
        d += 1
    return d
def T_stop(n):
    s = 0
    while n != 1:
        n = n // 2 if n % 2 == 0 else (3 * n + 1) // 2; s += 1
    return s
out = []
cmax = 200000
D = {}
for m in range(2, int(math.isqrt(cmax)) + 2):
    for n in range(1, m):
        if (m - n) % 2 == 1 and math.gcd(m, n) == 1 and m * m + n * n <= cmax:
            D[m * m + n * n, m, n] = depth(m, n)
out.append(f'{len(D)} primitive triples with c <= {cmax}')
for lo in (1000, 10000, 100000):
    v = [d for (c, m, n), d in D.items() if lo <= c < 2 * lo]
    L = math.log(1.5 * lo)
    mu = statistics.mean(v)
    out.append(f'c in [{lo}, {2 * lo}): {len(v)} triples, mean depth {mu:.2f} = {mu / L:.3f} log(1.5 lo), '
               f'(6/pi^2)(log sqrt c)^2 at c = 1.5 lo: {6 / math.pi ** 2 * (0.5 * L) ** 2:.2f}, ratio {mu / (6 / math.pi ** 2 * (0.5 * L) ** 2):.3f}')
for lo in (1000, 10000, 100000):
    v = [T_stop(n) for n in range(lo, 2 * lo)]
    out.append(f'Collatz n in [{lo}, {2 * lo}): mean T-stopping time {statistics.mean(v):.2f}')
print('\n'.join(out))
with open('audit_berggren.out', 'w') as fh: fh.write('\n'.join(out) + '\n')

# scaling beyond the note's range: random primitive triples with c in [lo, 2 lo), uniform over lattice points (m, n)
import random
rng = random.Random(2026)
out2 = []
for lo in (10 ** 6, 10 ** 8, 10 ** 10):
    v = []
    while len(v) < 40000:
        R = math.isqrt(2 * lo) + 1
        m = rng.randrange(1, R); n = rng.randrange(1, R)
        if n >= m: continue
        c = m * m + n * n
        if not (lo <= c < 2 * lo) or (m - n) % 2 == 0 or math.gcd(m, n) != 1: continue
        v.append(depth(m, n))
    L = math.log(1.5 * lo)
    mu = statistics.mean(v)
    out2.append(f'sampled c in [{lo:.0e}, {2 * lo:.0e}): 40000 triples, mean depth {mu:.2f}; (6/pi^2)(log sqrt c)^2 = '
                f'{6 / math.pi ** 2 * (0.5 * L) ** 2:.2f}; ratio {mu / (6 / math.pi ** 2 * (0.5 * L) ** 2):.3f}; mean/log c = {mu / L:.3f}; '
                f'mean/(log c)^2 = {mu / L ** 2:.4f}')
print('\n'.join(out2))
with open('audit_berggren.out', 'a') as fh: fh.write('\n'.join(out2) + '\n')
