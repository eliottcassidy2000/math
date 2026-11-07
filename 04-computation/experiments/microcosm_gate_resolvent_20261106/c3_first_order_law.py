# mac-mini-2026-11-06-microcosm. HYP-9490 evidence: over all isomorphism classes of orders 4..7, the relative deviation from the reference tower is kappa(T) eps_N with kappa = c3 excess.
import subprocess, itertools, math
from fractions import Fraction
from collections import Counter
import os; exec(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cover_tower_values.py')).read().split('cores = {')[0])   # reuse forward_counts, P_value, tower_H
def parse(line, n):
    A = [[False] * n for _ in range(n)]; k = 0
    for i in range(n):
        for j in range(i + 1, n):
            if line[k] == '1': A[i][j] = True
            else: A[j][i] = True
            k += 1
    return A
def c3(n, A):
    s = [sum(A[i]) for i in range(n)]
    return math.comb(n, 3) - sum(math.comb(x, 2) for x in s)
for n, N in ((4, 402), (5, 401), (6, 402), (7, 401)):
    out = subprocess.run(['/opt/homebrew/bin/gentourng', '-q', str(n)], capture_output=True, text=True).stdout.split()
    if n % 2:
        ref_cnt, ref_n = forward_counts(1, [[False]]), 1
    else:
        ref_cnt, ref_n = forward_counts(2, [[False, True], [False, False]]), 2
    ref = tower_H(ref_cnt, ref_n, (N - ref_n) // 2)
    # c3 of the reference tower at order n
    rr = (n - ref_n) // 2
    c3_ref = 0 + sum(ref_n + 2 * i for i in range(rr))
    rows = []
    for line in out:
        A = parse(line, n)
        dev = float(Fraction(tower_H(forward_counts(n, A), n, (N - n) // 2), ref) - 1)
        rows.append((c3(n, A) - c3_ref, dev))
    # universal factor estimate from the classes with kappa = +-1
    unit = [d / k for k, d in rows if k != 0]
    eps = sum(unit) / len(unit)
    spread = max(abs(d - k * eps) for k, d in rows)
    print('order n = %d (%d classes), N = %d: kappa = c3 - %d ; dev/kappa ranges [%.4e, %.4e] ; max |dev - kappa*eps| = %.2e (eps = %.4e)'
          % (n, len(rows), N, c3_ref, min(unit), max(unit), spread, eps), flush=True)
