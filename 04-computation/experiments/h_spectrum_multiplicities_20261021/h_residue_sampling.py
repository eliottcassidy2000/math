"""Pre-registered residue test for HYP-9355 (session opus-2026-10-21-S50; predictions committed in fad6a40581 before
this script was run).  Note: 05-knowledge/results/h_spectrum_multiplicities_20261021.md.

Uniform random labelled tournaments (each arc an independent fair coin; numpy default_rng with fixed seeds), exact
Hamiltonian-path counts by subset DP, at n = 10 (4,000,000 samples), n = 11 (2,000,000), n = 12 (500,000).
Reports, with binomial standard errors:
  P1  |p * freq(p | H) - 1| < 0.02 for p = 3, 5, 7, 11, 13 and n = 10, 11, 12;
  P2  residue 1 is the most frequent odd residue mod 8 at n = 10, 11, 12;
  P3  1 < 4 freq(H = 1 mod 8) < 1.0153 at n = 10;
  P4  (informative) the excess 4 freq(H = 1 mod 8) - 1 decreases from n = 10 to n = 12.
A prediction fails only if it is violated by more than 3 standard errors.  (A first run that drew 64-bit integer codes
crashed at n = 12, where C(12,2) = 66 bits overflow; its n = 10, 11 tables are kept in
h_residue_sampling_run1_crashed_at_n12.out.  This version draws one byte per arc for every n.)  The script's own checks are consistency
checks (every H odd by Redei; counts add up); prediction verdicts are printed, not asserted.
Prints ALL CHECKS PASSED.  Runtime about 2 minutes."""
import sys, os, time
import numpy as np
import numba as nb

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from h_multiplicities import ham_one                                     # noqa: E402

FAIL = []
T0 = time.time()


def check(c, m):
    print(('  ok   ' if c else '  FAIL ') + m, flush=True)
    if not c:
        FAIL.append(m)


@nb.njit(parallel=True, cache=True)
def sample_counts(n, bits):
    """bits: (m, C(n,2)) uint8 arc coins (pairs i<j in lexicographic order; 1 means i -> j).  Returns H for each."""
    m = bits.shape[0]
    out = np.zeros(m, dtype=np.int64)
    for t in nb.prange(m):
        adj = np.zeros(n, dtype=np.int64)
        e = 0
        for i in range(n):
            for j in range(i + 1, n):
                if bits[t, e]:
                    adj[i] |= 1 << j
                else:
                    adj[j] |= 1 << i
                e += 1
        out[t] = ham_one(n, adj)
    return out


SIZES = {10: 4_000_000, 11: 2_000_000, 12: 500_000}
SEEDS = {10: 20261021, 11: 20261022, 12: 20261023}
res = {}
for n, N in SIZES.items():
    rng = np.random.default_rng(SEEDS[n])
    E = n * (n - 1) // 2
    cnt_p = {p: 0 for p in (3, 5, 7, 11, 13)}
    c8 = np.zeros(8, dtype=np.int64); c16 = np.zeros(16, dtype=np.int64)
    odd = True; done = 0
    while done < N:
        k = min(250_000, N - done)
        bits = rng.integers(0, 2, size=(k, E), dtype=np.uint8)
        H = sample_counts(n, bits)
        odd &= bool(np.all(H % 2 == 1))
        for p in cnt_p:
            cnt_p[p] += int(np.sum(H % p == 0))
        c8 += np.bincount(H % 8, minlength=8); c16 += np.bincount(H % 16, minlength=16)
        done += k
    res[n] = (N, cnt_p, c8, c16)
    check(odd and int(c8.sum()) == N, f'n = {n}: {N} samples, every H odd (Redei) [{time.time() - T0:.0f}s]')
    print(f'     n = {n}: p | H: ' +
          ', '.join(f'{p}: {cnt_p[p] / N:.5f} (rel {p * cnt_p[p] / N:.4f} +- {p * np.sqrt((1 / p) * (1 - 1 / p) / N):.4f})' for p in cnt_p))
    print(f'     n = {n}: mod 8 relative masses (1,3,5,7): ' +
          ', '.join(f'{r}:{4 * c8[r] / N:.4f}' for r in (1, 3, 5, 7)) + f'  (s.e. {4 * np.sqrt(0.25 * 0.75 / N):.4f})')
    print(f'     n = {n}: mod 16 relative masses: ' + ', '.join(f'{r}:{8 * c16[r] / N:.3f}' for r in range(1, 16, 2)))

print('verdicts on the predictions committed in fad6a40581')
p1 = True
for n, (N, cnt_p, c8, c16) in res.items():
    for p, c in cnt_p.items():
        rel = p * c / N; se = p * np.sqrt((1 / p) * (1 - 1 / p) / N)
        if abs(rel - 1) >= 0.02 + 3 * se:
            p1 = False
print('  P1', 'PASSED' if p1 else 'FAILED', '(|p freq - 1| < 0.02 within 3 s.e.)')
p2 = True
for n, (N, cnt_p, c8, c16) in res.items():
    se = 4 * np.sqrt(0.25 * 0.75 / N)
    m1 = 4 * c8[1] / N
    for r in (3, 5, 7):
        if 4 * c8[r] / N > m1 + 3 * se:
            p2 = False
    print(f'     n = {n}: residue 1 mass {m1:.4f}; best other {max(4 * c8[r] / N for r in (3, 5, 7)):.4f}; strictly largest: {m1 > max(4 * c8[r] / N for r in (3, 5, 7))}')
print('  P2', 'PASSED' if p2 else 'FAILED', '(residue 1 mod 8 not beaten by more than 3 s.e.)')
N, cnt_p, c8, c16 = res[10]
m1 = 4 * c8[1] / N; se = 4 * np.sqrt(0.25 * 0.75 / N)
p3 = (m1 > 1 - 3 * se) and (m1 < 1.0153 + 3 * se)
print(f'  P3', 'PASSED' if p3 else 'FAILED', f'(n = 10 mass of 1 mod 8 = {m1:.4f} +- {se:.4f}, predicted in (1, 1.0153))',
      '; strictly inside without the tolerance:' , (m1 > 1) and (m1 < 1.0153))
ex = {n: 4 * res[n][2][1] / res[n][0] - 1 for n in res}
print('  P4 (informative): excess of residue 1 mod 8:', {n: round(float(v), 4) for n, v in ex.items()},
      '-> decreasing' if ex[10] > ex[11] > ex[12] else '-> not monotone')

print()
print('ALL CHECKS PASSED' if not FAIL else f'FAILED: {FAIL}', f' ({time.time() - T0:.0f}s)')
