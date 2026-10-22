"""Audit S50 - a7b (second replication, n = 11 with 4e6 and n = 12 with 1e6 samples, new seeds) - copy of a7: independent replication of the HYP-9355 sampling test (own generator, own seeds, own DP).
Uniform labelled tournaments drawn as random 64-bit words split into arc bits (numpy PCG64 via Generator.integers on
uint64 words, masked per arc), n = 10 (4e6 samples) and n = 11 (2e6).  Reports p * freq(p | H) and mod-8 masses with
binomial standard errors, and the s.e. of a difference of two mod-8 masses computed with the multinomial covariance."""
import time
import numpy as np
import numba as nb
import alib as L

T0 = time.time()


@nb.njit(parallel=True, cache=True)
def sample(n, words, out):
    m = words.shape[0]
    for t in nb.prange(m):
        adj = np.zeros(n, dtype=np.int64)
        e = 0
        for i in range(n):
            for j in range(i + 1, n):
                w = words[t, e >> 6]
                if (w >> np.uint64(e & 63)) & np.uint64(1):
                    adj[i] |= 1 << j
                else:
                    adj[j] |= 1 << i
                e += 1
        out[t] = L.ham_dp(n, adj)


for n, N, seed in ((11, 4_000_000, 777111), (12, 1_000_000, 777012)):
    rng = np.random.default_rng(seed)
    E = n * (n - 1) // 2
    nw = (E + 63) // 64
    cp = {p: 0 for p in (3, 5, 7, 11, 13)}
    c8 = np.zeros(8, dtype=np.int64)
    done = 0
    while done < N:
        k = min(250_000, N - done)
        words = rng.integers(0, np.iinfo(np.uint64).max, size=(k, nw), dtype=np.uint64, endpoint=True)
        out = np.zeros(k, dtype=np.int64)
        sample(n, words, out)
        assert np.all(out % 2 == 1)
        for p in cp:
            cp[p] += int(np.sum(out % p == 0))
        c8 += np.bincount(out % 8, minlength=8)
        done += k
    print(f'n = {n} ({N} samples, seed {seed}) [{time.time() - T0:.0f}s]')
    print('   p * freq(p | H): ' + ', '.join(f'{p}: {p * c / N:.4f} +- {p * np.sqrt((1 / p) * (1 - 1 / p) / N):.4f}' for p, c in cp.items()))
    se = 4 * np.sqrt(0.25 * 0.75 / N)
    sed = 4 * np.sqrt(0.5 / N)
    print('   mod-8 masses (1,3,5,7): ' + ', '.join(f'{r}: {4 * c8[r] / N:.4f}' for r in (1, 3, 5, 7)) +
          f'  (s.e. of one mass {se:.4f}; of a difference of two masses {sed:.4f})')
print(f'done ({time.time() - T0:.0f}s)')
