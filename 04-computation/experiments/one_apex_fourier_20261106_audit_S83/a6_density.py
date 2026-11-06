# Audit S83, step 6: the actual density of one-apex values among odd numbers <= 2^k (all word lengths),
# versus the note's "The one-apex values alone have density about 0.5, and at most ln 2 (Session 16)".
# Every x in X_n is x' 1^t with x' strong (or 1^n), H(x) = H(x'); so X_N contains every value of every word of length <= N.
# We check min{H(x) : x strong of length m} >= 2^(m-2) + 1 for m <= 26, so values <= 2^k come from words of length <= k + 2.
import numpy as np, numba as nb
from oa_engine import g_block

@nb.njit(cache=True)
def min_strong(vals):
    # vals = H on X_n in order v = 2^(n-1) + i; strong = last bit 0 = even v
    m = 1 << 62
    for i in range(0, vals.size, 2):
        if vals[i] < m:
            m = vals[i]
    return m

@nb.njit(cache=True)
def mark(vals, lim):
    seen = np.zeros(lim + 1, dtype=np.uint8)
    for i in range(vals.size):
        h = vals[i]
        if h <= lim:
            seen[h] = 1
    return seen

ok = True
for m in range(2, 27):
    mn = min_strong(g_block(m))
    ok &= mn >= (1 << (m - 2)) + 1
print(f'min H over strong words of length m >= 2^(m-2) + 1 for m = 2..26: {ok}')
print(' k | odd numbers <= 2^k attained by some one-apex word, as a fraction of 2^(k-1) | D_k/|X_k| (words of length <= k only)')
for k in range(8, 25):
    vals = g_block(k + 2)
    seen = mark(vals, 1 << k)
    att = int(seen.sum())
    vals_k = g_block(k)
    Dk = int(mark(vals_k, 1 << k).sum())
    print(f'{k:2d} | {att / (1 << (k - 1)):.4f} | {Dk / (1 << (k - 1)):.4f}')
