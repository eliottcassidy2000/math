"""leafcount_dp.py -- independent count of the normal-form-B search domain (numpy dynamic programming).

For d, k, a (b = k - a, delta = a - b), n = d - 1, the search domain is the set of pairs (A, B) with
  A in the representative file (a-subsets of Q_n), B any b-subset of Q_n, and
  |beta_c(A) + beta_c(B)| >= delta for all c = 0..n-1,   beta_c(X) = #{x in X: x_c = 0} - #{x: x_c = 1}.
This script counts it WITHOUT enumerating B:
  dp_b[beta] = number of b-subsets of Q_n with imbalance vector beta (built vertex by vertex),
  N(alpha)   = sum_beta dp_b[beta] * prod_c [ |alpha_c + beta_c| >= delta ],
  leaves     = sum over representatives A of N(beta(A)).
N is invariant under signed permutations of alpha (dp_b is invariant under Aut(Q_n)); this is used to
group representatives by the sorted vector |beta(A)|, and is itself checked numerically on random alpha.
No code is shared with the C search.

usage: python3 leafcount_dp.py d k [a1 a2 ...]   (default: all a = ceil(k/2)..k with an existing rep file)
       prints one line per (k, a):  DP d= k= a= b= reps= leaves=
"""
import sys, os, itertools, random
from collections import Counter
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, '..', 'data')

def repfile(n, a):
    return os.path.join(DATA, 'q%d' % n, 'reps_n%d_a%d.bin' % (n, a))

_DP = {}
def dp_subsets(n, b):
    """list dp[j] (j=0..b), dp[j] has shape (2j+1,)*n, index i <-> beta = i - j"""
    key = (n, b)
    if key in _DP: return _DP[key]
    dp = [np.zeros((2 * j + 1,) * n, dtype=np.int64) for j in range(b + 1)]
    dp[0][(0,) * n] = 1
    for w in range(1 << n):
        sig = [1 if not (w >> c) & 1 else -1 for c in range(n)]
        for j in range(b, 0, -1):
            # beta' = beta + sig ; index in dp[j] = (i_old - (j-1)) + sig + j = i_old + 1 + sig
            sl = tuple(slice(1 + s, 1 + s + 2 * j - 1) for s in sig)
            dp[j][sl] += dp[j - 1]
    # sanity: total number of b-subsets
    from math import comb
    for j in range(b + 1):
        assert int(dp[j].sum()) == comb(1 << n, j), (j, int(dp[j].sum()))
    _DP[key] = dp
    return dp

def N_of_alpha(dpb, b, n, delta, alpha):
    vals = np.arange(-b, b + 1)
    T = dpb
    for c in range(n):
        ind = (np.abs(alpha[c] + vals) >= delta).astype(np.int64)
        T = np.tensordot(ind, T, axes=([0], [0]))   # contract the leading axis
    return int(T)

_PAT = {}
def beta_patterns(n, a, chunk=1 << 20):
    """(number of reps, Counter of sorted |beta| vectors, up to 1000 sample beta rows); the file is read in
    chunks so that the large files (e.g. 16.9M reps for n = 6, a = 11) stay well below the RAM limit"""
    if (n, a) in _PAT: return _PAT[(n, a)]
    R = np.fromfile(repfile(n, a), dtype='<u8')
    pat = Counter(); samples = []
    for start in range(0, len(R), chunk):
        Rc = R[start:start + chunk]
        pc = np.zeros(len(Rc), dtype=np.int64)
        B = np.zeros((len(Rc), n), dtype=np.int64)
        for w in range(1 << n):
            bit = ((Rc >> np.uint64(w)) & np.uint64(1)).astype(np.int64)
            pc += bit
            for c in range(n):
                B[:, c] += bit * (1 if not (w >> c) & 1 else -1)
        assert (pc == a).all(), 'representative of wrong size'
        u, cnt = np.unique(np.sort(np.abs(B), axis=1), axis=0, return_counts=True)
        for row, cn in zip(u.tolist(), cnt.tolist()):
            pat[tuple(row)] += cn
        if start == 0:
            samples = B[:1000].tolist()
    _PAT[(n, a)] = (len(R), pat, samples)
    return _PAT[(n, a)]

def leafcount(d, k, a, check_symmetry=True):
    n = d - 1; b = k - a; delta = a - b
    assert b >= 0 and delta >= 0
    dp = dp_subsets(n, b); dpb = dp[b]
    nreps, pat, samples = beta_patterns(n, a)
    assert sum(pat.values()) == nreps
    tot = 0
    for p, cnt in pat.items():
        tot += cnt * N_of_alpha(dpb, b, n, delta, p)
    if check_symmetry and samples:
        rng = random.Random(1000 * k + a)
        for _ in range(3):
            row = samples[rng.randrange(len(samples))]
            perm = list(range(n)); rng.shuffle(perm)
            signed = [row[perm[c]] * rng.choice((-1, 1)) for c in range(n)]
            assert N_of_alpha(dpb, b, n, delta, row) == N_of_alpha(dpb, b, n, delta, signed) == \
                   N_of_alpha(dpb, b, n, delta, tuple(sorted(abs(x) for x in row)))
    return nreps, tot

if __name__ == '__main__':
    d, k = int(sys.argv[1]), int(sys.argv[2])
    n = d - 1
    avals = [int(x) for x in sys.argv[3:]] or [a for a in range((k + 1) // 2, k + 1) if os.path.exists(repfile(n, a))]
    total = 0
    for a in avals:
        nreps, lv = leafcount(d, k, a)
        total += lv
        print('DP d=%d k=%d a=%d b=%d reps=%d leaves=%d' % (d, k, a, k - a, nreps, lv), flush=True)
    print('DP-TOTAL d=%d k=%d a_values=%s leaves=%d' % (d, k, ','.join(map(str, avals)), total), flush=True)
