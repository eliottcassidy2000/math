"""
collatz_tau_runs_20261013.py -- runs of consecutive integers with equal total stopping time tau (Garner's clusters),
their mean length as a function of N, and the density of 'breaks' (n with tau(n) != tau(n+1)).
If pairs fail to coalesce with probability ~ c (log N)^(-1/2), breaks have density ~ c (log N)^(-1/2) and runs have
mean length ~ (log N)^(1/2)/c.  Also: the parity-flip density (t(n) != t(n+1)) and its ratio to the break density.
Run: python3 collatz_tau_runs_20261013.py <N_log2>
"""
import sys, time
import numpy as np
from numba import njit

NL = int(sys.argv[1]) if len(sys.argv) > 1 else 28

@njit(cache=True)
def tau_array(N):
    tau = np.zeros(N + 2, dtype=np.int32)
    for n in range(2, N + 2):
        m = n; steps = 0
        while m >= n:
            if m & 1:
                m = (3 * m + 1) >> 1
            else:
                m >>= 1
            steps += 1
        tau[n] = steps + tau[m]
    return tau

N = 1 << NL
t0 = time.time()
tau = tau_array(N)
print(f'tau to N=2^{NL} in {time.time()-t0:.1f}s')
print(' scale     break density   parity-flip density   flip/break   mean run length   longest run (start)')
for L in range(16, NL + 1, 2):
    M = 1 << L
    a = tau[1:M + 1]; b = tau[2:M + 2]
    breaks = (a != b)
    flips = ((a - b) % 2 == 1)
    bd = float(breaks.mean()); fd = float(flips.mean())
    # run lengths of equal tau
    idx = np.nonzero(breaks)[0]
    runs = np.diff(np.concatenate(([0], idx + 1)))
    longest = int(runs.max()); start = int(np.argmax(runs))
    print(f' 2^{L:2d}   {bd:.6f}        {fd:.6f}            {fd/bd:.4f}        {1/bd:8.4f}          {longest} (near index {start})')
# distribution of tau(n+1) - tau(n) at the breaks (top values)
d = (tau[2:N + 2] - tau[1:N + 1])
vals, cnt = np.unique(d[d != 0], return_counts=True)
order = np.argsort(-cnt)[:12]
print('most common nonzero differences tau(n+1)-tau(n):', [(int(vals[i]), int(cnt[i])) for i in order])
