"""mean and variance of tau per dyadic octave (drift per octave, for comparison with the parity mode)."""
import numpy as np
from numba import njit
@njit(cache=True)
def tau_array16(N):
    tau = np.zeros(N + 8, dtype=np.int16)
    for n in range(2, N + 8):
        m = n; steps = 0
        while m >= n:
            if m & 1:
                m = (3 * m + 1) >> 1
            else:
                m >>= 1
            steps += 1
        tau[n] = steps + tau[m]
    return tau
tau = tau_array16(1 << 29)
print('octave  mean tau  slope per octave  var tau')
prev = None
for j in range(18, 30):
    lo, hi = 1 << (j-1), 1 << j
    m = float(tau[lo:hi].mean()); v = float(tau[lo:hi].astype(np.float64).var())
    print(f' 2^{j:2d}  {m:8.3f}   {(m-prev) if prev is not None else float("nan"):+7.3f}        {v:8.2f}')
    prev = m
