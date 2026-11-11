# Auditor recompute-5128, extra: rate max_k rho(T_k)/2 for a large modulus by power iteration (no dense eigenvalues),
# validated on m = 79 and m = 127 against the dense eigenvalue computation of a1_rates.py.
# Usage: python3 a9b_rate_power.py m [iters]
import sys
import math
import numpy as np
from numba import njit


@njit(cache=True)
def power_rates(m, iters, burn):
    inv2 = (m + 1) // 2
    n = 2 * m
    j0 = np.zeros(n, dtype=np.int64)
    j1 = np.zeros(n, dtype=np.int64)
    c = np.zeros(n, dtype=np.int64)
    for eta in range(m):
        for beta in range(2):
            i = 2 * eta + beta
            j0[i] = 2 * ((2 * eta) % m)
            j1[i] = 2 * ((2 * eta + 1) % m) + 1
            c[i] = ((eta + beta) * inv2) % m
    zeta = np.zeros(m, dtype=np.complex128)
    for a in range(m):
        zeta[a] = np.exp(2j * np.pi * a / m)
    out = np.zeros((m - 1) // 2 + 1, dtype=np.float64)
    v = np.zeros(n, dtype=np.complex128)
    w = np.zeros(n, dtype=np.complex128)
    ph = np.zeros(n, dtype=np.complex128)
    for k in range(1, (m - 1) // 2 + 1):
        for i in range(n):
            ph[i] = zeta[(k * c[i]) % m]
            # deterministic pseudo-random start
            t = (i * 2654435761 + k * 40503) % 1000003
            v[i] = complex(math.cos(t * 0.001), math.sin(t * 0.0017 + 1.0))
        lg = 0.0
        for it in range(iters):
            s = 0.0
            for i in range(n):
                z = v[j0[i]] + ph[i] * v[j1[i]]
                w[i] = z
                s += z.real * z.real + z.imag * z.imag
            nv = math.sqrt(s)
            for i in range(n):
                v[i] = w[i] / nv
            if it >= burn:
                lg += math.log(nv)
        out[k] = math.exp(lg / (iters - burn)) / 2.0
    return out


def main():
    m = int(sys.argv[1])
    iters = int(sys.argv[2]) if len(sys.argv) > 2 else 6000
    out = power_rates(m, iters, iters // 3)
    k = int(np.argmax(out))
    print(f"m = {m}: power iteration ({iters} steps, first third discarded): max_k rho(T_k)/2 = {out[k]:.5f} at k = {k}; min over k = {out[1:].min():.5f}; five largest: {np.sort(out)[-5:][::-1].round(5).tolist()}")


if __name__ == "__main__":
    main()
