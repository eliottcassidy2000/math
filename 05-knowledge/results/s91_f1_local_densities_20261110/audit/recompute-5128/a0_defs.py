# Auditor recompute-5128: definitions only (written from the theorem text, not from the session scripts).
#   psi(0) = 0, psi(2h) = psi(h) + h, psi(2h+1) = psi(h) + floor(h/2) + 1
#   x(h) = 2 psi(h) + 1, y(h) = 2h + 1, q(h) = h - psi(h)
import numpy as np
from numba import njit


@njit(cache=True)
def psi_table(n):
    """psi(h) for 0 <= h < n, int64."""
    psi = np.zeros(n, dtype=np.int64)
    for h in range(1, n):
        k = h >> 1
        if h & 1:
            psi[h] = psi[k] + (k >> 1) + 1
        else:
            psi[h] = psi[k] + k
    return psi


@njit(cache=True)
def psi_table32(n):
    """psi(h) for 0 <= h < n, int32 (psi(h) <= h)."""
    psi = np.zeros(n, dtype=np.int32)
    for h in range(1, n):
        k = h >> 1
        if h & 1:
            psi[h] = psi[k] + (k >> 1) + 1
        else:
            psi[h] = psi[k] + k
    return psi


def psi_py(h):
    """Pure-python recursion (independent of the table), exact integers."""
    if h == 0:
        return 0
    k = h >> 1
    if h & 1:
        return psi_py(k) + (k >> 1) + 1
    return psi_py(k) + k


def q_py(h):
    """q(h) = sum over one-bits i of ceil(H_i/2), H_i = floor(h / 2^(i+1))  (closed form of (B3), used as a cross-check)."""
    tot = 0
    i = 0
    while (h >> i) > 0:
        if (h >> i) & 1:
            H = h >> (i + 1)
            tot += (H + 1) // 2
        i += 1
    return tot


if __name__ == "__main__":
    n = 1 << 16
    t = psi_table(n)
    t32 = psi_table32(n)
    for h in range(n):
        p = psi_py(h)
        assert t[h] == p and t32[h] == p
        assert h - p == q_py(h), h
        # recursion of q stated in the theorem
        if h >= 1:
            k = h >> 1
            if h & 1:
                assert (h - p) == (k - psi_py(k)) + (k + 1) // 2
            else:
                assert (h - p) == (k - psi_py(k))
        # 0 <= q <= h/2, y/2 < x <= y
        qq = h - p
        assert 0 <= 2 * qq <= h
        x = 2 * p + 1
        y = 2 * h + 1
        assert y < 2 * x <= 2 * y
    print("defs OK for h <", n)
    print("first states (h, psi, x, y):", [(h, int(t[h]), int(2 * t[h] + 1), 2 * h + 1) for h in range(12)])
