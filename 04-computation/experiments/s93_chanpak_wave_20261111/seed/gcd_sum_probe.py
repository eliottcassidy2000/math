import numpy as np
def rev(t, tau):
    r = np.zeros_like(t)
    for i in range(tau):
        r |= ((t >> i) & 1) << (tau - 1 - i)
    return r
for tau in range(4, 13):
    T = np.arange(1 << tau, dtype=np.int64)
    R = rev(T, tau)
    A = T[:, None] - T[None, :]
    B = R[:, None] - R[None, :]
    mask = ((T[:, None] ^ T[None, :]) & 1) == 1          # last bits differ
    al = np.abs(B[mask]); be = np.abs(A[mask])
    g = np.gcd(al, be)
    s1 = (g / al).sum(); s2 = (g / be).sum(); s3 = (g / np.maximum(al, be)).sum()
    print(f"tau = {tau:2d}: pairs {mask.sum():9d}; sum g/|alpha| = {s1 / 2**tau:8.3f} * 2^tau; sum g/|beta| = {s2 / 2**tau:8.3f} * 2^tau; "
          f"sum g/max = {s3 / 2**tau:7.3f} * 2^tau; max g = {g.max()}; pairs with g >= 2^(tau/2): {(g >= 2 ** (tau / 2)).sum()}")
