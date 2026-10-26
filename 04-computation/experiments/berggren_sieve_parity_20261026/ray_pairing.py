# Exact ray pairing for the branch character.  For a node w = (m, n) with sign s = (-1)^#A(w):
#   D-rays of its U- and A-children:  wU D^j = (2m-n+2jm, m),  wA D^j = (2m+n+2jm, m)   (same second coordinate),
#   U-rays of its D- and A-children:  wD U^j = (m+2n+j(m+n), n+j(m+n)),  wA U^j = (2m+n+j(m+n), m+j(m+n))
#   (same line x - y = m + n); heights interleave, signs are opposite.  Every ray except the root's two is paired,
#   and every node is the root, an A-child, or lies on exactly one ray.  Hence the identity
#     S(X) = 1 + N_U(X) + N_D(X) + E(X),   E(X) = sum_w s_w (delta_D(w,X) + delta_U(w,X) - [c(wA) <= X]),
#   with N_U, N_D the root rays and delta in {0, 1}.  This script checks the interleaving and the identity exactly.
import sys
from math import gcd, isqrt
X = int(sys.argv[1]) if len(sys.argv) > 1 else 10**6
def c(v): return v[0] ** 2 + v[1] ** 2
def cnt_ray(start, step, X):
    """number of j >= 1 with c(start + j*step) <= X (c increasing in j on these rays)"""
    j = 0
    # galloping then bisection
    lo, hi = 0, 1
    while c((start[0] + hi * step[0], start[1] + hi * step[1])) <= X: lo, hi = hi, hi * 2
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if c((start[0] + mid * step[0], start[1] + mid * step[1])) <= X: lo = mid
        else: hi = mid
    return lo
# enumerate all nodes with their signs via BFS from the root (heights <= X)
S = 0; E = 0; interleave_bad = 0
stack = [((2, 1), 1)]
while stack:
    (m, n), s = stack.pop()
    if m * m + n * n > X: continue
    S += s
    wU, wA, wD = (2 * m - n, m), (2 * m + n, m), (m + 2 * n, n)
    dD = cnt_ray(wU, (2 * m, 0), X) - cnt_ray(wA, (2 * m, 0), X)
    dU = cnt_ray(wD, (m + n, m + n), X) - cnt_ray(wA, (m + n, m + n), X)
    if dD not in (0, 1) or dU not in (0, 1): interleave_bad += 1
    E += s * (dD + dU - (1 if c(wA) <= X else 0))
    stack += [(wU, s), (wA, -s), (wD, s)]
NU = sum(1 for k in range(1, isqrt(X) + 2) if (k + 2) ** 2 + (k + 1) ** 2 <= X)
ND = sum(1 for k in range(1, isqrt(X) + 2) if (2 * k + 2) ** 2 + 1 <= X)
print(f"X = {X}: S = {S}, 1 + N_U + N_D + E = {1 + NU + ND + E} (N_U = {NU}, N_D = {ND}, E = {E}); "
      f"identity holds: {S == 1 + NU + ND + E}; interleaving violations: {interleave_bad}")
