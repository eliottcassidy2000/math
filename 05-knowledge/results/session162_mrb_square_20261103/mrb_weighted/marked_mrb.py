from fractions import Fraction

def floor_root_scaled(n, Q):
    target = n * Q ** n
    lo, hi = (0, Q)
    while lo < hi:
        m = (lo + hi + 1) // 2
        if (Q + m) ** n <= target:
            lo = m
        else:
            hi = m - 1
    if not (Q + lo) ** n <= target < (Q + lo + 1) ** n:
        raise AssertionError()
    return lo
N = 200
Q = 10 ** 9
p = {i: floor_root_scaled(i, Q) for i in range(2, N + 3)}
ne = sum((i % 2 == 0 for i in range(2, N + 1)))
no = N - 1 - ne
D = sum(((1 if i % 2 == 0 else -1) * p[i] for i in range(2, N + 1)))
lo = 2 * D - 2 * no - 2 * p[N + 1] - 2 + p[N + 2]
hi = 2 * D + 2 * ne - p[N + 1]
if not (D, p[N + 1], p[N + 2], lo, hi) == (201254860, 26735756, 26626878, 375664886, 375774164):
    raise AssertionError()
print('N Q even odd D p201 p202:', N, Q, ne, no, D, p[201], p[202])
print('rigorous bracket:', f'[{lo}/{2 * Q}, {hi}/{2 * Q}]')
K = 2 ** (N - 2) + 1
E = K * sum((p[i] for i in range(2, N + 1) if i % 2 == 0))
O = K * sum((p[i] for i in range(2, N + 1) if i % 2 == 1))
if not E - O == K * D:
    raise AssertionError()
print('normalized even-minus-odd marked count:', f'(E-O)/K={D}')
