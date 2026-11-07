"""Exact quadratic path readout of Session212's two-apex cycle comb."""
def require(test, message):
    if not test:
        raise RuntimeError(message)

def comb(m, q, tau=0, reversed_apex=False):
    require(m >= 1 and 0 <= q < 1 << m and tau in (0, 1), 'invalid comb parameters')
    n = m + 4
    a, b = n - 2, n - 1
    out = [0] * n
    def edge(i, j):
        out[i] |= 1 << j
    for i in range(m + 2):
        for j in range(i + 1, m + 2):
            edge(i, j)
        edge(a, i) if i == 1 else edge(i, a)
        incoming = 0 if i == 0 else tau if i == 1 else (q >> (i - 2)) & 1
        edge(i, b) if incoming else edge(b, i)
    edge(b, a) if reversed_apex else edge(a, b)
    return out

def coefficients(m, tau, reversed_apex):
    require(m >= 1 and tau in (0, 1), 'invalid coefficient parameters')
    v = int(reversed_apex)
    power = 1 << m
    base = power + 1 + (1 - v) * (4 * power - 2) + tau * (2 * power + 2) - tau * (1 - v) * power
    linear = []
    for i in range(m):
        u = (1 << (m - 1 - i)) + (1 << i)
        first, last = int(i == 0), int(i == m - 1)
        d = (1 << (m + 1 - i)) * 3 ** i + (1 << (i + 2)) - 2 * last
        linear.append(d - (1 - v) * (u - first + last) - tau * (u + first - last))
    pair = {d: (1 << (m - d)) * (3 ** (d - 1) + 1) + (1 << (d - 1)) + int(d == 1) for d in range(1, m)}
    return base, linear, pair

def path_formula(m, q, tau=0, reversed_apex=False):
    require(0 <= q < 1 << m, 'invalid word')
    base, linear, pair = coefficients(m, tau, reversed_apex)
    selected = [i for i in range(m) if q >> i & 1]
    return base + sum(linear[i] for i in selected) - sum(pair[j - i] for i in selected for j in selected if i < j)

def spectrum(limit):
    """All-length comb H values at most limit, certified by H >= 2**m+1."""
    require(limit >= 3, 'limit must be at least three')
    max_m = (limit - 1).bit_length() - 1
    found = {}
    scanned = 0
    for m in range(1, max_m + 1):
        for tau in (0, 1):
            for rev in (False, True):
                base, linear, pair = coefficients(m, tau, rev)
                # Gray traversal, maintaining every first discrete derivative.
                value, previous = base, 0
                marginal = linear.copy()
                for k in range(1 << m):
                    q = k ^ (k >> 1)
                    if k:
                        flip = q ^ previous
                        i = flip.bit_length() - 1
                        direction = 1 if q & flip else -1
                        value += direction * marginal[i]
                        for j in range(m):
                            if j != i:
                                marginal[j] -= direction * pair[abs(i - j)]
                    scanned += 1
                    require(value >= (1 << m) + 1 and value % 2 == 1, 'monotonicity/parity control failed')
                    if value <= limit and value not in found:
                        found[value] = (m, q, tau, rev)
                    previous = q
    return found, scanned, max_m
