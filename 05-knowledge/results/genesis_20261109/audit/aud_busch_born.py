"""AUDIT G1: values h in S(11) born at their Busch bound, n(h) = B(h) = max{m : f(m) <= h} (census T-lines)."""
S = {n: set(int(x) for x in open('S_%d.txt' % n)) for n in range(1, 12)}
def f(m):   # Busch floor for strong m-tournaments, m >= 3: min 3^a 5^b with 2a + 3b = m - 1
    return min(3 ** ((m - 1 - 3 * b) // 2) * 5 ** b for b in range(m) if m - 1 - 3 * b >= 0 and (m - 1 - 3 * b) % 2 == 0)
def B(h):
    m = 3
    while f(m + 1) <= h: m += 1
    return m
print('f(3..12) =', [f(m) for m in range(3, 13)])
print('h in S(11) with n(h) = B(h):', [h for h in sorted(S[11]) if h >= 3 and min(k for k in range(1, 12) if h in S[k]) == B(h)])
