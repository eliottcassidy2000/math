import math, json

def q(n):
    out = 0
    while n:
        a = n // 2
        if n % 2:
            out += a + a % 2
        n = a
    return out

def step(s, c):
    h, r, e, d = s
    if c == 'O':
        return (h, 2 * r, 0, d)
    if c == 'Z':
        return (h + r - e, 2 * r - d, d, d)
    return (h + r - e, h + r, h - e, h)
seeds = {(1, 2, 0, 1)}
for _ in range(3):
    seeds |= {step(s, c) for s in list(seeds) for c in 'OZJ'}
checks = 0
for M in range(2, 34, 2):
    k = 0
    u = M
    while u % 2 == 0:
        k += 1
        u //= 2
    b = max(k, 2)
    while pow(2, b, u) != 1 % u:
        b += 1
    x = 2 ** b - 1
    c = 2 ** b - 2
    for h, r, e, d in seeds:
        actual = set()
        for a in range(M):
            A = a.bit_length() + k
            m = a
            for t in range(1, M // 2 + 1):
                m += x << A + (t - 1) * (b + k)
                if not (m % M == a and q(m) % M == (q(a) + t * c) % M):
                    raise RuntimeError('Verification failed')
                H = h + r * m - d * q(m) - e * (m % 2)
                actual.add(H % M)
                checks += 1
        g = math.gcd(r, d, M // 2)
        expected = {n for n in range(M) if (n - h) % (2 * g) == 0 or (n - h - r + e) % (2 * g) == 0}
        if not actual == expected:
            raise RuntimeError((M, (h, r, e, d), actual, expected))
print(json.dumps({'status': 'PASS', 'physical_seeds': len(seeds), 'checks': checks, 'even_moduli_tested': 16}, indent=2))
