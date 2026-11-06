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
for M in range(2, 66, 2):
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
        if math.gcd(d, u) != 1:
            continue
        residues = set()
        m = 0
        for t in range(1, M // 2 + 1):
            m += x << k + (t - 1) * (b + k)
            Q = q(m)
            if not (m % M == 0 and Q % M == t * c % M):
                raise RuntimeError('Verification failed')
            for L in (0, 3):
                z = m << L
                H = h + r * z - d * q(z) - e * (z % 2)
                if not H % M == (h - d * t * c) % M:
                    raise RuntimeError('Verification failed')
                checks += 1
            residues.add(H % M)
        if not residues == set(range(1, M, 2)):
            raise RuntimeError('Verification failed')
print(json.dumps({'status': 'PASS', 'physical_seeds': len(seeds), 'checks': checks, 'even_moduli_tested': 32}, indent=2))
