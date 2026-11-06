import json

def q(n):
    if not n:
        return 0
    a = n // 2
    return q(a) + (a + a % 2 if n % 2 else 0)

def step(s, c):
    h, r, e, d = s
    if c == 'O':
        return (h, 2 * r, 0, d)
    if c == 'Z':
        return (h + r - e, 2 * r - d, d, d)
    return (h + r - e, h + r, h - e, h)
seeds = {(1, 2, 0, 1)}
for _ in range(5):
    seeds |= {step(s, c) for s in list(seeds) for c in 'OZJ'}
checks = 0
for k in range(1, 8):
    mod = 2 ** k
    for s in seeds:
        h, r, e, d = s
        residues = set()
        for t in range(1, mod // 2 + 1):
            m = 3 * sum((1 << k + j * (k + 2) for j in range(t)))
            if not q(m) % mod == 2 * t % mod:
                raise RuntimeError('Verification failed')
            for L in (0, 1, 5):
                z = m << L
                H = h + r * z - d * q(z) - e * (z % 2)
                if not H % mod == (h - 2 * d * t) % mod:
                    raise RuntimeError('Verification failed')
                checks += 1
            residues.add(H % mod)
        if not residues == set(range(1, mod, 2)):
            raise RuntimeError('Verification failed')
print(json.dumps({'status': 'PASS', 'physical_seeds': len(seeds), 'checks': checks, 'maximum_k': 7}, indent=2))
