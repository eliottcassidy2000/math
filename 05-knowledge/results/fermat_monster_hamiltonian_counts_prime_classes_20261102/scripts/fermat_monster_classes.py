import math
def primes_below(N):
    s = bytearray([1]) * N; s[0:2] = b'\x00\x00'
    for i in range(2, int(N ** 0.5) + 1):
        if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(N) if s[i]]
def isprime(q):
    if q < 2: return False
    if q % 2 == 0: return q == 2
    d = 3
    while d * d <= q:
        if q % d == 0: return False
        d += 2
    return True
PAIRS = {1: 35, 13: 23, 25: 11, 7: 29, 19: 17, 31: 5}
restricted = set(PAIRS); partner = set(PAIRS.values())
def kind(p):
    if p in (2, 3): return 'special'
    r = p % 36
    return 'restricted(%d->%d)' % (r, PAIRS[r]) if r in restricted else ('partner(%d)' % r if r in partner else 'free(%d)' % r)
print('--- Fermat primes ---')
for n in range(5):
    F = 2 ** (2 ** n) + 1
    s = math.isqrt(F) + 1; parts = []
    while len(parts) < 3:
        q = s * s - F
        if q > 1 and isprime(q): parts.append((q, s))
        s += 1
    print('F_%d = %d: mod 36 = %d, mod 3 = %d, mod 4 = %d, kind %s; first square-sum prime partners %s' % (n, F, F % 36, F % 3, F % 4, kind(F), parts))
P = primes_below(10 ** 5); Ps = set(P)
def degree(p, bound):
    return sum(1 for s in range(math.isqrt(p) + 1, math.isqrt(p + bound) + 1) if 0 < s * s - p < bound and (s * s - p) in Ps and s * s - p != p)
print('degrees below 10^5 of 3, 5, 17, 257, 65537:', [degree(F, 10 ** 5) for F in (3, 5, 17, 257, 65537)])
print('--- the 25 primes below 100 by mod-36 kind ---')
layers = {'Fi22': [2,3,5,7,11,13], 'B\\Fi22': [17,19,23,31,47], 'M\\B': [29,41,59,71], 'nonSS-1': [37,43,53,61,67], 'nonSS-2': [73,79,83,89,97]}
for name, L in layers.items():
    print('  %-8s' % name, [(p, kind(p)) for p in L])
SS = [2,3,5,7,11,13,17,19,23,29,31,41,47,59,71]; NS = [37,43,53,61,67,73,79,83,89,97]
print('restricted among the 15 supersingular primes: %d; among the 10 non-supersingular primes below 100: %d' % (
    sum(1 for p in SS if p % 36 in restricted), sum(1 for p in NS if p % 36 in restricted)))
print('isolated primes below 100 and their kinds:', [(p, kind(p)) for p in (37, 43, 67)], '; their first partners:', [(107, kind(107)), (101, kind(101)), (257, kind(257))])
# square classes mod 36 = {1, 13, 25}; 7 * squares = {7, 19, 31}
print('squares of units mod 36:', sorted({(u * u) % 36 for u in range(36) if math.gcd(u, 36) == 1}), '; 7*squares:', sorted({(7 * u * u) % 36 for u in range(36) if math.gcd(u, 36) == 1}))
# Chebyshev counts mod 36 among primes below 10^5 by class
c = {}
for p in P:
    if p > 3: c[p % 36] = c.get(p % 36, 0) + 1
print('counts below 10^5 by class mod 36 (restricted classes marked *):', {('%d*' % r if r in restricted else r): c[r] for r in sorted(c)})
