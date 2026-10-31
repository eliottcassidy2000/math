"""opus-2026-10-31-S70: PRIMEGAME over the first ten primes (THM-4852) - driver with assertions.

  A  the two published lists (Wikipedia/OEIS = owner's, Guy 1983): same powers of 2 for the first 40 primes;
     phase law 2^n 7^(D(n)-1) at every state-free configuration for n <= 300 (the C simulator goes to 2500)
  B  fraction anatomy over 2..29: exponent changes, F=U+S defects, lcm = 29#, product 55^2/(2*17^2), pure fractions
  C  the selector cellular automaton in base 29#: 20,000 steps digit by digit
  D  the register-flow graph versus the n = 5 merged tiling line metagraph (THM-4729 table): no adjacency bijection
  E  supersingular / groups-of-five arithmetic (numerology table)
"""
import time
from fractions import Fraction
from math import prod, lcm, isqrt

T0 = time.time()
P = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29]
OWNER = [(17, 91), (78, 85), (19, 51), (23, 38), (29, 33), (77, 29), (95, 23), (77, 19), (1, 17), (11, 13), (13, 11),
         (15, 2), (1, 7), (55, 1)]                     # Wikipedia/OEIS form (the owner's list)
GUY = [(17, 91), (78, 85), (19, 51), (23, 38), (29, 33), (77, 29), (95, 23), (77, 19), (1, 17), (11, 13), (13, 11),
       (15, 14), (15, 2), (55, 1)]                      # Guy 1983 form


def say(*a):
    print(*a, flush=True)


def vec(n):
    v = []
    for p in P:
        e = 0
        while n % p == 0:
            n //= p
            e += 1
        v.append(e)
    assert n == 1
    return v


def spf(n):
    p = 2
    while p * p <= n:
        if n % p == 0:
            return p
        p += 1
    return n


def run(fr, nmax):
    pr = [(vec(a), vec(b)) for a, b in fr]
    x = vec(2)
    pows, phases, inphase = [], [], False
    while True:
        nostate = all(x[i] == 0 for i in range(4, 10))
        if nostate and not inphase:
            phases.append(tuple(x[:4]))
            inphase = True
            if x[0] >= nmax:
                return pows, phases
        if not nostate:
            inphase = False
        for num, den in pr:
            if all(x[i] >= den[i] for i in range(10)):
                x = [x[i] - den[i] + num[i] for i in range(10)]
                break
        if all(x[i] == 0 for i in range(1, 10)):
            pows.append(x[0])


# ---------------------------------------------------------------- A
say('=' * 78)
say('A. Both published lists: powers of 2 and the phase law (n <= 300 here; C simulator to 2500)')
say('=' * 78)
po, pho = run(OWNER, 300)
pg, phg = run(GUY, 300)
primes = [p for p in range(2, 300) if spf(p) == p]
assert po[:len(primes)] == primes and pg[:len(primes)] == primes
for k, (r2, r3, r5, r7) in enumerate(pho):
    n = k + 1
    D = 1 if n == 1 else n // spf(n)
    assert (r2, r3, r5, r7) == (n, 0, 0, D - 1), (n, r2, r3, r5, r7)
assert pho == phg
say(f'  both lists emit 2^p for the {len(primes)} primes p < 300 in order; phase entries equal (n, 0, 0, D(n)-1) for n <= 300')
say('  first phase entries:', pho[:12])

# ---------------------------------------------------------------- B
say('')
say('=' * 78)
say('B. Fraction anatomy over the first ten primes')
say('=' * 78)


def defect(N):
    """F - U - S with F = tau - 2, U = omega - [prime], S = 2^omega - 1 - [squarefree]; None for N = 1."""
    if N == 1:
        return None
    v = vec(N)
    tau = prod(e + 1 for e in v)
    om = sum(1 for e in v if e)
    sqf = all(e <= 1 for e in v)
    prime = (om == 1 and sum(v) == 1)
    return (tau - 2) - (om - prime) - (2 ** om - 1 - sqf)


net = [0] * 10
for k, (a, b) in enumerate(OWNER):
    d = [x - y for x, y in zip(vec(a), vec(b))]
    net = [u + w for u, w in zip(net, d)]
    say(f'  {chr(65 + k)}  {a:>3}/{b:<3} change {d}  defect num/den {defect(a)} / {defect(b)}')
assert net == [-1, 0, 2, 0, 2, 0, -2, 0, 0, 0]
prodO = prod(Fraction(a, b) for a, b in OWNER)
prodG = prod(Fraction(a, b) for a, b in GUY)
assert prodO == Fraction(55 ** 2, 2 * 17 ** 2) and prodG == prodO * Fraction(15, 2)   # Guy has 15/14 AND 15/2
assert lcm(*[b for a, b in OWNER]) == lcm(*[b for a, b in GUY]) == prod(P) == 6469693230
parts = [a for a, b in OWNER] + [b for a, b in OWNER]
shapes = {defect(x) for x in parts if x > 1}
assert shapes == {0, -2, -3}
assert all((defect(x) == 0) == (spf(x) == x) for x in parts if x > 1)       # balanced parts = the primes
pure_owner = [(a, b) for a, b in OWNER if a == 1 or b == 1]
pure_guy = [(a, b) for a, b in GUY if a == 1 or b == 1]
assert pure_owner == [(1, 17), (1, 7), (55, 1)] and pure_guy == [(1, 17), (55, 1)]
say(f'  owner list: net change {net}, product = 55^2/(2*17^2) = {prodO}; Guy list product = {prodG} = owner * 15/2')
say(f'  lcm of denominators = 29# = {prod(P)} (both lists)')
say(f'  balanced (F=U+S) parts are exactly the primes; defects present {sorted(shapes)}; no p^3 or p^2qr')
say(f'  pure fractions: owner {pure_owner}  (the owner\'s {{7, 17, 55}}); Guy {pure_guy}')

# ---------------------------------------------------------------- C
say('')
say('=' * 78)
say('C. PRIMEGAME as a radius-1 selector cellular automaton in base 29#')
say('=' * 78)
B = prod(P)


def digits(x):
    d = []
    while x:
        d.append(x % B)
        x //= B
    return d or [0]


def ca_step(d, p, q):
    x = d + [0, 0]
    y = [x[i] // q + (x[i + 1] % q) * (B // q) for i in range(len(d) + 1)] + [0]
    z = [(p * y[i]) % B + ((p * y[i - 1]) // B if i > 0 else 0) for i in range(len(d) + 2)]
    for i in range(len(z)):
        assert z[i] < B
    while len(z) > 1 and z[-1] == 0:
        z.pop()
    return z


for name, fr in (('owner', OWNER), ('Guy', GUY)):
    x, d = 2, digits(2)
    for t in range(20000):
        k = next(k for k, (p, q) in enumerate(fr) if d[0] % q == 0)
        assert next(k2 for k2, (p, q) in enumerate(fr) if x % q == 0) == k   # selector = units digit
        p, q = fr[k]
        d = ca_step(d, p, q)
        x = x * p // q
        assert sum(v * B ** i for i, v in enumerate(d)) == x
    say(f'  {name} list: 20,000 steps reproduced digit by digit (radius 1; rule chosen by the units digit)')

# ---------------------------------------------------------------- D
say('')
say('=' * 78)
say('D. Register-flow graph versus the n = 5 merged tiling line metagraph (THM-4729)')
say('=' * 78)
arcs = set()
for a, b in OWNER:
    for p in P:
        if b % p == 0:
            for q in P:
                if a % q == 0:
                    arcs.add((p, q))
pairs = {frozenset(e) for e in arcs}
assert len(arcs) == 24 and len(pairs) == 19 and all(len(e) == 2 for e in pairs)
blue = [(0, 7), (3, 8), (4, 9), (5, 6), (5, 7), (5, 9), (7, 8), (8, 9)]
black = [(1, 6), (2, 4), (2, 5), (2, 6), (2, 7), (2, 8), (4, 7), (5, 5), (5, 8), (6, 7), (6, 8), (8, 8)]   # each x2
simple = {frozenset(e) for e in blue + black if e[0] != e[1]}
assert len(blue) == 8 and 2 * len(black) == 24 and len(simple) == 18
loops = [e for e in black if e[0] == e[1]]
pure_blue = [v for v in range(10) if not any(v in e for e in black)]
assert pure_blue == [0, 3, 9] and loops == [(5, 5), (8, 8)]
say(f'  flow graph: 24 arcs, 19 simple pairs, no loops; metagraph: 8 blue + 24 black lines, 18 simple pairs, black loops {loops}')
say(f'  pure-blue nodes {pure_blue} = TT5, 1->C3->1, T5; pure fractions 3 -- equal counts only')
say('  => no bijection primes -> nodes preserves adjacency (19 != 18 simple edges; loops on one side only)')
deg_flow = sorted(sum(1 for e in pairs if p in e) for p in P)
deg_meta = sorted(sum(1 for e in simple if v in e) for v in range(10))
say(f'  simple degree sequences: flow {deg_flow}; metagraph {deg_meta}')
assert deg_flow != deg_meta

# ---------------------------------------------------------------- E
say('')
say('=' * 78)
say('E. Supersingular primes and fives (CITED facts, NUMEROLOGY readings)')
say('=' * 78)
OGG = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 41, 47, 59, 71]
assert P == OGG[:10] and [p for p in OGG if p not in P] == [31, 41, 47, 59, 71]
assert [p for p in range(2, 72) if spf(p) == p and p not in OGG][0] == 37
FERMAT = [3, 5, 17, 257, 65537]
assert [p for p in FERMAT if p in OGG] == [3, 5, 17]
cyc, x = [], 17
while True:
    cyc.append(x)
    x = 3 * x - 1
    while x % 2 == 0:
        x //= 2
    if x == 17:
        break
assert sorted(cyc) == [17, 25, 37, 41, 55, 61, 91]
say(f'  registers = the ten smallest Ogg primes; unused Ogg five {[p for p in OGG if p not in P]}; first non-supersingular prime 37')
say(f'  Fermat cap Ogg = [3, 5, 17]; 3x-1 sporadic cycle (odd terms) {sorted(cyc)} contains 17, 55, 91 = 7*13 (NUMEROLOGY)')
say('')
say(f'ALL ASSERTIONS PASSED ({time.time() - T0:.1f}s)')
