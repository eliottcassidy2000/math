"""Exact cyclic clocks, integral recipes, and literal two-source repairs."""
from pathlib import Path
from fractions import Fraction
from functools import lru_cache
from itertools import product
import hashlib
import json
import math

HERE = Path(__file__).resolve().parent
T = (1, 0, 2, 1)

def need(ok, message):
    if not ok:
        raise RuntimeError(message)

def step(s, letter):
    D, E, U, d = s
    if letter == 'O':
        return D, 0, 2 * U, d
    if letter == 'Z':
        return D - E + U, d, 2 * U - d, d
    need(letter == 'J', 'known source letter')
    return D - E + U, D - E, D + U, D

def apply(s, word):
    for letter in word:
        s = step(s, letter)
    return s

@lru_cache(None)
def g(n):
    if n == 1:
        return 1
    if n % 2 == 0:
        return g(n // 2)
    h = (n + 1) // 2
    return g(h) + h - h % 2

def eval_bits(word):
    n, H = 0, 1
    for bit in word:
        if bit == '0':
            H += n + n % 2
            n *= 2
        else:
            n = 2 * n + 1
    return n, H

def modular_word(M, u, z, extra=0):
    """The inherited actual binary compiler, not an inverse on a cone."""
    q, k = M, 0
    while q % 2 == 0:
        q //= 2
        k += 1
    aa = 2
    while pow(2, aa, q) != 1 % q:
        aa += 1
    delta = ((1 << aa) - 1) * ((1 << k) - 1) + 1
    need(math.gcd(delta, M) == 2, 'joint compiler loop content')
    rr, m = (u - 1) % M, M // 2
    j = 0 if m == 1 else (((z - g(rr + 1)) // 2) * pow(delta // 2, -1, m)) % m
    j += (1 + extra) * m
    return ('1' * aa + '0' * k) * j + format(rr, f'0{(M - 1).bit_length()}b')

def egcd(a, b):
    if b == 0:
        return a, 1, 0
    dd, x, y = egcd(b, a % b)
    return dd, y, x - (a // b) * y

def coefficients(bits):
    k = len(bits)
    A = sum(e << i for i, e in enumerate(bits))
    R = sum(e << (k - 1 - i) for i, e in enumerate(bits))
    C = sum(bits[j] * bits[i] * (1 << (i - j - 1))
            for i in range(k) for j in range(i))
    C += sum(bits[i - 1] * bits[i] for i in range(k))
    return A, R, C

def clock(a, b):
    need(a > 0 and b % 2 == 1 and math.gcd(a, b) == 1 and b < 2 * a < 2 * b,
         'admissible slope')
    k = 1
    while (2 ** k - 1) % b:
        k += 1
    B, A = 2 ** k, a * (2 ** k - 1) // b
    bits = [(A >> i) & 1 for i in range(k)]
    AA, R, C = coefficients(bits)
    need(A == AA, 'binary coefficient')
    return dict(a=a, b=b, k=k, B=B, A=A, R=R, C=C,
                K=a * R - b * C, word=''.join('Z' if e else 'O' for e in bits))

def construct(a, b, target, extra=0):
    """Return a compressed actual source; do not expand its repeated graph."""
    need((target - a + b) % 2 == 0, 'permitted target parity')
    cc = clock(a, b)
    K = cc['K']
    need(K > 0 and K % 2 == 0, 'positive even clock')
    ell = 1
    while (1 << ell) * (2 * a - b) <= b:
        ell += 1
    c = (1 << ell) * (2 * a - b)
    gg, alpha, beta = egcd(c, b)
    need(gg == 1, 'coprime seed coefficients')
    z = (target - c + b + a) // 2
    D0, U0 = (1 - 2 * beta * z) % K, (1 + 2 * alpha * z) % K
    need(D0 % 2 == U0 % 2 == 1, 'odd seed residues')
    need((c * U0 - b * D0 - a - target) % K == 0, 'seed response residue')
    u = (U0 + 1) // 2
    hres = (D0 - u + u % 2) % K
    while True:
        seed_bits = modular_word(K, u, hres, extra)
        nn, gg = eval_bits(seed_bits)
        L = nn + 1
        D, U = gg + L - L % 2, 2 * L - 1
        need(L >= 2 and D % K == D0 and U % K == U0, 'actual modular seed')
        S = c * U - b * D - a
        need(U >= D >= L and S >= (c - b) * D - a, 'physical overshoot bound')
        if (c - b) * D - a >= target + K:
            break
        extra += 1
    h = (S - target) // K
    need(h >= 1 and S - h * K == target, 'integral exact response correction')
    return dict(a=a, b=b, target=target, clock=cc, seed_L=L,
                seed_bits=seed_bits, seed_profile=(D, 1, U, 1),
                prefix='O' * ell + 'Z', repetitions=h)

def N(v, w):
    D, E, U, d = v
    F, G, V, dd = w
    return (D - E + U) * F - (D - E + d) * G + (D + 3 * U - d) * V - (U + E) * dd

def p(s):
    return s[0] + s[3]

def rho(s):
    return s[0] - s[1] + s[2] - s[3]

def graph(word):
    rows, root = [2, 0], 0
    for letter in word:
        n = len(rows)
        rows.append(0)
        for j in range(n):
            if letter != 'O' and j == root:
                rows[n] |= 1 << j
            else:
                rows[j] |= 1 << n
        if letter == 'J':
            root = n
    return rows, root

def counts(rows, root=None):
    """Independent subset recurrence. Cycle DP fixes vertex 0."""
    n, full = len(rows), (1 << len(rows)) - 1
    end = [[0] * n for _ in range(full + 1)]
    fixed = [[0] * n for _ in range(full + 1)]
    H = [1] + [0] * full
    for j in range(n):
        end[1 << j][j] = 1
    fixed[1][0] = 1
    starts = [[0] * n for _ in range(full + 1)] if root is not None else None
    for S in range(1, full + 1):
        for j in range(n):
            if not (S >> j) & 1:
                continue
            Q = S ^ (1 << j)
            if Q:
                end[S][j] = sum(end[Q][i] for i in range(n)
                                if (Q >> i) & 1 and (rows[i] >> j) & 1)
                if S & 1:
                    fixed[S][j] = sum(fixed[Q][i] for i in range(n)
                                      if (Q >> i) & 1 and (rows[i] >> j) & 1)
            if starts is not None:
                starts[S][j] = 1 if not Q else sum(starts[Q][i] for i in range(n)
                                                  if (Q >> i) & 1 and (rows[j] >> i) & 1)
        H[S] = sum(end[S])
    cyc = sum(fixed[full][j] for j in range(1, n) if rows[j] & 1)
    if root is None:
        return H[full], cyc
    U = sum(H[S] * starts[full ^ S][root] for S in range(full + 1))
    return (H[full], end[full][root], U, H[full ^ (1 << root)]), starts[full][root], cyc

def graft_N(vword, wword):
    V, vr = graph(vword)
    W, wr = graph(wword)
    n = len(V)
    V = [row | (1 << n) for row in V] + [1 << vr]
    V[vr] ^= 1 << n
    entries = {vr, n}
    n = len(V)
    full = (1 << len(W)) - 1
    converse = [full ^ (1 << i) ^ row for i, row in enumerate(W)]
    rows = V + [row << n for row in converse]
    for i in range(n):
        for j in range(len(W)):
            if j == wr and i in entries:
                rows[n + j] |= 1 << i
            else:
                rows[i] |= 1 << (n + j)
    return rows

def main():
    clocks, recipes = [], []
    materialized, orbit_checks, formal_checks = 0, 0, 0
    basis = ((1, 0, 0, 0), (0, 0, 1, 0), (0, 1, 0, 1))
    for b in range(3, 32, 2):
        for a in range((b + 1) // 2, b):
            if math.gcd(a, b) > 1:
                continue
            cc = clock(a, b)
            B, A, R, C, K, W = (cc[x] for x in ('B', 'A', 'R', 'C', 'K', 'word'))
            need(W[-1] == 'Z' and (W[0] == 'Z') == (a % 2 == 1), 'clock boundary parity')
            need('O' in W and 'Z' in W and B // 2 <= A < B - 1, 'mixed interior word')
            for s in basis:
                D, E, U, d = s
                ss = apply(s, W)
                need(ss == (D + A * U - C * d, d, B * U - R * d, d), 'E=d basis transfer')
                need(a * ss[2] - b * ss[0] == a * U - b * D - K * d, 'E=d response transfer')
                formal_checks += 1
            eps = [int(ch == 'Z') for ch in W]
            u, states, terms = Fraction(R, B - 1), [], []
            for i, e in enumerate(eps):
                need(0 < u < 1 and e == (2 * u).__floor__(), 'periodic cover digit')
                states.append(u)
                terms.append(e * (u - eps[i - 1]))
                u = 2 * u - e
            need(u == states[0] and sum(terms) == Fraction(K, b), 'periodic increment')
            runs = []
            for i, e in enumerate(eps):
                if e and not eps[i - 1]:
                    j = 1
                    while eps[(i + j) % len(eps)]:
                        j += 1
                    v = states[(i + j) % len(eps)]
                    part = v + (1 - v) / 2 ** j
                    need(part > 0 and part == sum(terms[(i + h) % len(eps)] for h in range(j)), 'run contribution')
                    runs.append(part)
            need(sum(runs) == Fraction(K, b) and K > 0 and K % 2 == 0, 'positive even decrement')
            for shift in range(len(eps)):
                for complement in (False, True):
                    es = eps[shift:] + eps[:shift]
                    if complement:
                        es = [1 - e for e in es]
                    aa, rr, ccc = coefficients(es)
                    need((aa * b) % (B - 1) == 0, 'rotated rational numerator')
                    aaa = aa * b // (B - 1)
                    need(aaa * rr - b * ccc == K, 'cyclic/complement decrement')
                    need(apply((7, es[-1] * 3, 19, 3), ''.join('Z' if e else 'O' for e in es)) ==
                         (7 + aa * 19 - ccc * 3, es[-1] * 3, B * 19 - rr * 3, 3), 'formal rotated endpoint sidecar')
                    orbit_checks += 1
            clocks.append(cc)
            for target in range(-10, 11):
                if (target - a + b) % 2:
                    continue
                rec = construct(a, b, target)
                h = rec['repetitions']
                seed = rec['seed_profile']
                q = apply(seed, rec['prefix'])
                need(q[1] == q[3] == 1 and a * q[2] - b * q[0] - K * h == target, 'compressed physical recipe')
                if h * len(W) <= 10000:
                    for _ in range(h):
                        q = apply(q, W)
                    need(q[1] == q[3] == 1 and a * q[2] - b * q[0] == target, 'materialized response')
                    materialized += 1
                recipes.append(dict(a=a, b=b, target=target, prefix=rec['prefix'],
                                    seed_bits=len(rec['seed_bits']), seed_sha256=hashlib.sha256(rec['seed_bits'].encode()).hexdigest(),
                                    repetitions_bits=h.bit_length(), K=K))
    need(clock(9, 17)['K'] == 2 and clock(10, 17)['K'] == 30, 'period-only hostile')
    # Independently evaluate short literal O/Z graphs, including clocks starting Z.
    literals = []
    for W in ('Z' + clock(3, 5)['word'], 'Z' + clock(5, 7)['word'],
              'Z' + clock(5, 9)['word'], 'ZOOZ', 'ZOZOZ', 'ZJOZO'):
        rows, root = graph(W)
        measured, start, cyc = counts(rows, root)
        need(measured == apply(T, W) and start == 1, 'literal source response')
        literals.append(dict(word=W, profile=measured, start=start, cycles=cyc, order=len(rows)))
    # Algebraic J transfer, including arbitrary basis vectors, not just physical samples.
    raw_basis = [tuple(int(i == j) for i in range(4)) for j in range(4)]
    for v in raw_basis:
        for w in raw_basis:
            need(N(v, step(w, 'J')) - N(step(v, 'J'), w) == rho(v) * p(w) - p(v) * rho(w), 'full defect basis')
    grafts = []
    vword, wword = 'ZOOZ', 'ZOZOZ'
    v, w = apply(T, vword), apply(T, wword)
    need(Fraction(rho(v), p(v)) == Fraction(rho(w), p(w)) == Fraction(9, 4), 'same actual response fiber')
    need({v, step(w, 'J')} != {step(v, 'J'), w}, 'nontrivial unordered pair')
    common = N(v, step(w, 'J'))
    need(common == N(step(v, 'J'), w) == 7649, 'exact C2 collision')
    for words in ((vword, wword + 'J'), (vword + 'J', wword)):
        rows = graft_N(*words)
        measured = counts(rows)
        need(measured == (common, 2), 'literal C2 collision graft')
        grafts.append(dict(words=words, H=measured[0], C=measured[1], order=len(rows), rows=rows))
    need(N(T, apply(T, 'OJ')) == 41 and N(step(T, 'J'), step(T, 'O')) == 45, 'universal swap hostile')
    # Explicit new-parity companion, whose complete positive coordinates are manageable.
    companion = construct(5, 9, 24)
    q = apply(companion['seed_profile'], companion['prefix'])
    for _ in range(companion['repetitions']):
        q = apply(q, companion['clock']['word'])
    v = apply(T, 'ZJOZO')
    need(5 * q[2] - 9 * q[0] == 24 and q[1] == q[3] == 1, 'new odd-numerator companion')
    need(Fraction(rho(v), p(v)) == Fraction(rho(q), p(q)) == Fraction(14, 5), 'new physical fiber')
    need(N(v, step(q, 'J')) == N(step(v, 'J'), q), 'new target-preserving repair')
    companion['final_profile'] = q
    companion['left_profile'] = v
    companion['common_C2_count'] = N(v, step(q, 'J'))
    companion['claim'] = 'Compressed physical source recipe, not a literal graph enumeration at this order.'
    (HERE / 'COMPANION_5_9.json').write_text(json.dumps(companion, sort_keys=True, indent=2) + '\n')
    # Smaller exact companion discovered by bounded search, certified without it.
    short_seed = 'ZOZZOZZOOOZ'
    short_word = short_seed + 'ZZOOOZ' * 2
    small = apply(T, short_word)
    need(apply(T, short_seed) == (1463, 1, 2639, 1), 'small companion seed')
    need(small == (6003419, 1, 10806159, 1) and g(small[2]) == small[0], 'small companion exact arithmetic')
    need(5 * small[2] - 9 * small[0] == 24, 'small companion response')
    hsmall = (small[2] - 3) // 18
    need(hsmall == 600342 and 4952 * hsmall + 421 == 2972894005, 'simultaneous target arithmetic')
    need(N(v, step(small, 'J')) == N(step(v, 'J'), small) == 2972894005, 'small companion C2 target')
    need(len(short_word) + len('ZJOZO') + 6 == 34, 'large graph certificate order')
    short_record = dict(left_word='ZJOZO', right_word=short_word, right_profile=small,
                        h=hsmall, common_C2_count=2972894005, graph_order=34,
                        verification='Exact source word, independent g recurrence and bilinear graft law; no literal 34-vertex count.')
    (HERE / 'SHORT_COMPANION_5_9.json').write_text(json.dumps(short_record, sort_keys=True, indent=2) + '\n')
    # Verify the general linear readout and the exact h parameterization.
    for vv in raw_basis:
        V, E, HH, dd = vv
        AA, BB = V - E + 3 * HH - 2 * dd, 2 * V - E + 4 * HH - dd
        for DD, UU in ((0, 0), (1, 0), (0, 1)):
            need(N(vv, step((DD, 1, UU, 1), 'J')) == AA * DD + BB * UU + dd - HH, 'fixed-target readout identity')
    for cc in clocks:
        a, b = cc['a'], cc['b']
        for hh in range(1, 7):
            DD, UU = 2 * a * hh - 1, 2 * b * hh + 3
            need(a * UU - b * DD == 3 * a + b and DD % 2 == UU % 2 == 1, 'exact arithmetic fiber parameterization')
    for n in range(3, 5001, 2):
        need(5 * n - 9 * g(n) != 24, 'finite companion absence through U=5000')
        s = 3 * n - 4 * g(n)
        need(s <= 0 or s % 4 == 1, 'even-denominator hostile')
    out = dict(status='PASS', max_odd_denominator=31, clocks=clocks,
               formal_E_equals_d_checks=formal_checks, cyclic_complement_checks=orbit_checks,
               compressed_recipes=recipes, materialized_recipes=materialized,
               literal_sources=literals, literal_C2_grafts=grafts,
               short_companion_5_9=short_record,
               companion_5_9=dict(seed_L=companion['seed_L'], prefix=companion['prefix'],
                                  repetitions=companion['repetitions'], final_U_bits=q[2].bit_length(),
                                  common_C2_count_bits=companion['common_C2_count'].bit_length()),
               dependencies={f.name: hashlib.sha256(f.read_bytes()).hexdigest()
                             for f in sorted((HERE / 'dependencies').glob('*.md'))},
               scope='Exact signed parity-class coverage and physical fiber-level C2 reassociation; no universal positive-target coverage.')
    content = json.dumps(out, sort_keys=True, indent=2) + '\n'
    (HERE / 'VALIDATION.json').write_text(content)
    print(content, end='')

if __name__ == '__main__':
    main()
