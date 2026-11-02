#!/usr/bin/env python3
"""Independent audit of THM-4857 (the GF(2) automaton of a tournament).

Own code only (stdlib).  The similarity invariant used here is the list of INVARIANT FACTORS of A over GF(2),
obtained from a Smith normal form of xI - A over GF(2)[x] (the session used ranks of p(A)^k instead).
Tournament classes are computed with an own canonical form (colour refinement + permutations inside cells).

Checks:
  G0  sanity of the Smith-form invariant: brute force over ALL n x n GF(2) matrices, n <= 4, reproduces the classical
      number of similarity classes (series prod_i 1/(1 - 2x^i)); invariant unchanged under A -> A^T and under
      random conjugation.
  G1  classes / merged nodes / similarity types for n = 3..7, n = 6 collisions, the TT6 collision.
  G2  counting: exact number c(n) of similarity classes of n x n GF(2) matrices vs merged nodes; the session's crude
      bound p(n) 2^n; A000568 from Davis' formula.
  G3  circulant rank formula (all connection sets, odd q <= 27 incl. composite), all-ones kernel at q = 1 mod 4.
  G4  Paley: rank by q mod 8 (primes q = 3 mod 4 below 400), image code dimension / minimum distance at q = 7, 23,
      which QR side the image is, and the doubly-regular identities A A^T = I (q = 3 mod 8), A^2 = A (q = 7 mod 8),
      also on the non-circulant Paley tournament of GF(27).
  G5  the n = 5 table of THM-4857 (rank, nilpotent depth, cycle multiset, char poly) for the twelve classes of
      ../../groups_of_five_primegame_20261031/merged_metagraph_n5.json (read as data only).
usage (from the repo root): python3 04-computation/experiments/circulant_fractran_goldbach_20261101/audit_lenses/audit_gf2_automaton.py
"""
import itertools
import json
import math
import os
import random
import sys
import time
from collections import Counter, defaultdict

HERE = os.path.dirname(os.path.abspath(__file__))
FAILS = []


def report(name, ok, detail=''):
    print(f'  [{"OK" if ok else "MISMATCH"}] {name}' + (f' -- {detail}' if detail else ''))
    if not ok:
        FAILS.append(name)


# ------------------------------------------------------------------ GF(2)[x] arithmetic, polynomials as int bitmasks
def pdeg(a):
    return a.bit_length() - 1


def pmul(a, b):
    r = 0
    while b:
        if b & 1:
            r ^= a
        a <<= 1
        b >>= 1
    return r


def pdivmod(a, b):
    q = 0
    db = pdeg(b)
    while a and pdeg(a) >= db:
        s = pdeg(a) - db
        q ^= 1 << s
        a ^= b << s
    return q, a


def pgcd(a, b):
    while b:
        a, b = b, pdivmod(a, b)[1]
    return a


def pstr(p):
    if p == 0:
        return '0'
    t = []
    for k in range(pdeg(p), -1, -1):
        if p >> k & 1:
            t.append('1' if k == 0 else ('x' if k == 1 else f'x^{k}'))
    return '+'.join(t)


# ------------------------------------------------------------------ Smith normal form of xI - A over GF(2)[x]
def invariant_factors(A):
    """A: list of lists of 0/1 (n x n). Returns the tuple of non-unit invariant factors (each divides the next)."""
    n = len(A)
    M = [[(2 if i == j else 0) ^ (A[i][j] & 1) for j in range(n)] for i in range(n)]
    diag = []
    for t in range(n):
        while True:
            best = None
            for i in range(t, n):
                Mi = M[i]
                for j in range(t, n):
                    e = Mi[j]
                    if e and (best is None or e.bit_length() < best[0]):
                        best = (e.bit_length(), i, j)
            assert best is not None, 'xI - A is nonsingular, a zero block cannot occur'
            _, i0, j0 = best
            M[t], M[i0] = M[i0], M[t]
            if j0 != t:
                for row in M:
                    row[t], row[j0] = row[j0], row[t]
            p = M[t][t]
            clean = True
            for i in range(t + 1, n):
                if M[i][t]:
                    q, r = pdivmod(M[i][t], p)
                    Mi, Mt = M[i], M[t]
                    for j in range(t, n):
                        if Mt[j]:
                            Mi[j] ^= pmul(q, Mt[j])
                    if r:
                        clean = False
            for j in range(t + 1, n):
                if M[t][j]:
                    q, r = pdivmod(M[t][j], p)
                    for i in range(t, n):
                        if M[i][t]:
                            M[i][j] ^= pmul(q, M[i][t])
                    if r:
                        clean = False
            if not clean:
                continue
            bad = None
            for i in range(t + 1, n):
                for j in range(t + 1, n):
                    if M[i][j] and pdivmod(M[i][j], p)[1]:
                        bad = i
                        break
                if bad is not None:
                    break
            if bad is not None:
                Mt, Mb = M[t], M[bad]
                for j in range(t, n):
                    Mt[j] ^= Mb[j]
                continue
            break
        diag.append(M[t][t])
    inv = tuple(d for d in diag if pdeg(d) >= 1)
    assert sum(pdeg(d) for d in inv) == n
    for a, b in zip(inv, inv[1:]):
        assert pdivmod(b, a)[1] == 0
    return inv


# ------------------------------------------------------------------ GF(2) linear algebra on row bitmasks
def rows_of(A):
    n = len(A)
    return [sum(A[u][v] << v for v in range(n)) for u in range(n)]


def rank2(rows, ncols):
    basis = {}
    for r in rows:
        x = r
        while x:
            h = x.bit_length() - 1
            if h in basis:
                x ^= basis[h]
            else:
                basis[h] = x
                break
    return len(basis)


def matvec(rows, x):
    y = 0
    for u, r in enumerate(rows):
        if bin(r & x).count('1') & 1:
            y |= 1 << u
    return y


def dynamics(A):
    """rank, nilpotent depth (stabilisation index of the image chain), sorted cycle lengths on the periodic part"""
    n = len(A)
    rows = rows_of(A)
    f = [matvec(rows, x) for x in range(1 << n)]
    img = set(range(1 << n))
    depth = 0
    while True:
        nxt = {f[x] for x in img}
        if nxt == img:
            break
        img = nxt
        depth += 1
    seen = set()
    cyc = []
    for x in img:
        if x in seen:
            continue
        L, y = 0, x
        while True:
            seen.add(y)
            y = f[y]
            L += 1
            if y == x:
                break
        cyc.append(L)
    return rank2(rows, n), depth, tuple(sorted(cyc))


def charpoly(A):
    p = 1
    for d in invariant_factors(A):
        p = pmul(p, d)
    return p


def is_primitive(p):
    """irreducible of degree d with x of order 2^d - 1"""
    d = pdeg(p)
    # irreducible: no factor of degree <= d/2
    for q in range(2, 1 << (d // 2 + 1)):
        if pdeg(q) >= 1 and pdivmod(p, q)[1] == 0:
            return False
    order = (1 << d) - 1
    for r in range(1, order + 1):
        if order % r == 0:
            # x^r mod p
            acc, base, e = 1, 2, r
            while e:
                if e & 1:
                    acc = pdivmod(pmul(acc, base), p)[1]
                base = pdivmod(pmul(base, base), p)[1]
                e >>= 1
            if acc == 1:
                return r == order
    return False


# ------------------------------------------------------------------ tournaments: tilings and an own canonical form
def tilings(n):
    """path n-1 -> ... -> 0 plus an arbitrary orientation of the other pairs; every class appears (Redei)"""
    pairs = [(i, j) for i in range(n) for j in range(i + 2, n)]
    for bits in range(1 << len(pairs)):
        A = [[0] * n for _ in range(n)]
        for i in range(n - 1):
            A[i + 1][i] = 1
        for k, (i, j) in enumerate(pairs):
            if bits >> k & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
        yield A


def canon(A):
    n = len(A)
    col = [sum(A[v]) for v in range(n)]
    while True:
        sig = [(col[v], tuple(sorted(col[w] for w in range(n) if A[v][w])),
                tuple(sorted(col[w] for w in range(n) if A[w][v]))) for v in range(n)]
        keys = sorted(set(sig))
        new = [keys.index(s) for s in sig]
        if len(set(new)) == len(set(col)):
            col = new
            break
        col = new
    cells = defaultdict(list)
    for v in range(n):
        cells[col[v]].append(v)
    order = [cells[c] for c in sorted(cells)]
    best = None
    for choice in itertools.product(*[itertools.permutations(c) for c in order]):
        perm = [v for part in choice for v in part]
        key = tuple(A[perm[i]][perm[j]] for i in range(n) for j in range(n))
        if best is None or key < best:
            best = key
    return best


def transpose(A):
    n = len(A)
    return [[A[j][i] for j in range(n)] for i in range(n)]


def key_to_matrix(key, n):
    return [list(key[i * n:(i + 1) * n]) for i in range(n)]


# ------------------------------------------------------------------ counting
def similarity_class_counts(N):
    """coefficients of prod_{i>=1} 1/(1 - 2 x^i) up to x^N (all n x n matrices over GF(2); the series
    prod (1 - x^i)/(1 - 2x^i) counts only GL(n,2) classes), and an independent recount from irreducibles"""
    c = [1] + [0] * N
    for i in range(1, N + 1):
        # multiply by 1/(1 - 2x^i)
        for k in range(i, N + 1):
            c[k] += 2 * c[k - i]
    # second method: each monic irreducible f gets a partition; GF prod_f P(x^deg f)
    def mobius(m):
        r, k, mm = 1, 2, m
        while k * k <= mm:
            if mm % k == 0:
                mm //= k
                if mm % k == 0:
                    return 0
                r = -r
            k += 1
        if mm > 1:
            r = -r
        return r
    Nirr = [0] + [sum(mobius(d) * 2 ** (m // d) for d in range(1, m + 1) if m % d == 0) // m for m in range(1, N + 1)]
    part = [0] * (N + 1)
    part[0] = 1
    for k in range(1, N + 1):
        for s in range(k, N + 1):
            part[s] += part[s - k]
    g = [1] + [0] * N
    for d in range(1, N + 1):
        for _ in range(Nirr[d]):
            h = [0] * (N + 1)
            for s in range(N + 1):
                if g[s]:
                    for t in range(0, N // d + 1):
                        if s + d * t > N:
                            break
                        h[s + d * t] += g[s] * part[t]
            g = h
    assert g == c, (g, c)
    return c, part


def a000568(n):
    """Davis: sum over partitions of n into odd parts of 2^t / z_lambda"""
    from fractions import Fraction
    total = Fraction(0)

    def parts(m, maxp):
        if m == 0:
            yield []
            return
        for p in range(min(m, maxp), 0, -1):
            if p % 2 == 1:
                for rest in parts(m - p, p):
                    yield [p] + rest
    for lam in parts(n, n):
        z = 1
        for p, mlt in Counter(lam).items():
            z *= p ** mlt * math.factorial(mlt)
        t = sum((p - 1) // 2 for p in lam) + sum(math.gcd(lam[i], lam[j]) for i in range(len(lam))
                                                  for j in range(i + 1, len(lam)))
        total += Fraction(2 ** t, z)
    assert total.denominator == 1
    return int(total)


# ------------------------------------------------------------------ main
def main():
    t0 = time.time()
    random.seed(11012026)
    print('G0  sanity of the Smith-form similarity invariant')
    cser, part = similarity_class_counts(12)
    for n in (1, 2, 3, 4):
        types = set()
        for bits in range(1 << (n * n)):
            A = [[bits >> (i * n + j) & 1 for j in range(n)] for i in range(n)]
            types.add(invariant_factors(A))
        report(f'n={n}: brute-force similarity classes of all {1 << (n * n)} matrices = series coefficient',
               len(types) == cser[n], f'{len(types)} vs {cser[n]}')
    ok = True
    for _ in range(300):
        n = random.randint(3, 8)
        A = [[random.randint(0, 1) for _ in range(n)] for _ in range(n)]
        while True:
            P = [[random.randint(0, 1) for _ in range(n)] for _ in range(n)]
            if rank2(rows_of(P), n) == n:
                break
        # P A P^{-1}: compute P^{-1} by Gauss-Jordan
        aug = [P[i][:] + [1 if i == j else 0 for j in range(n)] for i in range(n)]
        for c in range(n):
            piv = next(r for r in range(c, n) if aug[r][c])
            aug[c], aug[piv] = aug[piv], aug[c]
            for r in range(n):
                if r != c and aug[r][c]:
                    aug[r] = [x ^ y for x, y in zip(aug[r], aug[c])]
        Pinv = [row[n:] for row in aug]
        mm = lambda X, Y: [[sum(X[i][k] & Y[k][j] for k in range(n)) & 1 for j in range(n)] for i in range(n)]
        B = mm(mm(P, A), Pinv)
        if invariant_factors(B) != invariant_factors(A) or invariant_factors(transpose(A)) != invariant_factors(A):
            ok = False
    report('invariant factors unchanged under random conjugation and under transposition (300 random A)', ok)

    print('G1  tournament classes, merged nodes, similarity types (n = 3..7)')
    claimed = {3: (2, 2, 2), 4: (4, 3, 3), 5: (12, 10, 10), 6: (56, 34, 24), 7: (456, 272, 51)}
    store = {}
    for n in (3, 4, 5, 6, 7):
        cls = {}
        type_of_tiling = Counter()
        for A in tilings(n):
            k = canon(A)
            if k not in cls:
                cls[k] = A
        merged = {}
        for k, A in cls.items():
            merged[k] = min(k, canon(transpose(A)))
        nodes = set(merged.values())
        selfconv = sum(1 for k in cls if merged[k] == k and canon(transpose(cls[k])) == k)
        types = defaultdict(set)
        for k, A in cls.items():
            types[invariant_factors(A)].add(merged[k])
        sizes = sorted(len(v) for v in types.values() if len(v) > 1)
        got = (len(cls), len(nodes), len(types))
        report(f'n={n}: classes, merged nodes, similarity types = {claimed[n]}', got == claimed[n],
               f'got {got}; self-converse classes {selfconv}; types holding >1 node: sizes {sizes}')
        store[n] = (cls, merged, types)
    # n = 6 collisions
    cls, merged, types = store[6]
    print('    n = 6 types holding more than one merged node (score sequences; invariant factors):')
    tt6 = None
    for inv, nds in sorted(types.items(), key=lambda kv: -len(kv[1])):
        if len(nds) > 1:
            desc = [tuple(sorted(sum(r) for r in key_to_matrix(k, 6))) for k in sorted(nds)]
            print(f'      {[pstr(d) for d in inv]}: {desc}')
            if any(d == (0, 1, 2, 3, 4, 5) for d in desc):
                tt6 = (inv, desc)
    ok = tt6 is not None and tt6[0] == (1 << 6,) and (2, 2, 2, 3, 3, 3) in tt6[1] and len(tt6[1]) == 2
    report('TT6 shares its type with exactly one other merged node, scores (2,2,2,3,3,3); type = single invariant '
           'factor x^6 (nilpotent of index 6)', ok, f'{tt6}')

    print('G2  counting (exact similarity-class numbers vs merged nodes)')
    for n in range(3, 13):
        an = a000568(n)
        crude = part[n] * 2 ** n
        lower = math.ceil(an / 2)
        print(f'    n={n:2d}: A000568={an:>14d}  merged >= {lower:>13d}  exact similarity classes c(n)={cser[n]:>6d}'
              f'  crude p(n)2^n={crude:>8d}  counting settles: exact {"YES" if cser[n] < lower else "no "}'
              f' crude {"YES" if crude < lower else "no"}')
    report('A000568(9) = 191536, so A000568(9)/2 = 95768; p(9) 2^9 = 15360', a000568(9) == 191536 and part[9] * 512 == 15360)
    report('crude bound fails at n = 8 (p(8) 2^8 = 5632 >= 3440)', part[8] * 256 == 5632)
    report('crude bound valid: c(n) <= p(n) 2^n for n <= 12', all(cser[n] <= part[n] * 2 ** n for n in range(1, 13)))
    report('EXACT count c(8) = 746 < 3440 <= merged nodes at n = 8, so n = 8 is NOT open (fails by counting); '
           'c(7) = 350 >= 272, so n = 7 needed the computation (S4 was a genuine test)',
           cser[8] == 746 and cser[8] < 3440 and cser[7] == 350 and cser[7] >= 272,
           f'c(6)={cser[6]}, c(7)={cser[7]}, c(8)={cser[8]}')
    # self-converse 8-tournaments are not needed: merged(8) >= 6880/2 = 3440 > 746

    print('G3  circulant rank formula and the all-ones kernel')
    ok, cnt = True, 0
    for q in (3, 5, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27):
        pairs = [(x, q - x) for x in range(1, (q + 1) // 2)]
        for choice in itertools.product((0, 1), repeat=len(pairs)):
            S = [p[c] for p, c in zip(pairs, choice)]
            A = [[1 if (v - u) % q in S else 0 for v in range(q)] for u in range(q)]
            a = sum(1 << s for s in S)
            g = pgcd((1 << q) | 1, a)
            if rank2(rows_of(A), q) != q - pdeg(g):
                ok = False
            if q % 4 == 1 and matvec(rows_of(A), (1 << q) - 1) != 0:
                ok = False
            cnt += 1
    report('rank = q - deg gcd(a, x^q - 1) on ALL connection sets (both halves), odd q <= 27; ones in kernel at q=1 mod 4',
           ok, f'{cnt} connection sets')
    ok_even = True
    for _ in range(2000):
        q = random.choice((2, 4, 6, 8, 10, 12, 14, 16))
        a = random.randrange(1, 1 << q)
        A = [[(a >> ((v - u) % q)) & 1 for v in range(q)] for u in range(q)]
        if rank2(rows_of(A), q) != q - pdeg(pgcd((1 << q) | 1, a)):
            ok_even = False
    report('the circulant rank formula needs no separability: it also holds for 2000 random 0/1 circulants of EVEN '
           'order q <= 16 (x^q - 1 inseparable over GF(2))', ok_even)

    print('G4  Paley tournaments')
    def isprime(m):
        return m > 1 and all(m % d for d in range(2, int(m ** 0.5) + 1))
    ok, rows4 = True, []
    for q in [p for p in range(3, 400) if isprime(p) and p % 4 == 3]:
        QR = {(x * x) % q for x in range(1, q)}
        A = [[1 if (v - u) % q in QR else 0 for v in range(q)] for u in range(q)]
        rk = rank2(rows_of(A), q)
        want = q if q % 8 == 3 else (q + 1) // 2
        ok &= rk == want
        rows4.append((q, rk))
    report('rank = q (q = 3 mod 8), (q+1)/2 (q = 7 mod 8) for every prime q = 3 mod 4 below 400', ok,
           f'{len(rows4)} primes')

    def image_code(A):
        n = len(A)
        cols = [sum(A[u][v] << u for u in range(n)) for v in range(n)]
        basis = []
        for c in cols:
            x = c
            for b in basis:
                x = min(x, x ^ b)
            if x:
                basis.append(x)
        return basis

    def min_wt(basis):
        best = None
        k = len(basis)
        # Gray-code enumeration
        w = 0
        for i in range(1, 1 << k):
            j = (i & -i).bit_length() - 1
            w ^= basis[j]
            c = bin(w).count('1')
            if best is None or c < best:
                best = c
        return best

    for q, want in ((7, (4, 3)), (23, (12, 7))):
        QR = {(x * x) % q for x in range(1, q)}
        A = [[1 if (v - u) % q in QR else 0 for v in range(q)] for u in range(q)]
        B = image_code(A)
        d = min_wt(B)
        report(f'P_{q}: image of x -> Ax has dimension {want[0]} and minimum distance {want[1]}', (len(B), d) == want,
               f'got [{q},{len(B)},{d}]')
        a = sum(1 << r for r in QR)
        b = sum(1 << r for r in range(1, q) if r not in QR)
        ga, gb = pgcd((1 << q) | 1, a), pgcd((1 << q) | 1, b)
        # which generator generates the image (column space)? test: every column polynomial divisible by g
        cols = [sum(A[u][v] << u for u in range(q)) for v in range(q)]
        div_a = all(pdivmod(c, ga)[1] == 0 for c in cols)
        div_b = all(pdivmod(c, gb)[1] == 0 for c in cols)
        print(f'    P_{q}: gcd(a, x^q-1) = {pstr(ga)}; gcd(b, x^q-1) = {pstr(gb)} (b = non-residue polynomial); '
              f'image (column space) generated by gcd(a,.): {div_a}, by gcd(b,.): {div_b}')
        report(f'P_{q}: the image is the code generated by gcd(b, x^q-1) (non-residue side), not by gcd(a, x^q-1)',
               div_b and not div_a)

    # doubly regular identities, including the non-circulant Paley tournament on GF(27)
    def mm2(A, B):
        n = len(A)
        return [[sum(A[i][k] & B[k][j] for k in range(n)) & 1 for j in range(n)] for i in range(n)]
    ok = True
    for q in (3, 7, 11, 19, 23, 31, 43, 47):
        QR = {(x * x) % q for x in range(1, q)}
        A = [[1 if (v - u) % q in QR else 0 for v in range(q)] for u in range(q)]
        I = [[int(i == j) for j in range(q)] for i in range(q)]
        if q % 8 == 3:
            ok &= mm2(A, transpose(A)) == I
        else:
            ok &= mm2(A, A) == A
    # GF(27) = GF(3)[t]/(t^3 - t - 1)
    els = [(a, b, c) for a in range(3) for b in range(3) for c in range(3)]
    def gmul(x, y):
        p = [0] * 5
        for i in range(3):
            for j in range(3):
                p[i + j] = (p[i + j] + x[i] * y[j]) % 3
        for k in (4, 3):   # t^3 = t + 1
            c = p[k]
            p[k] = 0
            p[k - 3] = (p[k - 3] + c) % 3
            p[k - 2] = (p[k - 2] + c) % 3
        return tuple(p[:3])
    sq = {gmul(x, x) for x in els if x != (0, 0, 0)}
    A27 = [[1 if tuple((y[i] - x[i]) % 3 for i in range(3)) in sq else 0 for y in els] for x in els]
    I27 = [[int(i == j) for j in range(27)] for i in range(27)]
    tour = all(A27[i][j] + A27[j][i] == (i != j) for i in range(27) for j in range(27))
    ok27 = tour and mm2(A27, transpose(A27)) == I27 and rank2(rows_of(A27), 27) == 27
    report('doubly regular identities: A A^T = I over GF(2) for q = 3 mod 8 (so A^-1 = A^T), A^2 = A for q = 7 mod 8 '
           '(the automaton is idempotent); also the non-circulant Paley tournament of GF(27) is reversible', ok and ok27)

    print('G5  the n = 5 table (classes read from merged_metagraph_n5.json as data)')
    M = json.load(open(os.path.join(HERE, '..', '..', 'groups_of_five_primegame_20261031', 'merged_metagraph_n5.json')))
    table = {0: (4, 5, (1,)), 1: (4, 2, (1, 1, 3, 3)), 2: (4, 1, (1, 1, 2, 4, 4, 4)), 3: (3, 1, (1, 1, 3, 3)),
             4: (5, 0, (1, 1, 2, 7, 7, 14)), 6: (5, 0, (1, 31)), 8: (5, 0, (1, 1, 5, 5, 5, 5, 5, 5)),
             9: (5, 0, (1, 3, 7, 21)), 10: (4, 1, (1, 1, 1, 1, 2, 2, 4, 4)), 11: (4, 1, (1, 15))}
    ok = True
    mats = {}
    for c in M['classes']:
        A = [[int(c['sig'][i * 5 + j]) for j in range(5)] for i in range(5)]
        mats[c['ci']] = A
        dyn = dynamics(A)
        cp = charpoly(A)
        inv = invariant_factors(A)
        print(f'    #{c["ci"]:2d} {c["type"]} scores {c["scores"]} H {c["H"]}: rank {dyn[0]}, depth {dyn[1]}, cycles '
              f'{dyn[2]}, invariant factors {[pstr(d) for d in inv]}')
        if c['ci'] in table and dyn != table[c['ci']]:
            ok = False
    for a, b in ((1, 7), (2, 5)):
        ok &= canon(transpose(mats[a])) == canon(mats[b]) and dynamics(mats[a]) == dynamics(mats[b])
    report('THM-4857 n = 5 table (rank, nilpotent depth, cycles) reproduced; #1/#7 and #2/#5 are converse pairs', ok)
    mixed = [c['ci'] for c in M['classes'] if c['type'] == 'MX']
    rev = [ci for ci in mixed if dynamics(mats[ci])[0] == 5]
    allrev = [c['ci'] for c in M['classes'] if dynamics(mats[c['ci']])[0] == 5]
    report('reversible (rank 5) classes are exactly four of the five mixed nodes', len(mixed) == 5 and len(rev) == 4
           and sorted(allrev) == sorted(rev), f'mixed {mixed}, reversible {allrev}')
    cp6 = charpoly(mats[6])
    report('#6: characteristic polynomial primitive of degree 5 (single 31-cycle)', is_primitive(cp6) and
           dynamics(mats[6])[2] == (1, 31), pstr(cp6))
    print()
    print(f'elapsed {time.time() - t0:.1f} s')
    print('ALL AUDIT CHECKS AGREE' if not FAILS else f'MISMATCHES: {FAILS}')
    return 0 if not FAILS else 1


if __name__ == '__main__':
    sys.exit(main())
