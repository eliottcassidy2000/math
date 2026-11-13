"""Standard-library verifier for kind-pasteur-2026-11-13-S26 (THM-5151).

Checks, from literal adjacency rows (outgoing-neighbour bit masks, labels from 0):
  1. the full slot profiles F_j (ordered covers by j nonempty paths) of every claimed exact pair;
  2. Q_X(t) = Q_Z(t) for t = 0..n (so the polynomials agree) and the SCC shape of each side;
  3. Session 310's B is isomorphic to the one-apex tournament O(2,5) = O(0011111): body 1 -> ... -> 7 (i -> j for
     i < j), apex beats the first two (top-ranked) body vertices and loses to the other five;
  4. the one-apex closed form of Q against brute force for every apex word of length <= 7;
  5. the C=1 counts F_(2n-6) for n = 4..9 by one-vertex extension (n = 10 with --full);
  6. the Pi-clock values of the prime pair, their recurrence and their primality (deterministic Miller-Rabin).
Run: python verify_s26.py [--full]
"""
import sys, itertools
from math import comb

def need(cond, msg=''):
    """explicit check that survives python -O"""
    if not cond:
        raise SystemExit('FAIL: %s' % (msg,))

def profile(out):
    n = len(out); N = 1 << n
    end = [[0] * n for _ in range(N)]; Hs = [0] * N
    for S in range(1, N):
        tot = 0
        for v in range(n):
            if not S >> v & 1:
                continue
            P = S & ~(1 << v)
            a = 1 if P == 0 else sum(end[P][w] for w in range(n) if P >> w & 1 and out[w] >> v & 1)
            end[S][v] = a; tot += a
        Hs[S] = tot
    g = Hs[:]; F = [g[N - 1]]
    for _ in range(2, n + 1):
        ng = [0] * N
        for S in range(N):
            s = 0; B = S
            while B:
                A = S & ~B
                if A:
                    s += g[A] * Hs[B]
                B = (B - 1) & S
            ng[S] = s
        g = ng; F.append(g[N - 1])
    return F

def Q(F, t):
    return sum(F[j - 1] * comb(t, j) for j in range(1, len(F) + 1))

def ordinal(X, Y):
    n, m = len(X), len(Y)
    return [X[i] | (((1 << m) - 1) << n) for i in range(n)] + [Y[k] << n for k in range(m)]

def strong(out):
    n = len(out); full = (1 << n) - 1
    def reach(adj):
        seen = fr = 1
        while fr:
            nf = 0
            for v in range(n):
                if fr >> v & 1:
                    nf |= adj[v]
            nf &= ~seen; seen |= nf; fr = nf
        return seen
    inn = [sum(1 << v for v in range(n) if out[v] >> w & 1) for w in range(n)]
    return reach(out) == full and reach(inn) == full

def ham_cycles(out):
    n = len(out); N = 1 << n
    st = {(1, 0): 1}
    layer = {(1, 0): 1}
    for _ in range(n - 1):
        nxt = {}
        for (S, v), a in layer.items():
            for w in range(n):
                if not S >> w & 1 and out[v] >> w & 1:
                    key = (S | 1 << w, w); nxt[key] = nxt.get(key, 0) + a
        layer = nxt
    return sum(a for (S, v), a in layer.items() if S == N - 1 and out[v] & 1)

def one_apex(word):
    """apex = vertex 0; body vertices 1..m transitive (i -> j for i < j); word[i] = '1' iff body vertex i+1 beats apex"""
    m = len(word); n = m + 1; out = [0] * n
    for i in range(1, n):
        for j in range(i + 1, n):
            out[i] |= 1 << j
    for i, c in enumerate(word, 1):
        if c == '1':
            out[i] |= 1
        else:
            out[0] |= 1 << i
    return out

def one_apex_Q_coeffs(word):
    """closed form: Q(t) = t^n + 2 * sum over (a<b, w_a=0, w_b=1) sum_l C(b-a-1, 2l-2) t^(n-2l)"""
    m = len(word); n = m + 1; c = [0] * (n + 1); c[n] = 1
    for a in range(m):
        for b in range(a + 1, m):
            if word[a] == '0' and word[b] == '1':
                g = b - a - 1
                for l in range(1, g // 2 + 2):
                    if 2 * l - 2 <= g:
                        c[n - 2 * l] += 2 * comb(g, 2 * l - 2)
    return c

def decode(code, n):
    out = [0] * n; k = 0
    for i in range(n):
        for j in range(i + 1, n):
            if code >> k & 1:
                out[i] |= 1 << j
            else:
                out[j] |= 1 << i
            k += 1
    return out

def is_prime(n):
    if n < 2:
        return False
    small = (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41)
    for p in small:
        if n % p == 0:
            return n == p
    need(n < 3317044064679887385961981)  # deterministic range for bases <= 41
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2; s += 1
    for a in small:
        x = pow(a, d, n)
        if x in (1, n - 1):
            continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1:
                break
        else:
            return False
    return True

C3 = [2, 4, 1]; TT1 = [0]
PAIRS = [
    # (label, strong X rows, decomposable partner rows, expected H, expected C(X))
    ('order 6, H=15, C=1 one-apex word 01001', [6, 4, 8, 19, 39, 15], ordinal(decode(152, 5), TT1), 15, 1),
    ('order 8, H=95 (Session 310 B) vs A->TT1', [2, 80, 3, 7, 173, 143, 189, 15], ordinal([2, 28, 25, 81, 97, 79, 7], TT1), 95, 1),
    ('order 8, H=45, C3 => Y5', [30, 12, 8, 16, 38, 79, 159, 63], ordinal(C3, decode(739, 5)), 45, 1),
    ('order 10, H=153, C3 => Y7', [126, 60, 56, 48, 32, 64, 158, 319, 639, 255], ordinal(C3, decode(1116481, 7)), 153, 1),
    ('order 10, H=219, C3 => Y7', [254, 508, 504, 240, 224, 192, 128, 256, 633, 255], ordinal(C3, decode(202882, 7)), 219, 1),
    ('order 10, H=317 (prime), Y9 -> TT1', [254, 508, 120, 112, 224, 192, 128, 268, 637, 255], ordinal([378, 380, 505, 496, 224, 448, 384, 259, 16], TT1), 317, 1),
    # pre-registered S3: refutes HYP-9678 (A) - a strong member with five Hamiltonian cycles
    ('order 10, H=333, C3 => Y7, C(X)=5', [1002, 1020, 1017, 1008, 481, 704, 896, 768, 544, 16], ordinal(C3, decode(43104, 7)), 333, 5),
]

def check_pairs():
    for lab, X, Z, H, CX in PAIRS:
        FX, FZ = profile(X), profile(Z)
        need(FX == FZ, lab)
        need(FX[0] == H, lab)
        need(strong(X) and not strong(Z), lab)
        need(all(Q(FX, t) == Q(FZ, t) for t in range(len(X) + 1)), lab)
        c = ham_cycles(X)
        print('  ok  %-44s F = %s  C(X) = %d' % (lab, FX, c))
        need(c == CX, lab)

def iso(A, B):
    n = len(A)
    if sorted(bin(r).count('1') for r in A) != sorted(bin(r).count('1') for r in B):
        return False
    for p in itertools.permutations(range(n)):
        if all((A[i] >> j & 1) == (B[p[i]] >> p[j] & 1) for i in range(n) for j in range(n) if i != j):
            return True
    return False

def check_B_one_apex():
    B = [2, 80, 3, 7, 173, 143, 189, 15]
    need(iso(B, one_apex('0011111')))  # apex beats body vertices 1, 2 (the two top-ranked) and loses to 3..7
    need(not iso(B, one_apex('0000011')))  # that is the converse O(5,2), which has the same Q
    need(profile(one_apex('0000011')) == profile(B))
    print('  ok  Session 310 B is the one-apex tournament O(0011111): apex beats the two top-ranked vertices of TT7 (converse O(0000011))')

def mono_from_F(F):
    n = len(F); pts = [Q(F, t) for t in range(n + 1)]
    from fractions import Fraction
    coef = [Fraction(0)] * (n + 1)
    for i in range(n + 1):
        num = [Fraction(1)]; den = 1
        for j in range(n + 1):
            if j == i:
                continue
            q = [Fraction(0)] * (len(num) + 1)
            for k, c in enumerate(num):
                q[k + 1] += c; q[k] -= j * c
            num = q; den *= (i - j)
        for k, c in enumerate(num):
            coef[k] += pts[i] * c / den
    need(all(c.denominator == 1 for c in coef))
    return [int(c) for c in coef]

def check_one_apex_formula():
    cnt = 0
    for m in range(1, 8):
        for bits in range(1 << m):
            w = ''.join('1' if bits >> i & 1 else '0' for i in range(m))
            X = one_apex(w)
            need(mono_from_F(profile(X)) == one_apex_Q_coeffs(w), w)
            if w[0] == '0' and w[-1] == '1':
                need(strong(X) and ham_cycles(X) == 1, w)
            cnt += 1
    print('  ok  one-apex closed form of Q on all %d apex words of length <= 7; strong ones have C = 1' % cnt)

def c1_extend(classes, m):
    """classes: list of row lists of order m (C=1 strong); return canonical C=1 classes of order m+1"""
    seen = {}
    for S in classes:
        N = 1 << m
        P = [[0] * m for _ in range(m)]
        for y in range(m):
            dp = {(1 << y, y): 1}; layer = dict(dp)
            for _ in range(m - 1):
                nxt = {}
                for (U, v), a in layer.items():
                    for w in range(m):
                        if not U >> w & 1 and S[v] >> w & 1:
                            key = (U | 1 << w, w); nxt[key] = nxt.get(key, 0) + a
                layer = nxt
            for (U, x), a in layer.items():
                if U == N - 1:
                    P[y][x] = a
        for O in range(1, N - 1):
            C = 0
            for x in range(m):
                if not O >> x & 1:
                    for y in range(m):
                        if O >> y & 1:
                            C += P[y][x]
            if C != 1:
                continue
            out = [S[i] | (0 if O >> i & 1 else 1 << m) for i in range(m)] + [O]
            n = m + 1
            # unique Hamiltonian cycle starting at 0
            cyc = None
            def dfs(path, used):
                nonlocal cyc
                if cyc is not None:
                    return
                if len(path) == n:
                    if out[path[-1]] & 1:
                        cyc = list(path)
                    return
                v = path[-1]
                for w in range(n):
                    if not used >> w & 1 and out[v] >> w & 1:
                        path.append(w); dfs(path, used | 1 << w); path.pop()
            dfs([0], 1)
            best = None
            for s in range(n):
                lab = [cyc[(s + i) % n] for i in range(n)]
                word = tuple(out[lab[i]] >> lab[j] & 1 for i in range(n) for j in range(i + 1, n))
                if best is None or word < best[0]:
                    best = (word, lab)
            lab = best[1]; pos = {v: i for i, v in enumerate(lab)}
            rows = [0] * n
            for v in range(n):
                for w in range(n):
                    if out[v] >> w & 1:
                        rows[pos[v]] |= 1 << pos[w]
            seen[best[0]] = rows
    return list(seen.values())

def check_c1_counts(full):
    fib = [0, 1]
    while len(fib) < 30:
        fib.append(fib[-1] + fib[-2])
    cls = [[2, 4, 1]]; m = 3
    top = 10 if full else 9
    while m < top:
        cls = c1_extend(cls, m); m += 1
        need(len(cls) == fib[2 * m - 6], (m, len(cls)))
        print('  ok  C=1 strong classes of order %d: %d = F_%d' % (m, len(cls), 2 * m - 6))

def check_clock():
    Y9 = [378, 380, 505, 496, 224, 448, 384, 259, 16]
    FY = profile(Y9); q2 = Q(FY, 2)
    need(FY[0] == 317 and q2 == 9280)
    def run(P3, h, r):
        return [60 * P3 * q2 ** r * 2 ** (r - k) * 3 ** k + h * 317 ** r for k in range(r + 1)]
    R = run(27, 1, 2)          # M = TT3: P_M(3) = 27, H(M) = 1
    need(R == [558047332489, 837070948489, 1255606372489])
    need(all((3 * R[k] - 317 ** 2) // 2 == R[k + 1] and (3 * R[k] - 317 ** 2) % 2 == 0 for k in range(2)))
    need(is_prime(R[0]) and not is_prime(R[1]))
    need(R[2] % 7 == 0 and is_prime(R[2] // 7) and R[2] // 7 == 179372338927)
    R1 = run(3, 1, 2)          # M = TT1
    need(R1[2] == 139511908489 and is_prime(R1[2]))
    R3 = run(3, 1, 3)
    need(R3[3] == 3884008766575013 and is_prime(R3[3]))
    # contrast with Session 310's pair: every Pi value with B^r is divisible by 5
    for P3, h in ((3, 1), (27, 1), (33, 3)):
        for r in (1, 2, 3):
            for k in range(r + 1):
                need((60 * P3 * 1276 ** r * 2 ** (r - k) * 3 ** k + h * 95 ** r) % 5 == 0)
    print('  ok  Pi clock with the 317 pair: prime 558047332489 -> 837070948489 -> 7 x 179372338927 (M = TT3, r = 2)')
    print('  ok  M = TT1: primes 139511908489 (r=2,k=2), 3884008766575013 (r=3,k=3); B^r clock values all divisible by 5')

if __name__ == '__main__':
    full = '--full' in sys.argv
    print('exact pairs:'); check_pairs()
    print('Session 310 B:'); check_B_one_apex()
    print('one-apex formula:'); check_one_apex_formula()
    print('C=1 enumeration:'); check_c1_counts(full)
    print('Pi clock:'); check_clock()
    print('ALL CHECKS PASSED')
