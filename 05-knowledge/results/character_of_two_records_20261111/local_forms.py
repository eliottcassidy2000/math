#!/usr/bin/env python3
"""Local statements around conjecture (bS) (H mod 4 on skew Hadamard switching classes), each tested; section 5.3
and section 10 of the note.  Usage: local_forms.py EXACT   (EXACT = compiled hampaths_exact.c)

 (1) det S = 4H - 3 (mod 16)?   (orders 2-10, random; exhaustive for orders 2-5)
 (2) one orthogonal row: if row v of S is orthogonal to all other rows, does switching at v keep H mod 4?
     (a) random tournaments of order 4 and 8 having such a row; (b) (source) + random regular tournament of
     order 3..15, switched at a random set, then at the source.
 (3) score parities only: <s, g> = 0?  s = score parities of R, g(x) = parity of the number of even-order paths of
     R with endpoint x (exhaustive for orders 3 and 5, random for 7 and 9).
 (4) the proved reduction: for |R| odd, g(x) = parity of the number of paths on |R| - 1 vertices with endpoint x;
     and (H(U) - H(U^v))/2 = sum_{x in N^-(v)} g_R(x) (mod 2).
"""
import sys, random, itertools, subprocess
from fractions import Fraction
from fold import ham

def det_int(M):
    n = len(M); A = [[Fraction(x) for x in r] for r in M]; d = Fraction(1)
    for i in range(n):
        p = next((r for r in range(i, n) if A[r][i] != 0), None)
        if p is None: return 0
        if p != i: A[i], A[p] = A[p], A[i]; d = -d
        d *= A[i][i]
        for r in range(i + 1, n):
            f = A[r][i] / A[i][i]
            if f:
                for c in range(i, n): A[r][c] -= f * A[i][c]
    assert d.denominator == 1; return int(d)
def rand_tournament(n, rng):
    out = [0] * n
    for i in range(n):
        for j in range(i + 1, n):
            if rng.getrandbits(1): out[i] |= 1 << j
            else: out[j] |= 1 << i
    return out
def all_tournaments(n):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    for bits in range(1 << len(pairs)):
        out = [0] * n
        for k, (i, j) in enumerate(pairs):
            if (bits >> k) & 1: out[i] |= 1 << j
            else: out[j] |= 1 << i
        yield out
def skew(out):
    n = len(out); return [[0 if i == j else (1 if (out[i] >> j) & 1 else -1) for j in range(n)] for i in range(n)]
def switch_v(out, v):
    n = len(out); o = out[:]
    for j in range(n):
        if j == v: continue
        if (out[v] >> j) & 1: o[v] &= ~(1 << j); o[j] |= (1 << v)
        else: o[v] |= (1 << j); o[j] &= ~(1 << v)
    return o
def exact_H(exe, tours):
    inp = "".join("%d\n%s\n" % (len(o), "\n".join(str(x) for x in o)) for o in tours)
    res = subprocess.run([exe], input=inp, capture_output=True, text=True, check=True).stdout.split("\n")
    return [int(l.split()[1]) for l in res if l.strip()]
def endpoint_counts(out, n):
    """E[k][x] = number of paths on k vertices ending at x."""
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n): dp[1 << v][v] = 1
    E = [[0] * n for _ in range(n + 1)]
    for S in range(1, 1 << n):
        k = bin(S).count('1'); row = dp[S]
        for v in range(n):
            c = row[v]
            if not c: continue
            E[k][v] += c; t = out[v] & ~S
            while t:
                w = (t & -t).bit_length() - 1; dp[S | (1 << w)][w] += c; t &= t - 1
    return E
def converse(out):
    n = len(out); return [((1 << n) - 1) & ~out[i] & ~(1 << i) for i in range(n)]
def g_vector(out):
    n = len(out); E = endpoint_counts(out, n); S = endpoint_counts(converse(out), n)
    g = [sum(E[k][x] + S[k][x] for k in range(2, n + 1, 2)) % 2 for x in range(n)]
    near = [(E[n - 1][x] + S[n - 1][x]) % 2 for x in range(n)] if n % 2 else None
    return g, near
def rand_regular(n, rng):
    k = (n - 1) // 2; A = [[1 if 0 < ((j - i) % n) <= k else 0 for j in range(n)] for i in range(n)]
    for _ in range(40 * n * n):
        a, b, c = rng.sample(range(n), 3)
        if A[a][b] and A[b][c] and A[c][a]:
            A[a][b] = A[b][c] = A[c][a] = 0; A[b][a] = A[c][b] = A[a][c] = 1
        elif A[b][a] and A[c][b] and A[a][c]:
            A[b][a] = A[c][b] = A[a][c] = 0; A[a][b] = A[b][c] = A[c][a] = 1
    return A

if __name__ == "__main__":
    exe = sys.argv[1]; rng = random.Random(20261111)
    print("(1) det S = 4H - 3 (mod 16) for tournaments of even order")
    for n in (2, 4):
        bad = sum(1 for o in all_tournaments(n) if (det_int(skew(o)) - 4 * ham(o, n) + 3) % 16)
        print("   order %d, all labelled tournaments: %d violations" % (n, bad))
    for n in (6, 8, 10):
        tr = 300 if n <= 8 else 80; bad = sum(1 for _ in range(tr) for o in [rand_tournament(n, rng)] if (det_int(skew(o)) - 4 * ham(o, n) + 3) % 16)
        print("   order %d, %d random tournaments: %d violations" % (n, tr, bad))
    print("(2a) one orthogonal row (row 0 of S orthogonal to all other rows): switching at 0")
    for n, trials in ((4, 4000), (8, 600000)):
        hits = same = nonsh = 0
        for _ in range(trials):
            o = rand_tournament(n, rng); S = skew(o)
            if any(sum(S[0][b] * S[b][c] for b in range(n)) for c in range(1, n)): continue
            hits += 1
            if any(sum(S[i][k] * S[j][k] for k in range(n)) for i in range(n) for j in range(i + 1, n)): nonsh += 1
            if ham(o, n) % 4 == ham(switch_v(o, 0), n) % 4: same += 1
        print("   order %d: %d tournaments with such a row (%d of them not skew Hadamard); H mod 4 unchanged in %d, changed in %d" % (n, hits, nonsh, same, hits - same))
    print("(2b) U = ((source) + regular T') switched at a random set; U against U switched at the source")
    for n in (3, 5, 7, 9, 11, 13, 15):
        tours = []
        for _ in range(40 if n <= 11 else 12):
            A = rand_regular(n, rng); assert all(sum(r) == (n - 1) // 2 for r in A)
            for _ in range(4):
                Y = set(x for x in range(n) if rng.getrandbits(1)); N = n + 1; out = [0] * N
                for i in range(n):
                    for j in range(n):
                        if i != j and (A[i][j] ^ ((i in Y) != (j in Y))): out[i] |= 1 << j
                    if i in Y: out[i] |= 1 << n
                    else: out[n] |= 1 << i
                tours += [out, switch_v(out, n)]
        hs = exact_H(exe, tours); same = sum(1 for i in range(0, len(hs), 2) if hs[i] % 4 == hs[i + 1] % 4)
        print("   regular T' of order %d (N = %d): H mod 4 equal in %d pairs, different in %d" % (n, n + 1, same, len(hs) // 2 - same))
    print("(3) score parities only: <s, g> for tournaments R of odd order")
    for n in (3, 5):
        bad = tot = 0
        for o in all_tournaments(n):
            g, near = g_vector(o); tot += 1; bad += sum(g[x] for x in range(n) if bin(o[x]).count('1') % 2) % 2
        print("   order %d, all %d labelled tournaments: <s, g> = 1 in %d" % (n, tot, bad))
    for n, tr in ((7, 400), (9, 80)):
        bad = 0
        for _ in range(tr):
            o = rand_tournament(n, rng); g, near = g_vector(o); bad += sum(g[x] for x in range(n) if bin(o[x]).count('1') % 2) % 2
        print("   order %d, %d random tournaments: <s, g> = 1 in %d" % (n, tr, bad))
    print("(4) the proved reduction")
    for n in (3, 4, 5, 6, 7, 8, 9):
        badnear = badred = tot = 0
        for _ in range(300 if n <= 7 else 100):
            R = rand_tournament(n, rng); g, near = g_vector(R); tot += 1
            if near is not None and near != g: badnear += 1
            X = rng.getrandbits(n); N = n + 1; U = [R[i] | ((1 << n) if not (X >> i) & 1 else 0) for i in range(n)] + [X]
            d = ((ham(U, N) - ham(switch_v(U, n), N)) // 2) % 2
            if d != sum(g[x] for x in range(n) if not (X >> x) & 1) % 2: badred += 1
        print("   |R| = %d, %d random (R, N^+(v)): reduction formula fails in %d; near-Hamiltonian form of g fails in %s" % (n, tot, badred, badnear if n % 2 else "n/a (|R| even)"))
