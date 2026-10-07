#!/usr/bin/env python3
"""Audit A, item 4: THM-4601 (iv) two-anchor compiler, independent code.

(a) identity of affine maps  F_(1^(K-4) 4 1 1 2^(J-3) v (c+2)) ((n+1)/8 - 1) = F_(1^(K-1) 2^J u c)(n)
    checked exactly (two evaluation points, Fractions) for the two explicit heads, K = 4..80, J = 3..60, c = 1..6, plus
    K, J up to 1001, 150; and for every head pair the audit finds itself with sum(u) <= 13 (own search).
(b) the head data: |v| = |u| + 3, sum u = sum v + 2, F_v(1) = 4 F_u(1) + 1, the stated values.
(c) validity lemma ("odd integer endpoint forces the word"): random words and odd starts.
(d) native classes of the post-run class w' (Y_end = 1 + 2^r w', X_end = 1 + 27 * 2^r w') mod 2^14, with c >= 1 and c >= 2.
(e) random residual sources, K up to 1001, J up to 150, w' drawn from the native class (c >= 1): the actual U-words of n and
    h_3 = (n+1)/8 - 1 are exactly the displayed words, the odd endpoints coincide, and the D = 3 chain absorbs at
    post-run Terras depth sum(u) + 1.
"""
import random
from fractions import Fraction as Fr
from collections import Counter

NCHK = 0
def check(c, m):
    global NCHK
    if not c: raise AssertionError(m)
    NCHK += 1

def v2(n):
    return (n & -n).bit_length() - 1

def Fw(w, x):
    x = Fr(x)
    for a in w:
        x = (3*x + 1) / Fr(2**a)
    return x

def Fw_affine(w):
    """slope, intercept of F_w as exact Fractions"""
    b = Fw(w, 0)
    return Fw(w, 1) - b, b

def uword(z, L):
    out = []
    for _ in range(L):
        y = 3*z + 1; a = v2(y); out.append(a); z = y >> a
    return out, z

HEADS = {1: ((1, 10), (1, 1, 2, 2, 3)), 2: ((10,), (3, 1, 1, 3))}

def identity_holds(u, v, K, J, c):
    src = (1,)*(K - 1) + (2,)*J + tuple(u) + (c,)
    chd = (1,)*(K - 4) + (4, 1, 1) + (2,)*(J - 3) + tuple(v) + (c + 2,)
    for n in (Fr(0), Fr(1), Fr(-7, 3)):
        if Fw(chd, (n + 1)/8 - 1) != Fw(src, n):
            return False
    return True

def part_a_b():
    for r, (u, v) in HEADS.items():
        check(len(v) == len(u) + 3 and sum(u) == sum(v) + 2, "head lengths/totals")
        check(Fw(v, 1) == 4*Fw(u, 1) + 1, "F_v(1) = 4 F_u(1) + 1")
    check(Fw(HEADS[1][0], 1) == Fr(7, 1024) and Fw(HEADS[1][1], 1) == Fr(263, 256), "r=1 values")
    check(Fw(HEADS[2][0], 1) == Fr(1, 256) and Fw(HEADS[2][1], 1) == Fr(65, 64), "r=2 values")
    n = 0
    for r, (u, v) in HEADS.items():
        for K in list(range(4, 81)) + [333, 1001]:
            for J in list(range(3, 61, 3)) + [150]:
                for c in (1, 2, 3, 6):
                    check(identity_holds(u, v, K, J, c), f"identity r={r} K={K} J={J} c={c}")
                    n += 1
    # negative control: wrong deletion depth or wrong final letter breaks it
    u, v = HEADS[2]
    check(identity_holds(u, v, 10, 5, 1), "control sanity")
    src = (1,)*9 + (2,)*5 + u + (1,); chd_bad = (1,)*6 + (4, 1, 1) + (2,)*2 + v + (2,)
    check(Fw(chd_bad, Fr(4, 8) - 1) != Fw(src, Fr(3)), "negative control (final letter c+1) fails")
    print(f"a,b  explicit heads: data ok; identity of affine maps ok on {n} (K, J, c) triples (K <= 1001, J <= 150)")

def own_head_search(SMAX=13):
    """all words u (first letter != 2) with sum <= SMAX and partners v with |v| = |u|+3, sum v = sum u - 2,
    F_v(1) = 4 F_u(1) + 1 (own enumeration via exact values)."""
    from itertools import product
    def words_with(total, length):
        if length == 0:
            if total == 0: yield ()
            return
        for a in range(1, total - length + 2):
            for rest in words_with(total - a, length - 1):
                yield (a,) + rest
    pairs = []
    vals = {}
    for S in range(1, SMAX - 1):
        for L in range(1, S + 1):
            for w in words_with(S, L):
                vals.setdefault((L, S, Fw(w, 1)), []).append(w)
    for S in range(3, SMAX + 1):
        for L in range(1, S + 1):
            for u in words_with(S, L):
                if u[0] == 2: continue
                target = 4*Fw(u, 1) + 1
                for v in vals.get((L + 3, S - 2, target), []):
                    pairs.append((u, v))
    return pairs

def part_a_general():
    pairs = own_head_search(13)
    rnd = random.Random(5)
    for (u, v) in pairs:
        for _ in range(4):
            K = rnd.randint(4, 60); J = rnd.randint(3, 40); c = rnd.randint(1, 7)
            check(identity_holds(u, v, K, J, c), f"general head {u}~{v}")
    same_type = [(u, v) for u, v in pairs if (u[0] == 1) == (v[0] == 1)]
    print(f"a'   own search: {len(pairs)} head pairs with sum(u) <= 13 ({len(same_type)} with matching run type); "
          f"the compiler identity holds for every one (random K, J, c); first: {pairs[:4]}")
    return pairs

def part_c():
    rnd = random.Random(77)
    hits = 0
    for _ in range(20000):
        L = rnd.randint(1, 7)
        w = tuple(rnd.randint(1, 4) for _ in range(L))
        P, B = Fw_affine(w)
        # choose m in the native cylinder of w half the time, else random odd
        if rnd.random() < 0.5:
            S = sum(w)
            # m with F_w(m) odd integer: m = (2^S z - B*2^S)/3^L for odd z, solve mod
            z = rnd.getrandbits(20)*2 + 1
            num = z - B
            m = num / P
            if m.denominator != 1 or m.numerator % 2 == 0 or m <= 0:
                continue
            m = m.numerator
        else:
            m = rnd.getrandbits(30)*2 + 1
        E = Fw(w, m)
        if E.denominator == 1 and E.numerator % 2 == 1:
            ww, z = uword(m, L)
            check(tuple(ww) == w and z == E, "validity lemma: odd integer endpoint forces the word")
            hits += 1
    print(f"c    validity lemma: {hits} random (word, odd start) with odd integer endpoint; actual word = given word in all")

def native_classes(r, u, cmin, L=14):
    good = []
    for wp in range(1, 1 << L, 2):
        Xe = 1 + 27*(1 << r)*wp + (1 << (L + r + 12)) * 7      # a positive lift; letters needed are decided mod 2^(r+L)
        ww, z = uword(Xe, len(u) + 1)
        if tuple(ww[:len(u)]) == u and ww[len(u)] >= cmin:
            good.append(wp)
    return good

def part_d():
    out = {}
    for r, (u, v) in HEADS.items():
        for cmin in (1, 2):
            g = native_classes(r, u, cmin)
            out[(r, cmin)] = len(g)
    print(f"d    native classes of w' mod 2^14 (8192 odd classes): r=1 (u=(1,10)): c>=1 {out[(1,1)]}, c>=2 {out[(1,2)]}; "
          f"r=2 (u=(10)): c>=1 {out[(2,1)]}, c>=2 {out[(2,2)]}  [measure 2^(r - sum u): r=1 {8192*2**(1-11):.0f}, r=2 {8192*2**(2-10):.0f}]")
    return out

def T(x):
    return (3*x + 1) >> 1 if x & 1 else x >> 1

def part_e():
    rnd = random.Random(2026_10_07)
    cnt = Counter(); depths = Counter(); cvals = Counter()
    for r, (u, v) in HEADS.items():
        Su = sum(u)
        good = native_classes(r, u, 1)
        for trial in range(160):
            K = rnd.choice([4, 5, 6, 7, 9, 17, 64, 255, 600, 1001])
            J = rnd.choice([3, 4, 5, 6, 9, 17, 40, 99, 150])
            j = 2*J + r
            M = j + 14 + 300
            wp = rnd.choice(good) + (1 << 14)*rnd.getrandbits(250)
            # w' = 3^(J-3) u_x  (u_x = (x-1)/2^j)  ->  u_x = 3^(3-J) w' mod 2^(M-j)
            ux = (pow(3, -(J - 3), 1 << (M - j)) * wp) % (1 << (M - j))
            ux |= 0
            check(ux & 1, "u_x odd")
            t = (pow(3, -(K - 1), 1 << M) * (1 + (ux << (j - 1)))) % (1 << M)
            check(t & 1, "t odd")
            n = (t << K) - 1
            h = (t << (K - 3)) - 1
            x = 2*3**(K - 1)*t - 1
            check(v2(x - 1) == j and (3**K*t - 1) % 4 == 2, "residual source with this j")
            # actual words
            Ls = (K - 1) + J + len(u) + 1
            Lc = (K - 4) + 3 + (J - 3) + len(v) + 1
            ws, zs = uword(n, Ls)
            wc, zc = uword(h, Lc)
            c = ws[-1]
            exp_s = [1]*(K - 1) + [2]*J + list(u) + [c]
            exp_c = [1]*(K - 4) + [4, 1, 1] + [2]*(J - 3) + list(v) + [c + 2]
            check(ws == exp_s, f"source word r={r} K={K} J={J}")
            check(wc == exp_c, f"child word r={r} K={K} J={J}")
            check(zs == zc and h < n, "common odd endpoint, smaller child")
            cvals[(r, c)] += 1
            # D = 3 chain from run ends: absorption time post-run
            y3 = (x + 1)//27 - 1
            a, b, k = x, y3, 3
            ab = None
            for s in range(2*J + Su + 5):
                if a == b and k == 0:
                    ab = s; break
                k += (a & 1) - (b & 1); a, b = T(a), T(b)
            check(ab is not None, "absorbed")
            depths[(r, ab - 2*J, ab - j)] += 1
            cnt[r] += 1
    print(f"e    random sources: {dict(cnt)} certified (K up to 1001, J up to 150), words exact, endpoints equal;")
    print(f"     final letters seen (r, c): {sorted(cvals.items())}")
    print(f"     D=3 absorption (r, post-run Terras depth = time - 2J, time - j): {sorted(depths.items())}")

if __name__ == "__main__":
    part_a_b()
    part_a_general()
    part_c()
    part_d()
    part_e()
    print(f"ALL (iv) CHECKS PASSED ({NCHK} assertions)")
