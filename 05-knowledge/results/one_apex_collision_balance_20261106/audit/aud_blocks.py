"""Independent audit check of THM-4890 (A) and (B) (auditor's own code; no session modules imported).

g is computed by the Session 193 state machine read MSB-first: start (H, A, p) = (1, 0, 0);
bit 1: (H, 2A+1, 1); bit 0: (H + A + p, 2A, 0).  g(v) = H after reading the binary digits of v.
This is equivalent to g(0)=1, g(2v+1)=g(v), g(2v)=g(v)+v+(v mod 2) and is checked against that recursion.

Checks:
  (0) state machine == recursion for v < 2^16; 1 <= g(x) <= x+1 for x < 2^20; g(x) = x+1 iff x = 0 or x = 2^k with k >= 1.
  (1) C_w involution, g(C_w x) = g(x) for all x < 2^w, w <= 14 (leading zeros allowed); v -> v* maps Y_n onto Y_n.
  (A) g(h 2^m + l) = g(h) + g(c) + h c + [h odd][c odd] - 1, c = C_m(l): exhaustive m <= 8, 1 <= h < 2^9;
      random big-integer cases (m, bitlen(h) <= 150).
  (B) g(v) = g(h)+g(x)+g(tau)+h tau 2^w + h C_w(x) + tau x + eps - 2, eps = [x odd][tau odd] + [h odd][x_(w-1)=0]:
      exhaustive over all words of length n <= 13 and all splits a >= 1, w >= 1, s >= 0 (s = 0: tau = 0);
      random big-integer cases with total length up to 300.
"""
import random, sys
from functools import lru_cache

sys.setrecursionlimit(10000)


def g(v):
    H, A, p = 1, 0, 0
    for ch in bin(v)[2:] if v > 0 else '':
        if ch == '1':
            A, p = 2 * A + 1, 1
        else:
            H, A, p = H + A + p, 2 * A, 0
    return H


@lru_cache(maxsize=None)
def g_rec(v):
    if v == 0:
        return 1
    if v & 1:
        return g_rec(v >> 1)
    u = v >> 1
    return g_rec(u) + u + (u & 1)


def Cw(x, w):
    # reversed complement of the w-bit word x: sum_(i<w) (1 - x_i) 2^(w-1-i)
    assert 0 <= x < (1 << w) or w == 0
    s = format(x, '0%db' % w) if w > 0 else ''
    # s[0] is bit w-1; reversed complement: read bits from bit 0 upward as the new MSB-first string, complemented
    comp_rev = ''.join('1' if ch == '0' else '0' for ch in reversed(s))
    return int(comp_rev, 2) if w > 0 else 0


def main():
    out = []
    # (0)
    bad = sum(1 for v in range(1 << 16) if g(v) != g_rec(v))
    out.append('(0) state machine vs recursion, v < 2^16: mismatches %d' % bad)
    mx = 0; eq = []; bad_lo = 0
    G = [0] * (1 << 20)
    G[0] = 1
    for v in range(1, 1 << 20):
        G[v] = G[v >> 1] if v & 1 else G[v >> 1] + (v >> 1) + ((v >> 1) & 1)
    for x in range(1 << 20):
        if G[x] < 1:
            bad_lo += 1
        d = G[x] - x - 1
        mx = max(mx, d)
        if d == 0:
            eq.append(x)
    pow2 = [0] + [1 << k for k in range(1, 20)]  # 2^0 = 1 is not an equality case: g(1) = 1
    out.append('(0) x < 2^20: min violations of g >= 1: %d; max g(x)-x-1 = %d; equality set == {0} u {2^k : k >= 1}: %s'
               % (bad_lo, mx, sorted(eq) == sorted(pow2)))
    # (1)
    bad_inv = bad_conv = 0
    for w in range(1, 15):
        for x in range(1 << w):
            c = Cw(x, w)
            if Cw(c, w) != x:
                bad_inv += 1
            if G[c] != G[x]:
                bad_conv += 1
    out.append('(1) w <= 14, all x < 2^w: C_w involution failures %d; g(C_w x) != g(x) failures %d' % (bad_inv, bad_conv))
    bad_Y = 0
    for n in range(2, 19):
        Y = set(range(1 << (n - 1), 1 << n, 2))
        img = set(Cw(v, n) for v in Y)
        if img != Y:
            bad_Y += 1
    out.append('(1) v -> v* maps Y_n onto Y_n for n = 2..18: failures %d' % bad_Y)
    # (A) exhaustive
    badA = 0; cnt = 0
    for m in range(1, 9):
        for h in range(1, 1 << 9):
            for l in range(1 << m):
                c = Cw(l, m)
                lhs = G[(h << m) + l] if (h << m) + l < (1 << 20) else g((h << m) + l)
                rhs = G[h] + G[c] + h * c + (h & 1) * (c & 1) - 1
                cnt += 1
                if lhs != rhs:
                    badA += 1
    out.append('(A) exhaustive m <= 8, 1 <= h < 2^9: %d cases, violations %d' % (cnt, badA))
    random.seed(20261106)
    badA2 = 0
    for _ in range(20000):
        m = random.randint(1, 150)
        hb = random.randint(1, 150)
        h = random.randrange(1 << (hb - 1), 1 << hb)
        l = random.randrange(0, 1 << m)
        c = Cw(l, m)
        if g((h << m) + l) != g(h) + g(c) + h * c + (h & 1) * (c & 1) - 1:
            badA2 += 1
    out.append('(A) random big-integer (m, bitlen h <= 150): 20000 cases, violations %d' % badA2)
    # (B) exhaustive
    badB = 0; cnt = 0
    for n in range(2, 14):
        for v in range(1 << (n - 1), 1 << n):
            for s in range(0, n - 1):
                for w in range(1, n - s):
                    a = n - w - s
                    if a < 1:
                        continue
                    h = v >> (w + s)
                    x = (v >> s) & ((1 << w) - 1)
                    sig = v & ((1 << s) - 1)
                    tau = Cw(sig, s) if s > 0 else 0
                    eps = (x & 1) * (tau & 1) + (h & 1) * (1 - ((x >> (w - 1)) & 1))
                    rhs = G[h] + G[x] + G[tau] + h * tau * (1 << w) + h * Cw(x, w) + tau * x + eps - 2
                    cnt += 1
                    if rhs != G[v]:
                        badB += 1
    out.append('(B) exhaustive, all words of length 2..13 and all splits (a>=1, w>=1, s>=0): %d cases, violations %d'
               % (cnt, badB))
    badB2 = 0
    for _ in range(20000):
        a = random.randint(1, 100); w = random.randint(1, 100); s = random.randint(0, 100)
        h = random.randrange(1 << (a - 1), 1 << a)
        x = random.randrange(0, 1 << w)
        sig = random.randrange(0, 1 << s) if s > 0 else 0
        v = (h << (w + s)) | (x << s) | sig
        tau = Cw(sig, s) if s > 0 else 0
        eps = (x & 1) * (tau & 1) + (h & 1) * (1 - ((x >> (w - 1)) & 1))
        rhs = g(h) + g(x) + g(tau) + h * tau * (1 << w) + h * Cw(x, w) + tau * x + eps - 2
        if rhs != g(v):
            badB2 += 1
    out.append('(B) random big-integer (a, w, s <= 100): 20000 cases, violations %d' % badB2)
    print('\n'.join(out))


if __name__ == '__main__':
    main()
