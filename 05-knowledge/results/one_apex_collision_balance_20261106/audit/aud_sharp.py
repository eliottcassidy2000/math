"""Audit UPGRADE check: the sharp balance lemma |A - S| <= 1 (no exceptional region), and X * Xi < 0.

Proof outline being checked (see AUDIT.md, upgrade U1):
  (L0) exact local formula: for every w-bit string z (leading zeros allowed),
         2^w g(z) = z C_w(z) + 2^w (1 + #10(z)) - 2^w R(z),   R(z) = sum_{i<j, b_i=0, b_j=1} 2^(i-j-1)
       (positions MSB-first).  Put L(z) = #10(z) - R(z).
  (L1) same-sign impossibility: x = 1..0 (even, top bit 1) has g(x) >= 1 + 2^(w-2); x' = 0..1 (odd, top bit 0) has
       g(x') <= 2^(w-2); so h Xi + tau X >= 2 > 1 >= g(x') - g(x) + [tau odd] + [h odd].
  (L2) opposite signs, oriented X > 0 > Xi, w >= 2: x = 1..1, x' = 0..0, and for EVERY such configuration
         g(v) - g(v') = T X - H |Xi| - delta,   T = tau + C_w(x)/2^w,  H = h + x'/2^w,
         delta = L(x') - L(x) + [h odd] - [tau odd];
       with 0 <= L(x) <= #10(x) <= 1 + (k-1)/2 (k = w - dX), L(x') > -1, X, |Xi| >= 3, C_w(x) < 2^(w-1), x' < 2^(w-1).
       If A >= S + 2 then H|Xi| - T X > 2^(a-1) + 2^s + 2^(k-1) - 1/2 and so g(v) < g(v')  (sign rule);
       if S >= A + 2 then g(v) > g(v') (by the converse).  Hence collisions have |A - S| <= 1.
  (L3) w = 1: collision iff ceil(h/2) = ceil(tau/2), so |a - s| <= 1, and A = a, S = s.
Checks: (L0) exhaustively for w <= 16; (L1) bounds exhaustively for w <= 22; (L2) identity, bounds and sign rule on
300,000 random configurations (biased to carries and to A - S >= 2), including non-collisions; all intermediate claims
on every collision pair of Y_n, n <= 20.
"""
import random
from fractions import Fraction

NG = 22
G = [0] * (1 << NG)
G[0] = 1
for v in range(1, 1 << NG):
    G[v] = G[v >> 1] if v & 1 else G[v >> 1] + (v >> 1) + ((v >> 1) & 1)


def g(v):
    if v < (1 << NG):
        return G[v]
    H, A, p = 1, 0, 0
    for ch in bin(v)[2:]:
        if ch == '1':
            A, p = 2 * A + 1, 1
        else:
            H, A, p = H + A + p, 2 * A, 0
    return H


def Cw(x, w):
    s = format(x, '0%db' % w)
    return int(''.join('1' if c == '0' else '0' for c in reversed(s)), 2)


def n10(z, w):
    s = format(z, '0%db' % w)
    return sum(1 for i in range(w - 1) if s[i] == '1' and s[i + 1] == '0')


def R2w(z, w):
    """2^w R(z) as an integer: sum over zero bits z_k (LSB index k) of 2^(w-1-k) (z mod 2^k)."""
    tot = 0
    for k in range(w):
        if not (z >> k) & 1:
            tot += (z & ((1 << k) - 1)) << (w - 1 - k)
    return tot


def R2w_slow(z, w):
    s = format(z, '0%db' % w)
    tot = 0
    for i in range(w):
        if s[i] == '0':
            for j in range(i + 1, w):
                if s[j] == '1':
                    tot += 1 << (w + i - j - 1)
    return tot


def L(z, w):
    return Fraction(n10(z, w)) - Fraction(R2w(z, w), 1 << w)


def bitlen(x):
    return abs(x).bit_length()


def main():
    out = []
    # (L0)
    bad = 0; cnt = 0
    for w in range(1, 17):
        for z in range(1 << w):
            cnt += 1
            if (G[z] << w) != z * Cw(z, w) + (1 << w) * (1 + n10(z, w)) - R2w(z, w):
                bad += 1
    out.append('(L0) exact local formula, all w-bit strings, w <= 16: %d strings, violations %d' % (cnt, bad))
    badR = sum(1 for w in range(1, 13) for z in range(1 << w) if R2w(z, w) != R2w_slow(z, w))
    out.append('(L0) fast vs direct double-sum R, w <= 12: mismatches %d' % badR)
    # bounds on L for x = 1..1 and x' = 0..0
    badL = 0
    for w in range(2, 15):
        for z in range(1 << w):
            top, bot = z >> (w - 1), z & 1
            Lz = L(z, w)
            if top == 1 and bot == 1 and not (0 <= Lz <= n10(z, w)):
                badL += 1
            if top == 0 and bot == 0 and not (-1 < Lz <= n10(z, w)):
                badL += 1
    out.append('(L0) bounds 0 <= L(x) <= #10(x) for x = 1..1 and -1 < L(x\') <= #10(x\') for x\' = 0..0, w <= 14: violations %d' % badL)
    # (L1)
    badL1 = 0
    for w in range(2, NG + 1):
        mo = max(G[x] for x in range(1, 1 << (w - 1), 2))
        me = min(G[x] for x in range(1 << (w - 1), 1 << w, 2))
        if not (mo <= (1 << (w - 2)) and me >= 1 + (1 << (w - 2))):
            badL1 += 1
    out.append('(L1) max g(odd x\' < 2^(w-1)) <= 2^(w-2) and min g(even x >= 2^(w-1)) >= 1 + 2^(w-2), w = 2..%d: violations %d' % (NG, badL1))
    # (L2) random configurations, identity + sign rule
    random.seed(4890)
    bad_id = bad_sign = bad_bounds = 0; tested = 0; n_AS2 = 0; n_SA2 = 0
    for it in range(300000):
        w = random.randint(2, 40)
        a = random.randint(1, 40)
        s = random.randint(1, 40)
        if random.random() < 0.5:
            a = random.randint(1, 6); s = random.randint(1, 6)
        h = random.randrange(1 << (a - 1), 1 << a)
        if random.random() < 0.3:
            h = (1 << (a - 1)) + random.randint(0, min(3, (1 << (a - 1)) - 1)) if a > 1 else 1
        sig = random.randrange(0, 1 << s) & ~1            # bottom bit 0: tau has top bit 1
        if random.random() < 0.3:
            sig = 0
        tau = Cw(sig, s)
        mode = random.random()
        if w == 2:
            x, xp = 3, 0
        else:
            mid = random.randrange(0, 1 << (w - 2)); midp = random.randrange(0, 1 << (w - 2))
            if mode < 0.25:      # carry in X: x = 1 0^p .., x' = 0 1^p ..
                p = random.randint(0, w - 2)
                mid &= ~(((1 << p) - 1) << (w - 2 - p)); midp |= (((1 << p) - 1) << (w - 2 - p))
            elif mode < 0.5:     # carry in Xi: x = .. 0^p 1, x' = .. 1^p 0
                p = random.randint(0, w - 2)
                mid &= ~((1 << p) - 1); midp |= (1 << p) - 1
            elif mode < 0.75:    # both carries
                p = random.randint(0, (w - 2) // 2); q = random.randint(0, (w - 2) // 2)
                mid &= ~(((1 << p) - 1) << (w - 2 - p)); midp |= (((1 << p) - 1) << (w - 2 - p))
                mid &= ~((1 << q) - 1); midp |= (1 << q) - 1
            x = (1 << (w - 1)) | (mid << 1) | 1
            xp = (midp << 1)
        n = a + w + s
        v = (h << (w + s)) | (x << s) | sig
        vp = (h << (w + s)) | (xp << s) | sig
        X = x - xp; y = Cw(x, w); yp = Cw(xp, w); Xi = y - yp
        if not (X >= 3 and Xi <= -3 and y < (1 << (w - 1)) and xp < (1 << (w - 1))):
            bad_bounds += 1
        dX = w - bitlen(X); dXi = w - bitlen(Xi); k = w - dX
        A = a - dXi; S = s - dX
        T = Fraction(tau * (1 << w) + y, 1 << w); Hh = Fraction(h * (1 << w) + xp, 1 << w)
        delta = L(xp, w) - L(x, w) + (h & 1) - (tau & 1)
        lhs = g(v) - g(vp)
        if lhs != T * X - Hh * (-Xi) - delta:
            bad_id += 1
        if not (0 <= L(x, w) <= n10(x, w) <= 1 + Fraction(k - 1, 2) and L(xp, w) > -1):
            bad_bounds += 1
        if A >= S + 2:
            n_AS2 += 1
            if not (Hh * (-Xi) - T * X > (1 << (a - 1)) + (1 << s) + (1 << (k - 1)) - Fraction(1, 2)):
                bad_bounds += 1
            if not (lhs < 0):
                bad_sign += 1
        if S >= A + 2:
            n_SA2 += 1
            if not (lhs > 0):
                bad_sign += 1
        tested += 1
    out.append('(L2) random configurations x = 1..1, x\' = 0..0 (w >= 2): %d tested (A >= S+2: %d, S >= A+2: %d); '
               'identity violations %d; bound violations %d; sign-rule violations %d' % (tested, n_AS2, n_SA2, bad_id, bad_bounds, bad_sign))
    # (L3) w = 1 rule on random cases
    bad3 = 0
    for it in range(200000):
        a = random.randint(1, 30); s = random.randint(1, 30)
        h = random.randrange(1 << (a - 1), 1 << a); sig = random.randrange(0, 1 << s) & ~1; tau = Cw(sig, s)
        v = (h << (1 + s)) | (1 << s) | sig; vp = (h << (1 + s)) | sig
        coll = (g(v) == g(vp))
        if coll != ((h + 1) // 2 == (tau + 1) // 2):
            bad3 += 1
        if coll and abs(a - s) > 1:
            bad3 += 1
    out.append('(L3) w = 1: collision <=> ceil(h/2) = ceil(tau/2), and then |a - s| <= 1: 200000 random cases, violations %d' % bad3)
    # all collision pairs n <= 20: every intermediate claim
    tot_pairs = 0; bad_pairs = 0; hist = {}
    for n in range(4, 21):
        vals = {}
        for v in range(1 << (n - 1), 1 << n, 2):
            vals.setdefault(G[v], []).append(v)
        for grp in vals.values():
            if len(grp) < 2:
                continue
            for p in grp:
                for q in grp:
                    if p == q:
                        continue
                    d = p ^ q; t = d.bit_length() - 1; s = (d & -d).bit_length() - 1
                    a = n - 1 - t; w = t - s + 1
                    x = (p >> s) & ((1 << w) - 1); xp = (q >> s) & ((1 << w) - 1)
                    X = x - xp; Xi = Cw(x, w) - Cw(xp, w)
                    if X < 0:   # orient
                        x, xp, X, Xi = xp, x, -X, -Xi
                        p, q = q, p
                    dX = w - bitlen(X); dXi = w - bitlen(Xi); A = a - dXi; S = s - dX
                    ok = Xi < 0 and abs(A - S) <= 1
                    h = p >> (w + s); sig = p & ((1 << s) - 1); tau = Cw(sig, s)
                    if w == 1:
                        ok = ok and x == 1 and xp == 0 and (h + 1) // 2 == (tau + 1) // 2 and A == a and S == s
                    else:
                        ok = ok and (x >> (w - 1)) == 1 and (x & 1) == 1 and (xp >> (w - 1)) == 0 and (xp & 1) == 0
                        ok = ok and X >= 3 and -Xi >= 3 and -Xi >= (1 << (w - 1 - dXi)) + 1 and X <= (1 << (w - dX)) - 1
                        T = Fraction(tau * (1 << w) + Cw(x, w), 1 << w); Hh = Fraction(h * (1 << w) + xp, 1 << w)
                        delta = L(xp, w) - L(x, w) + (h & 1) - (tau & 1)
                        ok = ok and (T * X - Hh * (-Xi) == delta)
                    tot_pairs += 1
                    hist[A - S] = hist.get(A - S, 0) + 1
                    if not ok:
                        bad_pairs += 1
    out.append('(all pairs) every collision pair of Y_n, n = 4..20 (converse included): %d pairs, failures of any step %d; '
               'distribution of A - S: %s' % (tot_pairs, bad_pairs, dict(sorted(hist.items()))))
    print('\n'.join(out))


if __name__ == '__main__':
    main()
