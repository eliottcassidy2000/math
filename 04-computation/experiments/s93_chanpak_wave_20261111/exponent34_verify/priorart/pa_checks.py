#!/usr/bin/env python3
"""Prior-art verifier (independent code): credit equivalences and the Hilbert-cube reading of (A).

Checks (exact integer arithmetic):
 P1. Q = 2 q where Q is Session 206's tail function (Q(0)=0, Q(2m)=Q(m), Q(2m+1)=Q(m)+m+(m mod 2)) and q is
     THM-5141's (q(0)=0, q(2h)=q(h), q(2h+1)=q(h)+ceil(h/2)).
 P2. Session 206's identity Q(Mv+r) = Q(v)+Q(r)+v rev_k(r)+topbit_k(r)(v mod 2) is THM-5141's block law
     2q(A 2^m + B) = 2q(A)+2q(B)+A rev_m(B)+(A mod 2) B_(m-1)  (exhaustive small range + random big integers).
 P3. Session 206's digit identity Q(m) = sum_j b_j floor(m/2^(j+1)) + sum_j b_j b_(j+1) equals 2 * (THM-5141's
     unrolled q(h) = sum_i h_i ceil(floor(h/2^(i+1))/2)).
 P4. THM-5093 / five-core-web pair-energy form: Q(m) = sum_(i<j) w_(j-i) x_i x_j with w_1 = 2, w_d = 2^(d-1).
 P5. x(h) = g(2h), g(v) = v + 1 - 2q(v) - (v mod 2) with THM-4890's g.
 P6. THM-4890 (B) window identity (g-normalisation, reversed complement C_w) versus THM-5141 (C)
     (q-normalisation, rev_w): for v = 2h, with tau = C_(s+1)(2Z): c = tau + 1, Xi_rev = -Xi_C, and
     c X - P Xi = gamma holds identically for all (P, x, x', Z) (not only collisions), where
     LHS := c X - P Xi - gamma + (x(h) - x(h')) must vanish... (we test x(h)-x(h') = cX - P Xi - gamma).
 P7. (A)1 exhaustively: no nonzero delta in {-2..2}^w with delta(2) = 0 = delta*(2) (w <= 9), and the base-b
     analogue (digits in {-2(b-1)..2(b-1)}) for b = 3 (w <= 6), b = 4 (w <= 5), b = 5 (w <= 4).
 P8. (A)3 by brute force over pairs: additive energy of {(x, rev_w x)} = 6^w (w <= 10), equal to the energy of
     the Boolean cube {0,1}^w in Z^w (computed independently as the number of quadruples), while the energy of
     the one-coordinate projection [0, 2^w) is (2*8^w + 2^w)/3 (reported), so the reversal coordinate is what
     makes the digit map a Freiman 2-isomorphism.
"""
import itertools, random, sys
from math import gcd

def q_thm(h):
    # THM-5141: q(0)=0, q(2h)=q(h), q(2h+1)=q(h)+ceil(h/2)
    r = 0
    while h:
        if h & 1:
            hh = h >> 1
            r += (hh + 1) // 2
        h >>= 1
    return r

def Q206(m):
    # Session 206: Q(0)=0, Q(2m)=Q(m), Q(2m+1)=Q(m)+m+(m mod 2)
    r = 0
    while m:
        if m & 1:
            mm = m >> 1
            r += mm + (mm & 1)
        m >>= 1
    return r

def g4890(v):
    # THM-4890: g(0)=1, g(2v+1)=g(v), g(2v)=g(v)+v+(v mod 2)
    r = 1
    # iterative: strip bits from the bottom
    stack = []
    while v:
        stack.append(v & 1)
        v >>= 1
    # rebuild from the top: start value 0 -> g = 1
    val = 0
    for b in reversed(stack):
        if b:
            val = 2 * val + 1          # g unchanged
        else:
            r += val + (val & 1)       # g(2v) = g(v) + v + (v mod 2)
            val = 2 * val
    return r

def rev(x, w):
    r = 0
    for i in range(w):
        if (x >> i) & 1:
            r |= 1 << (w - 1 - i)
    return r

def Cw(x, w):
    # reversed complement of a w-bit word
    r = 0
    for i in range(w):
        if not ((x >> i) & 1):
            r |= 1 << (w - 1 - i)
    return r

def xval(h):
    return 2 * h + 1 - 2 * q_thm(h)

out = []
def say(*a):
    s = " ".join(str(t) for t in a)
    print(s); out.append(s); sys.stdout.flush()

random.seed(20261111)

# P1
bad = sum(1 for m in range(1 << 18) if Q206(m) != 2 * q_thm(m))
say("P1  Q206(m) == 2 q(m) for m < 2^18: violations", bad)

# P2 exhaustive small + random big
bad = 0; n = 0
for k in range(1, 9):
    for v in range(0, 1 << 9):
        for r in range(1 << k):
            n += 1
            lhs = Q206((v << k) + r)
            rhs = Q206(v) + Q206(r) + v * rev(r, k) + ((r >> (k - 1)) & 1) * (v & 1)
            lhs2 = 2 * q_thm((v << k) + r)
            rhs2 = 2 * q_thm(v) + 2 * q_thm(r) + v * rev(r, k) + (v & 1) * ((r >> (k - 1)) & 1)
            if lhs != rhs or lhs2 != rhs2 or lhs != lhs2:
                bad += 1
say("P2  Session206 identity == THM-5141 block law, exhaustive k<=8, v<2^9:", n, "cases, violations", bad)
bad = 0
for _ in range(20000):
    k = random.randint(1, 80); v = random.getrandbits(random.randint(1, 200)); r = random.getrandbits(k)
    lhs = 2 * q_thm((v << k) + r)
    rhs = 2 * q_thm(v) + 2 * q_thm(r) + v * rev(r, k) + (v & 1) * ((r >> (k - 1)) & 1)
    if lhs != rhs: bad += 1
say("P2  block law, 20000 random big-integer cases (k<=80, v<2^200): violations", bad)

# P3, P4
bad3 = bad4 = 0
for m in range(1 << 14):
    bits = [(m >> j) & 1 for j in range(16)]
    d206 = sum(bits[j] * (m >> (j + 1)) for j in range(15)) + sum(bits[j] * bits[j + 1] for j in range(15))
    unroll = sum(bits[i] * (((m >> (i + 1)) + 1) // 2) for i in range(15))
    if d206 != Q206(m) or 2 * unroll != d206:
        bad3 += 1
    pe = 0
    for i in range(15):
        if bits[i]:
            for j in range(i + 1, 16):
                if bits[j]:
                    d = j - i
                    pe += 2 if d == 1 else (1 << (d - 1))
    if pe != Q206(m):
        bad4 += 1
say("P3  Session206 digit identity == 2 * (THM-5141 unrolled q), m < 2^14: violations", bad3)
say("P4  pair-energy form (w_1 = 2, w_d = 2^(d-1), THM-5093 / Session220) == Q, m < 2^14: violations", bad4)

# P5
bad = 0
for h in range(1 << 17):
    if xval(h) != g4890(2 * h): bad += 1
    v = h
    if g4890(v) != v + 1 - 2 * q_thm(v) - (v & 1): bad += 1
say("P5  x(h) == g(2h) and g(v) == v+1-2q(v)-(v mod 2), h < 2^17: violations", bad)
# sizes
bad = sum(1 for h in range(1, 1 << 17) if not (0 <= 2 * q_thm(h) <= h - 1 and h + 2 <= xval(h) <= 2 * h + 1))
say("P5b 0 <= q(h) <= (h-1)/2 and h+2 <= x(h) <= 2h+1 for 1 <= h < 2^17: violations", bad)

# P6: three-block identity difference, all (P, x, x', Z), small sizes, plus equivalence with THM-4890 (B)
bad = 0; n = 0; badB = 0
for tp in range(1, 5):
    for w in range(1, 6):
        for s in range(0, 5):
            for P in range(1 << (tp - 1), 1 << tp):
                for Z in range(1 << s):
                    c = (1 << (s + 1)) - rev(Z, s)
                    # THM-4890 normalisation for v = 2h: suffix sigma = Z|0 of length s+1
                    tau = Cw(2 * Z, s + 1)
                    if c != tau + 1: badB += 1
                    for x in range(1 << w):
                        hx = (P << (w + s)) + (x << s) + Z
                        for x2 in range(x):
                            n += 1
                            hx2 = (P << (w + s)) + (x2 << s) + Z
                            X = x - x2; Xi = rev(x, w) - rev(x2, w)
                            gam = 2 * (q_thm(x) - q_thm(x2)) + (P & 1) * (((x >> (w - 1)) & 1) - ((x2 >> (w - 1)) & 1)) \
                                  + (c & 1) * ((x & 1) - (x2 & 1))
                            if xval(hx) - xval(hx2) != c * X - P * Xi - gam:
                                bad += 1
                            # THM-4890 (B): g(v)-g(v') = (g(x)-g(x')) + h Xi_C + tau X + (eps-eps')
                            XiC = Cw(x, w) - Cw(x2, w)
                            eps = (x & 1) * (tau & 1) + (P & 1) * (1 - ((x >> (w - 1)) & 1))
                            eps2 = (x2 & 1) * (tau & 1) + (P & 1) * (1 - ((x2 >> (w - 1)) & 1))
                            if g4890(2 * hx) - g4890(2 * hx2) != (g4890(x) - g4890(x2)) + P * XiC + tau * X + (eps - eps2):
                                badB += 1
                            if XiC != -Xi: badB += 1
say("P6  x(h)-x(h') == c X - P Xi - gamma for all (P,x,x',Z), t'<=4, w<=5, s<=4:", n, "cases, violations", bad)
say("P6b same cases: THM-4890 (B) window identity holds with tau = c - 1, Xi_C = -Xi_rev: violations", badB)

# P7: (A)1 exhaustive, base 2 and other bases
def a1_check(b, w):
    D = 2 * (b - 1)
    cnt = 0; n = 0
    pw = [b ** i for i in range(w)]
    for delta in itertools.product(range(-D, D + 1), repeat=w):
        n += 1
        if sum(d * p for d, p in zip(delta, pw)) == 0 and sum(d * pw[w - 1 - i] for i, d in enumerate(delta)) == 0:
            if any(delta): cnt += 1
    return n, cnt
for b, wmax in ((2, 9), (3, 6), (4, 5), (5, 4)):
    for w in range(1, wmax + 1):
        n, cnt = a1_check(b, w)
        say("P7  base", b, "w", w, "strings", n, "nonzero common zeros of (D, D*):", cnt)

# P8: energy by brute force
for w in range(1, 11):
    N = 1 << w
    rv = [rev(x, w) for x in range(N)]
    cnt = {}
    for u in range(N):
        ru = rv[u]
        for u2 in range(N):
            key = (u - u2, ru - rv[u2])
            cnt[key] = cnt.get(key, 0) + 1
    E = sum(c * c for c in cnt.values())
    # projection energy
    c1 = {}
    for u in range(N):
        for u2 in range(N):
            c1[u - u2] = c1.get(u - u2, 0) + 1
    E1 = sum(c * c for c in c1.values())
    # every multiplicity a power of two, number of distinct differences 3^w
    pow2 = all((c & (c - 1)) == 0 for c in cnt.values())
    say("P8  w", w, "energy", E, "6^w", 6 ** w, "equal", E == 6 ** w, "| distinct differences", len(cnt), "3^w", 3 ** w,
        "| multiplicities powers of 2:", pow2, "| energy of [0,2^w) alone", E1, "(2*8^w+2^w)/3 =", (2 * 8 ** w + 2 ** w) // 3)

with open(__file__.replace(".py", ".out"), "w") as f:
    f.write("\n".join(out) + "\n")
