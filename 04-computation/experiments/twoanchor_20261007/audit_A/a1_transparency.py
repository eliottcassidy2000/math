#!/usr/bin/env python3
"""Audit A, item 1: THM-4600 (1)-(3), independent code (no import from the session's scripts).

All arithmetic is exact: 2-adic integers are represented by rationals with odd denominators (Fraction),
the U-letter of such x is v2(numerator(3x+1)), U(x) = (3x+1)/2^letter, the Terras step uses the parity of
the numerator. Checks:
  1a/1b  letter(v) = a  <=>  v2(v - c_a) >= a+1, the same for u = c_a + 3^k (v - c_a) (k in Z, negative too),
         and U(u) - c_a = 3^k (U(v) - c_a), v2 drops by exactly a;  exhaustive mod 2^M for odd integer residues too.
  1c     common run lengths: run(u) = run(v) = max(0, floor((m-1)/a)) where m = v2(v - c_a).
  2      periodic anchors: c_w has U-word w; v2(v - c_w) > sum(w) => u, v, c_w share w and F_w(u)-c = 3^k(F_w(v)-c);
         sharpness: v2(v - c_w) = sum(w) never gives the word w; c = another point of the cycle does NOT give w.
  3      conjugation tau A tau^-1 = 3^k x + c(1-3^k) + lambda (e - c(1-3^k)); fixed iff e = c(1-3^k);
         consistency with THM-4581's Terras pair-chain table over the a Terras steps of a common letter;
         exit law: a non-anchored state leaves a common a-run after exactly min(run(v), floor((v2(eps)-1)/a)) letters
         (when v2(eps) >= 1).
"""
import random
from fractions import Fraction as Fr

rnd = random.Random(7_000_001)
NCHK = 0
def check(cond, msg):
    global NCHK
    if not cond:
        raise AssertionError(msg)
    NCHK += 1

def v2(q):
    q = Fr(q)
    if q == 0:
        return 10**9
    n, d = q.numerator, q.denominator
    assert d % 2 == 1 or n % 2 == 1
    r = 0
    while n % 2 == 0:
        n //= 2; r += 1
    while d % 2 == 0:
        d //= 2; r -= 1
    return r

def letter(x):
    return v2(3*Fr(x) + 1)

def Ustep(x):
    a = letter(x)
    return (3*Fr(x) + 1) / Fr(2)**a, a

def Fa(a, x):
    return (3*Fr(x) + 1) / Fr(2)**a

def Fw(w, x):
    for a in w:
        x = Fa(a, x)
    return x

def run_length(x, a, cap=10**4):
    L = 0
    while L < cap and letter(x) == a:
        x, _ = Ustep(x); L += 1
    return L

def rand_odd_2adic(bits=40):
    """random 2-adic unit as a rational with odd denominator (sometimes an integer)."""
    num = rnd.getrandbits(bits) * 2 + 1
    if rnd.random() < 0.5:
        num = -num
    den = 1 if rnd.random() < 0.4 else (rnd.getrandbits(12) * 2 + 1)
    return Fr(num, den)

# ---------- 1: run transparency ----------
def part1():
    for a in range(1, 9):
        c = Fr(1, 2**a - 3)
        check(3*c + 1 == 2**a * c, "c_a fixed point of F_a")
        check(v2(c) == 0, "c_a is a 2-adic unit")
        for _ in range(400):
            m = rnd.randint(0, 8*a + 3)              # m = v2(v - c)
            s = rand_odd_2adic()
            v = c + Fr(2)**m * s
            if v2(v) != 0:                           # need v odd (a 2-adic unit)
                continue
            k = rnd.randint(-12, 12)
            u = c + Fr(3)**k * (v - c)
            check(v2(u - c) == v2(v - c) == m, "v2 preserved by 3^k")
            la, lu = letter(v), letter(u)
            check((la == a) == (m >= a + 1), "(1a) letter a iff v2(v-c) >= a+1")
            check((lu == a) == (la == a), "(1a) u has letter a iff v does")
            if la == a:
                Uv, _ = Ustep(v); Uu, _ = Ustep(u)
                check(Uu - c == Fr(3)**k * (Uv - c), "(1b) relation persists")
                check(v2(Uv - c) == m - a, "(1b) valuation drops by a")
            Lv, Lu = run_length(v, a), run_length(u, a)
            check(Lv == Lu == max(0, (m - 1)//a), "(1c) common run length max(0,floor((m-1)/a))")
            # after the run, the relation still holds
            uu, vv = u, v
            for _ in range(Lv):
                uu, _ = Ustep(uu); vv, _ = Ustep(vv)
            check(uu - c == Fr(3)**k * (vv - c), "(1c) relation after the run")
    # exhaustive version on integers mod 2^M: letter(v) == a <=> v == c_a mod 2^(a+1) (for a+1 <= M)
    M = 14
    for a in range(1, 11):
        cm = pow(2**a - 3, -1, 2**M) % 2**M
        for vres in range(1, 2**M, 2):
            vv = vres + 2**M * 12345          # any lift; letter < M is decided by the residue
            la = v2(3*vv + 1)
            if la < M - 1 and a < M - 1:
                check((la == a) == ((vres - cm) % 2**(a+1) == 0), "exhaustive iff mod 2^M")
    print(f"1  run transparency: a=1..8 random exact (k in [-12,12]), exhaustive a<=10 mod 2^14: ok")

# ---------- 2: periodic anchors ----------
def cyc_point(w):
    # F_w(x) = (3^|w| x + B)/2^S ; fixed point c = B/(2^S - 3^|w|)
    S = sum(w); P = 3**len(w)
    B = Fw(w, Fr(0)) * 2**S
    return B / (2**S - P)

def part2():
    ex = {(1,): -1, (2,): 1, (1, 2): -5, (1, 1, 1, 2, 1, 1, 4): -17}
    for w, cval in ex.items():
        check(cyc_point(w) == cval, f"integral cycle point {w}")
    nsharp = 0; nrot_fail = 0; nrot = 0
    for _ in range(3000):
        L = rnd.randint(1, 6)
        w = tuple(rnd.randint(1, 5) for _ in range(L))
        S = sum(w)
        if 2**S == 3**L:
            continue
        c = cyc_point(w)
        check(Fw(w, c) == c, "fixed point")
        # the 2-adic U-word of c is w (Banach fixed point on the cylinder)
        x = c; word = []
        for _ in range(L):
            x, a = Ustep(x); word.append(a)
        check(tuple(word) == w and x == c, "c_w has U-word w")
        m = rnd.randint(S - 2, S + 6)
        s = rand_odd_2adic()
        v = c + Fr(2)**m * s
        k = rnd.randint(-8, 8)
        u = c + Fr(3)**k * (v - c)
        def word_of(z, n):
            out = []
            for _ in range(n):
                z, a = Ustep(z); out.append(a)
            return tuple(out), z
        wv, Fv = word_of(v, L); wu, Fu = word_of(u, L)
        if m > S:
            check(wv == w and wu == w, "(2) share the word w when v2(v-c) > sum w")
            check(Fu - c == Fr(3)**k * (Fv - c), "(2) relation carried around the cycle")
        elif m == S:
            check(wv != w and wu != w, "(2) sharp: v2(v-c) = sum w never gives w")
            nsharp += 1
        # rotation points of the cycle: c' = F_{w[:i]}(c) has word w[i:]+w[:i], not w (unless w is a power)
        if L >= 2:
            i = rnd.randint(1, L - 1)
            c2 = Fw(w[:i], c)
            rotw = w[i:] + w[:i]
            if rotw != w:
                nrot += 1
                v2_ = c2 + Fr(2)**(S + 3) * s
                wv2, _ = word_of(v2_, L)
                check(wv2 == rotw, "rotation point reads the rotated word")
                nrot_fail += (wv2 != w)
    print(f"2  periodic anchors: integral examples -1,+1,-5,-17 ok; 3000 random words ok; boundary v2 = sum(w) "
          f"never reads w ({nsharp} cases: condition is sharp, no off-by-one); "
          f"{nrot_fail}/{nrot} rotation points read a rotated word (so 'c a point of that cycle' must be c = c_w)")

# ---------- 3: Borel-torus conjugation ----------
def chain_step(k, e, beta):
    """THM-4581 Terras pair-chain table, own implementation: u = 3^k v + e, beta = parity of v."""
    e = Fr(e)
    sigma = v2(e) == 0 if e != 0 else False      # parity of e (e in Z[1/3] -> 2-adic unit or even)
    sigma = 1 if (e != 0 and v2(e) == 0) else 0
    if sigma == 0 and beta == 0:
        return k, e / 2
    if sigma == 0 and beta == 1:
        return k, (3*e + 1 - Fr(3)**k) / 2
    if sigma == 1 and beta == 0:
        return k + 1, (3*e + 1) / 2
    return k - 1, (e - Fr(3)**(k - 1)) / 2

def terras(x):
    x = Fr(x)
    return (3*x + 1)/2 if v2(x) == 0 else x/2

def part3():
    for _ in range(3000):
        a = rnd.randint(1, 7)
        c = Fr(1, 2**a - 3); lam = Fr(3, 2**a)
        k = rnd.randint(-6, 6)
        e = Fr(rnd.randint(-10**6, 10**6), 3**max(0, -k))
        tau = lambda x: c + lam*(x - c)
        taui = lambda x: c + (x - c)/lam
        A = lambda x: Fr(3)**k * x + e
        # conjugate evaluated at two points gives slope and intercept
        f0 = tau(A(taui(Fr(0)))); f1 = tau(A(taui(Fr(1))))
        check(f1 - f0 == Fr(3)**k, "slope preserved")
        check(f0 == c*(1 - Fr(3)**k) + lam*(e - c*(1 - Fr(3)**k)), "(3) conjugation formula")
        eps = e - c*(1 - Fr(3)**k)
        check((f0 == e) == (eps == 0), "(3) fixed under conjugation iff anchored")
        # commuting with tau iff anchored
        x0 = Fr(rnd.randint(-999, 999), 7)
        check((A(tau(x0)) == tau(A(x0))) == (eps == 0), "commutes iff anchored")
    # consistency with the THM-4581 Terras chain on actual pairs sharing one letter a
    n = 0
    for _ in range(4000):
        a = rnd.randint(1, 6)
        c = Fr(1, 2**a - 3)
        k = rnd.randint(0, 6)
        # v with letter a, u = 3^k v + e with u also letter a (generic e with v2(eps) large)
        v = c + Fr(2)**(a + 1 + rnd.randint(0, 5)) * rand_odd_2adic()
        if v2(v) != 0 or letter(v) != a:
            continue
        eps = Fr(2)**rnd.randint(a + 1, a + 8) * Fr(rnd.randint(-999, 999) * 2 + 1)
        e = c*(1 - Fr(3)**k) + eps
        u = Fr(3)**k * v + e
        if v2(u) != 0 or letter(u) != a:
            continue
        kk, ee = k, e; vv = v
        for _s in range(a):
            beta = 1 if v2(vv) == 0 else 0
            kk, ee = chain_step(kk, ee, beta); vv = terras(vv)
        lam = Fr(3, 2**a)
        check(kk == k and ee == c*(1 - Fr(3)**k) + lam*eps, "THM-4581 chain over a common letter = conjugation")
        check(Fa(a, u) == Fr(3)**kk * Fa(a, v) + ee, "chain state is the actual relation")
        n += 1
    # exit law for a non-anchored state inside a long v-run of letter a
    m_ok = 0
    for _ in range(3000):
        a = rnd.randint(1, 5)
        c = Fr(1, 2**a - 3)
        k = rnd.randint(-5, 5)
        mv = rnd.randint(a + 1, 10*a)
        v = c + Fr(2)**mv * rand_odd_2adic()
        if v2(v) != 0:
            continue
        g = rnd.randint(1, 10*a)
        eps = Fr(2)**g * Fr(rnd.randint(-999, 999)*2 + 1, 3**max(0, -k))
        u = c + Fr(3)**k*(v - c) + eps
        if v2(u) != 0:
            continue
        # common run length
        L = 0; uu, vv = u, v
        while letter(uu) == a and letter(vv) == a:
            uu, _ = Ustep(uu); vv, _ = Ustep(vv); L += 1
        Lv = max(0, (mv - 1)//a)
        pred = min(Lv, max(0, (g - 1)//a))
        check(L == pred, f"exit law a={a} mv={mv} g={g} L={L} pred={pred}")
        m_ok += 1
    print(f"3  Borel-torus: conjugation formula / fixed iff anchored / commutes iff anchored ok (3000 random); "
          f"THM-4581 chain over a common letter = conjugation ({n} pairs); "
          f"exit law L = min(run(v), floor((v2(eps)-1)/a)) exact ({m_ok} cases)")

if __name__ == "__main__":
    part1(); part2(); part3()
    print(f"ALL THM-4600 CHECKS PASSED ({NCHK} assertions)")
