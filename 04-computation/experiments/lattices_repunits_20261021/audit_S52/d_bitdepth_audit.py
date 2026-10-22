"""Audit S52 (independent code) - part D: the e8.pdf 'Bit Depth trace lock' (a_p = 5 at q = 5/4), without SEA.

Curve (as printed in the PDF's listing): y^2 = x^3 + A x + B over F_p, p = 2^521 - 1, A = int(-zeta(5) * 10^160),
B = int((q * pi^2) * 10^160) reduced mod p, q a Python float (5/4, 1, 3/2, 2, 240/96 are exact binary floats), mpmath at
its default 15 digits (the listing sets no precision) and, as the 'intended' exact reading, at 200 digits.

D1  For every curve: a random point P (p = 3 mod 4, sqrt = f^((p+1)/4)), and [p + 1 - 5] P computed with my own Jacobian
    double-and-add.  [p + 1 - 5] P != O proves #E != p - 4, i.e. a_p != 5 (Lagrange).  Same for a_p = -5 (twist reading).
D2  The session's SEA values N = p + 1 - a_p (bitdepth_trace_test.out / scratch bitdepth_test.out, read as data):
    [N] P = O for 3 random points per curve, and a certificate that N = #E(F_p): trial division of N, Miller-Rabin on the
    cofactor; if N = m * c with c a probable prime > 4 sqrt(p) and [m] P != O then ord(P) >= c > 4 sqrt(p), so N is the
    only multiple of ord(P) in the Hasse interval.
D3  If cypari is importable, cross-check [N] P = O with PARI's ellmul for two curves.
"""
import sys, random, time, math
import mpmath

p = 2 ** 521 - 1
S = 10 ** 160
FAIL = []
T0 = time.time()


def check(c, msg):
    print(('  ok   ' if c else '  FAIL ') + msg, flush=True)
    if not c:
        FAIL.append(msg)


def jdouble(P, a):
    if P is None:
        return None
    X, Y, Z = P
    if Y == 0:
        return None
    YY = Y * Y % p
    Sv = 4 * X * YY % p
    Z2 = Z * Z % p
    M = (3 * X * X + a * Z2 * Z2) % p
    X3 = (M * M - 2 * Sv) % p
    Y3 = (M * (Sv - X3) - 8 * YY * YY) % p
    Z3 = 2 * Y * Z % p
    return (X3, Y3, Z3)


def jadd(P, Q, a):
    if P is None:
        return Q
    if Q is None:
        return P
    X1, Y1, Z1 = P
    X2, Y2, Z2 = Q
    Z1s, Z2s = Z1 * Z1 % p, Z2 * Z2 % p
    U1, U2 = X1 * Z2s % p, X2 * Z1s % p
    S1, S2 = Y1 * Z2s * Z2 % p, Y2 * Z1s * Z1 % p
    if U1 == U2:
        if S1 != S2:
            return None
        return jdouble(P, a)
    H = (U2 - U1) % p
    R = (S2 - S1) % p
    H2 = H * H % p
    H3 = H2 * H % p
    X3 = (R * R - H3 - 2 * U1 * H2) % p
    Y3 = (R * (U1 * H2 - X3) - S1 * H3) % p
    Z3 = H * Z1 * Z2 % p
    return (X3, Y3, Z3)


def jmul(k, P, a):
    R = None
    Q = P
    while k > 0:
        if k & 1:
            R = jadd(R, Q, a)
        Q = jdouble(Q, a)
        k >>= 1
    return R


def to_affine(P):
    if P is None:
        return None
    X, Y, Z = P
    zi = pow(Z, -1, p)
    return (X * zi * zi % p, Y * zi * zi * zi % p)


def random_point(A, B, rng):
    while True:
        x = rng.randrange(p)
        f = (x * x * x + A * x + B) % p
        if f == 0:
            continue
        if pow(f, (p - 1) // 2, p) == 1:
            y = pow(f, (p + 1) // 4, p)
            assert y * y % p == f
            return (x, y, 1)


def coeffs(dps, q):
    mpmath.mp.dps = dps
    A = int(-mpmath.zeta(5) * S) % p
    B = int((q * mpmath.pi ** 2) * S) % p
    return A, B


# session a_p values (data copied from bitdepth_trace_test.out lines 2-8 and the scratch bitdepth_test.out)
SESSION = {
    (15, '5/4'): -1302531715252292841663517700830180722205684373862232912635285944879446963906240,
    (15, '1/1'): 1930792719564996148492179509670253032096186738556489007322740281714242002349897,
    (15, '3/2'): -3085033469148451729370559745186191284513507956521804865967525131285591503349324,
    (15, '2/1'): -2340817031207933222034824492927345198928402869572894254097796679950123360939704,
    (15, '240/96'): 1828796663417158657210921705529160554650936802477122606085979220176000218332213,
    (200, '5/4'): 2545057371761942293574147350006463908150653644234194483686312706195515222162558,
    (200, '1/1'): 1565542834791342731980242635473913454319625572816367907667308199623514318593140,
    (200, '3/2'): 908907202839395486299605655358368125923951482499962276692243976095417374735272,
    (200, '2/1'): -4183214284888864027608389101485960671805942643423924515734636572284290273638648,
    (200, '240/96'): -860158236365072690035737218894445773260182725912836975844973021915092162789280,
}
QS = {'5/4': 5 / 4, '1/1': 1.0, '3/2': 3 / 2, '2/1': 2.0, '240/96': 240 / 96}


def miller_rabin(n, rounds=40, rng=random.Random(5)):
    if n < 2:
        return False
    for sp in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if n % sp == 0:
            return n == sp
    d, s = n - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for _ in range(rounds):
        a = rng.randrange(2, n - 1)
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


def small_factors(n, bound=2 * 10 ** 6):
    fs = []
    for q in [2, 3] + list(range(5, bound, 2)):
        while n % q == 0:
            fs.append(q)
            n //= q
        if q * q > n:
            break
    return fs, n


if __name__ == '__main__':
    rng = random.Random(20261022)
    four_sqrt_p = 4 * math.isqrt(p) + 4
    curves = {}
    print('D1  is a_p = 5 (or -5)?  [p + 1 -+ 5] P for a random point P')
    for dps in (15, 200):
        for qs, q in QS.items():
            A, B = coeffs(dps, q)
            disc = (4 * A ** 3 + 27 * B ** 2) % p
            P = random_point(A, B, rng)
            Q5 = jmul(p + 1 - 5, P, A)
            Qm5 = jmul(p + 1 + 5, P, A)
            curves[(dps, qs)] = (A, B)
            check(disc != 0 and Q5 is not None and Qm5 is not None,
                  f'dps {dps:3d} q = {qs:6s}: nonsingular; [p - 4] P != O and [p + 6] P != O, so a_p != 5 and a_p != -5 (PROVED)  [{time.time() - T0:.1f}s]')
    print('D2  the session SEA values: [N] P = O and an order certificate')
    for key, ap in SESSION.items():
        A, B = curves[key]
        N = p + 1 - ap
        ok = abs(ap) <= 2 * math.isqrt(p) + 2
        for _ in range(3):
            P = random_point(A, B, rng)
            ok &= jmul(N, P, A) is None
        fs, cof = small_factors(N)
        cert = ''
        if cof > four_sqrt_p and miller_rabin(cof):
            m = N // cof
            P = random_point(A, B, rng)
            if jmul(m, P, A) is not None:
                cert = f'CERTIFIED: N = {m} * c, c a {cof.bit_length()}-bit probable prime > 4 sqrt(p), [N/c] P != O'
        if not cert:
            cert = f'uncertified beyond [N]P = O (small factors {fs[:8]}{"..." if len(fs) > 8 else ""}, cofactor {cof.bit_length()} bits, prime: {miller_rabin(cof)})'
        check(ok, f'dps {key[0]:3d} q = {key[1]:6s}: a_p = session value ({ap.bit_length()} bits), [N] P = O for 3 random P; {cert}')
    print('D3  PARI cross-check (optional)')
    sys.path.insert(0, r'C:/Users/Eliott/AppData/Local/Temp/claude/C--Users-Eliott-Documents-GitHub-ephrepos-math/e1fb236f-ec78-4211-aa52-c187130b32b0/scratchpad/s52/pylib')
    try:
        from cypari import pari
        for key in ((15, '5/4'), (200, '5/4'), (200, '2/1')):
            A, B = curves[key]
            E = pari.ellinit([0, 0, 0, A, B], p)
            P = random_point(A, B, rng)
            Pa = to_affine(P)
            N = p + 1 - SESSION[key]
            Z = pari.ellmul(E, [pari.Mod(Pa[0], p), pari.Mod(Pa[1], p)], N)
            Z5 = pari.ellmul(E, [pari.Mod(Pa[0], p), pari.Mod(Pa[1], p)], p - 4)
            check(str(Z) == '[0]' and str(Z5) != '[0]', f'PARI ellmul: [N] P = O and [p - 4] P != O on curve {key}')
    except ImportError as e:
        print(f'     cypari not importable ({e}); skipped')
    print()
    print('ALL CHECKS PASSED' if not FAIL else f'FAILED: {FAIL}', f'({time.time() - T0:.0f}s)')
