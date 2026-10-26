#!/usr/bin/env python3
"""Referee: direct TM verification of the transcribed simplified rule systems A6, A11, A4, A2, A5, A10
(paper's second systems), branch by branch on exhaustive small and random moderate parameters.
Rules below are transcribed by the referee from the paper text, independently of the session code."""
import re, sys, random, time
sys.path.insert(0, __file__.rsplit('\\', 1)[0] if '\\' in __file__ else __file__.rsplit('/', 1)[0])
from tmcheck import check, HALT

# ---------------------------------------------------------------- A6
A6 = '1RB0LF_1RC1RA_0RD0RB_1LE---_0LA1LE_1RD1LA'


def a6_enc(p):
    u, v, w = p
    return 0, '11111' + '10' * u + '1', '11' + '01' * v + '1' * w + '01'


RA6L = re.compile(r'^11111((?:10)*)1$')
RA6R = re.compile(r'^11((?:01)*)(1*)01$')


def a6_parse(state, left, right):
    if state != 0:
        return None
    m1 = RA6L.match(left)
    if not m1:
        return None
    m2 = RA6R.match(right)
    if not m2:
        return None
    return (len(m1.group(1)) // 2, len(m2.group(1)) // 2, len(m2.group(2)))


def a6_rule(p):
    u, v, w = p
    if u == 0:
        if w == 0:
            return (0, 0, 2 * v + 12) if v % 2 == 0 else (2, 4, 2 * v + 7)
        if w == 1:
            return (2, 4, 2 * v + 9) if v % 2 == 0 else (0, 0, 2 * v + 14)
        if w == 2:
            return (2, v + 4, 1) if v % 2 == 0 else (0, v + 4, 1)
        return (v + 3, 0, w - 3)
    if w >= 2:
        return (u - 1, v + 2, w - 2)
    if w == 0:
        if u == 1:
            return (0, 0, 2 * v + 12)
        return (2, 4, 2 * u + 2 * v + 5) if u % 2 == 0 else (0, 0, 2 * u + 2 * v + 10)
    if u == 1:
        return HALT
    if u == 2:
        return (0, 0, 2 * v + 10)
    return (2, u, 2 * v + 5) if u % 2 == 0 else (0, u, 2 * v + 5)


# ---------------------------------------------------------------- A11
A11 = '1RB1LA_1LC0RE_1LF1LD_1RE0LC_1LB1RA_---0LA'


def a11_enc(p):
    a, b, c, d = p
    return 2, '1' * (a + 1), '00' + '1' * (b + 4) + '01' * c + '1' * d + '011'


RA11L = re.compile(r'^(1+)$')
RA11R = re.compile(r'^00(1{4,})((?:01)*)(1*)011$')


def a11_parse(state, left, right):
    if state != 2:
        return None
    m1 = RA11L.match(left)
    if not m1:
        return None
    m2 = RA11R.match(right)
    if not m2:
        return None
    a = len(m1.group(1)) - 1
    b = len(m2.group(1)) - 4
    c = len(m2.group(2)) // 2
    d = len(m2.group(3))
    if c == 0:
        # "00 1^m 011": normal form d = 3, b = m - 7 (C(a,b,0,6) = C(a,b+3,0,3))
        m = b + 4 + d
        if m < 7:
            return None
        return (a, m - 7, 0, 3)
    if d not in (3, 6):
        return None
    return (a, b, c, d)


def a11_rule(p):
    """Paper A.11 second system; the identity C(a,b,0,6) = C(a,b+3,0,3) is applied to the INPUT normal form."""
    a, b, c, d = p
    if a >= 2:
        return (a - 2, b + 3, c, d)
    if c >= 2:
        return (a + b + 4, 3, c - 2, d)
    if c == 1 and d == 6:
        return (a + b + 4, 0, 3, 3)
    if c == 0:
        if d == 6:
            b, d = b + 3, 3
        if b % 2 == 0:
            return (a, 0, b // 2 + 3, 3)
        return HALT if a == 0 else (6, 3, (b - 1) // 2, 3)
    # c == 1, d == 3
    if b % 2 == 1:
        return (a, 0, (b - 1) // 2 + 3, 6)
    if a == 0:
        return HALT
    return (6, 0, 3, 3) if b == 0 else (6, 3, (b - 2) // 2, 6)


def a11_norm(p):
    if p == HALT:
        return p
    a, b, c, d = p
    if c == 0 and d == 6:
        return (a, b + 3, 0, 3)
    return p


# ---------------------------------------------------------------- A4
A4 = '1RB0LF_1LC1RA_0RE0RD_---1LE_1LF1RC_1LC1LA'


def a4_enc(p):
    a, b = p
    return 4, '', '10' * a + '1' * b


RA4R = re.compile(r'^((?:10)*)(1*)$')


def a4_parse(state, left, right):
    if state != 4 or left != '':
        return None
    m = RA4R.match(right)
    if not m:
        return None
    return (len(m.group(1)) // 2, len(m.group(2)))


def a4_rule(p):
    a, b = p
    if b > 2 * a + 2:
        return (2 * a + 2, b - 2 * a - 2)
    if b == 1:
        return (2, 2 * a + 2)
    if b == 2:
        return HALT
    if b % 2 == 0:
        k = b // 2
        if a >= 1 and 2 <= k <= a + 1:
            x = 6 * a + 2 * k + 2
            return (x // 4, x % 4 + 1)
        return None
    k = (b - 1) // 2
    if a >= 1 and 1 <= k <= a:
        x = 6 * a + 2 * k + 6
        return (x // 4, x % 4 + 1)
    return None


# ---------------------------------------------------------------- A2 (second system, D(a,b,w))
A2 = '1RB1RD_0LC1RA_0RA1LC_1RA1LE_1RF0LD_---0RE'


def a2_enc(p):
    a, b, w = p
    return 3, '10' * a + '1' * b, '1' * (w + 1)


RA2L = re.compile(r'^((?:10)*)(1+)$')
RA2R = re.compile(r'^(1+)$')


def a2_parse(state, left, right):
    if state != 3:
        return None
    m1 = RA2L.match(left)
    m2 = RA2R.match(right)
    if not m1 or not m2:
        return None
    a, b, w = len(m1.group(1)) // 2, len(m1.group(2)), len(m2.group(1)) - 1
    if a >= 1 and b % 2 == 0:
        return None
    return (a, b, w)


def a2_rule(p):
    a, b, w = p
    if a == 0:
        if b % 2 == 1:
            k = (b - 1) // 2
            if w <= 1:
                return (0, 1, 2 * k + 2 * w + 2)
            return (0, 2 * k + 4, w - 2)
        k = (b - 2) // 2
        if w == 0:
            return HALT
        if w <= 3:
            return (k + 1, 1, 2 * w - 1)
        return (k + 1, 5, w - 4)
    k = (b - 1) // 2
    if w <= 1:
        return (a - 1, 1, 2 * k + 2 * w + 3)
    return (a - 1, 2 * k + 5, w - 2)


# ---------------------------------------------------------------- A5 (second system, D_1(a,z,w))
A5 = '1RB1RD_0LC1RA_1RA1LC_0RB1LE_1RF0LD_---0RE'


def a5_enc(p):
    a, z, w = p
    return 3, '10' * a + '1' * z, '1' + '1' * w


def a5_parse(state, left, right):
    if state != 3:
        return None
    m1 = RA2L.match(left)
    m2 = RA2R.match(right)
    if not m1 or not m2:
        return None
    a, z, w = len(m1.group(1)) // 2, len(m1.group(2)), len(m2.group(1)) - 1
    return (a, z, w)


def a5_rule(p):
    a, z, w = p
    if a == 0:
        if z % 2 == 1:
            k = (z - 1) // 2
            if w >= 2:
                return (0, 2 * k + 4, w - 2)
            return (0, 2, 2 * k + 2) if w == 0 else (0, 2, 2 * k + 4)
        k = (z - 2) // 2
        if w >= 4:
            return (k + 1, 5, w - 4)
        if w == 0:
            return HALT
        return (k, 3, 2 * w - 1)
    if z % 2 == 0:
        return None
    k = (z - 1) // 2
    if w >= 2:
        return (a - 1, 2 * k + 5, w - 2)
    if a == 1:
        return (0, 2, 2 * k + 3) if w == 0 else (0, 2, 2 * k + 5)
    return (a - 2, 3, 2 * k + 3) if w == 0 else (a - 2, 3, 2 * k + 5)


# ---------------------------------------------------------------- A10 (second system, E(a,b,c))
A10 = '1RB1LE_1RC1RA_0LD1RB_1LE1LD_1RF0LA_---0RE'


def a10_enc(p):
    a, b, c = p
    return 4, '10' * b + '1' * (2 * a + 2), '1' * (c + 1)


def a10_parse(state, left, right):
    if state != 4:
        return None
    m1 = RA2L.match(left)
    m2 = RA2R.match(right)
    if not m1 or not m2:
        return None
    n1 = len(m1.group(2))
    if n1 % 2 == 1 or n1 < 2:
        return None
    return ((n1 - 2) // 2, len(m1.group(1)) // 2, len(m2.group(1)) - 1)


def a10_rule(p):
    a, b, c = p
    if c == 0:
        if a >= 1:
            return (a - 1, b, 1)
        if b == 0:
            return (1, 0, 1)
        if b == 1:
            return (1, 0, 2)
        if b == 2:
            return (1, 0, 6)
        return (1, b - 3, 4)
    if b == 0:
        if c <= 2:
            return (1, 0, 2 * a + 2 * c + 1)
        if c == 3:
            return HALT
        if c <= 6:
            return (1, a + 1, 2 * (c - 4) + 2)
        return (1, a + 3, c - 6)
    if c >= 3:
        return (a + 2, b - 1, c - 2)
    if b == 1:
        return (1, 0, 2 * a + 2 * c + 2)
    if b == 2:
        return (1, 2, 2 * a + 2 * c)
    return (1, b - 3, 2 * a + 2 * c + 4)


FAM = {
    'A6': (A6, a6_enc, a6_parse, a6_rule, 3),
    'A11': (A11, a11_enc, a11_parse, a11_rule, 4),
    'A4': (A4, a4_enc, a4_parse, a4_rule, 2),
    'A2': (A2, a2_enc, a2_parse, a2_rule, 3),
    'A5': (A5, a5_enc, a5_parse, a5_rule, 3),
    'A10': (A10, a10_enc, a10_parse, a10_rule, 3),
}


def params_for(name, rng, small, nrand, big):
    dim = FAM[name][4]
    pts = []
    import itertools
    for t in itertools.product(range(small), repeat=dim):
        pts.append(t)
    for _ in range(nrand):
        pts.append(tuple(rng.randrange(0, big) for _ in range(dim)))
    out = []
    for t in pts:
        if name == 'A11':
            a, b, c, d = t
            t = (a, b, c, 3 if d % 2 == 0 else 6)
        if name == 'A2':
            a, b, w = t
            if b == 0:
                continue
            if a >= 1 and b % 2 == 0:
                continue
        if name == 'A5':
            a, z, w = t
            if z == 0 or (a >= 1 and z % 2 == 0):
                continue
        out.append(t)
    return out


if __name__ == '__main__':
    rng = random.Random(12345)
    names = sys.argv[1:] or list(FAM)
    for name in names:
        code, enc, parse, rule, dim = FAM[name]
        small = {2: 14, 3: 9, 4: 6}[dim]
        pts = params_for(name, rng, small, 400, 45)
        t0 = time.time()
        res = {}
        bad = []
        for p in pts:
            q = rule(p)
            if q is None:
                res['no-rule'] = res.get('no-rule', 0) + 1
                continue
            r, st, got = check(code, enc, parse, rule, p, max_steps=3 * 10 ** 6, norm=(a11_norm if name == 'A11' else (lambda x: x)))
            res[r] = res.get(r, 0) + 1
            if not r.startswith('ok'):
                bad.append((p, q, r, st, got))
        print(name, 'tested', len(pts), res, 'bad sample', bad[:6], f'({time.time() - t0:.0f}s)', flush=True)
