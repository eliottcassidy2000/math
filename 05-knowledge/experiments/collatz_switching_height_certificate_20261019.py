#!/usr/bin/env python3
"""Focused exact controls for source-guarded switched Collatz height ranks."""
from fractions import Fraction as Q
from itertools import product
import json

CHECKS = 0

def check(test, label):
    global CHECKS
    CHECKS += 1
    if not test:
        raise RuntimeError(label)


def vp(z, prime=2):
    z = Q(z)
    if not z:
        return None
    a, b, out = abs(z.numerator), z.denominator, 0
    while a % prime == 0:
        a //= prime
        out += 1
    while b % prime == 0:
        b //= prime
        out -= 1
    return out


def T(x):
    check(type(x) is int and x > 0, 'actual positive integer source')
    return (3*x+1)//2 if x % 2 else x//2


def walk(x, n):
    a = C = 0
    bits = []
    for t in range(n):
        bit = x % 2
        bits.append(str(bit))
        if bit:
            a += 1
            C = 3*C + 2**t
        x = T(x)
    return x, ''.join(bits), a, C


def coefficients(bits):
    a = C = 0
    for t, bit in enumerate(bits):
        if bit == '1':
            a += 1
            C = 3*C + 2**t
    n = len(bits)
    r = (-C*pow(3**a, -1, 2**n)) % 2**n
    return n, a, C, r


def word(name, bits):
    n, a, C, r = coefficients(bits)
    return dict(name=name, bits=bits, n=n, a=a, C=C, residue=r,
                slope=Q(3**a, 2**n), intercept=Q(C, 2**n),
                center=Q(C, 2**n-3**a))


F = word('F', '10')
H = word('H', '110101000')
J = word('J', '0111111111001101101000')
WORDS = {w['name']: w for w in (F, H, J)}


def apply(w, x):
    return w['slope']*x + w['intercept']


def legal(w, x):
    return type(x) is int and x > 0 and x % 2**w['n'] == w['residue']


# Symbolic coefficients establish each whole return cylinder; literal lifts
# are independent controls, not a residue census.
returns = []
for w in WORDS.values():
    r, n = w['residue'], w['n']
    lo, bits, a, C = walk(r, n)
    hi, ubits, ua, uC = walk(r+1, n)
    check((bits, a, C) == (w['bits'], w['a'], w['C']), 'lower word reconstruction')
    check(ua == a and 3**a + uC-C == 2**n, 'symbolic full-cylinder adjacent return')
    check(hi == lo+1, 'actual adjacent return anchor')
    check(coefficients(ubits)[3] == r+1, 'upper cylinder matches source plus one')
    for t in (0, 1, 17, 103):
        x = r + 2**n*t
        y, bb, aa, cc = walk(x, n)
        u, bu, au, cu = walk(x+1, n)
        check((bb, aa, cc) == (bits, a, C), 'lower cylinder lift')
        check((bu, au, cu) == (ubits, ua, uC), 'upper cylinder lift')
        check(y == apply(w, x) and u == y+1, 'lifted pair map')
    returns.append({**w, 'upper_bits': ubits, 'upper_C': uC})
check((J['n'], J['a'], J['C'], J['residue']) ==
      (22, 14, 10623730, 1817598), 'expanding return exact data')
check(J['slope'] > 1 and J['center'] == Q(-2124746, 117733), 'expanding negative center')

# Exact resonance formula on lawful source lifts.
resonance_controls = 0
for wi in WORDS.values():
    for t in (0, 1, 2, 7, 19, 61):
        x = wi['residue'] + 2**wi['n']*t
        z = apply(wi, x)
        check(z.denominator == 1 and z > 0, 'lawful integer resonance endpoint')
        ri = vp(x-wi['center'])
        if ri is None:
            continue
        for wj in WORDS.values():
            rj = vp(z-wj['center'])
            if wi['center'] == wj['center']:
                check(rj == ri-wi['n'], 'same-center exact budget')
            else:
                s = vp(wi['center']-wj['center'])
                if ri-wi['n'] != s:
                    check(rj == min(ri-wi['n'], s), 'nonresonant exact minimum')
                else:
                    check(rj is None or rj >= s, 'resonance minimum lower bound')
            resonance_controls += 1

# Universal forward recharge via independent 2-adic target and 3-adic inverse
# integrality. The initial rank is fixed even while the next rank grows.
recharge = []
for wi in WORDS.values():
    for wj in WORDS.values():
        if wi is wj:
            continue
        s = vp(wi['center']-wj['center'])
        for M in (wj['n'], 2*wj['n']+3):
            check(M > s, 'target rank exceeds center separation')
            m2, m3 = 2**(M+1), 3**wi['a']
            p = wj['center']
            r2 = (p.numerator*pow(p.denominator, -1, m2)+2**M) % m2
            r3 = wi['C']*pow(2**wi['n'], -1, m3) % m3
            z = r2 + m2*((r3-r2)*pow(m2, -1, m3) % m3)
            z += m2*m3
            x = Q(2**wi['n']*z-wi['C'], 3**wi['a'])
            check(x.denominator == 1 and x > 0, 'positive integer universal recharge source')
            x = int(x)
            check(legal(wi, x) and walk(x, wi['n'])[0] == z, 'actual recharge predecessor')
            check(vp(x-wi['center']) == wi['n']+s, 'fixed exact resonance gate')
            check(vp(z-wj['center']) == M, 'arbitrarily specified next rank')
            y, bb, _, _ = walk(z, wj['n']*(M//wj['n']))
            check(bb == wj['bits']*(M//wj['n']), 'complete next-word budget lawful')
            check(walk(z, wj['n']*(M//wj['n']+1))[1] !=
                  wj['bits']*(M//wj['n']+1), 'next-word budget exhausted exactly')
            recharge.append(dict(first=wi['name'], second=wj['name'], M=M,
                                 initial_source=x, endpoint=z,
                                 initial_rank=wi['n']+s, next_rank=M))

# The full fixed-center vector can recharge from a fixed vector (3,1,0).
full_vector = []
qpre = (H['center']-F['intercept'])/F['slope']
check(qpre == Q(-9, 431), 'off-bank preimage center')
check(qpre not in {w['center'] for w in WORDS.values()}, 'preimage outside center bank')
for M in (12, 24, 48, 96):
    mod = 2**(M+1)
    x = (qpre.numerator*pow(qpre.denominator, -1, mod)+2**M) % mod
    x += mod
    old = [vp(x-w['center']) for w in WORDS.values()]
    z = apply(F, x)
    new = [vp(z-w['center']) for w in WORDS.values()]
    check(old == [3, 1, 0], 'all initial center ranks fixed')
    check(z.denominator == 1 and legal(F, x), 'full-vector recharge legal')
    check(new[1] == M-2, 'full-vector output rank grows without bound')
    full_vector.append(dict(M=M, initial_source=x, endpoint=int(z), before=old, after=new))

# Exact two-vertex common-height certificate.
q = Q(3, 4)
omega = {'A': Q(1), 'B': Q(2)}
cap = {'A': Q(1), 'B': Q(4)}
edges = [('A', 'F', 'A'), ('A', 'H', 'A'), ('A', 'J', 'B'),
         ('B', 'F', 'B'), ('B', 'H', 'A')]
edge_receipts = []
for v, name, u in edges:
    w = WORDS[name]
    scaled = w['slope']*omega[v]/omega[u]
    image = apply(w, cap[v])
    check(scaled <= q, 'weighted slope contracts')
    check(image <= cap[u], 'affine cap preserved')
    edge_receipts.append(dict(source_vertex=v, word=name, target_vertex=u,
                              weighted_slope=scaled, cap_image=image,
                              target_cap=cap[u]))
delta = min((int(cap[v])+1-cap[v])/omega[v] for v in cap)
check(delta == Q(1, 2), 'integer rank floor')

# Explicit guarded descent and genuinely excluded JJ transition.
JHbits = J['bits']+H['bits']
n, a, C, residue = coefficients(JHbits)
check((n, a, C, residue) == (31, 18, 1284146834, 2119941118), 'composite return data')
center = Q(C, 2**n-3**a)
check(center == Q(1284146834, 1760063159) and center < 1, 'JH descent threshold')
t_class = ((H['residue']-int(apply(J, J['residue'])))*pow(3**J['a'], -1, 512)) % 512
check(t_class == 505 and J['residue']+2**22*t_class == residue, 'exact reset-cylinder lifting')
JH_after_J = walk(residue, 22)[0]
JH_after_H = walk(JH_after_J, 9)[0]
check((JH_after_J, JH_after_H) == (2417472043, 382451632) and JH_after_H < residue, 'actual JH smaller-source receipt')
check(not legal(H, int(apply(J, J['residue']))), 'J does not force H')
JJbits = J['bits']*2
_, _, _, JJsource = coefficients(JJbits)
check(JJsource == 4555502500862, 'exact lawful JJ source')
JJmiddle = walk(JJsource, 22)[0]
check(JJmiddle == 5194861231102 and legal(J, JJmiddle), 'actual next J rather than forced reset')
check(not legal(F, JJmiddle) and not legal(H, JJmiddle), 'JJ escapes the B graph')

# Long symbolic mixed paths are decoded to actual positive integer sources.
# All prefixes are independently replayed, with source guards and rank checks.
patterns = ['FHF', 'JH', 'JFH', 'JFJFH', 'FJFFHJFH', 'JH'*8, 'F'*20]
# JFJFH is deliberately invalid at the second J while in B.
mixed = []
for pattern in patterns:
    bits = ''.join(WORDS[c]['bits'] for c in pattern)
    nn, aa, cc, rr = coefficients(bits)
    x0 = rr + 2**nn
    x, vertex, rank0, count = x0, 'A', Q(x0-1), 0
    rows = []
    exit_word = None
    for name in pattern:
        options = [(v, c, u) for v, c, u in edges if v == vertex and c == name]
        if not options:
            exit_word = name
            break
        _, _, target = options[0]
        w = WORDS[name]
        check(legal(w, x), 'mixed path exact source guard')
        prior_rank = (x-cap[vertex])/omega[vertex]
        y, bb, _, _ = walk(x, w['n'])
        check(bb == w['bits'] and y == apply(w, x), 'mixed path literal replay')
        new_rank = (y-cap[target])/omega[target]
        check(new_rank <= q*prior_rank, 'mixed common-rank decrease')
        count += 1
        if y > cap[target]:
            check(new_rank >= delta and delta <= q**count*rank0, 'exact path-length budget')
        rows.append(dict(vertex=vertex, word=name, target=target,
                         source=x, endpoint=y, rank_before=prior_rank, rank_after=new_rank))
        x, vertex = y, target
    mixed.append(dict(pattern=pattern, initial_source=x0, admitted_edges=count,
                      exit_word=exit_word, transitions=rows))
check(next(r for r in mixed if r['pattern'] == 'JFJFH')['exit_word'] == 'J',
      'forbidden policy edge is an uncertified exit')

# Arbitrarily long J prefixes are controlled by an exact rational-center
# valuation, rather than an invented finite-state contraction.
expansive_prefixes = []
for L in (1, 2, 3, 5, 8, 13):
    bits = J['bits']*L
    n, a, C, r = coefficients(bits)
    x0 = r + 2**n
    x, bb, aa, cc = walk(x0, n)
    rank = vp(x0-J['center'])
    check(bb == bits and rank >= 22*L, 'long lawful positive J prefix')
    check(x > x0 and J['slope']**L > 1, 'expansive path cannot be contracting cycle')
    check((x-J['center']) == J['slope']**L*(x0-J['center']), 'exact expanding centered iterate')
    expansive_prefixes.append(dict(L=L, initial_source=x0, endpoint=x, initial_J_rank=rank))

# Actual tiny core receipts, stopped on first visit to 1.
core = []
for source in range(1, 5):
    x = source
    orbit = [x]
    while x != 1:
        x = T(x)
        orbit.append(x)
        check(len(orbit) <= 10, 'finite core receipt length')
    core.append(dict(source=source, orbit=orbit))


def encode(z):
    if isinstance(z, Q):
        return str(z)
    raise TypeError(type(z).__name__)


print(json.dumps(dict(
    status='PROVED conditional switching certificate; no universal guard coverage',
    checks=CHECKS, return_words=returns, resonance_controls=resonance_controls,
    universal_recharge_controls=recharge, full_center_vector_recharge=full_vector,
    graph_certificate=dict(q=q, weights=omega, caps=cap, delta=delta, edges=edge_receipts),
    guarded_JH_descent=dict(source_residue=residue, modulus=2**31,
                            J_lift_class_mod_512=t_class, after_J=JH_after_J, after_H=JH_after_H,
                            center=center),
    excluded_JJ=dict(source=JJsource, after_J=JJmiddle),
    mixed_path_controls=mixed, expanding_J_prefixes=expansive_prefixes,
    finite_core_receipts=core), indent=2, sort_keys=True, default=encode))
