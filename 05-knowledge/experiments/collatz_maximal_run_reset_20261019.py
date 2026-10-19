#!/usr/bin/env python3
"""Exact first-descent-safe expanding runs, maximal-run decoder, and reset tails."""
from fractions import Fraction as F
import json

P, Q, C = 4782969, 4194304, 5311865
d, B = 117733, 1062373
R_BITS = '1111111110011011010000'
CHECKS = 0


def check(test, label):
    global CHECKS
    CHECKS += 1
    if not test:
        raise RuntimeError(label)


def v2(n):
    if n == 0:
        raise ValueError('infinite valuation outside positive-source controls')
    n = abs(n)
    return (n & -n).bit_length()-1


def T(x):
    check(type(x) is int and x > 0, 'positive integer trajectory guard')
    return (3*x+1)//2 if x % 2 else x//2


def walk(x, n):
    start = x
    bits = []
    odds = carry = 0
    first_descent = None
    for t in range(n):
        bit = x % 2
        bits.append(str(bit))
        if bit:
            odds += 1
            carry = 3*carry+2**t
        x = T(x)
        if first_descent is None and x < start:
            first_descent = t+1
    return dict(endpoint=x, bits=''.join(bits), odds=odds, carry=carry,
                first_descent=first_descent)


def coefficients(bits):
    odds = carry = 0
    for t, bit in enumerate(bits):
        if bit == '1':
            odds += 1
            carry = 3*carry+2**t
    n = len(bits)
    residue = (-carry*pow(3**odds, -1, 2**n)) % 2**n
    return n, odds, carry, residue


RESETS = {}
for name, bits, shell in [('E', '0', 0), ('F', '10', 1), ('H', '110101000', 2)]:
    n, odds, carry, residue = coefficients(bits)
    a, q = 3**odds, 2**n
    RESETS[name] = dict(name=name, bits=bits, n=n, odds=odds, a=a, q=q,
                        b=carry, residue=residue, center=F(carry, q-a), shell=shell)


def run_data(L):
    PL, QL = P**L, Q**L
    CL = C*(PL-QL)//(P-Q)
    rL = (-CL*pow(PL, -1, QL)) % QL
    yL = (PL*rL+CL)//QL
    return PL, QL, CL, rL, yL


def compress(x):
    if type(x) is not int or x <= 0:
        raise ValueError('positive integer source required')
    rank = v2(d*x+B)
    L, s = divmod(rank, 22)
    m = (d*x+B) >> rank
    num = P**L*2**s*m-B
    check(num % d == 0, 'compressed exit integrality')
    y = num//d
    check(y > 0 and m % 2 == 1 and v2(d*y+B) == s, 'positive exit and exact residual shell')
    return dict(source=x, rank=rank, L=L, shell=s, odd_part=m, exit=y)


def required(w, x, y):
    j = 1
    center = w['center']
    while F(w['a'], w['q'])**j*(y-center) >= x-center:
        j += 1
    return j


def iterate(w, y, j):
    return w['center']+F(w['a'], w['q'])**j*(y-w['center'])


def decode(x):
    rec = compress(x)
    L, s, m, y = rec['L'], rec['shell'], rec['odd_part'], rec['exit']
    if L == 0:
        return {**rec, 'certified': False, 'reason': 'outside R cylinder'}
    if s == 0:
        name = 'E'
        available = v2(y)
        check(available == v2(P**L*m-B), 'even reset capacity from odd part')
    elif s == 1:
        name = 'F'
        rank = v2(y-1)
        check(rank == 1+v2(P**L*m-590053) and rank >= 2, 'F capacity from odd part')
        available = rank//2
    elif s == 2:
        name = 'H'
        rank = v2(431*y-101)
        check(rank == 2+v2(431*P**L*m-117443449) and rank >= 3, 'H capacity from odd part')
        available = rank//9
    else:
        check(y % 8 == 7, 'unsupported exit shell modulo eight')
        return {**rec, 'certified': False, 'reason': 'no E, F, or H reset in this exit shell'}
    w = RESETS[name]
    need = required(w, x, y)
    if available < need:
        return {**rec, 'certified': False, 'reset': name, 'available_blocks': available,
                'required_blocks': need, 'reason': 'insufficient lawful reset capacity'}
    z = iterate(w, y, need)
    check(z.denominator == 1 and 0 < z < x, 'actual positive smaller reset endpoint')
    return {**rec, 'certified': True, 'reset': name, 'available_blocks': available,
            'reset_blocks': need, 'endpoint': int(z), 'steps': 22*L+w['n']*need,
            'exact_first_stopping_time': name in ('E', 'F')}


def composite(L, w, j):
    PL, QL, CL, _, _ = run_data(L)
    a, q, b = w['a'], w['q'], w['b']
    reset_carry = b*(q**j-a**j)//(q-a)
    slope_num, modulus = a**j*PL, q**j*QL
    carry = a**j*CL+QL*reset_carry
    residue = (-carry*pow(slope_num, -1, modulus)) % modulus
    return slope_num, modulus, carry, residue


# Calibrate the original J: its odd upper source already descends at step 20.
J_upper = walk(1817599, 20)
check(J_upper['odds'] == 12 and J_upper['carry'] == 532465, 'old upper-prefix coefficients')
check(F(J_upper['carry'], 2**20-3**12) == F(106493, 103427) < 2, 'old upper prefix is descending')
check(T(1817598) < 1817598, 'old lower J source descends immediately')

# New R survives every prefix, with a separately guarded phase companion.
check(coefficients(R_BITS) == (22, 14, C, 908799), 'first-descent-safe R word')
prefixes = []
for n in range(1, 23):
    _, odds, carry, _ = coefficients(R_BITS[:n])
    check(3**odds > 2**n and carry >= 0, 'every nonempty R prefix expands')
    prefixes.append(dict(length=n, odds=odds, multiplier_numerator=3**odds,
                         multiplier_denominator=2**n))
phase_anchor = walk(3*908799+2, 22)
check(phase_anchor['odds'] == 14 and phase_anchor['carry'] == 3*C+2*Q-2*P,
      'companion affine relation identity')
for lift in (0, 1, 17):
    x = 908799+Q*lift
    lo, hi = walk(x, 22), walk(3*x+2, 22)
    check(lo['bits'] == R_BITS and hi['bits'] == phase_anchor['bits'], 'separate phase source guards')
    check(hi['endpoint'] == 3*lo['endpoint']+2, 'actual relation returns to (3u+2,u)')
    check(lo['first_descent'] is None, 'odd driver has no prior descent')
check(P-Q == 5*d and C == 5*B, 'reduced center constants')
check(431*B+101*d == 4*117443449, 'H shell capacity constant')
for w in RESETS.values():
    gap = w['center']+F(B, d)
    check(v2(gap.numerator) == w['shell'], 'exact reset-center separation')

# All-source automatic low-L comparisons, proved with exact rationals.
check(8*P > 9*Q and 7*P < 8*Q, 'simple multiplier bracket')
check(100000*B < d*908799, 'small correction at least R source')
check(100001*8**5 < 200000*7**5 and 9**6 > 2*8**6, 'sharp single-halving bounds')
_, q2, c2, _, _ = run_data(2)
y2 = F(P**2*908799+c2, q2)
check(iterate(RESETS['F'], y2, 1) < 908799, 'one F reset sufficient through L two')
check(3*P**3 > 4*Q**3, 'one F reset insufficient from L three')

families, sources, whole_cylinders = [], [], []
for L in (1, 2, 3, 5, 6, 10, 14, 15, 16, 32, 64):
    PL, QL, CL, rL, yL = run_data(L)
    check(CL == sum(C*P**(L-1-i)*Q**i for i in range(L)), 'geometric carry identity')
    check(v2(CL) == 0 and 0 < rL < QL and rL % Q == 908799, 'least expanding cylinder source')
    check(coefficients(R_BITS*L) == (22*L, 14*L, CL, rL), 'whole expanding word cylinder')
    for name, w in RESETS.items():
        j = 1
        while w['a']**j*PL >= w['q']**j*QL:
            j += 1
        slope_num, modulus, carry, residue = composite(L, w, j)
        delta = modulus-slope_num
        first_lift = max(0, (carry-delta*residue)//(delta*modulus)+1)
        cq = F(B, d)+w['center']
        head_bound = -(-cq.numerator//cq.denominator)
        check(head_bound == {'E': 10, 'F': 11, 'H': 10}[name], 'uniform exceptional-head constant')
        check(0 < carry < cq*modulus and delta >= 1, 'carry versus positive integer gap')
        check(first_lift <= head_bound, 'uniform tail cutoff')
        check(delta*(residue+modulus*first_lift) > carry, 'first lift passes strict descent threshold')
        if first_lift:
            check(delta*(residue+modulus*(first_lift-1)) <= carry, 'previous lift fails exact threshold')
        bits = R_BITS*L+w['bits']*j
        n, odds, cc, rr = coefficients(bits)
        check((3**odds, 2**n, cc, rr) == (slope_num, modulus, carry, residue), 'complete optimal reset cylinder')
        j_all = required(w, rL, yL)
        check(j <= j_all <= j+1, 'whole-cylinder reset within one block of optimal')
        families.append(dict(L=L, reset=name, optimal_blocks=j, modulus=modulus,
                             residue=residue, first_certified_lift=first_lift,
                             uniform_head_bound=head_bound, full_cylinder_blocks=j_all,
                             source_rank=22*L+w['shell'], receipt_steps=n,
                             exact_first_stopping_time=name in ('E', 'F')))
        for lift in sorted({first_lift, first_lift+1, head_bound+7}):
            x = residue+modulus*lift
            rec = compress(x)
            check((rec['L'], rec['shell']) == (L, w['shell']), 'source guard gives exact maximal run')
            rising = walk(x, 22*L)
            check(rising['bits'] == R_BITS*L and rising['endpoint'] == rec['exit'], 'actual compressed rising run')
            check(rising['first_descent'] is None and rising['endpoint'] > x, 'no prefix descent before reset')
            check(QL <= d*x+B and F(d*rec['exit']+B, d*x+B) == F(PL, QL), 'source-wise endpoint-height identities')
            for i in sorted({0, 1, L-1}):
                if i >= L:
                    continue
                ci = C*(P**i-Q**i)//(P-Q)
                xx = (P**i*x+ci)//Q**i
                rem = compress(xx)
                check((rem['L'], rem['shell'], rem['exit']) ==
                      (L-i, w['shell'], rec['exit']), 'exit invariant along maximal run')
            actual = walk(x, n)
            check(actual['bits'] == bits and actual['endpoint'] < x, 'literal optimal reset receipt')
            check(actual['endpoint']*modulus == slope_num*x+carry, 'exact composite endpoint')
            if name in ('E', 'F'):
                check(actual['first_descent'] == n, 'exact first stopping time at optimal reset')
            else:
                check(actual['first_descent'] is not None and actual['first_descent'] <= n, 'H receipt upper bound on stopping time')
            dec = decode(x)
            check(dec['certified'] and dec['reset'] == name and dec['reset_blocks'] == j, 'source decoder realizes optimal reset')
            check(dec['endpoint'] == actual['endpoint'] and dec['steps'] == n, 'decoded clock and endpoint')
            sources.append(dict(L=L, reset=name, lift=lift, source=x, exit=rec['exit'],
                                endpoint=actual['endpoint'], receipt_steps=n,
                                first_stopping_time=actual['first_descent']))
        aa, mm, ca, ra = composite(L, w, j_all)
        for lift in (0, 1):
            x = ra+mm*lift
            rec = compress(x)
            check((rec['L'], rec['shell']) == (L, w['shell']), 'full-cylinder maximal-run guard')
            z = iterate(w, rec['exit'], j_all)
            check(z.denominator == 1 and 0 < z < x, 'all-positive-cylinder example descends')
        whole_cylinders.append(dict(L=L, reset=name, blocks=j_all, residue=ra, modulus=mm))


def source_for_exit(L, target, modulus):
    PL, QL, CL, rL, yL = run_data(L)
    lift = ((target-yL)*pow(PL, -1, modulus)) % modulus
    return rL+QL*lift


# Exact source-selected success and failure controls, not a residue sweep.
exit_controls = []
for L in (1, 2, 5, 6, 10):
    x = source_for_exit(L, 2, 4)  # exactly one even step
    dec = decode(x)
    check(dec['L'] == L and dec['shell'] == 0 and dec['available_blocks'] == 1,
          'prescribed one-halving exit capacity')
    check(dec['certified'] == (L <= 5), 'sharp automatic one-halving dichotomy')
    exit_controls.append(dec)
for L in (1, 2, 3, 8):
    x = source_for_exit(L, 5, 8)  # F rank exactly two
    dec = decode(x)
    check(dec['L'] == L and dec['shell'] == 1 and dec['available_blocks'] == 1,
          'prescribed one-F exit capacity')
    check(dec['certified'] == (L <= 2), 'sharp automatic one-F dichotomy')
    exit_controls.append(dec)
for L in (1, 4, 16):
    for target, modulus, shell in ((3, 16, 2), (7, 16, 3)):
        x = source_for_exit(L, target, modulus)
        dec = decode(x)
        check(dec['L'] == L and dec['shell'] == shell and not dec['certified'],
              'insufficient or absent reset reported explicitly')
        exit_controls.append(dec)


def encode(value):
    if isinstance(value, F):
        return str(value)
    raise TypeError(type(value).__name__)


print(json.dumps(dict(
    status='PROVED first-descent-safe maximal-run resets; universal source entry OPEN',
    checks=CHECKS,
    constants=dict(P=P, Q=Q, C=C, d=d, B=B, center=str(-F(B, d))),
    R_word=R_BITS, companion_word=phase_anchor['bits'],
    positive_prefix_certificate=prefixes,
    optimal_reset_families=families,
    literal_source_receipts=sources,
    all_positive_lift_cylinders=whole_cylinders,
    exact_exit_success_and_failure_controls=exit_controls,
), indent=2, sort_keys=True, default=encode))
