#!/usr/bin/env python3
"""Exact odd-lower expanding adjacent return and uniformly bounded-head resets.

No finite search for return words is included or needed. The all-L statements
are proved in the companion paper; these are independent exact replay controls.
"""
from fractions import Fraction as F
import json

CHECKS = 0
LOWER = '110110110110110101011011101010'
UPPER = '001111111111111101010110000100'
P, Q, C = 1162261467, 1073741824, 4912087879
D, RHO = P-Q, 815464443
CENTER_SIZE = F(C, D)


def check(test, label):
    global CHECKS
    CHECKS += 1
    if not test:
        raise RuntimeError(label)


def v2(n):
    if not n:
        raise ValueError('infinite valuation')
    n = abs(n)
    return (n & -n).bit_length()-1


def ceil(r):
    return -(-r.numerator//r.denominator)


def coefficients(bits):
    a = c = 0
    for t, bit in enumerate(bits):
        if bit == '1':
            a += 1
            c = 3*c+2**t
        elif bit != '0':
            raise ValueError('binary parity word required')
    p, q = 3**a, 2**len(bits)
    return p, q, c, (-c*pow(p, -1, q)) % q


def branch(bits, upper=None):
    a, q, b, r = coefficients(bits)
    out = dict(bits=bits, a=a, q=q, b=b, residue=r,
               center=F(b, q-a), length=len(bits))
    if upper is not None:
        au, qu, bu, ru = coefficients(upper)
        check((au, qu, bu, ru) == (a, q, b+q-a, (r+1) % q),
              'separately guarded adjacent companion')
        out['upper_bits'] = upper
    return out


RESETS = {
    'F': branch('10', '01'),
    'H': branch('110101000', '001101001'),
}
RESETS['F']['shell'] = 1
RESETS['H']['shell'] = 4


def walk(x, bits=None, length=None):
    check(type(x) is int and x > 0, 'positive integer trajectory source')
    if bits is not None:
        length = len(bits)
    source, z = x, x
    actual = []
    first = None
    for i in range(length):
        bit = z % 2
        actual.append(str(bit))
        if bits is not None:
            check(bit == int(bits[i]), 'literal source parity guard at every step')
        z = (3*z+1)//2 if bit else z//2
        check(type(z) is int and z > 0, 'literal positive integer trajectory')
        if first is None and z < source:
            first = i+1
    return dict(endpoint=z, first_stopping_time=first, bits=''.join(actual))


def repeated(p, q, c, L):
    if L == 0:
        return 1, 1, 0
    return p**L, q**L, c*(p**L-q**L)//(p-q)


def run_data(L):
    pl, ql, cl = repeated(P, Q, C, L)
    rho = (-cl*pow(pl, -1, ql)) % ql
    return pl, ql, cl, rho


def minimal_reset(p, q, L, w):
    j = 1
    while w['a']**j*p**L >= w['q']**j*q**L:
        j += 1
    return j


def composite(p, q, c, L, w, j):
    pl, ql, cl = repeated(p, q, c, L)
    aj, qj = w['a']**j, w['q']**j
    bj = w['b']*(qj-aj)//(w['q']-w['a'])
    aa, mm, cc = aj*pl, qj*ql, aj*cl+ql*bj
    rho = (-cc*pow(aa, -1, mm)) % mm
    delta = mm-aa
    if delta <= 0:
        raise ValueError('strictly contracting composite required')
    first = max(0, (cc-delta*rho)//(delta*mm)+1)
    center_sum = F(c, p-q)+w['center']
    check(cc == F(c, p-q)*aa+w['center']*mm-center_sum*aj*ql,
          'general centered composite carry identity')
    check(0 < cc < center_sum*mm, 'general strict carry bound')
    check(first <= ceil(center_sum), 'general uniform initial-lift bound')
    check(delta*(rho+first*mm) > cc, 'first successful lift satisfies strict descent')
    if first:
        check(delta*(rho+(first-1)*mm) <= cc, 'preceding lift is excluded')
    return dict(A=aa, M=mm, C=cc, residue=rho, delta=delta,
                first_certified_lift=first, uniform_head_bound=ceil(center_sum))


def compress(x):
    check(type(x) is int and x > 0, 'positive maximal-run source')
    rank = v2(D*x+C)
    L, s = divmod(rank, 30)
    m = (D*x+C) >> rank
    numerator = P**L*2**s*m-C
    check(numerator % D == 0, 'compressed exit integrality')
    y = numerator//D
    check(y > 0 and m % 2 == 1 and v2(D*y+C) == s, 'exact positive exit and residual shell')
    return dict(L=L, shell=s, rank=rank, odd_part=m, exit=y)


def decode(x):
    rec = compress(x)
    L, s, m, y = rec['L'], rec['shell'], rec['odd_part'], rec['exit']
    if not L:
        return {**rec, 'certified': False, 'reason': 'no initial K block'}
    if s == 1:
        name = 'F'
        rank = v2(y-1)
        check(rank == 1+v2(P**L*m-2500303761), 'F source-sensitive capacity identity')
    elif s == 4:
        name = 'H'
        rank = v2(431*y-101)
        check(rank == 4+v2(431*P**L*m-132878147487), 'H source-sensitive capacity identity')
    else:
        return {**rec, 'certified': False, 'reason': 'neither adjacent reset has this exit shell'}
    w = RESETS[name]
    available = rank//w['length']
    need = 1
    p = w['center']
    while F(w['a'], w['q'])**need*(y-p) >= x-p:
        need += 1
    if available < need:
        return {**rec, 'reset': name, 'available_blocks': available,
                'required_blocks': need, 'certified': False,
                'reason': 'insufficient lawful adjacent-reset capacity'}
    z = p+F(w['a'], w['q'])**need*(y-p)
    check(z.denominator == 1 and 0 < z < x, 'decoder smaller positive integer endpoint')
    return {**rec, 'reset': name, 'available_blocks': available,
            'reset_blocks': need, 'endpoint': int(z), 'certified': True,
            'steps': 30*L+w['length']*need,
            'exact_first_stopping_time': name == 'F'}


# Fixed word, symbolic guards, and all 30 prefix multiplier comparisons.
check((P, Q) == (3**19, 2**30), 'K powers')
check(D == 88519643, 'K center denominator')
check(coefficients(LOWER) == (P, Q, C, RHO), 'K lower word coefficients and exact cylinder')
check(coefficients(UPPER) == (P, Q, C+Q-P, RHO+1), 'K upper word coefficients and independent cylinder')
check(C+Q-P == 4823568236, 'K upper carry')
prefixes = []
for t in range(1, 31):
    aa, mm, cc, _ = coefficients(LOWER[:t])
    check(aa > mm and cc >= 0, 'every nonempty lower prefix is strictly expanding')
    prefixes.append(dict(length=t, odds=LOWER[:t].count('1'), A=aa, M=mm, C=cc))
check(min(F(r['A'], r['M']) for r in prefixes) == F(531441, 524288),
      'smallest prefix multiplier')
anchor_controls = []
for lift in (0, 1, 17, 104729):
    x = RHO+Q*lift
    lo, hi = walk(x, LOWER), walk(x+1, UPPER)
    check(lo['endpoint']*Q == P*x+C, 'lower affine anchor identity')
    check(hi['endpoint'] == lo['endpoint']+1, 'oriented adjacent return')
    check(lo['endpoint'] > x and lo['first_stopping_time'] is None, 'odd lower screen counterexample')
    check(hi['first_stopping_time'] == 1, 'upper endpoint is separately screened')
    anchor_controls.append(dict(source=x, endpoint=lo['endpoint'],
                                lower_first_stopping_time=lo['first_stopping_time'],
                                upper_first_stopping_time=hi['first_stopping_time']))
check(anchor_controls[0]['endpoint'] == 882691615, 'K anchor endpoint')
for name, w in RESETS.items():
    gap = CENTER_SIZE+w['center']
    check(v2(gap.numerator) == w['shell'] and gap.denominator % 2 == 1,
          'exact center separation for maximal K run')
check(C+D == 2*2500303761, 'F capacity constant')
check(431*C+101*D == 16*132878147487, 'H capacity constant')

# The general ceiling bound is attained, not merely approached.
odd_step = coefficients('1')
even_step = branch('0')
sharp = composite(*odd_step[:3], 1, even_step, 1)
check((sharp['A'], sharp['M'], sharp['C'], sharp['residue'],
       sharp['first_certified_lift'], sharp['uniform_head_bound']) == (3, 4, 1, 1, 1, 1),
      'at-most ceiling bound is sharp at U=odd step, W=even step')
check(walk(1, '10')['endpoint'] == 1 and walk(5, '10')['endpoint'] == 4,
      'strict inequality excludes the first lift in sharpness example')

families, receipts, whole_cylinders = [], [], []
for L in (1, 2, 3, 4, 7, 10, 23, 24, 29, 64, 128):
    pl, ql, cl, rl = run_data(L)
    check(coefficients(LOWER*L) == (pl, ql, cl, rl), 'complete K run word and cylinder')
    check(coefficients(UPPER*L) == (pl, ql, cl+ql-pl, (rl+1) % ql),
          'complete upper K run independently guarded')
    check(0 < rl < ql and rl % Q == RHO, 'least positive K run source')
    yl = F(pl*rl+cl, ql)
    for name, w in RESETS.items():
        j = minimal_reset(P, Q, L, w)
        comp = composite(P, Q, C, L, w, j)
        check(comp['uniform_head_bound'] == {'F': 57, 'H': 56}[name], 'K-specific uniform head bound')
        lower_bits, upper_bits = LOWER*L+w['bits']*j, UPPER*L+w['upper_bits']*j
        aa, mm, cc, rr = coefficients(lower_bits)
        check((aa, mm, cc, rr) == (comp['A'], comp['M'], comp['C'], comp['residue']),
              'lower composite complete source cylinder')
        check(coefficients(upper_bits) == (aa, mm, cc+mm-aa, (rr+1) % mm),
              'upper composite complete source cylinder independently checked')
        check(w['a']**j*pl < w['q']**j*ql, 'minimal reset contracts')
        if j > 1:
            check(w['a']**(j-1)*pl >= w['q']**(j-1)*ql, 'previous full reset cannot contract')
        for lift in (comp['first_certified_lift'], comp['first_certified_lift']+1,
                     comp['uniform_head_bound']+7):
            x = rr+mm*lift
            rec = compress(x)
            check((rec['L'], rec['shell']) == (L, w['shell']), 'exact maximal K run and exit shell')
            check((pl*x+cl)//ql == rec['exit'], 'compressed run endpoint matches affine run')
            for i in sorted({0, 1, L-1}):
                if i >= L:
                    continue
                pi, qi, ci = repeated(P, Q, C, i)
                xx = (pi*x+ci)//qi
                rem = compress(xx)
                check((rem['L'], rem['shell'], rem['exit']) == (L-i, w['shell'], rec['exit']),
                      'maximal run exit and shell invariant')
            lower, upper = walk(x, lower_bits), walk(x+1, upper_bits)
            check(0 < lower['endpoint'] < x, 'literal positive smaller-source receipt')
            check(upper['endpoint'] == lower['endpoint']+1 < x+1, 'paired receipt retains orientation and contracts both endpoints')
            check(lower['endpoint']*mm == aa*x+cc, 'literal receipt agrees with composite affine endpoint')
            if name == 'F':
                check(lower['first_stopping_time'] == len(lower_bits), 'exact F first stopping time')
            else:
                check(30*L < lower['first_stopping_time'] <= len(lower_bits), 'H stopping time receipt bounds')
            check(upper['first_stopping_time'] == 1, 'no false upper survival claim')
            dec = decode(x)
            check(dec['certified'] and dec['reset'] == name and dec['reset_blocks'] == j,
                  'source-sensitive decoder obtains the minimally contracting reset')
            check(dec['endpoint'] == lower['endpoint'] and dec['steps'] == len(lower_bits),
                  'source decoder endpoint and receipt clock')
            receipts.append(dict(L=L, reset=name, reset_blocks=j, lift=lift, source=x,
                                 expanding_exit=rec['exit'], endpoint=lower['endpoint'],
                                 upper_source=x+1, upper_endpoint=upper['endpoint'],
                                 receipt_steps=len(lower_bits), first_stopping_time=lower['first_stopping_time']))
        # An all-lift refinement: evaluate the reset requirement at the least K-run source.
        j_all = 1
        while F(w['a'], w['q'])**j_all*(yl-w['center']) >= rl-w['center']:
            j_all += 1
        check(rl > 3*CENTER_SIZE+4*w['center'], 'uniform source lower bound for all-lift refinement')
        check(j <= j_all <= j+1, 'at most one extra reset block certifies every positive lift')
        all_comp = composite(P, Q, C, L, w, j_all)
        check(all_comp['first_certified_lift'] == 0, 'all-lift composite has no excluded positive source')
        whole_cylinders.append(dict(L=L, reset=name, blocks=j_all,
                                    residue=all_comp['residue'], modulus=all_comp['M']))
        families.append(dict(L=L, reset=name, minimal_blocks=j,
                             receipt_steps=len(lower_bits), source_rank=30*L+w['shell'],
                             exact_first_stopping_time=name == 'F', **comp))


def source_for_exit(L, target, modulus):
    pl, ql, cl, rl = run_data(L)
    yl = (pl*rl+cl)//ql
    lift = ((target-yl)*pow(pl, -1, modulus)) % modulus
    return rl+ql*lift


# Source-derived failures: an exact residual shell is not a universal reset policy.
exit_controls = []
for L in (1, 3, 4, 10):
    x = source_for_exit(L, 5, 8)  # exactly one F block is lawful
    rec = decode(x)
    check((rec['L'], rec['shell'], rec['available_blocks']) == (L, 1, 1), 'one-F exit capacity guard')
    check(rec['certified'] == (L <= 3), 'one-F success and honest capacity failure')
    exit_controls.append(dict(source=x, **rec))
for L in (1, 23, 24, 64):
    # H center numerator chosen modulo 2^10 to have v2(431y-101)=9.
    y = ((101+512)*pow(431, -1, 1024)) % 1024
    x = source_for_exit(L, y, 1024)
    rec = decode(x)
    check((rec['L'], rec['shell'], rec['available_blocks']) == (L, 4, 1), 'one-H exit capacity guard')
    check(rec['certified'] == (L <= 23), 'one-H success and honest capacity failure')
    exit_controls.append(dict(source=x, **rec))
for L in (1, 7, 32):
    for target, modulus, shell in ((2, 4, 0), (7, 8, 2), (43+32, 64, 4)):
        x = source_for_exit(L, target, modulus)
        rec = decode(x)
        check((rec['L'], rec['shell']) == (L, shell), 'prescribed unsupported or short-capacity exit shell')
        check(not rec['certified'], 'uncertified source is not declared descending')
        exit_controls.append(dict(source=x, **rec))

print(json.dumps(dict(
    status='PROVED odd-lower expanding adjacent return and all-L paired reset cylinders; universal entry OPEN',
    checks=CHECKS,
    constants=dict(P=P, Q=Q, C=C, D=D, residue=RHO, center=str(-CENTER_SIZE)),
    lower_word=LOWER, upper_word=UPPER,
    lower_prefix_certificate=prefixes,
    odd_lower_screen_counterexamples=anchor_controls,
    sharp_general_ceiling_example=sharp,
    optimal_adjacent_reset_families=families,
    literal_paired_receipts=receipts,
    all_positive_lift_cylinders=whole_cylinders,
    exact_exit_controls=exit_controls,
), indent=2, sort_keys=True))
