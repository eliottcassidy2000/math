#!/usr/bin/env python3
"""Typed Zagier recurrence parameters, singular-location carriers and golden data.

Standard library only. The finite replay is not a classification of integral
recurrences. All public receipt validators retain the original integer triple.
Run with --out PATH to write the same deterministic JSON printed on stdout.
"""
import argparse
from dataclasses import dataclass
from fractions import Fraction as Q
from itertools import product
from math import comb
from pathlib import Path
import json


def integer(x):
    if type(x) is not int:
        raise TypeError('an exact integer is required')
    return x


def rational(x):
    if type(x) not in (int, Q):
        raise TypeError('an exact rational is required')
    return Q(x)


@dataclass(frozen=True)
class Recurrence:
    a: int
    b: int
    c: int

    def __post_init__(self):
        for x in (self.a, self.b, self.c):
            integer(x)


def source(t):
    if type(t) is not Recurrence:
        raise TypeError('a source Recurrence is required')
    t.__post_init__()
    return t


CANON = {
    'A': Recurrence(7, 2, -8), 'B': Recurrence(9, 3, 27),
    'C': Recurrence(10, 3, 9), 'D': Recurrence(11, 3, -1),
    'E': Recurrence(12, 4, 32), 'F': Recurrence(17, 6, 72),
}
SUPPLIED = (CANON['B'], CANON['D'], CANON['C'], Recurrence(11, 5, 125),
            CANON['E'], Recurrence(17, 5, 72))


def coefficients(t, last):
    source(t)
    integer(last)
    if last < 0:
        raise ValueError('nonnegative final index required')
    u = [Q(1)]
    for n in range(last):
        old = t.c*n*n*u[n-1] if n else 0
        u.append(((t.a*n*(n+1)+t.b)*u[n]-old)/(n+1)**2)
    return tuple(u)


@dataclass(frozen=True)
class Prefix:
    triple: Recurrence
    values: tuple


def audit_prefix(t, packet):
    source(t)
    if type(packet) is not Prefix or type(packet.values) is not tuple:
        raise TypeError('an exact Prefix packet is required')
    source(packet.triple)
    if packet.triple != t or not packet.values:
        raise ValueError('source mismatch or empty prefix')
    for x in packet.values:
        rational(x)
    if tuple(map(Q, packet.values)) != coefficients(t, len(packet.values)-1):
        raise ValueError('incorrect recurrence prefix')
    return True


@dataclass(frozen=True)
class BranchCarrier:
    triple: Recurrence
    root_discriminant: int
    elliptic_discriminant: int
    c4: int
    j: Q


def branch_carrier(t):
    source(t)
    d = t.a*t.a-4*t.c
    if not t.c or not d:
        raise ValueError('the singular-location double cover is not smooth')
    return BranchCarrier(t, d, 16*t.c*t.c*d, 16*(t.a*t.a-3*t.c),
                         Q(256*(t.a*t.a-3*t.c)**3, t.c*t.c*d))


def audit_carrier(t, packet):
    source(t)
    if type(packet) is not BranchCarrier:
        raise TypeError('BranchCarrier required')
    source(packet.triple)
    for v in (packet.root_discriminant, packet.elliptic_discriminant, packet.c4):
        integer(v)
    if type(packet.j) is not Q:
        raise TypeError('j must be an exact Fraction')
    if packet != branch_carrier(t):
        raise ValueError('incorrect carrier or wrong original source')
    return True


def scale(t, s):
    source(t)
    integer(s)
    if not s:
        raise ValueError('nonzero integer scale required')
    return Recurrence(s*t.a, s*t.b, s*s*t.c)


def decode_jet(b, u2, u3):
    b, u2, u3 = map(rational, (b, u2, u3))
    den = 6*u2-8*b*b
    if not den:
        raise ValueError('degenerate jet: the three terms do not determine a')
    a = (9*u3-17*b*u2+4*b**3)/den
    c = (2*a+b)*b-4*u2
    if any(x.denominator != 1 for x in (a, b, c)):
        raise ValueError('jet has no integer triple in this nondegenerate chart')
    return Recurrence(int(a), int(b), int(c))


# Polynomial arithmetic, ascending coefficients, exact rational coefficient field.
def poly(p):
    if type(p) is not tuple or not p:
        raise TypeError('nonempty exact coefficient tuple required')
    p = tuple(map(rational, p))
    while len(p) > 1 and p[-1] == 0:
        p = p[:-1]
    return p


def add(*ps):
    ps = tuple(map(poly, ps))
    return poly(tuple(sum(p[i] if i < len(p) else 0 for p in ps)
                      for i in range(max(map(len, ps)))))


def mul(p, q):
    p, q = poly(p), poly(q)
    out = [Q(0)]*(len(p)+len(q)-1)
    for i, x in enumerate(p):
        for j, y in enumerate(q):
            out[i+j] += x*y
    return poly(tuple(out))


def smul(s, p):
    return poly(tuple(rational(s)*v for v in poly(p)))


def diff(p):
    p = poly(p)
    return poly(tuple(i*p[i] for i in range(1, len(p))) or (0,))


def evaluate(p, x):
    x = rational(x)
    v = Q(0)
    for c in reversed(poly(p)):
        v = v*x+c
    return v


def gauge_numerators(t, new_b, numerator, denominator):
    """Coefficients of f',f after clearing Q^2,Q^3 respectively.

    No claim about a selected single solution: this is the entire differential
    operator intertwining test at the unchanged independent variable.
    """
    source(t)
    integer(new_b)
    P, D = poly(numerator), poly(denominator)
    if P == (0,) or D == (0,):
        raise ValueError('nonzero rational scalar required')
    p = (0, 1, -t.a, t.c)
    W = add(mul(diff(P), D), smul(-1, mul(P, diff(D))))
    R2 = add(mul(diff(diff(P)), mul(D, D)),
             smul(-1, mul(P, mul(D, diff(diff(D))))),
             smul(-2, mul(diff(P), mul(D, diff(D)))),
             smul(2, mul(P, mul(diff(D), diff(D)))))
    first = smul(2, mul(p, W))
    second = add(mul(p, R2), mul(diff(p), mul(W, D)),
                 smul(t.b-new_b, mul(P, mul(D, D))))
    return first, second


def fib_lucas(n):
    integer(n)
    if n < 0:
        raise ValueError('nonnegative index required')
    f, g = 0, 1
    for _ in range(n):
        f, g = g, f+g
    return f, 2*g-f


def signature(n):
    integer(n)
    if n < 1:
        raise ValueError('positive power index required')
    f, l = fib_lucas(n)
    return (l, f, 5*f*f)


def recover_c_from_signature(a, root_coefficient, discriminant):
    for x in (a, root_coefficient, discriminant):
        integer(x)
    if root_coefficient <= 0 or discriminant != 5*root_coefficient**2:
        raise ValueError('not the stated positive golden-root signature')
    c = Q(a*a-discriminant, 4)
    if c.denominator != 1:
        raise ValueError('signature is not integral in this chart')
    return int(c)  # Deliberately cannot return an accessory b.


def divisor_count(n):
    integer(n)
    if n < 1:
        raise ValueError('positive integer required')
    return sum(n % d == 0 for d in range(1, n+1))


def record_numbers(bound):
    integer(bound)
    if bound < 1:
        raise ValueError('positive bound required')
    best, out = 0, []
    for n in range(1, bound+1):
        d = divisor_count(n)
        if d > best:
            best = d
            out.append((n, d))
    return out


ATLAS = [
    (-1, 'D: signed c; odd golden norm', 'abc_difference_families_20261011c.md',
     'L5^2-5F5^2=-4 determines c=-1, but not accessory b'),
    (3, 'B,C,D: accessory b; B also has discriminant -27', 'lucas_block_cm_detector_sporadic_ad_hauptmoduln_20261018.md',
     'accessory coefficient and bad-prime/splitting data are different roles'),
    (4, 'E: accessory b and one reciprocal singular slope', 'divisor_equation_recursions_20261011f.md',
     'also a divisor-record input; no map from divisor records to an operator'),
    (5, 'D: square-root field and F5; supplied accessory mistakes', 'abc_difference_families_20261011c.md',
     'the square-root coefficient is not the recurrence accessory'),
    (9, 'B: a; C: c and reciprocal singular slope', 'digit_idoneal_difference_20261018b.md',
     'the C9 digit tournament is a different marked object'),
    (10, 'C: a', 'paley_next_cores_20261014i.md',
     'ten-edge core count does not supply a recurrence transformation'),
    (11, 'D: a=L5', 'jacobian_elliptic_bridge_20261016c.md',
     'golden singular splitting is Q(sqrt5), not the HN field Q(sqrt(-11))'),
    (12, 'E: a', 'divisor_equation_recursions_20261011f.md',
     'strict divisor record12 is distinct from twelve-valued torsion markings'),
    (17, 'F: a', 'fermat_wieferich_order_20261014e.md',
     'Fermat prime17 gives no accessory b; canonical F requires6'),
    (27, 'B: c=3^3, discriminant=-27', 'wheel_apery_extension_20261014h.md',
     'prime-power field size and recurrence coefficient require separate inputs'),
    (32, 'E: c=2^5', 'divisor_equation_recursions_20261011f.md',
     'a prime power need not be a strict divisor record'),
    (72, 'F: c=8*9', 'lucas_block_cm_detector_sporadic_ad_hauptmoduln_20261018.md',
     'the two slopes8,9 survive accessory perturbation b6 to b5'),
    (125, 'D: discriminant=5^3; supplied alternative c', 'abc_difference_families_20261011c.md',
     'putting a discriminant in the c slot changes the singular set'),
]


CHECKS = 0


def check(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise AssertionError(message)


def rejected(f):
    try:
        f()
    except (ValueError, TypeError):
        check(True, 'rejected')
        return
    raise AssertionError('malformed input was accepted')


def verify():
    global CHECKS
    CHECKS = 0
    rows = []
    for label, t in CANON.items():
        u = coefficients(t, 24)
        check(all(x.denominator == 1 for x in u), ('canonical finite integrality', label))
        check(audit_prefix(t, Prefix(t, u)), 'prefix')
        check(decode_jet(*u[1:4]) == t, 'jet inverse')
        check(audit_carrier(t, branch_carrier(t)), 'carrier audit')
        for s in (-3, -1, 1, 2, 4):
            target = scale(t, s)
            v = coefficients(target, 24)
            for n in range(25):
                check(v[n] == s**n*u[n], 'coefficient scaling')
            check(branch_carrier(target).j == branch_carrier(t).j, 'j scaling')
            check(branch_carrier(target).root_discriminant == s*s*branch_carrier(t).root_discriminant,
                  'root discriminant scaling')
        bc = branch_carrier(t)
        rows.append({'name': label, 'triple': [t.a,t.b,t.c],
                     'coefficients_0_to_8': [str(x) for x in u[:9]],
                     'root_discriminant': bc.root_discriminant,
                     'branch_carrier_j': str(bc.j)})

    # Independent classical finite binomial sums, not the recurrence producer.
    for n in range(17):
        check(coefficients(CANON['A'], n)[n] == sum(comb(n,k)**3 for k in range(n+1)), 'Franel')
        check(coefficients(CANON['C'], n)[n] == sum(comb(n,k)**2*comb(2*k,k) for k in range(n+1)), 'C sum')
        check(coefficients(CANON['D'], n)[n] == sum(comb(n,k)**2*comb(n+k,k) for k in range(n+1)), 'Apery zeta2')
        check(coefficients(CANON['E'], n)[n] == sum(comb(n,k)*comb(2*k,k)*comb(2*n-2*k,n-k) for k in range(n+1)), 'E sum')

    # Independent coefficient checks of the operator and the gauge commutator.
    jets_decoded = 0
    for a,b,c in product(range(-3,6), repeat=3):
        t = Recurrence(a,b,c)
        u = coefficients(t, 8)
        p = (0,1,-a,c)
        f = u
        Lf = add(diff(mul(p,diff(f))), mul((-b,c),f))
        for n in range(8):
            check((Lf[n] if n < len(Lf) else 0) == 0, 'formal differential coefficients')
        if 6*u[2]-8*b*b:
            check(decode_jet(*u[1:4]) == t, 'complete nondegenerate jet inverse')
            jets_decoded += 1
        else:
            rejected(lambda u=u: decode_jet(*u[1:4]))
        if c*(a*a-4*c):
            bc = branch_carrier(t)
            check(bc.c4**3 == bc.j*bc.elliptic_discriminant, 'Weierstrass invariant')
            check((bc.j == 0) == (a*a == 3*c), 'equianharmonic iff')
            check((bc.j == 1728) == (a == 0 or 2*a*a == 9*c), 'harmonic iff')

    for t in (*CANON.values(), SUPPLIED[3], SUPPLIED[5]):
        for new_b in (t.b-1,t.b,t.b+1):
            for P,D in (((1,),(1,)), ((2,2),(1,1)), ((1,1),(1,)),
                        ((1,),(1,1)), ((2,3,1),(3,-1,1))):
                first,second = gauge_numerators(t,new_b,P,D)
                W = add(mul(diff(P),D),smul(-1,mul(P,diff(D))))
                check((first == (0,) and second == (0,)) == (W == (0,) and new_b == t.b),
                      'rational scalar rigidity')
                # Check two polynomials f directly without symbolic differentiation of R.
                for f in ((1,), (0,1), (1,2,3)):
                    p = (0,1,-t.a,t.c)
                    Pf = mul(P,f)
                    w = add(mul(diff(Pf),D),smul(-1,mul(Pf,diff(D))))
                    r2 = add(mul(diff(w),D),smul(-2,mul(w,diff(D))))
                    direct = add(mul(p,r2),mul(diff(p),mul(w,D)),
                                 mul((-new_b,t.c),mul(Pf,mul(D,D))),
                                 smul(-1,mul(P,mul(add(diff(mul(p,diff(f))),mul((-t.b,t.c),f)),mul(D,D)))))
                    rhs = add(mul(first,mul(D,diff(f))),mul(second,f))
                    check(direct == rhs, 'independent cleared commutator')

    golden = []
    for n in range(1,25):
        F,L = fib_lucas(n)
        check(L*L-5*F*F == 4*(-1)**n, 'signed golden norm')
        sig = signature(n)
        check(recover_c_from_signature(*sig) == (-1)**n, 'signature decoder')
        check(branch_carrier(Recurrence(L,0,(-1)**n)).root_discriminant == 5*F*F, 'golden splitting')
        # phi^n=(L+F sqrt5)/2, proved independently by multiplying (1+sqrt5)/2.
        aa,bb = Q(1),Q(0)
        for _ in range(n):
            aa,bb = (aa+5*bb)/2,(aa+bb)/2
        check((aa,bb) == (Q(L,2),Q(F,2)), 'power slopes')
        if n <= 6:
            golden.append({'power':n,'signature':list(sig),'c':(-1)**n})
    check(signature(5) == (11,5,125), 'literal signature')
    check(coefficients(SUPPLIED[3],2)[2] == Q(5,2), 'wrong semantic slots')
    check(coefficients(SUPPLIED[5],2)[2] == Q(123,4), 'accessory perturbation')
    check(branch_carrier(SUPPLIED[5]).j == branch_carrier(CANON['F']).j, 'same carrier not integrality')
    rejected(lambda: audit_carrier(SUPPLIED[5],branch_carrier(CANON['F'])))
    check(branch_carrier(CANON['A']).j == branch_carrier(CANON['C']).j == branch_carrier(CANON['F']).j,
          'inherited reanchoring branch orbit')
    check(branch_carrier(CANON['A']).j == Q(1556068,81), 'ACF invariant')
    check(branch_carrier(CANON['D']).j == Q(488095744,125), 'D invariant')
    check(branch_carrier(SUPPLIED[3]).j == Q(4195088384,5921875), 'literal alternative invariant')
    check(coefficients(CANON['C'],4)[4] == 639, 'actual C coefficient, not a digit-object map')

    recs = record_numbers(125)
    check([n for n,d in recs] == [1,2,4,6,12,24,36,48,60,120], 'bounded strict divisor records')
    integers = sorted({v for t in SUPPLIED for v in (t.a,t.b,t.c)})
    check(integers == [-1,3,4,5,9,10,11,12,17,27,32,72,125], 'complete supplied alphabet')
    check([n for n,d in recs if n in integers] == [4,12], 'record input role')
    check(sorted({d for n,d in recs if d in integers}) == [3,4,9,10,12], 'record output role')
    check(divisor_count(72) == divisor_count(60) == 12, '72 is not a new record')

    # Exact-type and source-binding hostiles, including valid-value aliases.
    for bad in (True,1.0,Q(1)):
        rejected(lambda bad=bad: Recurrence(bad,3,9))
        rejected(lambda bad=bad: coefficients(CANON['C'],bad))
    for bad in (True,1.0):
        rejected(lambda bad=bad: decode_jet(bad,15,93))
        rejected(lambda bad=bad: gauge_numerators(CANON['C'],3,(bad,),(1,)))
    rejected(lambda: branch_carrier(Recurrence(2,1,1)))
    rejected(lambda: branch_carrier(Recurrence(1,1,0)))
    rejected(lambda: decode_jet(0,0,0))
    check(coefficients(Recurrence(1,0,0),5) == coefficients(Recurrence(100,0,0),5), 'degenerate jet ambiguity')
    check(decode_jet(1,1,1) == Recurrence(2,1,1), 'singular operator still has a valid jet')
    rejected(lambda: decode_jet(1,0,1))
    rejected(lambda: scale(CANON['D'],0))
    rejected(lambda: audit_prefix(CANON['C'],Prefix(CANON['D'],coefficients(CANON['D'],3))))
    rejected(lambda: audit_prefix(CANON['C'],Prefix(CANON['C'],(1,3,15,94))))
    rejected(lambda: gauge_numerators(CANON['C'],3,(0,),(1,)))
    rejected(lambda: recover_c_from_signature(11,5,124))
    return {
        'status':'FINITE-EXACT controls of the proved typed identities; not global integrality classification',
        'checks':CHECKS,
        'universe':{'canonical_last_index':24,'binomial_sum_last_index':16,
                    'integer_parameter_cube':[-3,5],'jet_decoded_count':jets_decoded,
                    'golden_power_indices':[1,24],'divisor_record_bound':125},
        'canonical':rows,
        'literal_supplied':[{'triple':[t.a,t.b,t.c],
                             'prefix':[str(x) for x in coefficients(t,5)],
                             'canonical_names':[k for k,v in CANON.items() if v == t],
                             'branch_carrier_j':str(branch_carrier(t).j)} for t in SUPPLIED],
        'golden_signatures':golden,'divisor_records':recs,
        'integer_atlas':[{'integer':n,'recurrence_role':role,'context_route':route,
                          'retained_boundary':boundary} for n,role,route,boundary in ATLAS],
        'not_claimed':['all integral triples classified','accessory recovered from j',
                       'auxiliary branch cover identified with a Picard-Fuchs fibre',
                       'all golden powers produce integral sporadic recurrences'],
    }


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--out')
    args = parser.parse_args()
    result = json.dumps(verify(),indent=2,sort_keys=True)+'\n'
    if args.out:
        Path(args.out).write_text(result,encoding='utf-8',newline='\n')
    print(result,end='')
