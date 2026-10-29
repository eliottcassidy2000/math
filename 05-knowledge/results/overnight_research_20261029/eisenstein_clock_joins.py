"""Exact matched C3 cross-block classification and minimal-carrier checks."""
from collections import Counter
from fractions import Fraction as F
from itertools import product
from math import gcd, isqrt
from pathlib import Path
import json
from port_queries import mm, tr, invdet

GATES = 0
A = [[0, -1], [1, -1]]


def need(ok, message):
    global GATES
    GATES += 1
    if not ok:
        raise RuntimeError(message)


def norm(a, b):
    return a*a + a*b + b*b


def block(a, b):
    return [[a, b], [-a-b, a]]


def source(d, e, a, b):
    return [[0,d,a,b],[-d,0,-a-b,a],[-a,a+b,0,e],[-b,-a,-e,0]]


def factors(n):
    out = {}
    p = 2
    while p*p <= n:
        while n % p == 0:
            out[p] = out.get(p, 0) + 1
            n //= p
        p += 1
    if n > 1:
        out[n] = out.get(n, 0) + 1
    return out


def divisors(n):
    return [d for d in range(1,n+1) if n % d == 0]


def mobius(n):
    ff = factors(n)
    return 0 if any(v > 1 for v in ff.values()) else (-1)**len(ff)


def r_factor(n):
    if n < 0:
        return 0
    if n == 0:
        return 1
    ans = 6
    for p, v in factors(n).items():
        if p % 3 == 2 and v % 2:
            return 0
        if p % 3 == 1:
            ans *= v+1
    return ans


def chi(n):
    return (0,1,-1)[n % 3]


def count_group(d, e, s, n):
    if d % s or e % s:
        return 0
    dd, ee = d//s, e//s
    ans = 0
    for sign in (-1,1):
        target = dd*ee-sign*n
        if target < 0:
            continue
        for h in divisors(gcd(dd,ee)):
            if target % (h*h) == 0:
                ans += mobius(h)*r_factor(target//(h*h))
    return ans


def shell_and_group_counts():
    limit = 4096
    bound = isqrt(4*limit//3)
    pairs = [(a,b,norm(a,b)) for a,b in product(range(-bound,bound+1), repeat=2)
             if norm(a,b) <= limit]
    histogram = Counter(n for a,b,n in pairs)
    for n in range(limit+1):
        need(histogram[n] == r_factor(n), 'full norm shells versus factor formula')
        if n:
            # Independent divisor sieve, not factoring again.
            sm = sum(chi(d) + (chi(n//d) if d*d != n else 0)
                     for d in range(1,isqrt(n)+1) if n % d == 0)
            need(histogram[n] == 6*sm, 'divisor character formula')
    labelled = nonsingular = queries = 0
    for d,e in product(range(1,7), repeat=2):
        direct = Counter()
        for a,b,nn in pairs:
            t = d*e-nn
            labelled += 1
            if not t:
                continue
            nonsingular += 1
            s = gcd(gcd(d,e),gcd(a,b))
            need(t % (s*s) == 0, 'paired invariant chain divisibility')
            n = abs(t)//(s*s)
            if n <= 50:
                direct[s,n] += 1
        for s in range(1,7):
            for n in range(1,51):
                need(direct[s,n] == count_group(d,e,s,n), 'all-order content count versus complete norm shells')
                queries += 1
    primes = [p for p in range(2,2001) if factors(p) == {p:1}]
    admitted = []
    for p in primes:
        number = count_group(3,3,1,p)
        need(number == r_factor(9-p)+r_factor(9+p), 'prime-plane formula')
        if number:
            admitted.append(p)
            need(p <= 9 or p % 12 == 7, 'necessary prime residue')
    need(admitted[:10] == [2,3,5,7,19,43,67,103,139,163], 'named prime examples')
    need(count_group(3,3,1,31) == 0, '7 modulo12 is not sufficient')
    special = {str(n): count_group(3,3,1,n) for n in [5,9,11,25,27,57]}
    need([special[str(n)] for n in (25,27,57)] == [0,0,0], 'requested cyclic paired planes excluded')
    need(count_group(3,3,3,3) == 6, 'norm36 belongs to H3 plus H9')
    need(count_group(3,3,3,1) == 1, 'zero coupling content case')
    return {'norm_limit':limit,'coordinate_box':[-bound,bound],
            'norm_pairs':len(pairs),'source_d_e':[1,6],
            'labelled_pairs_examined':labelled,'nonsingular_pairs_examined':nonsingular,
            'complete_group_queries':queries,'group_parameter_n':[1,50],
            'primes_through':2000,'admitted_primes':admitted,'special_Hn_counts':special,
            'H3_plus_H9_count':6}


def chart_and_carrier_checks():
    # Exhaust arbitrary matrices in a box before accepting the 2-coordinate chart.
    legal = 0
    for x,y,z,w in product(range(-2,3),repeat=4):
        b = [[x,y],[z,w]]
        actual = mm(tr(A),mm(b,A)) == b
        need(actual == (z == -x-y and w == x), 'complete linear clock equations')
        legal += int(actual)
    v = [[F(t)**k for k in range(5)] for t in range(5)]
    vi,_ = invdet(v)
    decoded = terminal_samples = 0
    for d,e in product(range(1,7),repeat=2):
        for a,b in product(range(-2,3),repeat=2):
            if d*e == norm(a,b):
                continue
            coeffs = []
            for axis in range(2):
                vals = []
                for w in range(5):
                    aa,bb = (a+w,b) if axis == 0 else (a,b+w)
                    determinant = invdet(source(d,e,aa,bb))[1]
                    need(determinant == (d*e-norm(aa,bb))**2, 'direct determinant response')
                    terminal_samples += int(not determinant)
                    vals.append([determinant])
                coeffs.append([row[0] for row in mm(vi, vals)])
            l1,l2 = coeffs[0][3]/2,coeffs[1][3]/2
            t = (l1*l1-coeffs[0][2])/2
            recovered = ((2*l1-l2)/3,(2*l2-l1)/3)
            need(recovered == (a,b) and t+norm(*recovered) == d*e,
                 'nine-query minimal determinant carrier decoder')
            cancelled = source(d,e,0,0)
            first_factor = gcd(*(abs(x) for row in cancelled for x in row))
            need(first_factor == gcd(d,e), 'cancelling query recovers content coordinate')
            decoded += 1
    # Clock-centralizing unimodular re-markings are precisely these six units.
    at = tr(A)
    units,cur = [],[[1,0],[0,1]]
    for _ in range(3):
        units += [cur, [[-x for x in row] for row in cur]]
        cur = mm(cur,at)
    remarks = 0
    for r1,r2 in product(units,repeat=2):
        i1,d1 = invdet(r1);i2,d2 = invdet(r2)
        need(d1 == d2 == 1 and mm(r1,at) == mm(at,r1), 'clock unit')
        for a,b in product(range(-2,3),repeat=2):
            bb = mm(i1,mm(block(a,b),tr(i2)))
            aa,bbb = bb[0]
            need(bb == block(aa,bbb) and norm(aa,bbb) == norm(a,b), 'closed unit chart and norm')
            need(aa.denominator == bbb.denominator == 1, 'remark retains integral coordinates')
            need(gcd(3,gcd(int(aa),int(bbb))) == gcd(3,gcd(a,b)), 'content survives remark')
            remarks += 1
    # Phase and source-content hostiles are symbolic all-future claims in the note.
    need(norm(1,0) == norm(0,1) and norm(2,0) != norm(1,1), 'norm alone is not closed')
    need(invdet(source(3,3,2,0))[1] == 25 and invdet(source(3,3,1,1))[1] == 36,
         'same attachment separates equal-norm states')
    for a,b in product(range(-3,4),repeat=2):
        need(invdet(source(1,9,a,b))[1] == invdet(source(3,3,a,b))[1],
             'same product distinct content determinant hostile')
    need(gcd(1,9) != gcd(3,3), 'paired group distinguishes source factors')
    need(1*36 == 4*9 and gcd(1,36) == gcd(4,9) == 1,
         'same unmarked carrier can have different named port orders')
    for d,e in product(range(1,7),repeat=2):
        z = d*e
        dp = invdet(source(d+1,e,0,0))[1]
        dm = invdet(source(d-1,e,0,0))[1]
        need((dp-dm)/(4*z) == e, 'one internal direction restores source factor')
    bad = [[1,0],[0,-2]]
    need(mm(tr(A),mm(bad,A)) != bad, 'H11 witness is illegal with clock')
    need(gcd(3,3,1,-2) == 1 and 9-invdet(bad)[1] == 11, 'H11 witness without clock')
    need(norm(1,7) == 57 and r_factor(66) == 0, 'old norm57 is not shifted join norm66')
    return {'arbitrary_cross_matrices':625,'legal_in_box':legal,
            'nonsingular_states_decoded':decoded,'singular_query_values':terminal_samples,
            'clock_unit_remark_checks':remarks,'phase_and_content_hostiles':True}


def main():
    out = {'counts':shell_and_group_counts(),'carriers':chart_and_carrier_checks(),
           'gates':GATES,'scope':'matched irreducible C3 clocks; labelled integer cross couplings',
           'classical_norm_source':'https://categorified.net/23Winter3032/HW6.pdf'}
    Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(out,indent=2))


if __name__ == '__main__':
    main()
