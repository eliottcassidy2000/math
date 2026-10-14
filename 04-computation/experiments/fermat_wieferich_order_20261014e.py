"""Exact doubling-order towers and source-bound inverse-port guards.

Theorems have no finite bound. This elementary implementation uses trial
division for primes and a linear discrete-log search only at the first level.
Run python -B [-O] -X utf8; stdout is the reproducible JSON artifact.
"""
from dataclasses import dataclass, replace, asdict
from functools import lru_cache
from math import gcd, isqrt
from itertools import permutations
import json


def integer(n, low=None):
    if type(n) is not int or (low is not None and n < low):
        raise ValueError('exact integer outside domain')
    return n


def is_prime(n):
    integer(n, 0)
    return n >= 2 and all(n % d for d in range(2, isqrt(n) + 1))


def factors(n):
    integer(n, 1)
    result = []
    d = 2
    while d*d <= n:
        if n % d == 0:
            result.append(d)
            while n % d == 0:
                n //= d
        d += 1
    if n > 1:
        result.append(n)
    return tuple(result)


def valuation(n, p):
    integer(n)
    integer(p, 2)
    if n == 0:
        raise ValueError('zero has no finite valuation')
    v = 0
    while n % p == 0:
        n //= p
        v += 1
    return v


@dataclass(frozen=True)
class Profile:
    prime: int
    base_order: int
    delay: int
    first_defect: int


@lru_cache(None)
def _profile(p):
    if not is_prime(p) or p == 2:
        raise ValueError('odd prime required')
    d = p - 1
    for ell in factors(d):
        while d % ell == 0 and pow(2, d//ell, p) == 1:
            d //= ell
    s = 1
    while pow(2, d, p**(s+1)) == 1:
        s += 1
    c = (pow(2, d, p**(s+1)) - 1) // p**s
    return Profile(p, d, s, c)


def profile(p):
    integer(p, 3)  # before cache lookup: bool/float aliases are not accepted
    return _profile(p)


def audit_profile(z):
    if type(z) is not Profile:
        raise ValueError('exact Profile required')
    for a in (z.prime, z.base_order, z.delay, z.first_defect):
        integer(a, 1)
    if z != profile(z.prime):
        raise ValueError('incorrect order/defect profile')
    return z


def level(p, k):
    z = profile(p)
    integer(k, 1)
    order = z.base_order * p**max(0, k-z.delay)
    cycles = ((p-1)//z.base_order) * p**min(k-1, z.delay-1)
    return order, cycles


@dataclass(frozen=True)
class Phase:
    prime: int
    depth: int
    target: int
    exponent: int
    period: int


def audit_phase(z):
    if type(z) is not Phase:
        raise ValueError('exact Phase required')
    for a in (z.prime, z.depth, z.period):
        integer(a, 1)
    integer(z.target, 0)
    integer(z.exponent, 0)
    order, _ = level(z.prime, z.depth)
    modulus = z.prime**z.depth
    if not (0 < z.target < modulus and z.target % z.prime and
            z.period == order and z.exponent < order and
            pow(2, z.exponent, modulus) == z.target):
        raise ValueError('incorrect native phase')
    return z


def lift_phase(z, target):
    """One more target digit: either an exact phase or certified no lift."""
    audit_phase(z)
    integer(target, 1)
    p, k, a, L = z.prime, z.depth, z.exponent, z.period
    old, new = p**k, p**(k+1)
    if target >= new or target % old != z.target:
        raise ValueError('target must be a canonical lift of this source target')
    value = pow(2, a, new)
    if k < profile(p).delay:
        if target != value:
            return None
        out = Phase(p, k+1, target, a, L)
    else:
        c = (pow(2, L, new)-1)//old
        digit = ((target-value)//old * pow((value % p)*c, -1, p)) % p
        out = Phase(p, k+1, target, a+digit*L, p*L)
    return audit_phase(out)


def compile_phase(p, k, target):
    """Complete phase of 2**a == target mod p**k, a an integer >=0."""
    profile(p)
    integer(k, 1)
    integer(target, 1)
    if target >= p**k or target % p == 0:
        raise ValueError('canonical unit target required')
    residue = target % p
    L, _ = level(p, 1)
    value, a = 1, None
    for e in range(L):
        if value == residue:
            a = e
            break
        value = 2*value % p
    if a is None:
        return None
    result = Phase(p, 1, residue, a, L)
    for depth in range(2, k+1):
        result = lift_phase(result, target % p**depth)
        if result is None:
            return None
    return audit_phase(result)


@dataclass(frozen=True)
class PortGuard:
    prime: int
    depth: int
    endpoint: int
    sheet: int
    source_class: int
    phase: Phase


def compile_port_guard(p, k, q, b, c):
    """Source n=(q*2**a-b)/p, with n == c mod p**(k-1)."""
    profile(p)
    integer(k, 1)
    integer(q, 1)
    integer(b)
    integer(c, 0)
    if q % 2 == 0 or q % p == 0 or b not in (-1, 1) or c >= p**(k-1):
        raise ValueError('odd unit endpoint, signed sheet and native source class required')
    modulus = p**k
    target = ((b+p*c)*pow(q, -1, modulus)) % modulus
    phase = compile_phase(p, k, target)
    return None if phase is None else PortGuard(p, k, q, b, c, phase)


def audit_guard(z):
    if type(z) is not PortGuard:
        raise ValueError('exact PortGuard required')
    audit_phase(z.phase)
    for a in (z.prime, z.depth, z.endpoint):
        integer(a, 1)
    integer(z.sheet)
    integer(z.source_class, 0)
    if compile_port_guard(z.prime, z.depth, z.endpoint, z.sheet, z.source_class) != z:
        raise ValueError('wrong source-bound guard')
    return z


@dataclass(frozen=True)
class Port:
    guard: PortGuard
    exponent: int
    source: int


def audit_port(z):
    if type(z) is not Port:
        raise ValueError('exact Port required')
    g = audit_guard(z.guard)
    integer(z.exponent, 1)
    integer(z.source, 2)
    if (z.source % 2 != 1 or z.source % g.prime**(g.depth-1) != g.source_class or
        (z.exponent-g.phase.exponent) % g.phase.period or
        g.prime*z.source+g.sheet != g.endpoint*(1 << z.exponent)):
        raise ValueError('port is not an actual source-owned odd edge')
    return z


def make_port(g, a):
    audit_guard(g)
    integer(a, 1)
    if (a-g.phase.exponent) % g.phase.period:
        raise ValueError('exponent misses the retained phase')
    numerator = g.endpoint*(1 << a)-g.sheet
    if numerator % g.prime:
        raise ValueError('nonintegral source')
    return audit_port(Port(g, a, numerator//g.prime))


def direct_cycles(p, k):
    """Independent finite orbit walk, used only by controls."""
    modulus = p**k
    unused = {a for a in range(1, modulus) if a % p}
    result = []
    while unused:
        first = min(unused)
        cycle, x = [], first
        while x in unused:
            unused.remove(x)
            cycle.append(x)
            x = 2*x % modulus
        if x != first:
            raise RuntimeError('not a disjoint permutation orbit')
        result.append(tuple(cycle))
    return result


def permutation_order(z):
    seen, order = set(), 1
    for i in range(len(z)):
        if i in seen:
            continue
        j, length = i, 0
        while j not in seen:
            seen.add(j)
            length += 1
            j = z[j]
        order = order*length//gcd(order, length)
    return order


def main():
    checks = 0
    def check(condition, description):
        nonlocal checks
        checks += 1
        if not condition:
            raise RuntimeError(description)
    def rejects(f):
        try:
            f()
        except (ValueError, TypeError):
            check(True, 'rejected malformed input')
        else:
            check(False, 'accepted malformed input')

    ps = [p for p in range(3, 200, 2) if is_prime(p)]
    for p in ps + [257, 1093, 3511, 65537]:
        z = audit_profile(profile(p))
        for k in range(1, 5):
            L, count = level(p, k)
            check(pow(2, L, p**k) == 1, 'claimed order returns')
            for ell in factors(L):
                check(pow(2, L//ell, p**k) != 1, 'order is minimal')
            check(L*count == (p-1)*p**(k-1), 'cycle count')
        check((pow(2,p-1,p*p) == 1) == (z.delay >= 2), 'Fermat quotient equivalence')

    cycle_rows = []
    for p in (3,5,7,11,13,17,19):
        for k in range(1,4):
            cycles = direct_cycles(p,k)
            L, count = level(p,k)
            check(len(cycles) == count, 'direct cycle count')
            check({len(c) for c in cycles} == {L}, 'direct cycle lengths')
            if k < 3:
                next_cycles = direct_cycles(p,k+1)
                for cy in cycles:
                    lifted = [nc for nc in next_cycles if set(x % p**k for x in nc) == set(cy)]
                    check(len(lifted) == (p if k < profile(p).delay else 1), 'cycle splitting')
            cycle_rows.append([p,k,L,count])

    for p in (3,5,7,11,13):
        for k in range(1,4):
            L,_ = level(p,k)
            orbit = {pow(2,a,p**k):a for a in range(L)}
            for target in range(1,p**k):
                if target % p == 0:
                    continue
                z = compile_phase(p,k,target)
                check((z is not None) == (target in orbit), 'complete phase membership')
                if z is not None:
                    check(z.exponent == orbit[target] and z.period == L, 'exact exponent phase')

    special = []
    for p,d in ((1093,364),(3511,1755)):
        z = profile(p)
        check((z.base_order,z.delay) == (d,2), 'Wieferich exact depth')
        good = compile_phase(p,1,1)
        for t in range(p):
            lifted = lift_phase(good,1+p*t)
            check((lifted is not None) == (t == 0), 'delay rejects other digits')
        good = compile_phase(p,2,1)
        for t in range(p):
            lifted = lift_phase(good,1+p*p*t)
            check(lifted is not None, 'growth accepts every digit')
            check(pow(2,lifted.exponent,p**3) == 1+p*p*t, 'growth digit exact')
        for t in (1,2,3,p,p*p):
            # Evaluate only modulo the first nonzero digit, not the enormous source.
            v = z.delay-1+valuation(t,p)
            numerator = (pow(2,d*t,p**(v+2))-1) % p**(v+2)
            check(numerator % p**(v+1) == 0 and numerator % p**(v+2) != 0,
                  'endpoint-one source valuation')
        special.append(dict(profile=asdict(z), levels=[level(p,k) for k in range(1,5)],
                            endpoint_one_source_multiple=p))

    # Retain the exponent address: normalized source digits are an isometry,
    # even when the unnormalized sources all collapse modulo p.
    for p in (3,5,7,17,1093,3511):
        z = profile(p)
        for t in range(1,13):
            for u in range(1,t):
                v = z.delay-1+valuation(t-u,p)
                mod = p**(v+2)
                difference = (pow(2,z.base_order*t,mod)-pow(2,z.base_order*u,mod)) % mod
                check(difference % p**(v+1) == 0 and difference % p**(v+2) != 0,
                      'source-address valuation isometry')
        depth = 2 if p < 20 else 1
        modulus = p**depth
        values = {((pow(2,z.base_order*t,p**(z.delay+depth))-1)//p**z.delay)
                  for t in range(modulus)}
        check(values == set(range(modulus)), 'normalized source-digit bijection')

    fermat = []
    for m in range(5):
        p = (1 << (1 << m))+1
        d = 1 << (m+1)
        z = profile(p)
        check((z.base_order,z.delay) == (d,1), 'Fermat profile')
        check(((1 << d)-1)//p == p-2, 'exact Fermat defect')
        quotient = (pow(2,p-1,p*p)-1)//p
        check(quotient == pow(2,-m,p), 'Fermat quotient residue')
        fermat.append(dict(index=m,prime=p,order=d,delay=1,cycles=level(p,1)[1],
                           fermat_quotient_mod_prime=quotient,
                           endpoint_one_first_positive_source=p-2,
                           first_source_is_strict=p>3))
    check((1 << 32)+1 == 641*6700417, 'composite next Fermat number')
    F5 = (1 << 32)+1
    check((pow(2,F5-1,F5*F5)-1)//F5 == pow(2,-5,F5),
          'Fermat-shape quotient congruence does not certify primality')
    rejects(lambda: profile(F5))
    repunits = []
    for spacing in range(1,9):
        for count in range(2,13):
            if spacing*(count-1) > 20:
                continue
            p = sum(1 << (spacing*i) for i in range(count))
            if is_prime(p):
                z = profile(p)
                check((z.base_order,z.delay) == (spacing*count,1), 'prime binary repunit')
                check(((1 << (spacing*count))-1)//p == (1 << spacing)-1, 'repunit defect')
                repunits.append([spacing,count,p])

    ports = 0
    for p in (3,5,7,17):
        for k in (1,2,3):
            for q in (1,3,5,7,11):
                if q % p == 0:
                    continue
                for b in (-1,1):
                    for c in range(p**(k-1)):
                        g = compile_port_guard(p,k,q,b,c)
                        # Independent literal powers at this bounded depth.
                        target = ((b+p*c)*pow(q,-1,p**k)) % p**k
                        reachable = target in {pow(2,a,p**k) for a in range(level(p,k)[0])}
                        check((g is not None) == reachable, 'source guard complete')
                        if g is None:
                            continue
                        a = g.phase.exponent
                        while a < 1 or q*(1 << a) <= p+b:
                            a += g.phase.period
                        z = make_port(g,a)
                        x = p*z.source+b
                        actual = 0
                        while x % 2 == 0:
                            x //= 2
                            actual += 1
                        check((x,actual) == (q,a), 'actual odd endpoint and valuation')
                        check(z.source % p**(k-1) == c, 'retained source residue')
                        ports += 1

    check(compile_port_guard(17,1,3,1,0) is None, 'Fermat does not imply transitivity')
    check(compile_port_guard(17,1,3,-1,0) is None, 'both sheets can miss an endpoint')
    for p in (1093,3511):
        check(compile_port_guard(p,2,1,1,1) is None, 'Wieferich source residue blocked')
        check(compile_port_guard(p,2,1,1,0) is not None, 'Wieferich surviving source residue')
    g = compile_port_guard(3,3,1,1,1)
    check(g.phase.exponent == 2 and g.phase.period == 18, 'classical source-cell phase')
    rejects(lambda: make_port(g,2))  # n=1: no strict first edge from ROOT
    check(make_port(g,20).source == 349525, 'same native phase above ROOT')
    good = make_port(compile_port_guard(5,2,1,1,3),4)
    for bad in (replace(good,source=float(good.source)),replace(good,exponent=True),
                replace(good,source=good.source+2),replace(good,guard=replace(good.guard,sheet=-1))):
        rejects(lambda bad=bad:audit_port(bad))
    for p in (True,3.0,2,9):
        rejects(lambda p=p:profile(p))
    for bad in (replace(profile(3),prime=3.0),replace(profile(3),delay=True),
                replace(profile(1093),delay=1)):
        rejects(lambda bad=bad:audit_profile(bad))
    for bad in (replace(g.phase,target=float(g.phase.target)),replace(g.phase,period=True),
                replace(g.phase,depth=2)):
        rejects(lambda bad=bad:audit_phase(bad))
    rejects(lambda:lift_phase(g.phase,g.phase.target+1))
    rejects(lambda:compile_port_guard(3,2,1,1,True))

    schlafli = [(a,b) for a in range(3,10) for b in range(3,10) if (a-2)*(b-2)<4]
    check(schlafli == [(3,3),(3,4),(3,5),(4,3),(5,3)], 'five curvature solutions')
    group_rows = []
    for n,even_only,name in ((4,True,'A4'),(4,False,'S4'),(5,True,'A5')):
        group = []
        for z in permutations(range(n)):
            inversions = sum(z[i]>z[j] for i in range(n) for j in range(i+1,n))
            if not even_only or inversions % 2 == 0:
                group.append(z)
        orders = sorted(set(map(permutation_order,group)))
        check(17 not in orders, 'no Platonic element of order17')
        group_rows.append(dict(name=name,size=len(group),element_orders=orders))
    check([r['size'] for r in group_rows] == [12,24,60], 'three rotation-group sizes')

    print(json.dumps(dict(status='PROVED elementary order/guard identities; FINITE-EXACT controls',
        checks=checks, fermat=fermat, wieferich=special, binary_repunits=repunits,
        direct_cycle_rows=cycle_rows, actual_ports_checked=ports,
        platonic_groups=group_rows, schlafli=schlafli,
        limits='No complete list of Fermat or Wieferich primes; no Collatz convergence or graph-colouring theorem'),
        indent=2,sort_keys=True))


if __name__ == '__main__':
    main()
