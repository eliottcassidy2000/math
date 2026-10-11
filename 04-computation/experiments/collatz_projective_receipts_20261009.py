"""Source-owned run ports and marked projective receipt rays.

No giant Mersenne expansion and no production ROOT discovery.  The exact
inherited fan is consumed only as a set of exponent labels in main().
"""
from dataclasses import dataclass, replace
from fractions import Fraction as F
import collatz_fibre_integration_20261009 as fibres
import collatz_child_normal_forms_20261007 as inverse
from collatz_debt_observer_integration_20261009 import replay_u

CHECKS = 0


def check(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise ValueError(message)


def integer(n, minimum=0):
    if type(n) is not int or n < minimum:
        raise ValueError('exact integer in the stated domain required')
    return n


def v3(n):
    integer(n, 1)
    r = 0
    while n % 3 == 0:
        n //= 3
        r += 1
    return r


def available(source):
    integer(source, 1)
    return v3(2*source+1) if source % 4 == 3 else 0


def carrier(word):
    if type(word) is not tuple or any(type(a) is not int or a < 1 for a in word):
        raise ValueError('exact positive valuation tuple required')
    p, q, b = 1, 1, 0
    for a in word:
        p, q, b = 3*p, (1 << a)*q, 3*b+q
    return p, q, b


@dataclass(frozen=True)
class Row:
    run: int
    t0: int
    base: tuple
    period: tuple


def row(run):
    integer(run, 4)
    t0 = 7 if run % 2 == 0 else 5
    p = 3**run
    q = 1 << (run+1)
    return Row(run, t0, ((p*t0-1)//2, q*t0-1, (3*p*t0-1)//4),
               (4*p, 8*q, 6*p))


def ray(run, parameter):
    integer(parameter)
    spec = row(run)
    return tuple(a+parameter*b for a, b in zip(spec.base, spec.period))


@dataclass(frozen=True)
class Join:
    source: int
    child: int
    endpoint: int
    run: int
    source_word: tuple
    child_word: tuple


def join(source, run=None):
    integer(source, 1)
    fuel = available(source)
    if run is None:
        run = fuel
    integer(run, 4)
    if run > fuel:
        raise ValueError('requested run exceeds source-owned ternary fuel')
    t = (2*source+1)//3**run
    child = (1 << (run+1))*t-1
    endpoint = (3*source+1)//2
    u, v = (1,), (1,)*run+(3,)
    if not (1 < child < source < endpoint):
        raise ValueError('strict source/child/payment boundary failed')
    if replay_u(source, u) != endpoint or replay_u(child, v) != endpoint:
        raise ValueError('actual native words failed')
    return Join(source, child, endpoint, run, u, v)


def audit(spec):
    if type(spec) is not Join:
        raise ValueError('typed Join required')
    for n in (spec.source, spec.child, spec.endpoint, spec.run):
        integer(n, 1)
    carrier(spec.source_word); carrier(spec.child_word)
    if spec != join(spec.source, spec.run):
        raise ValueError('packet differs from its actual source')
    return spec


def clocks(spec):
    audit(spec)
    return (sum(spec.source_word)-sum(spec.child_word),
            len(spec.source_word)-len(spec.child_word))


def discharge(spec, supplied_child_root):
    audit(spec)
    replay_u(spec.child, supplied_child_root, True)
    k = len(spec.child_word)
    if supplied_child_root[:k] != spec.child_word:
        raise ValueError('child ROOT prefix does not match the common future')
    result = spec.source_word+supplied_child_root[k:]
    replay_u(spec.source, result, True)
    return result


def chain(source, runs):
    """Return child, common endpoint, and actual words for chosen legal runs.

    The unchanged original source has the word (1).  No step substitutes a
    different source for that original obligation or assumes ROOT.
    """
    integer(source, 1)
    if type(runs) is not tuple or not runs:
        raise ValueError('nonempty exact run tuple required')
    current, child_word = source, (1,)
    endpoint = (3*source+1)//2
    for r in runs:
        spec = join(current, r)
        child_word = spec.child_word+child_word[1:]
        current = spec.child
        if replay_u(current, child_word) != endpoint:
            raise ValueError('same-source clock-alignment seam failed')
    return current, endpoint, (1,), child_word


def simplex_point(run, parameter):
    n, m, j = ray(run, parameter)
    z = (n, m, j, 1)
    return tuple(F(x, sum(z)) for x in z)


def projective_components(run, parameter):
    """Exact convex decomposition; it does not drop the integer parameter."""
    integer(parameter)
    spec = row(run)
    z0, d = spec.base+(1,), spec.period+(0,)
    a, b = sum(z0), sum(d)
    weight = F(a, a+parameter*b)
    return weight, tuple(F(x, a) for x in z0), tuple(F(x, b) for x in d)


@dataclass(frozen=True)
class MersennePort:
    exponent: int
    run: int


def mersenne_fuel(exponent):
    integer(exponent, 2)
    return 1+v3(exponent+1) if exponent % 2 else 0


def symbolic_port(exponent, run=None):
    fuel = mersenne_fuel(exponent)
    if run is None:
        run = fuel
    integer(run, 4)
    if run > fuel:
        raise ValueError('Mersenne exponent fails the native divisibility')
    return MersennePort(exponent, run)


def audit_port(port):
    if type(port) is not MersennePort or port != symbolic_port(port.exponent, port.run):
        raise ValueError('authenticated symbolic Mersenne port required')
    return port


def child_mod(port, modulus):
    audit_port(port); integer(modulus, 1)
    p, q = 3**port.run, 1 << (port.run+1)
    # Retain the full denominator, also when 3 divides the requested modulus.
    numerator = q*(pow(2, port.exponent+1, p*modulus)-1) % (p*modulus)
    if numerator % p:
        raise ValueError('symbolic divisibility failed')
    return (numerator//p-1) % modulus


def relative_factor(port):
    audit_port(port)
    return F(1 << (port.run+2), 3**port.run)


@dataclass(frozen=True)
class ExtendedPort:
    port: MersennePort
    inverse_word: tuple


def extend(port, word):
    """Authenticate a selected inherited inverse word at the actual child."""
    audit_port(port)
    record = inverse.encode(word)
    residue, modulus = inverse.native_cell(record)
    if child_mod(port, modulus) != residue:
        raise ValueError('inverse word is not native at the retained child')
    return ExtendedPort(port, word)


def extension_mod(spec, modulus):
    if type(spec) is not ExtendedPort or spec != extend(spec.port, spec.inverse_word):
        raise ValueError('authenticated ExtendedPort required')
    integer(modulus, 1)
    record = inverse.encode(spec.inverse_word)
    den = 3**record.B
    numerator = (1 << record.A)*child_mod(spec.port, den*modulus)-record.C
    if numerator % den:
        raise ValueError('full inverse denominator was lost')
    return (numerator//den) % modulus


def extension_words(spec):
    extension_mod(spec, 1)
    return (1,), inverse.forward_word(spec.inverse_word)+(1,)*spec.port.run+(3,)


def follow_inverse(port, budget=12):
    integer(budget)
    spec = extend(port, ())
    for _ in range(budget):
        z = extension_mod(spec, 2187)
        choices = tuple(g for g in inverse.GENERATORS
                        if z % inverse.native_cell(inverse.encode((g,)))[1]
                        == inverse.native_cell(inverse.encode((g,)))[0])
        if not choices:
            return spec, 'STALLED'
        if len(choices) != 1:
            raise ValueError('inherited disjoint guard invariant failed')
        spec = extend(port, spec.inverse_word+choices)
    return spec, 'BUDGET'


def above_lower_mersenne(port, lower_exponent):
    """Exact sufficient comparison using only the small exponent gap.

    True proves child > M_lower; False makes no comparison claim.
    """
    audit_port(port); integer(lower_exponent, 2)
    if lower_exponent > port.exponent:
        raise ValueError('ordered exponent comparison required')
    gap = port.exponent-lower_exponent
    return (1 << (port.run+1))*((1 << (gap+1))-1) > 3**port.run


def bad(call):
    try:
        call()
    except (TypeError, ValueError):
        check(True, 'hostile rejected')
    else:
        check(False, 'hostile accepted')


ROOT223 = (1,1,1,1,3,2,1,1,1,1,1,3,1,1,1,4,2,2,4,3,1,1,5,4)


def main():
    global CHECKS
    CHECKS = 0
    count = 0
    for r in range(4, 17):
        p, q, b = carrier((1,)*r+(3,))
        check((p, q, b) == (3**(r+1), 1 << (r+3), 3**(r+1)-(1 << (r+1))), 'closed carry identity')
        for k in range(65):
            n, m, j = ray(r, k)
            spec = join(n, r)
            check((spec.child, spec.endpoint) == (m, j), 'native ray versus source-owned reader')
            check(clocks(spec) == (-r-2, -r), 'true clock defect')
            check(m < (8*n-5)//9 < n and v3(n+5) == 2, 'strict improvement over the stopped G5 port')
            x = m
            for a in spec.child_word:
                x = replay_u(x, (a,))
                check(x != n, 'child prefix does not pass through parent')
            point = simplex_point(r, k)
            w, start, infinity = projective_components(r, k)
            check(point == tuple(w*a+(1-w)*d for a, d in zip(start, infinity)), 'projective simplex identity')
            check(point[0] > point[1] and infinity[-1] == 0, 'payment and boundary are marked')
            count += 1
    print('All-height formula controls:', count, 'rows; r=4..16, parameter0..64; actual strict word replay.')

    # A complete bounded source-domain comparison, independent congruence side.
    for n in range(1, 20001, 2):
        expected = n % 324 == 283
        check((available(n) >= 4) == expected, 'whole ordinary source domain iff')
        if expected:
            selected = join(n)
            for r in range(4, available(n)+1):
                check(selected.child <= join(n, r).child, 'maximal fuel is the smallest offered child')
                if r < available(n):
                    current = join(n,r).child
                    check(inverse.apply(inverse.encode((1,)),current) == join(n,r+1).child,
                          'successive run choices normalize by an actual G1 predecessor')
    first = join(283, 4)
    root283 = discharge(first, ROOT223)
    check(len(root283) == 20 and sum(root283) == 40, 'supplied223 proof grounds283 with exact clock saving')
    print('Least member:283 ->425 <-223; G5 instead selects251. Supplied ROOT223 gives ROOT283, ranks24->20, costs46->40.')

    # Real convexity preserves incidence, but not the odd integer native lattice.
    n0, m0, j0 = ray(4, 0)
    n1, m1, j1 = ray(4, 1)
    n, m, j = tuple(F(a+b, 2) for a, b in zip((n0,m0,j0), (n1,m1,j1)))
    check((n,m,j) == (445,351,668), 'integral midpoint hostile')
    check(3*n+1 == 2*j and 243*m+211 == 128*j and m < n < j, 'midpoint keeps incidence and payment')
    bad(lambda: replay_u(int(n), (1,)))
    check((3*n0+1)/2 != F(243*m1+211,128), 'independent box coordinates destroy incidence')
    print('Convex midpoint(445,351,668) retains both affine equations/payment but loses exact odd-word guards.')

    for s in range(100):
        n, m, z = 13891+26244*s, 10975+20736*s, 8671+16384*s
        check(join(n,4).child == m and join(m,4).child == z, 'all-height two-return ray')
        out, endpoint, u, v = chain(n, (4,4))
        check((out,u,v) == (z,(1,),(1,1,1,1,3,1,1,1,3)), 'nonzero-clock seam composition')
        check((sum(u)-sum(v),len(u)-len(v)) == (-12,-8), 'composed clock defects')
    check(2*64**3 < 81**3, 'three run reductions more than halve shifted source')
    print('Repeated family:13891+26244s ->10975+20736s ->8671+16384s, runs(4,4), defect(-12,-8).')

    small_ports = 0
    for e in range(2, 1001):
        check((mersenne_fuel(e) >= 4) == (e % 54 == 53), 'Mersenne phase iff')
        if e % 54 == 53:
            port = symbolic_port(e)
            ordinary = join((1 << e)-1)
            check(ordinary.run == port.run, 'symbolic maximal fuel')
            for modulus in (1, 2, 16, 27, 81, 55, 1536):
                check(child_mod(port, modulus) == ordinary.child % modulus, 'full-denominator modular reader')
            extended, status = follow_inverse(port)
            actual = inverse.apply(inverse.encode(extended.inverse_word), ordinary.child)
            for modulus in (4,27,324,2187):
                check(extension_mod(extended,modulus) == actual % modulus, 'symbolic extension versus exact ordinary source')
            u,v = extension_words(extended)
            check(replay_u(actual,v) == ordinary.endpoint and u == (1,), 'native inverse extension keeps original parent')
            check((ordinary.child+1) & ordinary.child != 0, 'ordinary child is not Mersenne')
            small_ports += 1
    print('Small symbolic/literal Mersenne ports:', small_ports, '; exponents2..1000, seven moduli each.')

    fan = fibres.fan_members()
    ports = tuple(symbolic_port(e) for e in sorted(fan) if mersenne_fuel(e) >= 4)
    hist = {r:sum(p.run == r for p in ports) for r in sorted({p.run for p in ports})}
    check(len(fan) == 1949 and len(ports) == 36 and hist == {4:24,5:8,6:3,7:1}, 'exact inherited label set, not interval fill')
    best = min(ports, key=relative_factor)
    check((best.exponent, best.run, relative_factor(best)) == (99708992495,7,F(512,2187)), 'strongest multiplicative saving')
    check(ports[0].exponent == 99708991793, 'least eligible inherited exponent')
    for port in ports:
        check(above_lower_mersenne(port, min(fan)), 'new ports do not beat the old class minimum')
        check(child_mod(port, 1 << (port.run+2)) == (1 << (port.run+1))-1, 'non-Mersenne odd cofactor marker')
    expected = {99708992225:(5,),99708992711:(5,5),99708992819:(5,),
                99708993197:(5,1),99708993683:(5,)}
    found = {}
    for port in ports:
        extended,status = follow_inverse(port)
        check(status == 'STALLED', 'finite guard route actually stops, not a budget exhaustion')
        for k in range(len(extended.inverse_word)+1):
            stage = extend(port,extended.inverse_word[:k])
            check(extension_mod(stage,324) != 283, 'no hidden run-return alternative in this retained path')
        if extended.inverse_word:
            found[port.exponent] = extended.inverse_word
            check(extension_mod(extended,3) == 0, 'five terminal children have no ordinary odd predecessor')
    check(found == expected, 'exact five independent guarded continuations')
    full_count = 0
    for port in ports:
        maximal, _ = follow_inverse(port)
        for r in range(4,port.run+1):
            smaller_run = symbolic_port(port.exponent,r)
            extended,status = follow_inverse(smaller_run)
            check(status == 'STALLED' and extended.inverse_word == (1,)*(port.run-r)+maximal.inverse_word,
                  'all run choices are the same exact normalized terminal port')
            for modulus in (324,2187,65536):
                check(extension_mod(extended,modulus) == extension_mod(maximal,modulus),
                      'normalization keeps marked arithmetic source residues')
            full_count += 1
    check(full_count == 53, 'complete finite run-choice universe')
    print('Inherited fan:1949 exact labels;36 new non-Mersenne ports; maximal-fuel histogram',hist)
    print('Eligible exponent/run pairs:',tuple((p.exponent,p.run) for p in ports))
    print('Best leading factor512/2187 at E99708992495,R7. Every port exceeds existing minimum M_99708991741.')
    print('All36 port ROOT obligations remain OPEN; no giant source or certificate was materialized.')
    print('Further inherited inverse paths:',tuple(sorted(found.items())))
    print('All53 initial run choices normalize by G1 to the same36 maximal ports and their recorded continuations.')
    print('All36 stop under this declared run/G1/G5/G17 bank; five improved terminals are divisible3. Other rule banks remain untested.')

    bad(lambda: join(True))
    bad(lambda: join(283,4.0))
    bad(lambda: join(445,4))
    bad(lambda: join(283,5))
    bad(lambda: audit(replace(first,child=223.0)))
    bad(lambda: discharge(first,(4,)))
    bad(lambda: symbolic_port(1459))
    bad(lambda: symbolic_port(53,True))
    bad(lambda: child_mod(replace(symbolic_port(53),exponent=53.0),9))
    bad(lambda: child_mod(symbolic_port(53),True))
    bad(lambda: ray(4,-1))
    bad(lambda: chain(13891,(4,True)))
    bad(lambda: chain(283,(4,4)))
    bad(lambda: extend(symbolic_port(53),(True,)))
    bad(lambda: extension_mod(ExtendedPort(symbolic_port(53),(1,)),27))
    print('Exact checks:',CHECKS)


if __name__ == '__main__':
    main()
