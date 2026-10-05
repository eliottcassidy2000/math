"""Sharp branch tolls and a proper integral rank on the signed Collatz graph.

No universal convergence assumption enters the rank reduction or local minimum
classification. Frozen-source joins remain obligations unless a child is grounded.
"""
from collections import deque
from fractions import Fraction
from pathlib import Path
import json

CHECKS = 0


def check(ok, why):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise ArithmeticError(why)


def v2(n):
    if n == 0:
        raise ValueError('valuation at zero')
    n = abs(n)
    return (n & -n).bit_length()-1


def step(n, sign=1):
    z = 3*n+sign
    a = v2(z)
    return z//2**a, a


def replay(n, word, sign=1):
    nodes = [n]
    for a in word:
        z, count = 3*n+sign, 0
        while z % 2 == 0:
            z //= 2
            count += 1
        check(count == a, 'literal exact valuation')
        n = z
        nodes.append(n)
    return nodes


def rank(n):
    if n == 1:
        return (0, 0)
    k = v2(n-1)
    t = (n-1)//2**k
    return (3**k*t*t, k)


def ell(k):
    d = 1
    while 3**d <= 2**(d+2*k):
        d += 1
    return d


def cbrt(n):
    lo, hi = 0, 1 << ((n.bit_length()+2)//3)
    while lo+1 < hi:
        mid = (lo+hi)//2
        if mid**3 <= n:
            lo = mid
        else:
            hi = mid
    return lo


def critical(n):
    return n in (1, -1) or (n % 4 == 3 and n % 3 != 2 and n % 32 != 11)


def lower_neighbors(n):
    """Complete, using the proved properness bound; no inverse exponent cap."""
    R = rank(n)
    result = set()
    y, _ = step(n)
    if rank(y) < R:
        result.add(y)
    if n % 3 == 0:
        return sorted(result)
    bound = cbrt(R[0]**2)
    a = 1 if n % 3 == 2 else 2
    while True:
        m = (2**a*n-1)//3
        if abs(m-1) > bound:
            break
        if rank(m) < R:
            result.add(m)
        a += 2
    return sorted(result)


def branch_audit():
    sources = joins = 0
    for sign in (1, -1):
        for n in range(3, 502, 2):
            q = v2(n+sign)
            if q < 2:
                continue
            t = (n+sign)//2**q
            if (3**q*t-sign) % 4 != 2:
                continue
            J = (3**q*t-sign)//2
            sources += 1
            canonical = (2,)+(1,)*(q-1)
            for m in range(1, n, 2):
                x, word = m, ()
                for d in range(1, q+ell(6)+1):
                    x, a = step(x, sign)
                    word += (a,)
                    if x != J:
                        continue
                    reverse = tuple(reversed(word))
                    differing = [j for j in range(min(q, d))
                                 if reverse[j] != canonical[j]]
                    if not differing:
                        check(d > q, 'canonical smaller ancestor comes after source')
                        continue
                    j = differing[0]
                    delta = reverse[j]-canonical[j]
                    check(delta > 0 and delta % 2 == 0, 'first branch excess is positive even')
                    k = delta//2
                    check(d >= q+ell(k), 'general sharp branch toll')
                    joins += 1
    families = []
    for q in range(2, 9):
        B = 2**(q+2)
        b = 2**q*(3 if q % 2 == 0 else 1)-1
        for k in range(1, 13):
            d = ell(k)
            T = 3**d
            r = -((4**k+2)//3)*pow(4**k, -1, T) % T
            n0 = b+B*((r-b)*pow(B, -1, T) % T)
            S = 4**k*n0+(4**k-1)//3
            m0 = 2**d*((S+1)//T)-1
            P, Q = B*T, 2**(q+d+2*k+2)
            w = (1,)*(q-1)+(2,)
            v = (1,)*d+(1+2*k,)+(1,)*(q-2)+(2,)
            for h in (0, 1, 17, 10**6):
                n, m = n0+P*h, m0+Q*h
                left, right = replay(n, w), replay(m, v)
                check(0 < m < n and left[-1] == right[-1], 'sharp all-height family')
                check(n not in right, 'genuine bypass of source')
            families.append(dict(q=q, k=k, delay=d, n=n0, m=m0, P=P, Q=Q,
                                 source_word=w, child_word=v))
    return dict(signed_sources=sources, noncanonical_smaller_hits=joins, sharp_families=families)


def address_commutator():
    cases = 0
    for k in range(1, 25):
        for d in range(1, 9):
            defect = Fraction(2*(4**k-1)*(3**d-2**d), 3**(d+1))
            for n in (1, 3, 7, 27, 223, 233):
                S = lambda x: 4**k*x+Fraction(4**k-1, 3)
                E = lambda x: Fraction(2**d, 3**d)*(x+1)-1
                check(E(S(n))-S(E(n)) == defect, 'exact commutation defect')
                cases += 1
            check((defect.denominator == 1) == (k % 3**d == 0), 'integral square obstruction')
    return cases


def rank_audit():
    audited = 0
    for n in range(-2001, 2002, 2):
        E, k = rank(n)
        check(E**2 >= abs(n-1)**3, 'proper integral energy')
        check((not lower_neighbors(n)) == critical(n), 'complete signed critical classification')
        if n != 1 and step(n)[1] == 2:
            y, _ = step(n)
            check(rank(y) == (E, k-2), 'valuation-two invariant and fuel drop')
        audited += 1
    stops, larger, max_steps = set(), 0, (0, None)
    for n in range(-10001, 10002, 2):
        x, path = n, [n]
        while not critical(x):
            y, _ = step(x)
            if rank(y) >= rank(x):
                check(x % 3 == 2, 'remaining legal inverse-one edge')
                y = (2*x-1)//3
            check(rank(y) < rank(x), 'well-founded reduction')
            check(step(x)[0] == y or step(y)[0] == x, 'basin-preserving actual edge')
            larger += abs(y) > abs(x)
            path.append(y)
            x = y
        stops.add(x)
        if len(path)-1 > max_steps[0]:
            max_steps = (len(path)-1, path)
    return dict(complete_neighbor_sources=audited, reduction_sources=10002,
                distinct_critical_stops=len(stops), magnitude_increasing_edges=larger,
                longest_reduction=max_steps,
                critical_residues=[n for n in range(1, 96, 2) if critical(n) and n != 1])


def reflection_audit():
    for n in range(-2001, 2002, 2):
        y, a = step(n)
        z, b = step(2-n)
        check(rank(n) == rank(2-n), 'root-centered sign fold')
        if a == 1:
            check(b == 1 and 2-z == y-2, 'oriented carry at valuation one')
        elif a == 2:
            check(b == 2 and 2-z == y, 'reflection at valuation two')
        elif a == 3:
            check(b >= 4, 'deep reset exchanged with valuation three')
        else:
            check(b == 3, 'valuation three exchanged with deep reset')
    for modulus in (18, 36):
        for h in range(1, modulus, 6):
            check(h*(2-h) % modulus == 1, 'reflection becomes inversion on the six-spaced subgroup')
    return dict(folded_cycle_minima=[(3,-1),(7,-5),(19,-17)],
                mod18_involution=[(h,(2-h) % 18) for h in (1,7,13)])


def rank_lift(row):
    n, m, w, v = row['n'], row['m'], tuple(row['w']), tuple(row['v'])
    check(m != 1 and rank(m) < rank(n), 'nonterminal ranked child')
    r, s, A, D = len(w), len(v), sum(w), sum(v)
    Kn, Km = rank(n)[1], rank(m)[1]
    P, Q = 2**(A+1)*3**max(s-r, 0), 2**(D+1)*3**max(r-s, 0)
    h = max(0, Kn+1-v2(P), Km+1-v2(Q))
    P, Q = P*2**h, Q*2**h
    check(Q*Q*3**Km*4**Kn <= P*P*3**Kn*4**Km, 'weighted family slope')
    Wn, Wm = Fraction(3**Kn,4**Kn), Fraction(3**Km,4**Km)
    polynomial = (Wn*(n-1)**2-Wm*(m-1)**2,
                  2*(Wn*(n-1)*P-Wm*(m-1)*Q), Wn*P*P-Wm*Q*Q)
    check(all(c >= 0 for c in polynomial) and (polynomial[0] > 0 or Km < Kn),
          'independent all-height energy polynomial certificate')
    for t in (0, 1, 2, 17, 10**6):
        x, y = n+P*t, m+Q*t
        check((rank(x)[1], rank(y)[1]) == (Kn, Km), 'retained root-distance strata')
        check(replay(x, w)[-1] == replay(y, v)[-1] and rank(y) < rank(x),
              'ranked all-height certificate transport')
    return dict(n=n, m=m, w=w, v=v, P=P, Q=Q, source_k=Kn, child_k=Km,
                integer_slope=str(Fraction(Q, P)), refined_by=2**h,
                energy_difference_coefficients=list(map(str,polynomial)))


def ranked_inverse(J, N, maximum, r, A):
    original = rank(N)
    M = 1+cbrt(original[0]**2)
    threshold = M+1
    queue = deque([(J, (), 0, (J,))])
    visits = 0
    while queue:
        z, reverse, D, ancestry = queue.popleft()
        s = len(reverse)
        visits += 1
        Ez, Kz = rank(z)
        if (Ez, Kz) < original:
            if z == 1 or 2**(2*D)*3**(2*r+Kz)*4**original[1] <= \
                    2**(2*A)*3**(2*s+original[1])*4**Kz:
                return z, tuple(reversed(reverse)), visits
        remaining = maximum-s
        if remaining == 0 or z % 3 == 0:
            continue
        if 2**remaining*(z+1) >= 3**remaining*(threshold+1):
            continue
        a = 1 if z % 3 == 2 else 2
        while True:
            x = (2**a*z-1)//3
            if 2**(remaining-1)*(x+1) >= 3**(remaining-1)*(threshold+1):
                break
            if x not in ancestry:
                queue.append((x, reverse+(a,), D+a, ancestry+(x,)))
            a += 2
    return None, (), visits


def source_join(n, depth, forward=6):
    x, w, visits = n, (), 0
    for j in range(1, forward+1):
        x, a = step(x)
        w += (a,)
        m, v, nodes = ranked_inverse(x, n, depth, j, sum(w))
        visits += nodes
        if m is not None:
            return dict(n=n, m=m, w=w, v=v, join=x, visits=visits)
    return dict(n=n, m=None, visits=visits)


def independent_inverse_audit():
    cases = 0
    for N in range(3, 32, 2):
        R = rank(N)
        M = 1+cbrt(R[0]**2)
        for J in range(1, 32, 2):
            for depth in range(4):
                for r, A in ((0, 0), (1, 1), (2, 3)):
                    expected = False
                    for m in range(1, M+1, 2):
                        if rank(m) >= R:
                            continue
                        x, D = m, 0
                        for s in range(depth+1):
                            K = rank(m)[1]
                            good_slope = m == 1 or 2**(2*D)*3**(2*r+K)*4**R[1] <= \
                                2**(2*A)*3**(2*s+R[1])*4**K
                            if x == J and good_slope:
                                expected = True
                            x, a = step(x)
                            D += a
                    actual, _, _ = ranked_inverse(J, N, depth, r, A)
                    check((actual is not None) == expected, 'independent complete ranked inverse search')
                    cases += 1
    return cases


def normalize_child(N, m):
    """At a critical positive source, convert rank debt to ordinary size debt."""
    check(N > 1 and rank(N)[1] == 1 and 0 < m and rank(m) < rank(N),
          'normalization input')
    x, word = m, ()
    while x != 1 and rank(x)[1] >= 3:
        x, a = step(x)
        check(a == 2, 'root-centered normalization letter')
        word += (a,)
    if x != 1 and rank(x)[1] == 2:
        x, a = step(x)
        check(a >= 3, 'even stratum exit')
        word += (a,)
    check(x < N, 'rank descent pays original integer after explicit normalization')
    return dict(original=N, rank_child=m, smaller_successor=x, normalization_word=word)


def main():
    root = Path(__file__).resolve().parents[2]
    report = dict(status='PROVED scoped rank and branch laws; FINITE-EXACT; Collatz OPEN',
                  branch=branch_audit(), commutator_cases=address_commutator(), rank=rank_audit(),
                  independent_inverse_cases=independent_inverse_audit(), reflection=reflection_audit())
    prime = dict(n=223, m=233, w=(1,1,1,1,3), v=(2,1,1,1,2,3,1,1,2,1))
    report['larger_prime_child'] = rank_lift(prime)
    check((report['larger_prime_child']['P'], report['larger_prime_child']['Q']) == (62208,65536),
          'larger-child exact periods')
    inherited = json.loads((root/'05-knowledge/results/checked_switch_phase19_20261004.json').read_text())
    seeds = inherited['compiler'][-1]['seed_sources']
    check(len(seeds) == 239, 'frozen inherited universe')
    benchmarks, templates, normalizations = [], {}, []
    for depth in (2, 4, 6, 8):
        rows = [source_join(n, depth) for n in seeds]
        hits = [r for r in rows if r['m'] is not None]
        summary = dict(forward_cap=6, inverse_cap=depth, hits=len(hits),
                       larger_children=sum(r['m'] > r['n'] for r in hits),
                       inverse_visits=sum(r['visits'] for r in rows))
        benchmarks.append(summary)
        print('BENCHMARK '+json.dumps(summary), flush=True)
        for row in hits:
            if row['m'] == 1:
                continue
            if row['m'] > row['n'] and depth == 2:
                normalizations.append(normalize_child(row['n'], row['m']))
            key = (row['n'], row['m'], row['w'], row['v'])
            if key not in templates:
                templates[key] = rank_lift(row)
    report['benchmark'] = benchmarks
    report['ranked_families'] = list(templates.values())
    report['larger_child_normalizations'] = normalizations+[normalize_child(223,233)]
    report['checks'] = CHECKS
    stem = root/'05-knowledge/results/collatz_branch_toll_rank_20261004'
    stem.with_suffix('.json').write_text(json.dumps(report, indent=2)+'\n')
    lines = ['BRANCH TOLL AND ROOT-CENTERED RANK: universal coverage OPEN',
             'DELAYS '+str([(k,ell(k)) for k in range(1,13)]),
             'RANK '+json.dumps(report['rank']),
             'PRIME '+json.dumps(report['larger_prime_child'])]
    lines += ['BENCHMARK '+json.dumps(r) for r in benchmarks]
    lines += [f'PASS: {CHECKS} explicit checks.']
    output = '\n'.join(lines)+'\n'
    stem.with_suffix('.out').write_text(output)
    print(output, end='')


if __name__ == '__main__':
    main()
