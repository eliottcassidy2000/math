"""Exact depth-complete inverse joins, reset shields, and family transport.

Standard library, explicit checks survive -O. Frozen-source experiments and
independent all-smaller-source saturation audits; no universal coverage claim.
"""
from collections import deque
from fractions import Fraction
from hashlib import sha256
from pathlib import Path
import json

CHECKS = 0


def check(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise ArithmeticError(message)


def step(n):
    z = 3*n+1
    a = (z & -z).bit_length()-1
    return z >> a, a


def inverse_join(J, threshold, maximum, forward_length=0, forward_cost=0):
    """Exhaust all depths <= maximum that can yield 0<m<threshold.

    No arbitrary valuation cap. Even ideal exponent-one continuations bound
    every omitted branch. Return the first nonexpanding family witness.
    """
    queue = deque([(J, (), 0, (J,))])
    visits = 0
    while queue:
        z, reverse, A, ancestry = queue.popleft()
        d = len(reverse)
        visits += 1
        if 0 < z < threshold and 2**A*3**forward_length <= 2**forward_cost*3**d:
            return (z, tuple(reversed(reverse)), visits)
        remaining = maximum-d
        if remaining == 0 or z % 3 == 0:
            continue
        # Every further ancestor is >= (2/3)^remaining*(z+1)-1.
        if 2**remaining*(z+1) >= 3**remaining*(threshold+1):
            continue
        a = 1 if z % 3 == 2 else 2
        while True:
            x = (2**a*z-1)//3
            left = remaining-1
            if 2**left*(x+1) >= 3**left*(threshold+1):
                break
            if x not in ancestry:
                queue.append((x, reverse+(a,), A+a, ancestry+(x,)))
            a += 2
    return (None, (), visits)


def source_join(n, depth, forward_cap=6):
    x, w, A, visits = n, (), 0, 0
    for j in range(1, forward_cap+1):
        x, a = step(x)
        w += (a,)
        A += a
        m, v, count = inverse_join(x, n, depth, j, A)
        visits += count
        if m is not None:
            return {"n": n, "m": m, "w": w, "v": v, "join": x, "visits": visits}
    return {"n": n, "m": None, "visits": visits}


def compose(word):
    A = B = 0
    for a in word:
        B = 3*B+2**A
        A += a
    return len(word), A, B


def literal(n, word, sign=1):
    nodes = [n]
    for a in word:
        z, actual = 3*n+sign, 0
        while z % 2 == 0:
            z //= 2
            actual += 1
        check(actual == a, 'independent valuation replay')
        n = z
        nodes.append(n)
    return nodes


def family(row):
    n, m, w, v = row['n'], row['m'], tuple(row['w']), tuple(row['v'])
    r, A, B = compose(w)
    s, D, E = compose(v)
    P, Q = 2**(A+1)*3**max(s-r, 0), 2**(D+1)*3**max(r-s, 0)
    check(0 < m < n and P >= Q, 'all-height rank')
    for t in (0, 1, 2, 17, 10**6):
        x, y = n+P*t, m+Q*t
        J = literal(x, w)[-1]
        K = literal(y, v)[-1]
        check(J == K and 0 < y < x, 'lifted join and rank')
        check(2**A*J == 3**r*x+B and 2**D*K == 3**s*y+E,
              'independent affine identity')
    return {'source': n, 'child': m, 'source_period': P, 'child_period': Q,
            'w': w, 'v': v, 'slope': str(Fraction(Q, P)),
            'intercept': str(Fraction(P*m-Q*n, P)),
            'source_residue': n % P}


def independent_search_audit():
    """Forward-enumerate every possible smaller odd source, no inverse pruning."""
    cases = 0
    for n in range(3, 34, 2):
        for J in range(1, 66, 2):
            for depth in range(5):
                for r, A in ((0, 0), (1, 1), (2, 3), (3, 6)):
                    expected = False
                    for m in range(1, n, 2):
                        x, cost = m, 0
                        for s in range(depth+1):
                            if x == J and 2**cost*3**r <= 2**A*3**s:
                                expected = True
                            if s < depth:
                                x, a = step(x)
                                cost += a
                    actual, word, _ = inverse_join(J, n, depth, r, A)
                    check((actual is not None) == expected, 'complete inverse search audit')
                    if actual is not None:
                        check(literal(actual, word)[-1] == J, 'returned inverse witness')
                    cases += 1
    return cases


def shield_audit():
    sources = 0
    for sign in (1, -1):
        for q in range(2, 6):
            for t in range(1, 16, 2):
                n = 2**q*t-sign
                z = 3**q*t-sign
                if z % 4 != 2:
                    continue
                J = z//2
                w = (1,)*(q-1)+(2,)
                check(literal(n, w, sign)[-1] == J, 'reset-two source')
                # Complete finite audit of smaller possible ancestors; no valuation cap.
                for m in range(1, n, 2):
                    x, path = m, [m]
                    for d in range(1, q+4):
                        z = 3*x+sign
                        while z % 2 == 0:
                            z //= 2
                        x = z
                        path.append(x)
                        if x == J:
                            check(d > q and path[d-q] == n,
                                  'short smaller ancestor factors through original source')
                sources += 1
    sharp = {'n': 283, 'm': 223, 'w': (1, 2), 'v': (1, 1, 1, 1, 3, 2)}
    sharp_family = family(sharp)
    check(literal(223, sharp['v'])[-3] != 283, 'noncanonical sharp witness')
    check(sharp_family['source_period'] == 1296 and sharp_family['child_period'] == 1024,
          'sharp shield family')
    return {'signed_sources': sources, 'sharp_family': sharp_family}


def sibling_run_audit():
    rows = []
    for k in range(9):
        ell = 1
        while 3**ell <= 2**(ell+2*k):
            ell += 1
        modulus = 3**ell
        c = (4**k+2)//3
        residue = -c*pow(4**k, -1, modulus) % modulus
        n = residue if residue % 2 else residue+modulus
        m = 2**ell*(4**k*n+c)//modulus-1
        sibling = 4**k*n+(4**k-1)//3
        check(literal(m, (1,)*ell)[-1] == sibling and step(sibling)[0] == step(n)[0],
              'sibling-run common future')
        check(0 < m < n, 'sharp sibling-run rank')
        rows.append({'k': k, 'ell': ell, 'modulus': modulus, 'residue': residue,
                     'least_odd_source': n, 'child': m, 'gap': modulus-2**(ell+2*k)})
    # The rank sign cannot be discarded on the minus sheet.
    for n, w in ((5, (1, 2)), (17, (1, 1, 1, 2, 1, 1, 4))):
        r, A, B = compose(w)
        check(literal(n, w, -1)[-1] == n and 2**A < 3**r,
              'minus-cycle inverse contraction is not strict rank')
    return rows


def full_paths(limit):
    paths = {}
    for n in range(1, limit+1, 2):
        x, nodes, word = n, [n], []
        while x != 1:
            x, a = step(x)
            nodes.append(x)
            word.append(a)
            check(len(word) <= 1000, 'declared finite trajectory cap')
        paths[n] = (nodes, word)
    return paths


def saturation_audit(seeds, limit=10000):
    """Full routes are audit inputs, never supplied to bounded inverse search."""
    paths = full_paths(limit)
    records, owner = {}, {}
    all_depth, clocks, excluded = [], [], []
    for n, (nodes, word) in paths.items():
        if n in seeds:
            kappa = next(j for j, J in enumerate(nodes[1:], 1) if J in owner)
            tau = next(j for j, J in enumerate(nodes[1:], 1) if J < n)
            clocks.append({'n': n, 'coalescence': kappa, 'first_descent': tau,
                           'owner': owner[nodes[kappa]]})
            A, candidates, any_smaller = 0, [], False
            for r, (J, a) in enumerate(zip(nodes[1:7], word[:6]), 1):
                A += a
                for m, s, B in records.get(J, ()):
                    any_smaller = True
                    if 2**B*3**r <= 2**A*3**s:
                        candidates.append((s, r, m, J))
            if candidates:
                s, r, m, J = min(candidates)
                all_depth.append({'n': n, 'm': m, 'w': tuple(word[:r]),
                                  'v': tuple(paths[m][1][:s]), 'join': J})
            else:
                check(not any_smaller, 'remaining cases have no smaller ancestor, even without slope filter')
                excluded.append(n)
        A = 0
        for s, J in enumerate(nodes):
            if s:
                A += word[s-1]
            owner.setdefault(J, n)
            records.setdefault(J, []).append((n, s, A))
    check(len(all_depth) == 17 and len(excluded) == 222, 'all-depth saturation')
    profile = [(F, sum(row['coalescence'] <= F for row in clocks))
               for F in (1, 2, 4, 6, 8, 12, 16, 24, 32, 48, 64)]
    # Audit the expanding-lift hostile: a real smaller source alone is insufficient.
    nodes, word = paths[231]
    s = nodes.index(233)
    check(s == 17 and sum(word[:s]) == 27, '231 to 233 hostile')
    check(2**27 > 3**17, 'expanding inverse coefficient hostile')
    return {'full_route_inputs': len(paths), 'distinct_vertices': len(owner),
            'all_depth_hits': all_depth, 'all_depth_excluded': excluded,
            'forward_profile': profile, 'coalescence_rows': clocks,
            'expanding_hostile_word': word[:s]}


def first_merge_probe(limit=100000):
    """A conjecture probe, not a claimed theorem or input to the compiler."""
    best, owner = {}, {}
    failures, shortened, unseen_count, already_seen = [], 0, 0, 0
    max_saving = (0, None)
    for n in range(1, limit+1, 2):
        unseen = n not in owner
        x, nodes, word, A = n, [n], [], 0
        kappa = tau = data = None
        if n > 1 and not unseen:
            q, m, inverse_depth = best[n]
            kappa, data = 0, (q, m, inverse_depth, Fraction(1))
        while x != 1:
            x, a = step(x)
            nodes.append(x)
            word.append(a)
            A += a
            check(len(word) <= 1000, 'first-merge finite trajectory cap')
            if tau is None and x < n:
                tau = len(word)
            if kappa is None and x in owner:
                kappa = len(word)
                q, m, inverse_depth = best[x]
                data = (q, m, inverse_depth, Fraction(2**A, 3**kappa))
        if n > 1:
            if data[0] > data[3]:
                failures.append(n)
            if unseen:
                unseen_count += 1
                shortened += kappa < tau
                if tau-kappa > max_saving[0]:
                    max_saving = (tau-kappa, {'n': n, 'coalescence': kappa,
                                            'first_descent': tau, 'child': data[1],
                                            'inverse_depth': data[2]})
            else:
                already_seen += 1
        A = 0
        for r, J in enumerate(nodes):
            if r:
                A += word[r-1]
            q = Fraction(2**A, 3**r)
            owner.setdefault(J, n)
            if J not in best or q < best[J][0]:
                best[J] = (q, n, r)
    check(not failures and (unseen_count, shortened, already_seen) == (26599, 2419, 23400),
          'finite first-merge signal')
    return {'odd_nonroot_sources': (limit-1)//2, 'failed_nonexpanding_lift': failures,
            'previously_unseen_sources': unseen_count, 'already_in_smaller_route_union': already_seen,
            'unseen_with_earlier_than_actual_descent': shortened, 'maximum_unseen_saving': max_saving}


def beyond_short_bank():
    """Strict extension of a complete depth-six bank, with exact family guards."""
    for h in (0, 1, 2, 17, 10**6):
        n, m = 4647+69984*h, 4351+65536*h
        literal(n, (1, 1, 2))
        check(literal(n, (1,))[-1] == literal(m, (1,)*7+(5,))[-1],
              'hard-reset sibling subfamily')
    missed = source_join(144615, 6, 6)
    found = source_join(144615, 8, 1)
    check(missed['m'] is None and found['m'] == 101567,
          'strict extension beyond complete forward-six inverse-six bank')
    check(found['w'] == (1,) and found['v'] == (1,)*5+(2, 3),
          'seven-edge alternative')
    lifted = family({'n': 1731, 'm': 1215, 'w': (1,), 'v': (1,)*5+(2, 3)})
    check((lifted['source_period'], lifted['child_period']) == (2916, 2048),
          'strict-extension periods')
    return {'missed': missed, 'found': found, 'lift': lifted}


def main():
    root = Path(__file__).resolve().parents[2]
    inherited = json.loads((root/'05-knowledge/results/checked_switch_phase19_20261004.json').read_text())
    seeds = inherited['compiler'][-1]['seed_sources']
    check(len(seeds) == 239 and len(set(seeds)) == 239, 'frozen source universe')
    report = {'status': 'PROVED scoped mechanisms; FINITE-EXACT; universal coverage OPEN',
              'seed_sha256': sha256(json.dumps(seeds, separators=(',', ':')).encode()).hexdigest(),
              'independent_search_cases': independent_search_audit(),
              'shield': shield_audit(), 'sibling_runs': sibling_run_audit(),
              'strict_extension': beyond_short_bank()}
    lines, matrix, templates = [], [], {}
    for F in (6, 8, 12, 16):
        for D in (6, 8, 12, 16):
            rows = [source_join(n, D, F) for n in seeds]
            hits = [row for row in rows if row['m'] is not None]
            summary = {'forward_cap': F, 'inverse_cap': D, 'resolved': len(hits),
                       'visited_inverse_nodes': sum(row['visits'] for row in rows)}
            matrix.append(summary)
            lines.append('SEARCH '+json.dumps(summary, sort_keys=True))
            for row in hits:
                key = (tuple(row['w']), tuple(row['v']), row['n'], row['m'])
                if key not in templates:
                    templates[key] = family(row)
    report['search_matrix'] = matrix
    report['family_witnesses'] = list(templates.values())
    deeper = []
    for D in (20, 24, 28):
        rows = [source_join(n, D) for n in seeds]
        row = {'forward_cap': 6, 'inverse_cap': D,
               'resolved': sum(r['m'] is not None for r in rows),
               'visited_inverse_nodes': sum(r['visits'] for r in rows)}
        deeper.append(row)
        lines.append('DEEPER '+json.dumps(row, sort_keys=True))
    report['deeper_inverse'] = deeper
    report['saturation'] = saturation_audit(set(seeds))
    report['first_merge_probe'] = first_merge_probe()
    for row in report['saturation']['all_depth_hits']:
        family(row)
    lines.append('ALL_DEPTH '+json.dumps({'first_six_resolved': 17, 'first_six_excluded': 222,
                                         'forward_profile': report['saturation']['forward_profile']}))
    lines.append('FIRST_MERGE '+json.dumps(report['first_merge_probe'], sort_keys=True))
    report['checks'] = CHECKS
    lines.append(f'PASS: {CHECKS} explicit checks; no assertions removed by -O.')
    directory = root/'05-knowledge/results'
    stem = 'collatz_join_shields_and_lifts_20261004'
    (directory/(stem+'.json')).write_text(json.dumps(report, indent=2)+'\n')
    out = '\n'.join(lines)+'\n'
    (directory/(stem+'.out')).write_text(out)
    print(out, end='')


if __name__ == '__main__':
    main()
