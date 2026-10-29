"""Exact finite-clock recovery and infinite-clock hostile; standard library.

Run normally and with -O. All checks remain active. Imports only the adjacent
previous exact matrix helpers; no network, floating arithmetic or randomness.
"""
from fractions import Fraction as Q
from itertools import combinations, permutations, product
from math import comb, gcd
from pathlib import Path
import json
from port_queries import mm, tr, add, invdet

GATES = 0


def need(ok, label):
    global GATES
    GATES += 1
    if not ok:
        raise RuntimeError(label)


def zero(n):
    return [[0] * n for _ in range(n)]


def eye(n):
    return [[int(i == j) for j in range(n)] for i in range(n)]


def scale(a, q):
    return [[q * x for x in row] for row in a]


def trace(a):
    return sum(a[i][i] for i in range(len(a)))


def pf(a):
    if not a:
        return Q(1)
    return sum((-1)**(j+1) * a[0][j] * pf([
        [a[r][s] for s in range(1, len(a)) if s != j]
        for r in range(1, len(a)) if r != j])
        for j in range(1, len(a)))


def rank(a):
    if not a:
        return 0
    z = [[Q(x) for x in row] for row in a]
    r = 0
    for j in range(len(z[0])):
        k = next((k for k in range(r, len(z)) if z[k][j]), None)
        if k is None:
            continue
        z[r], z[k] = z[k], z[r]
        q = z[r][j]
        z[r] = [v / q for v in z[r]]
        for k in range(r + 1, len(z)):
            q = z[k][j]
            z[k] = [x - q*y for x, y in zip(z[k], z[r])]
        r += 1
        if r == len(z):
            break
    return r


def elementary(n, i, j):
    a = zero(n)
    a[i][j], a[j][i] = 1, -1
    return a


def block_diag(*blocks):
    a = zero(sum(len(b) for b in blocks))
    start = 0
    for b in blocks:
        for i, row in enumerate(b):
            a[start+i][start:start+len(b)] = row
        start += len(b)
    return a


def perm_matrix(perm):
    return [[int(i == perm[j]) for j in range(len(perm))]
            for i in range(len(perm))]


def cyclic(c):
    ans, a = [], eye(len(c))
    while a not in ans:
        ans.append(a)
        a = mm(a, c)
        if len(ans) > 100:
            raise RuntimeError('finite group constructor exceeded declared ceiling')
    need(a == ans[0], 'cyclic closure')
    return ans


def derivative_at_zero(values):
    n = len(values) - 1
    return (sum(Q((-1)**(t+1) * comb(n, t), t) * values[t]
                for t in range(1, n+1))
            - sum(Q(1, t) for t in range(1, n+1)) * values[0])


def orbit_probes(group):
    n = len(group[0])
    ans = []
    for i, j in combinations(range(n), 2):
        f, e = elementary(n, i, j), zero(n)
        for c in group:
            e = add(e, mm(c, mm(f, tr(c))))
        ans.append((i, j, e))
    return ans


def group_dimension(group):
    n = len(group[0])
    pairs = list(combinations(range(n), 2))
    basis = [elementary(n, i, j) for i, j in pairs]
    # Independent fixed-space equations, rather than just the orbit span.
    constraints = []
    for c in group:
        columns = [add(mm(c, mm(f, tr(c))), scale(f, -1)) for f in basis]
        constraints += [[a[i][j] for a in columns] for i, j in pairs]
    fixed_dim = len(pairs) - rank(constraints)
    orbit_dim = rank([[e[i][j] for i, j in pairs] for _, _, e in orbit_probes(group)])
    character_dim = Q(sum(trace(c)**2 - trace(mm(c, c)) for c in group), 2*len(group))
    need(fixed_dim == orbit_dim == character_dim, 'three invariant dimension paths')
    return fixed_dim


def check_state(name, m, ports, group):
    mi, delta = invdet(m)
    need(delta != 0 and delta == pf(m)**2, 'nondegenerate alternating source')
    g = mm(tr(ports), mm(mi, ports))
    p = len(g)
    for c in group:
        need(mm(tr(c), mm(g, c)) == g, 'actual port Gram invariance')
    singular = queries = 0
    for i, j, e in orbit_probes(group):
        for c in group:
            need(mm(c, mm(e, tr(c))) == e, 'integer orbit probe is legal')
        need(trace(mm(g, e)) == -2*len(group)*g[i][j], 'invariant trace separation')
        values = []
        for t in range(p + 1):
            k = add(m, scale(mm(ports, mm(e, tr(ports))), t))
            ki, kd = invdet(k)
            need(kd == pf(k)**2, 'Gaussian determinant versus matching Pfaffian')
            values.append(kd/delta)
            queries += 1
            singular += int(not kd)
            if kd:
                xi, xd = invdet(add(eye(p), scale(mm(g, e), t)))
                need(xd == kd/delta, 'determinant lemma')
                next_g = mm(tr(ports), mm(ki, ports))
                need(mm(xi, g) == next_g, 'clocked resolvent')
                for c in group:
                    need(mm(tr(c), mm(next_g, c)) == next_g, 'clock survives update')
        recovered = -derivative_at_zero(values) / (2*len(group))
        need(recovered == g[i][j], 'actual determinant values recover entry')
    return {'name': name, 'ports': p, 'clock_order': len(group),
            'invariant_dimension': group_dimension(group),
            'determinant_queries': queries, 'singular_queries': singular}


def partitions(n, least=1):
    if n == 0:
        yield ()
    else:
        for k in range(least, n+1):
            for tail in partitions(n-k, k):
                yield (k,) + tail


def cycle_dimension(lengths):
    return (sum((n-1)//2 for n in lengths)
            + sum(gcd(a, b) for a, b in combinations(lengths, 2)))


def pair_orbit_dimension(lengths):
    perm, start = [], 0
    for n in lengths:
        perm += list(range(start+1, start+n)) + [start]
        start += n
    unseen = set(combinations(range(start), 2))
    dim = 0
    while unseen:
        first = min(unseen)
        pair, sign = first, 1
        while True:
            unseen.discard(pair)
            a, b = perm[pair[0]], perm[pair[1]]
            if a > b:
                a, b, sign = b, a, -sign
            pair = (a, b)
            if pair == first:
                break
        dim += int(sign == 1)
    return dim


def unipotent_hostile():
    u = [[1, 1], [0, 1]]
    c = block_diag(u, u)
    samples = 0
    for d in range(-2, 3):
        g = block_diag([[0, 1], [-1, 0]], [[0, 1], [-1, 0]])
        g[1][3], g[3][1] = d, -d
        m, delta = invdet(g)
        need(delta == 1 and mm(tr(c), mm(g, c)) == g, 'unipotent invariant states')
        need(mm(c, mm(m, tr(c))) == m, 'unipotent source isometry')
        for a, b, uu, v in product(range(-2, 3), repeat=4):
            e = [[0, a, uu, v], [-a, 0, -v, 0],
                 [-uu, v, 0, b], [-v, 0, -b, 0]]
            need(mm(c, mm(e, tr(c))) == e, 'unipotent legal coupling')
            for sign in (-1, 1):
                k = add(scale(m, sign), e)
                need(pf(k) == (a-sign)*(b-sign)-v*v, 'symbolic hostile specialization')
                need(invdet(k)[1] == ((a-sign)*(b-sign)-v*v)**2, 'direct hostile determinant')
                samples += 1
    return samples


def word_and_four_query_checks(m4, cg):
    # Literal port re-marking (M fixed), followed by legal attachments and
    # a sign reversal. This is separate from changing ambient coordinates.
    m, ports, group = m4, eye(4), cg
    mi, delta = invdet(m)
    g = mi
    r = eye(4); r[0][1] = 2; r[2][3] = -1
    ri, _ = invdet(r)
    ports = mm(ports, r)
    group = [mm(ri, mm(c, r)) for c in group]
    g = mm(tr(r), mm(g, r))
    for c, cc in zip(cg, group):
        need(mm(c, ports) == mm(ports, cc), 'literal dual port transport')
    for step in range(4):
        if step == 2:
            m, g = scale(m, -1), scale(g, -1)
        e = orbit_probes(group)[step][2]
        # Avoid only an already detected singular branch; this is a
        # deterministic choice for testing nonterminal composition.
        t = 1
        while not invdet(add(m, scale(mm(ports, mm(e, tr(ports))), t)))[1]:
            t += 1
        m = add(m, scale(mm(ports, mm(e, tr(ports))), t))
        inv_step, ratio = invdet(add(eye(4), scale(mm(g, e), t)))
        delta, g = delta*ratio, mm(inv_step, g)
        mi, direct_delta = invdet(m)
        need(delta == direct_delta and g == mm(tr(ports), mm(mi, ports)),
             'full clocked word matches direct state')
    # Four signed queries reconstruct the two C3 coordinates.
    mi, delta = invdet(m4)
    for i, j in [(0,1), (1,2)]:
        e = next(e for ii,jj,e in orbit_probes(cg) if (ii,jj) == (i,j))
        dp = invdet(add(m4, e))[1]
        dm = invdet(add(m4, scale(e, -1)))[1]
        need((dm-dp)/(12*delta) == mi[i][j], 'four-query C3 recovery')
    # The rank-one averaging trap is not an invariant unknown.
    f = elementary(4, 0, 1)
    averaged = next(e for i,j,e in orbit_probes(cg) if (i,j) == (0,1))
    need(mm(cg[1], mm(f, tr(cg[1]))) != f, 'ordinary pair probe is illegal')
    need(trace(mm(f, averaged)) != -2*len(cg)*f[0][1],
         'dropping invariant-unknown hypothesis is rejected')


def main():
    c4 = perm_matrix([0, 2, 3, 1])
    cg = cyclic(c4)
    m4 = [[0,1,1,1],[-1,0,-1,1],[-1,1,0,-1],[-1,-1,1,0]]
    states = []
    for a, b in product((-2, -1, 1, 2), repeat=2):
        m = [[0,a,a,a],[-a,0,b,-b],[-a,-b,0,b],[-a,b,-b,0]]
        states.append(check_state(f'C3_star_{a}_triangle_{b}', m, eye(4), cg))
    states.append(check_state('two_C3_blocks', block_diag(m4, m4), eye(8),
                              [block_diag(c,c) for c in cg]))
    # Nonorthogonal integral clock obtained by a genuine unimodular re-mark.
    r = eye(4); r[0][1] = 2; r[2][3] = -1
    ri, rd = invdet(r)
    need(rd == 1, 'unimodular remark')
    transformed_m = mm(ri, mm(m4, tr(ri)))
    transformed_clock = [mm(ri, mm(c, r)) for c in cg]
    states.append(check_state('nonorthogonal_C3', transformed_m, eye(4), transformed_clock))
    # Dependent port declaration uses the same sufficient legality rule.
    pdep = [row + [int(i == 0)] for i, row in enumerate(eye(4))]
    states.append(check_state('dependent_fifth_port', m4, pdep,
                              [block_diag(c, [[1]]) for c in cg]))
    symplectic8 = [[0]*4 + row for row in eye(4)] + [row + [0]*4 for row in scale(eye(4), -1)]
    isotropic_ports = eye(4) + [[0]*4 for _ in range(4)]
    states.append(check_state('independent_isotropic_ports', symplectic8, isotropic_ports, cg))
    for order, c in [(3, [[0,-1],[1,-1]]), (4, [[0,-1],[1,0]]),
                     (6, [[0,-1],[1,1]])]:
        group = cyclic(c)
        need(len(group) == order, 'specified clock order')
        states.append(check_state(f'plane_order_{order}', [[0,1],[-1,0]], eye(2), group))
    s3 = [block_diag(perm_matrix(p), perm_matrix(p)) for p in permutations(range(3))]
    source = [[0]*3 + row for row in eye(3)] + [row + [0]*3 for row in scale(eye(3), -1)]
    states.append(check_state('noncommutative_S3_weighted_form', source, eye(6), s3))
    partition_tests = 0
    for n in range(1, 15):
        for lengths in partitions(n):
            need(cycle_dimension(lengths) == pair_orbit_dimension(lengths), 'cycle formula versus signed pair orbits')
            partition_tests += 1
    large_rows = []
    for n in (5,7,11,13,25,27,57):
        for lengths in [(n,), (1,n)]:
            dim = cycle_dimension(lengths)
            need(dim == pair_orbit_dimension(lengths), 'named long cycle exact orbit test')
            large_rows.append({'cycles': lengths, 'dimension': dim})
    hostile_cases = unipotent_hostile()
    word_and_four_query_checks(m4, cg)
    out = {'status': 'FINITE-EXACT supporting scoped proof', 'gates': GATES,
           'finite_states': states, 'cycle_partitions_through_14': partition_tests,
           'named_cycle_dimensions': large_rows, 'unipotent_specializations': hostile_cases,
           'unipotent_box': {'d': [-2,2], 'a_b_u_v': [-2,2], 'signs': [-1,1]},
           'scope': 'all integer invariant couplings; not native sign-cut closure or integral cokernel'}
    Path(__file__).with_suffix('.json').write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps(out, indent=2))


if __name__ == '__main__':
    main()
