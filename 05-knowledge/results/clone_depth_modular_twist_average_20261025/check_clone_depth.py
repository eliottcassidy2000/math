"""Clone-fibre cancellation depth over Q and modulo primes: the two 50-vertex tournaments and the p-group law.

Standard library only; every gate is an explicit exception (active under -O).
Run:  python3 check_clone_depth.py > clone-normal.json ; python3 -O check_clone_depth.py > clone-optimized.json
"""
import itertools, json, sys
from collections import Counter
from math import comb


def need(ok, msg):
    if not ok:
        raise RuntimeError(msg)


# ---------------------------------------------------------------- tournaments: automorphisms, pairs, inflation
def refine(n, out, colours):
    while True:
        sig = []
        for v in range(n):
            o = sorted(colours[w] for w in out[v])
            i = sorted(colours[w] for w in range(n) if v in out[w])
            sig.append((colours[v], tuple(o), tuple(i)))
        keys = sorted(set(sig))
        new = [keys.index(x) for x in sig]
        if len(set(new)) == len(set(colours)):
            return new
        colours = new


def automorphisms(n, out):
    """All automorphisms of a tournament (out-neighbourhood sets), by colour refinement and individualisation.
    Nothing about the group is assumed; every candidate map is checked on every ordered pair."""
    adj = [[1 if w in out[v] else 0 for w in range(n)] for v in range(n)]
    base = refine(n, out, [0] * n)
    results = []

    def search(cs, ct):
        classes, tclasses = {}, {}
        for v in range(n):
            classes.setdefault(cs[v], []).append(v)
            tclasses.setdefault(ct[v], []).append(v)
        if sorted(classes) != sorted(tclasses) or any(len(classes[c]) != len(tclasses[c]) for c in classes):
            return
        if all(len(c) == 1 for c in classes.values()):
            perm = [None] * n
            for c, vs in classes.items():
                perm[vs[0]] = tclasses[c][0]
            if all(adj[v][w] == adj[perm[v]][perm[w]] for v in range(n) for w in range(n) if v != w):
                results.append(tuple(perm))
            return
        c = min((len(vs), c) for c, vs in classes.items() if len(vs) > 1)[1]
        v = classes[c][0]
        m = max(cs) + 1
        cs2 = list(cs); cs2[v] = m; cs2 = refine(n, out, cs2)
        for w in tclasses[c]:
            ct2 = list(ct); ct2[w] = m; ct2 = refine(n, out, ct2)
            search(cs2, ct2)

    search(base, base)
    return sorted(set(results))


def homogeneous_pairs(n, out):
    return [(a, b) for a in range(n) for b in range(a + 1, n)
            if all(((c in out[a]) == (c in out[b])) and ((a in out[c]) == (b in out[c])) for c in range(n) if c not in (a, b))]


def inflate_tt2(n, out):
    """Inflate every vertex to the transitive pair TT2 (first -> second)."""
    out2 = [set() for _ in range(2 * n)]
    for v in range(n):
        for w in out[v]:
            for a in (2 * v, 2 * v + 1):
                for b in (2 * w, 2 * w + 1):
                    out2[a].add(b)
        out2[2 * v].add(2 * v + 1)
    return 2 * n, out2


def check_tournament(n, out):
    need(all(v not in out[v] for v in range(n)), 'no loops')
    need(all((w in out[v]) != (v in out[w]) for v in range(n) for w in range(n) if v != w), 'exactly one arc per pair')


# ---------------------------------------------------------------- polynomials and depths
def polymul(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            out[i + j] += x * y
    return out


def cycle_lengths(perm):
    seen, cs = set(), []
    for i in range(len(perm)):
        if i in seen:
            continue
        j, c = i, 0
        while j not in seen:
            seen.add(j); c += 1; j = perm[j]
        cs.append(c)
    return sorted(cs)


def burnside_subset_polynomial(perms):
    """P_A(t) = (1/|A|) sum_g prod_cycles (1 + t^len), from explicit permutations (the group must be closed)."""
    n = len(perms[0]); total = [0] * (n + 1)
    for g in perms:
        f = [1]
        for L in cycle_lengths(g):
            fac = [0] * (L + 1); fac[0] = fac[-1] = 1
            f = polymul(f, fac)
        total = [x + y for x, y in zip(total, f)]
    need(all(c % len(perms) == 0 for c in total), 'Burnside average is integral')
    return [c // len(perms) for c in total]


def shifted(P):
    """Coefficients of P(-1+x) in x."""
    return [sum(P[k] * comb(k, j) * (-1) ** (k - j) for k in range(j, len(P))) for j in range(len(P))]


def depth(Q, ell=None):
    for j, c in enumerate(Q):
        if (c != 0) if ell is None else (c % ell != 0):
            return j
    return None


def valuation(c, ell):
    if c == 0:
        return None
    v = 0
    while c % ell == 0:
        c //= ell; v += 1
    return v


PRIMES = [2, 3, 5, 7, 11, 13]


def analyze_action(perms, name):
    n = len(perms[0]); G = set(perms)
    ident = tuple(range(n))
    need(ident in G, name + ': identity present')
    need(all(tuple(g[h[i]] for i in range(n)) in G for g in perms for h in perms), name + ': closure')
    need(len(G) % 2 == 1, name + ': odd order')
    P = burnside_subset_polynomial(perms)
    Q = shifted(P)
    d = depth(Q)
    kappa = Q[d]
    cyc = Counter(len(cycle_lengths(g)) for g in perms)
    dmin = min(cyc)
    need(d == dmin, name + ': depth over Q equals minimum cycle count (Session62 theorem)')
    lead = sum(1 for _ in ())  # placeholder
    pred_kappa = sum(eval('*'.join(map(str, cycle_lengths(g)))) for g in perms if len(cycle_lengths(g)) == d)
    need(pred_kappa % len(perms) == 0 and pred_kappa // len(perms) == kappa, name + ': kappa formula')
    orbits = {frozenset(g[i] for g in perms) for i in range(n)}
    r = len(orbits)
    # signed binomial moments (7): vanish below d, equal kappa at d  (n ≡ d mod 2 so sign is +1)
    for j in range(d + 1):
        moment = sum(((-1) ** (n - k)) * comb(k, j) * P[k] for k in range(n + 1))
        need(moment == (0 if j < d else kappa), name + ': signed moment %d' % j)
    depths = {ell: depth(Q, ell) for ell in PRIMES}
    for ell in PRIMES:
        need(depths[ell] >= d, name + ': modular depth at least d')
        need((depths[ell] > d) == (kappa % ell == 0), name + ': Lemma B0 (jump iff ell | kappa) at %d' % ell)
    return {'name': name, 'points': n, 'group_order': len(perms), 'orbits_r': r, 'depth_Q': d, 'kappa': kappa,
            'cycle_count_census': dict(sorted(cyc.items())), 'P': P, 'orbit_count': sum(P),
            'depths_mod': depths, 'valuations': {ell: [valuation(c, ell) for c in Q] for ell in (2, 3, 5)} if n <= 30 else None}


# ---------------------------------------------------------------- the two 50-vertex tournaments
def z5sq_cayley():
    p = 5
    elts = [(a, b) for a in range(p) for b in range(p)]
    idx = {e: i for i, e in enumerate(elts)}
    S = [(0, 1), (0, 2), (1, 0), (1, 1), (1, 2), (1, 3), (1, 4), (2, 0), (2, 1), (2, 2), (2, 3), (3, 1)]
    need(all(((-s[0]) % p, (-s[1]) % p) not in S for s in S) and len(S) == 12, 'S is a half-system of Z5^2')
    out = [set(idx[((v[0] + s[0]) % p, (v[1] + s[1]) % p)] for s in S) for v in elts]
    translations = [tuple(idx[((e[0] + c) % p, (e[1] + d) % p)] for e in elts) for c in range(p) for d in range(p)]
    return out, translations, {'group': 'Z5xZ5', 'connection_set': S}


def z25_circulant():
    n = 25; S = list(range(1, 13))
    out = [set((v + s) % n for s in S) for v in range(n)]
    rotations = [tuple((i + c) % n for i in range(n)) for c in range(n)]
    return out, rotations, {'group': 'Z25', 'connection_set': S}


def tournament_record(builder):
    out, expected_group, info = builder()
    n = len(out)
    check_tournament(n, out)
    autR = automorphisms(n, out)
    need(len(autR) == 25 and set(autR) == set(expected_group), info['group'] + ': Aut of the 25-vertex base is exactly the regular group')
    need(homogeneous_pairs(n, out) == [], info['group'] + ': base is pair-free')
    N, T = inflate_tt2(n, out)
    check_tournament(N, T)
    autT = automorphisms(N, T)
    need(len(autT) == 25, info['group'] + ': |Aut| of the inflated 50-vertex tournament is 25')
    pairs = homogeneous_pairs(N, T)
    need(pairs == [(2 * v, 2 * v + 1) for v in range(n)], info['group'] + ': homogeneous pairs are exactly the 25 clone blocks')
    # gap action: block v -> block g(2v)//2, and g must preserve the inner order
    gap_perms = []
    for g in autT:
        need(all(g[2 * v + 1] == g[2 * v] + 1 and g[2 * v] % 2 == 0 for v in range(n)), 'automorphisms preserve blocks and their order')
        gap_perms.append(tuple(g[2 * v] // 2 for v in range(n)))
    gap_perms = sorted(set(gap_perms))
    need(len(gap_perms) == 25 and set(gap_perms) == set(expected_group), info['group'] + ': the gap action is the regular group')
    rec = analyze_action(gap_perms, info['group'] + ' gap action (50-vertex tournament)')
    rec.update(info); rec['vertices'] = N; rec['aut_order'] = len(autT); rec['gaps'] = n
    need(rec['orbits_r'] == 1, 'one invariant toggle direction')
    return rec


# ---------------------------------------------------------------- abstract regular actions: p-group law and coprime law
def regular_perms(elts, mul):
    idx = {e: i for i, e in enumerate(elts)}
    return [tuple(idx[mul(g, h)] for h in elts) for g in elts]


def abelian(moduli):
    elts = list(itertools.product(*[range(m) for m in moduli]))
    return elts, (lambda g, h: tuple((a + b) % m for a, b, m in zip(g, h, moduli)))


def semidirect(m, k, r):
    """C_m x| C_k with generator acting by multiplication by r (r^k = 1 mod m)."""
    need(pow(r, k, m) == 1, 'semidirect action well defined')
    elts = list(itertools.product(range(m), range(k)))
    return elts, (lambda g, h: ((g[0] + pow(r, g[1], m) * h[0]) % m, (g[1] + h[1]) % k))


def heisenberg(p):
    elts = list(itertools.product(range(p), repeat=3))
    return elts, (lambda g, h: ((g[0] + h[0]) % p, (g[1] + h[1]) % p, (g[2] + h[2] + g[0] * h[1]) % p))


def element_orders(elts, mul):
    idx = {e: i for i, e in enumerate(elts)}
    ident = next(e for e in elts if all(mul(e, h) == h for h in elts[:3]))
    orders = Counter()
    for g in elts:
        x, o = g, 1
        while x != ident:
            x = mul(g, x); o += 1
        orders[o] += 1
    return orders


def p_group_prediction(orders, p, N):
    e = 0
    while p ** e < N:
        e += 1
    need(p ** e == N, 'p-group')
    c = {j: orders.get(p ** j, 0) // ((p - 1) * p ** (j - 1)) for j in range(1, e + 1)}
    need(c[1] % p == 1, 'number of subgroups of order p is 1 mod p')
    J = max(j for j in c if c[j] % p != 0)
    return {'e': e, 'cyclic_subgroup_counts': c, 'J': J, 'predicted_depth_mod_p': N - p ** J + p ** (J - 1)}


def coprime_prediction(orders, N, ell):
    """Theorem B' for regular actions: depth mod ell = min{N/o : ell does not divide e_o}."""
    return min(N // o for o, cnt in orders.items() if cnt % ell != 0)


def mixed_candidates(orders, p, N):
    """Candidate orders for depth mod p of a regular action (Theorem B'' lower bound; equality is HYP-9374).
    E_{o'} = elements whose p'-part has order o'; e_{j,o'} = elements of order p^j o'."""
    e = 0
    while N % (p ** (e + 1)) == 0:
        e += 1
    m = N // p ** e
    cands = []
    for o_prime in [d for d in range(1, m + 1) if m % d == 0]:
        E = 0
        for o, cnt in orders.items():
            pp = o
            while pp % p == 0:
                pp //= p
            if pp == o_prime:
                E += cnt
        need(E % p ** e == 0, 'E_{o\'} divisible by p^e (Frobenius)')
        if (E // p ** e) % p != 0:
            cands.append(N // o_prime)
        for j in range(1, e + 1):
            cnt = orders.get(p ** j * o_prime, 0)
            if cnt:
                need(cnt % p ** (j - 1) == 0, 'e_{j,o\'} divisible by p^(j-1)')
                if (cnt // p ** (j - 1)) % p != 0:
                    cands.append(N // o_prime - p ** j + p ** (j - 1))
    return min(cands)


GROUPS = [
    ('C3xC3', 3, lambda: abelian((3, 3))), ('C9', 3, lambda: abelian((9,))),
    ('C5xC5', 5, lambda: abelian((5, 5))), ('C25', 5, lambda: abelian((25,))),
    ('C7xC7', 7, lambda: abelian((7, 7))),
    ('C27', 3, lambda: abelian((27,))), ('C9xC3', 3, lambda: abelian((9, 3))), ('C3xC3xC3', 3, lambda: abelian((3, 3, 3))),
    ('Heis27', 3, lambda: heisenberg(3)), ('M27=C9x|C3', 3, lambda: semidirect(9, 3, 4)),
    ('M81=C27x|C3', 3, lambda: semidirect(27, 3, 10)), ('C9x|C9', 3, lambda: semidirect(9, 9, 4)),
    ('C125', 5, lambda: abelian((125,))), ('Heis125', 5, lambda: heisenberg(5)), ('M125=C25x|C5', 5, lambda: semidirect(25, 5, 6)),
]
MIXED = [('C15', lambda: abelian((15,))), ('C21', lambda: abelian((21,))), ('C45', lambda: abelian((45,))),
         ('C7x|C3', lambda: semidirect(7, 3, 2)), ('C3xC5xC5', lambda: abelian((3, 5, 5)))]


def group_record(name, p, builder):
    elts, mul = builder()
    perms = regular_perms(elts, mul)
    rec = analyze_action(perms, name + ' regular')
    orders = element_orders(elts, mul)
    N = len(elts)
    rec['element_orders'] = dict(sorted(orders.items()))
    if p is not None:
        pred = p_group_prediction(orders, p, N)
        need(rec['depths_mod'][p] == pred['predicted_depth_mod_p'], name + ': Theorem B (p-group depth law)')
        rec['theorem_B'] = pred
    for ell in PRIMES:
        if N % ell != 0:
            need(rec['depths_mod'][ell] == coprime_prediction(orders, N, ell), name + ': Theorem B-prime at %d' % ell)
    rec['theorem_Bprime_checked_primes'] = [ell for ell in PRIMES if N % ell != 0]
    rec['candidate_lower_bound'] = {}
    for ell in PRIMES:
        if N % ell == 0:
            cand = mixed_candidates(orders, ell, N)
            need(rec['depths_mod'][ell] >= cand, name + ': depth mod %d is at least the candidate minimum' % ell)
            rec['candidate_lower_bound'][ell] = {'candidate_min': cand, 'equal': rec['depths_mod'][ell] == cand}
            need(rec['depths_mod'][ell] == cand, name + ': HYP-9374 equality holds in this finite universe at %d' % ell)
    return rec


def c21_cancellation_witness():
    """Mixed order: for C21 mod 3 the s^3 coefficients of the order-21 and order-7 terms cancel (760 + 98 = 858 = 0 mod 3)."""
    # exact contributions to [s^3] of (1/21) * e_o * (1 + (s-1)^o)^(21/o)
    from fractions import Fraction
    def term(o, e_o):
        poly = [0, 0, 0, 0]
        base = [0] * (o + 1)
        for i in range(0, o + 1):
            base[i] = comb(o, i) * (-1) ** (o - i)
        base[0] += 1  # 1 + (s-1)^o
        f = [1]
        for _ in range(21 // o):
            f = polymul(f, base)
        return Fraction(e_o * f[3], 21)
    contributions = {o: term(o, e) for o, e in [(1, 1), (3, 2), (7, 6), (21, 12)]}
    total = sum(contributions.values())
    need(total.denominator == 1 and total.numerator % 3 == 0, 'C21: s^3 coefficient vanishes mod 3 by cancellation')
    return {str(o): [c.numerator, c.denominator] for o, c in contributions.items()}


def main():
    out = {'tournaments': [tournament_record(z5sq_cayley), tournament_record(z25_circulant)]}
    rz = out['tournaments'][0]
    frozen = [1, 1, 12, 92, 506, 2130, 7084, 19228, 43263, 81719, 130760, 178296, 208012, 208012, 178296, 130760, 81719,
              43263, 19228, 7084, 2130, 506, 92, 12, 1, 1]
    need(rz['P'] == frozen, 'Z5xZ5 polynomial equals the frozen exact 2^25 census (census_orbits.c)')
    need((rz['depth_Q'], rz['depths_mod'][5], rz['depths_mod'][2], rz['depths_mod'][3], rz['depths_mod'][7]) == (5, 21, 25, 25, 5),
         'Z5xZ5: depths (Q,5,2,3,7) = (5,21,25,25,5)')
    need(rz['kappa'] == 3000, 'Z5xZ5: kappa = 3000')
    rc = out['tournaments'][1]
    need((rc['depth_Q'], rc['depths_mod'][5], rc['depths_mod'][2], rc['depths_mod'][3]) == (1, 5, 25, 1), 'Z25: depths (Q,5,2,3) = (1,5,25,1)')
    need(rc['kappa'] == 20, 'Z25: kappa = 20')
    out['groups'] = [group_record(name, p, b) for name, p, b in GROUPS]
    out['mixed_order_data'] = [group_record(name, None, b) for name, b in MIXED]
    out['c21_mod3_cancellation'] = c21_cancellation_witness()
    # even-action hostile: a transposition has subset polynomial 1+t+t^2 with no zero at -1
    Pt = burnside_subset_polynomial([(0, 1), (1, 0)])
    need(Pt == [1, 1, 1] and depth(shifted(Pt)) == 0, 'transposition hostile')
    out['hostile_transposition'] = Pt
    out['gates'] = 'all explicit exceptions; identical under -O'
    json.dump(out, sys.stdout, indent=1, sort_keys=True)
    print()


if __name__ == '__main__':
    main()
