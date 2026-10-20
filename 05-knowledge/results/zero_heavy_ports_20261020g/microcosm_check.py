"""Independent exact controls for the pointed-cube path product.

No producer imports, floating arithmetic, random sampling, or disabled asserts.
"""
from collections import Counter
from itertools import combinations, permutations, product
from math import gcd
from pathlib import Path
import hashlib
import json


def require(ok, reason):
    if not ok:
        raise RuntimeError(reason)


def tournament(a, b):
    n = a + b + 1
    rows = [sum(1 << v for v in range(u + 1, n)) for u in range(n)]
    for u, v in ((0, a), (a, a + b)):
        rows[u] ^= 1 << v
        rows[v] ^= 1 << u
    return tuple(rows)


def paths(rows):
    full = (1 << len(rows)) - 1

    def visit(word, unseen):
        if not unseen:
            yield word
            return
        choices = unseen if not word else unseen & rows[word[-1]]
        while choices:
            bit = choices & -choices
            choices ^= bit
            yield from visit(word + (bit.bit_length() - 1,), unseen ^ bit)

    return set(visit((), full))


def blue(rows, path):
    reversal = {v: path[-1 - i] for i, v in enumerate(path)}
    return all(bool(rows[u] & (1 << v)) !=
               bool(rows[reversal[u]] & (1 << reversal[v]))
               for u, v in combinations(range(len(rows)), 2))


def choice(k, x):
    top = 1 << (k - 1)
    return set(range(k)) if x == top else {
        i + 1 for i in range(k - 1) if x & (1 << i)}


def word(k, x, y):
    a = set(range(k))
    b = set(range(k + 1, 2 * k + 1))
    i, right = choice(k, x), choice(k, y)
    j = {2 * k - v for v in a - right}
    return tuple(sorted(i) + sorted(j) + [k] + sorted(a - i) + sorted(b - j))


def kappa(k, x):
    top = 1 << (k - 1)
    return x if x == top else top - 1 - x


def j_action(k, state):
    return tuple(kappa(k, x) for x in state)


def j_orbit(k, state):
    return min(state, j_action(k, state))


def homogeneous_chains(rows):
    n = len(rows)
    graph = [set() for _ in rows]
    for u, v in combinations(range(n), 2):
        if all(bool(rows[u] & (1 << w)) == bool(rows[v] & (1 << w))
               for w in range(n) if w not in (u, v)):
            graph[u].add(v)
            graph[v].add(u)
    unseen, answer = set(range(n)), []
    while unseen:
        todo, component = [min(unseen)], set()
        while todo:
            u = todo.pop()
            if u in component:
                continue
            component.add(u)
            todo.extend(graph[u] - component)
        unseen -= component
        answer.append(tuple(sorted(component)))
    return tuple(sorted(answer))


def weighted_feedback(rows, weights):
    return min(sum(weights[p[j]] * weights[p[i]]
                   for j in range(len(rows)) for i in range(j)
                   if rows[p[j]] & (1 << p[i]))
               for p in permutations(range(len(rows))))


def endpoint_formula(k):
    # Complete formula, including all empty-choice boundary terms.
    out = Counter()
    mass = [1] + [1 << (k - 1 - i) for i in range(1, k)]
    for i, j in product(range(k), repeat=2):
        out[i, 2 * k - j] += mass[i] * mass[j]
    for x in range(1, 1 << (k - 1)):
        subset = choice(k, x)
        alpha, beta = min(subset), max(set(range(k)) - subset)
        out[alpha, beta] += 1
        out[2 * k - beta, 2 * k - alpha] += 1
    out[0, k] += 1
    out[k, 2 * k] += 1
    out[k + 1, k - 1] += 1
    return out


def inspect(k):
    rows = tournament(k, k)
    top = 1 << (k - 1)
    q = top + 1
    states = list(product(range(q), repeat=2))
    by_state = {s: word(k, *s) for s in states}
    literal = paths(rows)
    require(len(set(by_state.values())) == q * q, "injective subset decoder")
    require(set(by_state.values()) == literal, "direct DFS versus decoder")
    require(all(tuple(2 * k - v for v in reversed(p)) == by_state[s[::-1]]
                for s, p in by_state.items()), "physical reflection gives swap")
    require(all(blue(rows, p) == (s[0] == s[1])
                for s, p in by_state.items()), "direct positional blue predicate")
    require(all(j_action(k, j_action(k, s)) == s and
                j_action(k, s[::-1]) == j_action(k, s)[::-1]
                for s in states), "commuting involutions")
    fixed_j = [s for s in states if j_action(k, s) == s]
    require(fixed_j == [(top, top)], "unique J fixed path")
    j_orbits = {j_orbit(k, s) for s in states}
    fixed_r_on_j = {s for s in j_orbits if j_orbit(k, s[::-1]) == s}
    joint = {min(j_orbit(k, s), j_orbit(k, s[::-1])) for s in states}
    m, n = (q + 1) // 2, (q - 1) // 2
    require(len(j_orbits) == m * m + n * n, "J orbit hypotenuse")
    require(len(j_orbits) - 1 == 2 * m * n, "nontrivial J orbit even leg")
    require(len(fixed_r_on_j) == q, "fixed R on J-orbits")
    require(len(joint) == m * m, "Klein-four joint orbits")
    require((len(j_orbits) - len(fixed_r_on_j)) // 2 == n * n,
            "paired R-orbits on J quotient")
    endpoints = Counter((p[0], p[-1]) for p in literal)
    require(endpoints == endpoint_formula(k), "complete endpoint response")
    diagonal = Counter((p[0], p[-1]) for p in literal if blue(rows, p))
    expected_diagonal = Counter({(0, 2 * k): 1, (k + 1, k - 1): 1})
    for i in range(1, k):
        expected_diagonal[i, 2 * k - i] = 1 << (k - 1 - i)
    require(diagonal == expected_diagonal, "blue endpoint response")
    chains = ((0,), tuple(range(1, k)), (k,),
              tuple(range(k + 1, 2 * k)), (2 * k,))
    require(homogeneous_chains(rows) == chains, "constant five-vertex chain core")
    reps = [c[0] for c in chains]
    core = tuple(sum(1 << j for j, v in enumerate(reps) if rows[u] & (1 << v))
                 for u in reps)
    require(core == tournament(2, 2), "literal weighted core identification")
    return dict(k=k, vertices=len(rows), H=len(literal), blue=q,
                J_orbits=len(j_orbits), fixed_R_on_J=len(fixed_r_on_j),
                Klein_orbits=len(joint), Euclid_parameters=[m, n],
                triple=[q, len(j_orbits)-1, len(j_orbits)],
                chains=[list(c) for c in chains],
                capped_weights=[1, min(k-1, 3), 1, min(k-1, 3), 1])


def growth(k):
    top = 1 << (k - 1)
    old_states = list(product(range(top + 1), repeat=2))
    fibres = Counter()
    colours = Counter()

    def project_choice(x):
        return top if x == 2 * top else x % top

    def project_path(p):
        kept = [v for v in p if v not in (k, k + 2)]
        return tuple(v if v < k else (k if v == k + 1 else v - 2) for v in kept)

    for state in product(range(2 * top + 1), repeat=2):
        parent = tuple(project_choice(x) for x in state)
        fibres[parent] += 1
        colours[parent[0] == parent[1], state[0] == state[1]] += 1
        require(project_path(word(k + 1, *state)) == word(k, *parent),
                "delete two actual vertices equals forget two bits")
        require(tuple(project_choice(x) for x in j_action(k + 1, state)) ==
                j_action(k, parent), "growth projection commutes with J")
    require(all(fibres[s] == 2 ** sum(x != top for x in s) for s in old_states),
            "four/two/one exact lift multiplicities")
    q = top + 1
    require(colours[False, True] == 0, "black parent has no blue lift")
    require(colours[True, False] == 2 * (q - 1), "blue-to-black birth count")
    require(sum(fibres.values()) == 4 * q * q - 4 * q + 1, "square growth")
    return dict(k=k, next_k=k+1, fibre_histogram=dict(sorted(Counter(fibres.values()).items())),
                blue_to_black=colours[True, False])


def main():
    levels = [inspect(k) for k in range(2, 8)]
    transitions = [growth(k) for k in range(2, 7)]
    core = tournament(2, 2)
    weighted = 0
    for left, right in product(range(1, 9), repeat=2):
        require(weighted_feedback(core, [1, left, 1, right, 1]) == 2,
                "weighted feedback remains two")
        weighted += 1
    rectangular = []
    for a, b in product(range(2, 6), repeat=2):
        rows = tournament(a, b)
        ps = paths(rows)
        count_blue = sum(blue(rows, p) for p in ps)
        require(len(ps) == (2**(a-1)+1)*(2**(b-1)+1), "rectangular product count")
        require(count_blue == (2**(a-1)+1 if a == b else 0), "blue requires balanced chain lengths")
        rectangular.append(dict(a=a, b=b, H=len(ps), blue=count_blue))
    balanced, unbalanced = tournament(5, 5), tournament(4, 6)
    require(tuple(r.bit_count() for r in balanced) == tuple(r.bit_count() for r in unbalanced),
            "same-order same-score hostile")
    require(sum(blue(unbalanced, p) for p in paths(unbalanced)) == 0,
            "same-order capped core does not preserve existence of blue paths")
    hostile = (16, 1, 3, 7, 14)
    hostile_paths = paths(hostile)
    hostile_blue = sorted(p for p in hostile_paths if blue(hostile, p))
    require(len(hostile_paths) == 9 and len(hostile_blue) == 1,
            "square H does not determine blue square root")
    same_ports = sorted(p for p in paths(tournament(3, 3)) if (p[0], p[-1]) == (1, 5))
    require(len(same_ports) == 4 and sum(blue(tournament(3, 3), p) for p in same_ports) == 2,
            "endpoint-only blue hostile")
    matrix = ((2, 1, -1), (-2, 2, 2), (-2, 1, 3))
    signature = (1, 1, -1)
    require(all(sum(matrix[t][i] * signature[t] * matrix[t][j] for t in range(3)) ==
                (4 * signature[i] if i == j else 0)
                for i, j in product(range(3), repeat=2)), "Lorentz similitude")
    for k in range(2, 101):
        n, m = 1 << (k - 2), (1 << (k - 2)) + 1
        q = 2 * n + 1
        t = (q, 2 * m * n, m * m + n * n)
        require(t[0]**2 + t[1]**2 == t[2]**2 and gcd(t[0], t[1]) == 1,
                "primitive Euclid triple")
        u, v, z = n, n + 1, n * (n + 1)
        require(v*z == u*(v+z) and (u+v, 2*z, v+2*z-u) == t,
                "inherited harmonic K4 face decoder")
        expected = (2*q-1, 4*n*(2*n+1), (2*n+1)**2 + (2*n)**2)
        require(tuple(sum(row[i]*t[i] for i in range(3)) for row in matrix) == expected,
                "doubling spinor gives exact triple transport")
        require((abs(n+1-n*n) == 1) == (k in (2, 3)), "only first two unit-Cassini windows")
    report = dict(status="FINITE-EXACT controls; all-order proofs in microcosm.md",
                  source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  levels=levels, growth=transitions, weighted_controls=weighted,
                  rectangular_controls=rectangular,
                  same_cap_same_order_hostile=dict(parameters=[[5,5],[4,6]],
                      vertices=11, capped_weights=[1,3,1,3,1], blue=[17,0], H=[289,297]),
                  same_H_hostile=dict(rows=hostile, H=9, blue=hostile_blue),
                  same_endpoints_hostile=dict(k=3, endpoints=[1, 5], paths=same_ports),
                  arithmetic_control="all k=2..100; symbolic matrix entries checked exactly")
    Path(__file__).with_suffix('.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    for row in levels:
        print('k={k}: H={H}, blue={blue}, J-orbits={J_orbits}, fixed R/J={fixed_R_on_J}, Klein-orbits={Klein_orbits}, triple={triple}'.format(**row))
    print('All five actual two-vertex deletions have exact 4/2/1 path fibres.')
    print('All 64 weighted five-core controls have feedback minimum 2.')
    print('T_(5,5) and T_(4,6): same order, capped core and scores; blue counts 17 and 0.')
    print('Hostiles: H=9 but blue=1; endpoints (1,5) contain two blue and two black paths.')
    print('PASS exact involutions, endpoint response, chain cores, growth and arithmetic.')


if __name__ == '__main__':
    main()
