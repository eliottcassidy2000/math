#!/usr/bin/env python3
"""Standalone finite controls for the fixed-graft-menu boundary theorem.

No repository imports, third-party packages, randomness, floating arithmetic,
or optimization-disabled assertions. All-order conclusions belong to PROOF.md.
"""
from collections import Counter
from fractions import Fraction
from itertools import combinations, permutations
import json
from math import gcd


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


RAYS = ((1, 1, 1, 1), (1, 0, 2, 1), (1, 0, 1, 0), (0, 0, 1, 0))
O = ((1, 0, 0, 0), (0, 0, 0, 0), (0, 0, 2, 0), (0, 0, 0, 1))
Z = ((1, -1, 1, 0), (0, 0, 0, 1), (0, 0, 2, -1), (0, 0, 0, 1))
J = ((1, -1, 1, 0), (1, -1, 0, 0), (1, 0, 1, 0), (1, 0, 0, 0))


def mv(matrix, x):
    return tuple(sum(a * b for a, b in zip(row, x)) for row in matrix)


def facets(x):
    D, E, U, d = x
    return (E, d - E, D - d, U - D - d + E)


def physical_controls(x):
    D, E, U, d = x
    return (D >= 1 and d >= 1 and min(facets(x)) >= 0
            and D % 2 == d % 2 == 1 and U % 2 == E % 2)


def tournaments(n):
    pairs = tuple(combinations(range(n), 2))
    for mask in range(1 << len(pairs)):
        rows = [0] * n
        for bit, (a, b) in enumerate(pairs):
            u, v = (a, b) if mask >> bit & 1 else (b, a)
            rows[u] |= 1 << v
        yield tuple(rows)


def transitive(n):
    return tuple(sum(1 << j for j in range(i + 1, n)) for i in range(n))


def cyclic_block(m):
    u, v = m, m + 1
    return tuple(sum(1 << j for j in range(i + 1, m)) | (1 << u)
                 for i in range(m)) + (1 << v, (1 << m) - 1)


def interval_paths(rows, p):
    n = len(p)
    good = [[True] * (n + 1) for _ in range(n + 1)]
    for i in range(n):
        for j in range(i + 2, n + 1):
            good[i][j] = good[i][j - 1] and bool(rows[p[j - 2]] >> p[j - 1] & 1)
    return good


def literal_unrooted(rows, vertices):
    """Count ordered words and all split positions, including empty bins."""
    H = Q = 0
    n = len(vertices)
    for p in permutations(vertices):
        good = interval_paths(rows, p)
        H += good[0][n]
        Q += sum(good[0][i] and good[i][n] for i in range(n + 1))
    return H, Q


def literal_statistics(rows):
    """All marked roots at once, by permutations and 2/3-bin cut positions."""
    n = len(rows)
    require(n >= 1, "gadgets must be nonempty")
    H = Q = 0
    s, e, alpha, delta, T = ([0] * n for _ in range(5))
    for p in permutations(range(n)):
        good = interval_paths(rows, p)
        if good[0][n]:
            H += 1
            s[p[0]] += 1
            e[p[-1]] += 1
        Q += sum(good[0][i] and good[i][n] for i in range(n + 1))
        for i in range(1, n + 1):
            if not good[0][i]:
                continue
            if good[i][n]:
                alpha[p[0]] += 1
                delta[p[i - 1]] += 1
            T[p[0]] += sum(good[i][j] and good[j][n] for j in range(i, n + 1))
    out = []
    for root in range(n):
        q, b = literal_unrooted(rows, tuple(i for i in range(n) if i != root))
        out.append(dict(H=H, Q=Q, s=s[root], e=e[root], alpha=alpha[root],
                        delta=delta[root], T=T[root], q=q, b=b))
    return out


def operator(stats, order):
    H, Q, s, e, alpha, delta, T, q, b = (
        stats[k] for k in ('H', 'Q', 's', 'e', 'alpha', 'delta', 'T', 'q', 'b'))
    return ((H, -H, Q - H, 2 * H - Q),
            (e, -e, delta - e, e - delta),
            (alpha, s - alpha, 2 * T - alpha, 3 * alpha - s - 2 * T),
            (q, -q, b - q, 2 * q - b) if order >= 2 else (1, 0, 0, 0))


def ray_ratios(stats, matrix):
    H, Q, alpha, T, s = (stats[k] for k in ('H', 'Q', 'alpha', 'T', 's'))
    numerators = (2 * alpha, 2 * alpha + 2 * T - s, 2 * T, 2 * T - alpha)
    denominators = (H, H + Q, Q, Q - H)
    require(H >= 1 and Q >= 2 * H, "positive denominators")
    out = []
    for r, num, den in zip(RAYS, numerators, denominators):
        y = mv(matrix, r)
        require(y[0] == den and y[2] == num, "literal-statistics ray formula")
        require(den > 0 and min(facets(y)) >= 0, "graft positive ray image")
        require(y[0] >= r[0], "graft D monotonicity on cone")
        out.append(Fraction(num, den))
    return tuple(out)


def dp_profile(rows, root):
    """Independent subset recurrences for start/end counts and cover split."""
    n = len(rows)
    size = 1 << n
    ends = [[0] * n for _ in range(size)]
    starts = [[0] * n for _ in range(size)]
    H = [0] * size
    H[0] = 1
    for mask in range(1, size):
        for v in range(n):
            if not (mask >> v & 1):
                continue
            rem = mask ^ (1 << v)
            if rem == 0:
                ends[mask][v] = starts[mask][v] = 1
            else:
                ends[mask][v] = sum(ends[rem][u] for u in range(n)
                                    if rem >> u & 1 and rows[u] >> v & 1)
                starts[mask][v] = sum(starts[rem][u] for u in range(n)
                                      if rem >> u & 1 and rows[v] >> u & 1)
        H[mask] = sum(ends[mask])
    full = size - 1
    U = sum(starts[mask][root] * H[full ^ mask] for mask in range(size)
            if mask >> root & 1)
    return (H[full], ends[full][root], U, H[full ^ (1 << root)])


def graft(a, ar, b, br):
    n, m = len(a), len(b)
    out = list(a) + [row << n for row in b]
    for i in range(n):
        for j in range(m):
            if i == ar:
                out[n + j] |= 1 << i
            else:
                out[i] |= 1 << (n + j)
    return tuple(out), n + br


def literal_O(a, root):
    n = len(a)
    return tuple(row | (1 << n) for row in a) + (0,), root


def rational_record(values):
    return [str(v) for v in values]


def main():
    counts = Counter()
    small_bank = []
    unique_states = set()
    maximum_by_order = {}
    for n in range(1, 6):
        max_kappa = Fraction(0)
        for rows in tournaments(n):
            all_stats = literal_statistics(rows)
            for root, stats in enumerate(all_stats):
                state = (stats['H'], stats['e'], stats['alpha'], stats['q'])
                require(dp_profile(rows, root) == state, "literal / subset profile agreement")
                require(physical_controls(state), "physical cone and parity")
                coeff = facets(state)
                require(tuple(sum(coeff[i] * RAYS[i][j] for i in range(4))
                              for j in range(4)) == state, "ray decomposition")
                matrix = operator(stats, n)
                ratios = ray_ratios(stats, matrix)
                max_kappa = max(max_kappa, *ratios)
                if n == 1:
                    require(matrix == J, "singleton full-matrix exception")
                if n <= 3:
                    small_bank.append((rows, root, state, matrix))
                unique_states.add(state)
                counts['literal_rooted_gadgets'] += 1
        maximum_by_order[str(n)] = str(max_kappa)
    require(counts['literal_rooted_gadgets'] == 5405, "unfiltered gadget universe")

    # Literal cross-edge construction, independently counted without a graft formula.
    for a, ar, x, _ in small_bank:
        for b, br, _, matrix in small_bank:
            rows, root = graft(a, ar, b, br)
            y = dp_profile(rows, root)
            require(y == mv(matrix, x), "literal graph / graft matrix agreement")
            require(y[0] >= x[0], "literal induced-source monotonicity")
            counts['literal_rooted_graft_pairs'] += 1
        rows, root = literal_O(a, ar)
        require(dp_profile(rows, root) == mv(O, x), "literal O")
        rows, _ = graft(a, ar, (0,), 0)
        require(dp_profile(rows, ar) == mv(Z, x), "literal Z")
        counts['literal_O_Z_pairs'] += 1

    gadgets = {
        'J': ((0,), 0),
        'L': ((2, 0), 0),
        'L_sink': ((2, 0), 1),
        'C3': ((2, 4, 1), 0),
        'B4': ((8, 1, 3, 6), 0),
        'T5_source': (transitive(5), 0),
    }
    menu = {'O': O, 'Z': Z}
    ratios_by_name = {}
    kappa_by_name = {'Z': Fraction(2)}
    gadget_records = {}
    for name, (rows, root) in gadgets.items():
        stats = literal_statistics(rows)[root]
        matrix = operator(stats, len(rows))
        menu[name] = matrix
        ratios = ray_ratios(stats, matrix)
        ratios_by_name[name] = rational_record(ratios)
        kappa_by_name[name] = max(ratios)
        gadget_records[name] = dict(rows=list(rows), root=root, statistics=stats)
    require(kappa_by_name['J'] == kappa_by_name['C3'] == kappa_by_name['B4'] == 2,
            "named kappa two controls")
    require(kappa_by_name['L'] == 4, "source edge sharp kappa")
    require(kappa_by_name['T5_source'] == 32, "no menu-independent kappa")
    require([mv(Z, r)[0] for r in RAYS] == [1, 3, 2, 1], "Z D rays")
    require([mv(Z, r)[2] for r in RAYS] == [1, 3, 2, 2], "Z U rays")
    require(mv(O, RAYS[3])[0] == 0 and mv(O, RAYS[3])[2] > 0, "O singular hostile")
    kappa = max(kappa_by_name.values())

    # Exhaustive finite words, with no positivity or size pruning.
    seed_states = {x for _, _, x, _ in small_bank}
    seed_states.update(dp_profile(gadgets['B4'][0], r) for r in range(4))
    big_U_seed = dp_profile(transitive(10), 0)
    require(big_U_seed == (1, 0, 512, 1), "large-U O-only seed")
    seed_states.add(big_U_seed)
    sink_rows = tuple(row | 8 for row in (2, 4, 1)) + (0,)
    sink_state = dp_profile(sink_rows, 3)
    require(sink_state == (3, 3, 3, 3), "nonunit sink seed")
    seed_states.add(sink_state)
    seeds = sorted(seed_states)
    depth = 5
    caps = (1, 5, 15, 63)
    X_cap = 100000
    observed = {C: set() for C in caps}
    for seed in seeds:
        # state, length, core pair, tail length, has a non-O operation
        stack = [(seed, 0, (seed[0], seed[2]), 0, False)]
        while stack:
            x, length, core, tail, has_nonO = stack.pop()
            D, E, U, d = x
            require(physical_controls(x), "finite-word physical constraints")
            require(D == core[0] and U == core[1] << tail, "last-nonO exact factorization")
            h, g = D + 2 * U, 2 * U
            endpoint = mv(J, mv(O, x))
            require(endpoint == (h, D, h, D), "final OJ state")
            closed = mv(J, endpoint)
            N = D + 4 * U
            require(closed[0] == N == h + g, "extra closing J scalar")
            require(h > g > 0 and 3 * g >= 2 * h and g % 2 == 0, "endpoint chamber")
            for C in caps:
                if D <= C:
                    if has_nonO:
                        require(1 <= core[1] <= (kappa * C).__floor__(), "finite core bound")
                    else:
                        require(core == (seed[0], seed[2]), "O-only core retained")
                    if N <= X_cap:
                        observed[C].add(N)
                    counts['bounded_intercept_word_controls'] += 1
            counts['finite_words'] += 1
            if length == depth:
                continue
            for name, matrix in menu.items():
                y = mv(matrix, x)
                require(y[0] >= D, "finite-word D monotonicity")
                if name == 'O':
                    stack.append((y, length + 1, core, tail + 1, has_nonO))
                else:
                    require(Fraction(y[2], y[0]) <= kappa_by_name[name], "individual operation ratio")
                    stack.append((y, length + 1, (y[0], y[2]), 0, True))
    require(counts['finite_words'] == len(seeds) * sum(len(menu) ** j for j in range(depth + 1)),
            "exhaustive word universe")
    boundary_records = {}
    for C in caps:
        pairs = {(D, u) for D in range(1, C + 1)
                 for u in range(1, (kappa * C).__floor__() + 1)}
        pairs.update((s[0], s[2]) for s in seeds if s[0] <= C)
        values = set()
        for D, u in pairs:
            value = D + 4 * u
            while value <= X_cap:
                values.add(value)
                u *= 2
                value = D + 4 * u
        B_C = C * (kappa * C).__floor__() + sum(s[0] <= C for s in seeds)
        bound = B_C * X_cap.bit_length()
        require(len(pairs) <= B_C and len(values) <= bound, "explicit envelope count")
        require(observed[C] <= values, "all finite-word counts in dyadic envelope")
        boundary_records[str(C)] = dict(core_bound=B_C, distinct_envelope_pairs=len(pairs),
                                        distinct_envelope_counts=len(values),
                                        distinct_word_counts=len(observed[C]),
                                        explicit_log_count_bound=bound)

    # Unbounded-history mechanisms are checked at large finite values separately.
    for x in ((1, 1, 1, 1), sink_state):
        y = x
        for _ in range(1000):
            y = mv(Z, y)
            require(y == x, "stationary sink-root Z run")
            counts['stationary_Z_steps'] += 1
    for seed in seeds:
        x = seed
        for k in range(257):
            require(x[0] == seed[0] and x[2] == seed[2] << k, "large O-tail exactness")
            require(mv(J, mv(J, mv(O, x)))[0] == seed[0] + (4 * seed[2] << k),
                    "large O-tail closing count")
            x = mv(O, x)
            counts['large_O_tail_controls'] += 1
    require(big_U_seed[2] > kappa * 1 and big_U_seed[0] == 1,
            "seed-tail exception cannot be removed")
    for n in range(1, 10):
        x = dp_profile(transitive(n), 0)
        expected = (1, 1, 1, 1) if n == 1 else (1, 0, 1 << (n - 1), 1)
        require(x == expected, "physical r4 sequence")
        require(mv(O, x)[2] == 1 << n, "O has no fixed ratio bound")
        y = dp_profile(cyclic_block(n), n)
        require(y == ((1 << n) + 1, 1, (1 << n) + 1, 1), "physical r3 sequence")
        counts['physical_ray_sequence_pairs'] += 1

    # Exact even-g congruence, discarded-line threshold, and missing-gap hostile.
    for b in range(2, 13):
        for a in range(1, b):
            if gcd(a, b) != 1:
                continue
            for c in range(-20, 21):
                for N in range(-30, 101):
                    numerator = a * N - c
                    arithmetic = numerator % (a + b) == 0 and (numerator // (a + b)) % 2 == 0
                    require(arithmetic == (numerator % (2 * (a + b)) == 0), "exact even-g congruence")
                    if abs(c) < a:
                        require((a - c) % (2 * (a + b)) != 0, "r=1 small-intercept gap")
                    if arithmetic:
                        g = numerator // (a + b)
                        h = N - g
                        if h > g > 0 and 3 * g >= 2 * h and 3 * a < 2 * b:
                            H_c = (-3 * c) // (2 * b - 3 * a)
                            require(h <= H_c and N <= 2 * H_c - 1, "discarded-line finite threshold")
                            counts['discarded_line_candidates'] += 1
                    counts['congruence_candidates'] += 1
    hostile = (0, -2, -4, -6, -8)
    for N in range(17, 1001, 2):
        matches = [c for c in hostile if (2 * N - c) % 10 == 0]
        require(len(matches) == 1, "no-common-residue-gap hostile")
        c = matches[0]
        h, g = (3 * N + c) // 5, (2 * N - c) // 5
        require(h > g > 0 and g % 2 == 0 and 3 * g >= 2 * h and 2 * h - 3 * g == c,
                "hostile positive even-g chamber candidates")
        counts['positive_hostile_chambers'] += 1

    report = dict(
        status='PASS', all_order_claims_source='PROOF.md', finite_counts=dict(sorted(counts.items())),
        labelled_tournament_orders=[1, 2, 3, 4, 5], distinct_small_physical_states=len(unique_states),
        literal_graft_source_and_gadget_orders=[1, 2, 3],
        max_kappa_in_literal_universe_by_order=maximum_by_order,
        mixed_menu_gadgets=gadget_records, ray_ratios=ratios_by_name,
        operation_kappa={k: str(v) for k, v in sorted(kappa_by_name.items())},
        menu_kappa=str(kappa), menu_word_alphabet=list(menu), maximum_word_length=depth,
        physical_seed_states=seeds, boundary_checks=boundary_records, scalar_ceiling=X_cap,
        no_common_residue_gap_hostile=list(hostile), hostile_uniform_positive_threshold=17,
        scope='Finite controls only; fixed-menu all-order proof is in PROOF.md.')
    print(json.dumps(report, sort_keys=True, indent=2))


if __name__ == '__main__':
    main()
