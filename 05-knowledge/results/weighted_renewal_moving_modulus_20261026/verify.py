#!/usr/bin/env python3
"""Exact, deterministic hostile controls; no removable assertions or dependencies."""
from fractions import Fraction as F
from itertools import combinations, permutations, product
import json
import random


def require(ok, message):
    if not ok:
        raise RuntimeError(message)


def forest(n, edges):
    parents = list(range(n))
    def root(x):
        while parents[x] != x:
            x = parents[x]
        return x
    for x, y in edges:
        a, b = root(x), root(y)
        if a == b:
            return False
        parents[a] = b
    return True


def edge(rank, defects, x, y):
    return (rank[x] < rank[y]) != (tuple(sorted((x, y))) in defects)


def literal(rank, defects):
    n = len(rank)
    return tuple(tuple(False if x == y else edge(rank, defects, x, y)
                       for y in range(n)) for x in range(n))


def hp(table, path):
    return all(table[x][y] for x, y in zip(path, path[1:]))


def insertion_path(table):
    path = []
    for x in range(len(table)):
        j = next((j for j, y in enumerate(path) if table[x][y]), len(path))
        path.insert(j, x)
    require(hp(table, path), "Redei insertion")
    return tuple(path)


def cut_insert_guard(table, path, i, j):
    # Remove position i, insert into slot j of the shortened list.
    v = path[i]
    q = path[:i] + path[i + 1:]
    if i == j:
        return True
    guards = []
    if 0 < i < len(path) - 1:
        guards.append((path[i - 1], path[i + 1]))
    if j:
        guards.append((q[j - 1], v))
    if j < len(q):
        guards.append((v, q[j]))
    # If the removed gap is also the insertion slot, we returned above.
    return all(table[a][b] for a, b in guards)


def triangle_or_cycle(edges, n):
    return not forest(n, edges)


def tournament_checks():
    counts = {}
    paths = moves = 0
    flip_hostile = None
    examples = []
    for n in range(1, 7):
        pairs = list(combinations(range(n), 2))
        count = 0
        rank = tuple(range(n))
        for mask in range(1 << len(pairs)):
            ds = frozenset(e for i, e in enumerate(pairs) if mask >> i & 1)
            if not forest(n, ds):
                continue
            count += 1
            table = literal(rank, ds)
            if n <= 3:
                examples.append((rank, ds, insertion_path(table)))
            for p in permutations(range(n)):
                if not hp(table, p):
                    continue
                paths += 1
                # Literal re-indexing versus retained rank/defect frame.
                new_rank = tuple(rank[x] for x in p)
                inverse = {x: i for i, x in enumerate(p)}
                new_ds = frozenset(tuple(sorted((inverse[x], inverse[y]))) for x, y in ds)
                require(literal(new_rank, new_ds) == tuple(tuple(table[x][y] for y in p) for x in p),
                        "relocation lost an arc")
                weights = tuple(F(x+1, n+1) for x in range(n))
                transported = tuple(weights[x] for x in p)
                require(all(transported[inverse[x]] == weights[x] for x in range(n)),
                        "named source weight changed under reindexing")
                require(p[inverse[0]] == 0 and p[inverse[n-1]] == n-1,
                        "named ports changed under reindexing")
                path_edges = frozenset(tuple(sorted(e)) for e in zip(p, p[1:]))
                flipped_ds = ds ^ path_edges
                flipped_rank = tuple(n - 1 - r for r in rank)
                ft = literal(flipped_rank, flipped_ds)
                require(all(ft[x][y] == (table[x][y] if tuple(sorted((x, y))) in path_edges
                                        else not table[x][y])
                            for x, y in pairs), "tau XOR update")
                if flip_hostile is None and triangle_or_cycle(flipped_ds, n):
                    flip_hostile = dict(n=n, rank=rank, defects=sorted(ds), path=p,
                                        new_defects=sorted(flipped_ds))
                for i, j in product(range(n), repeat=2):
                    q = p[:i] + p[i + 1:]
                    q = q[:j] + (p[i],) + q[j:]
                    require(cut_insert_guard(table, p, i, j) == hp(table, q),
                            "three guards versus full Hamilton check")
                    moves += 1
        counts[n] = count
    require(flip_hostile is not None, "missing forest flip hostile")
    # Both binary gluings, disjoint labels and exactly the same named hinge.
    joins = 0
    for (ra, da, pa), (rb, db, pb) in product(examples, repeat=2):
        a, b = len(ra), len(rb)
        pb = tuple(x + a for x in pb)
        db = frozenset((x + a, y + a) for x, y in db)
        hinge = tuple(sorted((pa[-1], pb[0])))
        for skew in (False, True):
            ranks = (tuple(r + b for r in ra) + rb) if skew else (ra + tuple(r + a for r in rb))
            ds = da | db | (frozenset((hinge,)) if skew else frozenset())
            require(forest(a + b, ds), "binary composition introduced a cycle")
            table = literal(ranks, ds)
            require(hp(table, pa + pb), "binary composition lost active path")
            require(all(table[x][y] == (not skew or (x, y) == hinge)
                        for x in range(a) for y in range(a, a + b)), "cross arcs")
            joins += 1
    return dict(forests_by_order=counts, legal_paths=paths, guarded_moves=moves,
                binary_joins=joins, first_frozen_forest_tau_hostile=flip_hostile)


def positive_definite(matrix):
    a = [row[:] for row in matrix]
    n = len(a)
    for k in range(n):
        if a[k][k] <= 0:
            return False
        for i in range(k + 1, n):
            for j in range(i, n):
                a[i][j] -= a[i][k] * a[k][j] / a[k][k]
                a[j][i] = a[i][j]
    return True


def hardy_check(u, v, gains):
    n = len(u)
    p = [F(1)]
    for a in gains:
        p.append(p[-1] * a)
    up = [u[i] * p[i] ** 2 for i in range(n)]
    vp = [v[i] * p[i] ** 2 for i in range(n)]
    cuts = [sum(up[k:]) * sum(F(1) / x for x in vp[:k + 1]) for k in range(n)]
    B = max(cuts)
    # Gram of the actual gain kernel, independently built from interval products.
    h = [[F(0) for _ in range(n)] for _ in range(n)]
    for i in range(n):
        h[i][i] = F(1)
        for j in range(i + 1, n):
            h[j][i] = h[j - 1][i] * gains[j - 1]
            require(h[j][i] == p[j] / p[i], "gain conjugation")
    gram = [[sum(u[j] * h[j][i] * h[j][k] for j in range(n))
             for k in range(n)] for i in range(n)]
    majorant = [[(4 * B * v[i] if i == k else 0) - gram[i][k]
                 for k in range(n)] for i in range(n)]
    require(positive_definite(majorant), "4B upper energy bound")
    for k in range(n):
        f = [F(1) / (v[i] * p[i]) if i <= k else F(0) for i in range(n)]
        denom = sum(v[i] * f[i] ** 2 for i in range(n))
        num = sum(u[j] * sum(h[j][i] * f[i] for i in range(n)) ** 2 for j in range(n))
        require(num >= cuts[k] * denom, "cut lower-bound witness")
    return B


def renewal_checks():
    count = 0
    for n in range(1, 5):
        for weights in product((F(1, 2), F(1), F(2)), repeat=2 * n):
            hardy_check(weights[:n], weights[n:], [F(1)] * (n - 1))
            count += 1
    rng = random.Random(48074808)
    for n in range(1, 13):
        for _ in range(20):
            u = [F(rng.randint(1, 7), rng.randint(1, 5)) for _ in range(n)]
            v = [F(rng.randint(1, 7), rng.randint(1, 5)) for _ in range(n)]
            g = [rng.choice((F(-2), F(-1, 2), F(1, 2), F(3, 2))) for _ in range(n - 1)]
            hardy_check(u, v, g)
            count += 1
    # Refinement uses a conditional mean, not a new physical transition.
    owners = (0, 0, 1, 1, 1, 2)
    vf = tuple(map(F, (1, 2, 2, 1, 3, 4)))
    uf = tuple(map(F, (3, 1, 1, 4, 2, 2)))
    vc = [sum(vf[a] for a, i in enumerate(owners) if i == j) for j in range(3)]
    uc = [sum(uf[a] for a, i in enumerate(owners) if i == j) for j in range(3)]
    def lift(x):
        return [x[i] for i in owners]
    def expectation(x):
        return [sum(vf[a] * x[a] for a, i in enumerate(owners) if i == j) / vc[j]
                for j in range(3)]
    for x in product((-1, 0, 1), repeat=3):
        require(expectation(lift(x)) == list(x), "refinement left inverse")
        require(sum(vf[a] * lift(x)[a] ** 2 for a in range(6)) ==
                sum(vc[i] * x[i] ** 2 for i in range(3)), "source mass isometry")
    for f in product((-1, 0, 1), repeat=6):
        e = expectation(f)
        coarse = [sum(e[:i + 1]) for i in range(3)]
        fine = lift(coarse)
        require(sum(uf[a] * fine[a] ** 2 for a in range(6)) ==
                sum(uc[i] * coarse[i] ** 2 for i in range(3)), "target mass isometry")
    # All a_i=1 are literal finite Collatz words; unit-weight gains amplify.
    collatz = []
    for L in range(1, 33):
        n = (1 << (L + 1)) - 1
        x, clock = n, 0
        for _ in range(L):
            z = 3 * x + 1
            require(z % 4 == 2, "all-one valuation cylinder")
            x = z // 2
            clock += 2
        require(x == 3 ** L * 2 - 1, "all-one endpoint")
        collatz.append(dict(L=L, source=n, endpoint=x, integer_clock=clock))
    return dict(weighted_gram_cases=count, refinement_vectors=27 + 729,
                artificial_refinement_rayleigh={r: str(F((r + 1) * (2 * r + 1), 6))
                                                for r in (2, 4, 8, 16, 32)},
                collatz_all_one_controls=collatz)


def apply_word(word):
    m, n, alpha, gamma = 2, 1, 1, 0
    chi = 1
    for letter in word:
        if letter == 'U':
            m, n = 2 * m - n, m
            alpha, gamma = 2 * alpha - gamma, alpha
        elif letter == 'A':
            m, n = 2 * m + n, m
            alpha, gamma = 2 * alpha + gamma, alpha
            chi = -chi
        else:
            m += 2 * n
            alpha += 2 * gamma
        require(alpha * n - gamma * m == chi, "pole determinant")
        require((alpha * n - gamma * m) % 3 == chi % 3, "fixed-modulus frame")
    return m, n, chi


def modulus_checks():
    for M in range(1, 257):
        m, n, chi = apply_word('A' + 'D' * (M - 1) + 'UD')
        require((m, n, chi) == (16 * M + 2, 4 * M + 1, -1), "every modulus hostile")
        require((m % M, n % M) == (2 % M, 1 % M), "root residue collision")
        require(m*m + n*n == 272 * M*M + 72 * M + 5, "height threshold")
    count = 0
    for depth in range(9):
        for word in product('UAD', repeat=depth):
            m, n, _ = apply_word(word)
            M = 1 << m.bit_length()
            require(M > m > n > 0 and (m % M, n % M) == (m, n), "dyadic exact lift")
            count += 1
    return dict(root_alias_moduli=256, all_words_through_depth=8,
                height_lift_and_pole_controls=count)


def passive_checks():
    branches = 0
    for n in range(1, 4096, 2):
        z = 3 * n + 1
        k = (z & -z).bit_length() - 1
        target = z >> k
        b = F(3 * n, z)
        d = 1 - b*b
        require(F(3, 1 << k) * n / target == b, "relative Collatz multiplier")
        for e, f in product((F(-2), F(-1, 3), F(0), F(3, 2)), repeat=2):
            out = b * e + f
            port = d * e - b * f
            require(out*out + port*port/d == e*e + f*f/d, "passive port identity")
            branches += 1
    trajectories = 0
    for n0 in range(1, 512, 2):
        n, e, initial = n0, F(2, 3), F(4, 9)
        loss = supplied = F(0)
        clock = 0
        for j in range(32):
            z = 3*n + 1
            k = (z & -z).bit_length() - 1
            b, f = F(3*n, z), F((-1)**j, j + 2)
            d = 1-b*b
            port = d*e - b*f
            loss += port*port/d
            supplied += f*f/d
            e = b*e + f
            n = z >> k
            clock += k+1
        require(e*e + loss == initial + supplied, "source-specific telescoping")
        require(clock >= 64, "integer clock is not odd-step count")
        trajectories += 1
    return dict(single_branch_rational_identities=branches,
                actual_32_odd_step_prefixes=trajectories)


def adaptive_source_checks():
    def mass(r, m):
        return F(1, (m+1)**r) - F(1, (m+2)**r)
    checked = 0
    schedules = {}
    for r in (1, 3, 6):
        mu = [mass(r, m) for m in range(128)]
        f = [F((-1)**(m % 3) * (m+1), m % 7 + 1) for m in range(128)]
        energy = sum(w*x*x for w, x in zip(mu, f))
        for m in range(16):
            for q in (16, 32, 64, 128):
                alias = sum(mu[s]*f[s] for s in range(m+q, 128, q))
                bound = F(r+1, q**(r+1))
                require(alias*alias <= bound*energy, "marked-source alias bound")
                checked += 1
        q, m = 1, 7
        schedule = []
        for t in range(17):
            E, epsilon = 9**t, F(1, 2**t)
            while q <= m or F(r+1, q**(r+1))*E > epsilon**2*mass(r, m)**2:
                q *= 2
            require(q & (q-1) == 0, "dyadic modulus")
            require(F(r+1, q**(r+1))*E <= epsilon**2*mass(r, m)**2,
                    "adaptive energy-to-precision budget")
            schedule.append(q)
        schedules[r] = schedule
    return dict(finite_alias_controls=checked, energy_9_to_t_precision_schedules=schedules)


def paley_escape_checks():
    n = 7
    t = tuple(tuple((y-x) % 7 in (1, 2, 4) for y in range(n)) for x in range(n))
    path = tuple(range(n))
    pe = {frozenset((i, i+1)) for i in range(n-1)}
    s = tuple(tuple(False if x == y else t[x][y] if frozenset((x, y)) in pe
                    else not t[x][y] for y in range(n)) for x in range(n))
    order = (0, 3, 1, 6, 4, 2, 5)
    rank = tuple(order.index(x) for x in range(n))
    defects = frozenset((x, y) for x, y in combinations(range(n), 2)
                        if s[x][y] != (rank[x] < rank[y]))
    require(forest(n, defects) and len(defects) == 6, "Paley predecessor forest")
    require(literal(rank, defects) == s and hp(s, path), "Paley predecessor ports")
    triangles = [q for q in combinations(range(n), 3)
                 if all(sum(t[x][y] for y in q if y != x) == 1 for x in q)]
    require(len(triangles) == 14, "Paley fourteen cycles")
    require(all(sum(x in q and y in q for q in triangles) == 2
                for x, y in combinations(range(n), 2)), "two triangles per edge")
    minimum = 21
    for p in permutations(range(n)):
        ds = [(x, y) for i, x in enumerate(p) for y in p[i+1:] if t[y][x]]
        minimum = min(minimum, len(ds))
        require(not forest(n, ds), "Paley cannot have any forest frame")
    require(minimum == 7, "Paley minimum feedback size")
    natural = [(x, y) for x, y in combinations(range(n), 2) if t[y][x]]
    # The natural frame has triangle 0-3-6-0 and trees attached to it.
    require(set(natural) == {(0, 3), (1, 4), (2, 5), (3, 6), (0, 5), (1, 6), (0, 6)},
            "explicit unicyclic repair")
    require(forest(n, [e for e in natural if e != (0, 6)]), "one cycle suffices")
    for copies in range(1, 13):
        N = 7*copies
        ri = tuple(7*(x//7) + rank[x % 7] for x in range(N))
        di = frozenset((7*b+x, 7*b+y) for b in range(copies) for x, y in defects)
        require(forest(N, di), "disjoint joined predecessors remain forests")
        before = literal(ri, di)
        require(hp(before, tuple(range(N))), "joined active path")
        ro = tuple(7*(copies-1-x//7) + x % 7 for x in range(N))
        do = frozenset((7*b+x, 7*b+y) for b in range(copies) for x, y in natural)
        do |= frozenset((7*b+6, 7*b+7) for b in range(copies-1))
        after = literal(ro, do)
        require(all(after[x][y] == (before[x][y] if y == x+1 else not before[x][y])
                    for x, y in combinations(range(N), 2)), "joined one-flip cycle debt")
        # Removing the one known cycle edge per block leaves a connected tree.
        removed = {(7*b, 7*b+6) for b in range(copies)}
        require(len(do) == N-1+copies and forest(N, do-removed), "sharp r-cycle repair")
    return dict(predecessor_rank_order=order, predecessor_defects=sorted(defects),
                all_rank_orders=5040, minimum_paley_feedback_edges=minimum,
                minimum_exception_cycle_rank=1, natural_repair=natural,
                joined_sharp_cycle_debt_copies=12)


if __name__ == '__main__':
    print(json.dumps(dict(status='PASS', tournaments=tournament_checks(),
                          renewal=renewal_checks(), moduli=modulus_checks(),
                          passive=passive_checks(), adaptive_sources=adaptive_source_checks(),
                          paley_escape=paley_escape_checks()),
                     indent=2, sort_keys=True))
