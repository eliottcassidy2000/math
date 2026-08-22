#!/usr/bin/env python3
"""Exact companion for THM-3339; dependency-free and -O safe."""

from fractions import Fraction
from itertools import combinations, permutations
from math import gcd


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def mmul(a, b):
    return tuple(
        tuple(sum(a[i][k] * b[k][j] for k in range(len(b))) for j in range(len(b[0])))
        for i in range(len(a))
    )


def mvec(a, v):
    return tuple(sum(a[i][j] * v[j] for j in range(len(v))) for i in range(len(a)))


def mpow(a, n):
    out = tuple(tuple(int(i == j) for j in range(len(a))) for i in range(len(a)))
    for _ in range(n):
        out = mmul(out, a)
    return out


def inv2(a):
    d = a[0][0] * a[1][1] - a[0][1] * a[1][0]
    return (
        (a[1][1] / d, -a[0][1] / d),
        (-a[1][0] / d, a[0][0] / d),
    )


def fibs(n):
    f = [0, 1]
    while len(f) <= n:
        f.append(f[-1] + f[-2])
    return f


def det(u, v):
    return u[0] * v[1] - u[1] * v[0]


def euclid(u):
    m, n = u
    return (n * n - m * m, 2 * m * n, n * n + m * m)


def primitive_triple(u):
    raw = euclid(u)
    d = gcd(gcd(raw[0], raw[1]), raw[2])
    out = tuple(x // d for x in raw)
    if out[0] % 2 == 0:
        out = (out[1], out[0], out[2])
    return out, d


def apply_word(word, u, branches):
    for letter in word:
        u = mvec(branches[letter], u)
    return u


def compose(p, q):
    return tuple(p[q[i]] for i in range(len(p)))


def parity(p):
    inversions = sum(p[i] > p[j] for i in range(len(p)) for j in range(i + 1, len(p)))
    return -1 if inversions % 2 else 1


def generated_group(generators):
    identity = tuple(range(len(generators[0])))
    group = {identity}
    frontier = [identity]
    while frontier:
        x = frontier.pop()
        for g in generators:
            y = compose(g, x)
            if y not in group:
                group.add(y)
                frontier.append(y)
    return group


def edge_image(p, edge):
    return tuple(sorted((p[edge[0]], p[edge[1]])))


def vadd(u, v):
    return tuple((u[i] + v[i]) % 2 for i in range(2))


def mod2mat(a):
    return tuple(tuple(x % 2 for x in row) for row in a)


def mod2vec(a, v):
    return tuple(x % 2 for x in mvec(a, v))


def affine_compose(f, g):
    """Return f o g for affine maps (linear, translation) over F_2^2."""
    lf, zf = f
    lg, zg = g
    return mod2mat(mmul(lf, lg)), vadd(mod2vec(lf, zg), zf)


def affine_group(generators):
    identity = (((1, 0), (0, 1)), (0, 0))
    group = {identity}
    frontier = [identity]
    while frontier:
        x = frontier.pop()
        for g in generators:
            y = affine_compose(g, x)
            if y not in group:
                group.add(y)
                frontier.append(y)
    return group


def signed_current(e):
    return (
        e[0] * e[3],
        2 * e[1] * e[2],
        e[0] * e[2] + e[1] * e[3],
        e[2] * e[3] - e[0] * e[1],
    )


def permutation_order(p):
    x = tuple(range(len(p)))
    for n in range(1, 25):
        x = compose(p, x)
        if x == tuple(range(len(p))):
            return n
    raise RuntimeError("permutation order overflow")


def main():
    A = ((0, 1), (-1, 2))
    B = ((0, 1), (1, 2))
    C = ((1, 0), (2, 1))
    G = ((0, 1), (1, 1))
    T = (
        (Fraction(-1, 2), Fraction(1, 2)),
        (Fraction(1, 2), Fraction(1, 2)),
    )
    branches = {"A": A, "B": B, "C": C}
    require(mpow(G, 3) == mmul(A, B), "G^3=AB failed")
    require(mmul(mmul(T, mpow(G, 3)), inv2(T)) == mmul(C, B), "T G^3 T^-1=CB failed")

    # Parameter A,B,C induce the standard three Berggren child matrices.
    U3 = ((1, -2, 2), (2, -1, 2), (2, -2, 3))
    M3 = ((1, 2, 2), (2, 1, 2), (2, 2, 3))
    D3 = ((-1, 2, 2), (-2, 1, 2), (-2, 2, 3))
    triple_branches = {"A": U3, "B": M3, "C": D3}
    for m in range(1, 40):
        for n in range(m + 1, 45):
            if gcd(m, n) != 1 or (m - n) % 2 == 0:
                continue
            u = (m, n)
            for letter in "ABC":
                require(
                    euclid(mvec(branches[letter], u)) == mvec(triple_branches[letter], euclid(u)),
                    f"symmetric-square branch failed at {letter},{u}",
                )

    f = fibs(220)
    u2 = (1, 2)
    ray_counts = {"(BA)^r": 0, "A(BA)^r": 0, "C(BC)^r": 0}
    content_two = 0
    for k in range(2, 202):
        uk = (f[k], f[k + 1])
        pell = uk[1] * uk[1] - uk[0] * uk[1] - uk[0] * uk[0]
        require(pell == (-1) ** k, f"golden Pell identity failed at k={k}")
        if k % 3 == 2:
            r = (k - 2) // 3
            word = "BA" * r
            target = uk
            lane = "(BA)^r"
        elif k % 3 == 0:
            r = (k - 3) // 3
            word = "A" + "BA" * r
            target = uk
            lane = "A(BA)^r"
        else:
            r = (k - 4) // 3
            word = "C" + "BC" * r
            target = mvec(T, uk)
            require(all(x.denominator == 1 for x in target), f"T not integral at k={k}")
            target = tuple(int(x) for x in target)
            lane = "C(BC)^r"
        got = apply_word(word, u2, branches)
        require(got == target, f"three-ray address failed at k={k}: {word}")
        triple, content = primitive_triple(uk)
        require(content == (2 if k % 3 == 1 else 1), f"content law failed at k={k}")
        if content == 2:
            content_two += 1
        require(euclid(got) == triple, f"primitive normalization transplant failed at k={k}")
        triple_from_word = (3, 4, 5)
        for letter in word:
            triple_from_word = mvec(triple_branches[letter], triple_from_word)
        require(triple_from_word == triple, f"Berggren triple address failed at k={k}")
        ray_counts[lane] += 1

    # Hostile finite search for the converse Pell classification.
    pell_solutions = []
    for n in range(1, 1000):
        for m in range(1, n + 1):
            d = n * n - m * n - m * m
            if abs(d) == 1:
                pell_solutions.append((m, n))
    fib_pairs = {(f[k], f[k + 1]) for k in range(1, len(f) - 1) if f[k + 1] < 1000}
    require(set(pell_solutions) == fib_pairs, "bounded Pell converse hostile failed")

    # The internal K4 current of a generalized Euclid window.
    for m in range(1, 80):
        for n in range(m + 1, 90):
            w = (n - m, m, n, n + m)
            e = {(i, j): w[i] * w[j] for i, j in combinations(range(4), 2)}
            a, b, c = euclid((m, n))
            require(e[(0, 3)] == a, "K4 odd-leg current failed")
            require(2 * e[(1, 2)] == b, "K4 even-leg current failed")
            require(e[(0, 2)] + e[(1, 3)] == c, "K4 sum current failed")
            require(e[(2, 3)] - e[(0, 1)] == c, "K4 difference current failed")

    a2, b2, c2 = (1, 0), (0, 1), (1, 1)
    residue_cycle = []
    for k in range(2, 8):
        left = (f[k], f[k + 1])
        right = (f[k + 1], f[k + 2])
        if det(left, right) < 0:
            left, right = right, left
        require(det(left, right) == 1, f"oriented Farey flank failed at k={k}")
        x = tuple(z % 2 for z in left)
        z = tuple(z % 2 for z in right)
        y = ((x[0] + z[0]) % 2, (x[1] + z[1]) % 2)
        pi = (x, y, z)
        residue_cycle.append(pi)
        reference = (a2, b2, c2)
        perm = tuple(reference.index(q) for q in pi)
        require(parity(perm) == (-1) ** k, f"channel-order sign failed at k={k}")
        require(
            f[k - 1] * f[k + 2] - f[k] * f[k + 1] == (-1) ** k,
            f"Cassini edge identity failed at k={k}",
        )
    expected_cycle = [
        (b2, c2, a2),
        (b2, a2, c2),
        (a2, b2, c2),
        (a2, c2, b2),
        (c2, a2, b2),
        (c2, b2, a2),
    ]
    require(residue_cycle == expected_cycle, "six-state residue hexagon failed")

    edge_orders = []
    for k in range(3, 100):
        window = (f[k - 1], f[k], f[k + 1], f[k + 2])
        require(window == (f[k + 1] - f[k], f[k], f[k + 1], f[k + 1] + f[k]), "window failed")
        weights = {(i, j): window[i] * window[j] for i, j in combinations(range(4), 2)}
        order = tuple(sorted(weights, key=weights.get))
        middle = ((0, 3), (1, 2)) if k % 2 else ((1, 2), (0, 3))
        expected = ((0, 1), (0, 2), *middle, (1, 3), (2, 3))
        require(order == expected, f"six-edge order failed at k={k}")
        require(abs(weights[(0, 3)] - weights[(1, 2)]) == 1, f"Cassini gap failed at k={k}")
        edge_orders.append(order)
    require(len(set(edge_orders)) == 2, "edge-order parity should have two states")
    require((1, 1, 2, 3)[0] == (1, 1, 2, 3)[1], "root tie hostile failed")
    require(tuple(sorted(range(4), key=(1, 2, 3, 5).__getitem__)) == (0, 1, 2, 3), "T4 control 1")
    require(tuple(sorted(range(4), key=(2, 3, 5, 8).__getitem__)) == (0, 1, 2, 3), "T4 control 2")

    # S4 on the six K4 edges: parity no-go and exact owner-lift census.
    edges = tuple(combinations(range(4), 2))
    edge_index = {e: i for i, e in enumerate(edges)}
    matchings = (
        frozenset(((0, 3), (1, 2))),
        frozenset(((0, 2), (1, 3))),
        frozenset(((0, 1), (2, 3))),
    )
    matching_index = {m: i for i, m in enumerate(matchings)}
    s4 = list(permutations(range(4)))
    matching_actions = {}
    for p in s4:
        ep = tuple(edge_index[edge_image(p, e)] for e in edges)
        require(parity(ep) == 1, f"S4 edge action is not even for {p}")
        image = []
        for matching in matchings:
            moved = frozenset(edge_image(p, e) for e in matching)
            image.append(matching_index[moved])
        matching_actions[p] = tuple(image)
    require(len(set(matching_actions.values())) == 6, "matching image is not S3")
    require(sum(action == (0, 1, 2) for action in matching_actions.values()) == 4, "matching kernel is not V4")
    for action in set(matching_actions.values()):
        require(sum(got == action for got in matching_actions.values()) == 4, "matching fibre is not fourfold")

    isolated = list(range(6))
    i03, i12 = edge_index[(0, 3)], edge_index[(1, 2)]
    isolated[i03], isolated[i12] = isolated[i12], isolated[i03]
    isolated = tuple(isolated)
    require(parity(isolated) == -1, "isolated Cassini swap should be odd")
    require(
        all(tuple(edge_index[edge_image(p, e)] for e in edges) != isolated for p in s4),
        "isolated Cassini swap unexpectedly induced by S4",
    )

    sigma = (1, 0, 2)
    tau = (0, 2, 1)
    lifts_sigma = [p for p in s4 if matching_actions[p] == sigma]
    lifts_tau = [p for p in s4 if matching_actions[p] == tau]
    require(len(lifts_sigma) == len(lifts_tau) == 4, "reflection lift count failed")
    owner_splittings = 0
    transitive_lifts = 0
    lift_profile = {}
    for ps in lifts_sigma:
        for pt in lifts_tau:
            group = generated_group((ps, pt))
            fixed = tuple(v for v in range(4) if all(g[v] == v for g in group))
            key = (len(group), len(fixed), permutation_order(ps), permutation_order(pt))
            lift_profile[key] = lift_profile.get(key, 0) + 1
            if len(group) == 6 and len(fixed) == 1:
                owner_splittings += 1
            elif len(group) == 24 and len(fixed) == 0:
                transitive_lifts += 1
            else:
                raise RuntimeError(f"unexpected lift pair {key}")
    require(owner_splittings == 4 and transitive_lifts == 12, "owner lift dichotomy failed")

    # A V4 relabeling fixes all three matchings but is not a recurrence symmetry.
    v4_move = (1, 0, 3, 2)
    require(matching_actions[v4_move] == (0, 1, 2), "V4 matching hostile failed")
    require(tuple((1, 2, 3, 5)[v4_move[i]] for i in range(4)) == (2, 1, 5, 3), "window hostile failed")

    # Branch-functorial affine owner cocycle and its signed-current obstruction.
    zero2, p2, q2, r2 = (0, 0), (1, 0), (0, 1), (1, 1)
    I2 = ((1, 0), (0, 1))
    J2 = ((0, 1), (1, 0))
    linear = {letter: mod2mat(branches[letter]) for letter in "ABC"}
    require(linear == {"A": J2, "B": J2, "C": I2}, "branch mod-two shadow failed")
    affine_branches = {letter: (linear[letter], p2) for letter in "ABC"}

    def cocycle(word):
        z = zero2
        for letter in word:
            z = vadd(mod2vec(linear[letter], z), p2)
        return z

    owners = []
    for k in range(2, 202):
        if k % 3 == 2:
            rr = (k - 2) // 3
            word = "BA" * rr
            expected_owner = r2 if rr % 2 else zero2
        elif k % 3 == 0:
            rr = (k - 3) // 3
            word = "A" + "BA" * rr
            expected_owner = vadd(p2, r2 if rr % 2 else zero2)
        else:
            rr = (k - 4) // 3
            word = "C" + "BC" * rr
            expected_owner = p2 if rr % 2 == 0 else q2
        got_owner = cocycle(word)
        require(got_owner == expected_owner, f"owner cocycle formula failed at k={k}")
        owners.append(got_owner)
    require(owners[:7] == [zero2, p2, p2, r2, q2, q2, zero2], "owner period failed")
    drift = [vadd(owners[i], owners[i + 1]) for i in range(6)]
    require(drift == [p2, zero2, q2, p2, zero2, q2], "owner drift failed")
    drift_sum = zero2
    for step in drift:
        drift_sum = vadd(drift_sum, step)
    require(drift_sum == zero2, "hexagon owner holonomy failed")

    for t in (zero2, p2, q2, r2):
        changed_c = vadd(vadd(p2, t), mod2vec(I2, t))
        require(changed_c == p2, "static origin unexpectedly killed C cocycle")
    ba_affine = affine_compose(affine_branches["A"], affine_branches["B"])
    require(ba_affine == (I2, r2), "BA translation failed")
    require(affine_branches["C"] == (I2, p2), "C translation failed")
    affine_image = affine_group(tuple(affine_branches.values()))
    require(len(affine_image) == 8, "affine branch image is not D4")
    require(len({g[0] for g in affine_image}) == 2, "affine linear image failed")
    require(sum(g[0] == I2 for g in affine_image) == 4, "translation kernel failed")

    points = (zero2, p2, q2, r2)
    point_index = {x: i for i, x in enumerate(points)}

    def translate(e, t):
        return tuple(e[point_index[vadd(x, t)]] for x in points)

    for m in range(1, 40):
        for n in range(m + 1, 45):
            e = (n - m, m, n, n + m)
            a, b, c, c_again = signed_current(e)
            require(c == c_again == n * n + m * m, "base signed current failed")
            expected_orbit = {
                zero2: (a, b, c, c),
                p2: (b // 2, 2 * a, c, c),
                q2: (b // 2, 2 * a, c, -c),
                r2: (a, b, c, -c),
            }
            orbit = {t: signed_current(translate(e, t)) for t in points}
            require(orbit == expected_orbit, f"signed-current orbit failed at {(m, n)}")
            require(len(set(orbit.values())) == 4, f"signed-current stabilizer failed at {(m, n)}")
            for alpha in range(2):
                for beta in range(2):
                    t = ((alpha + 0) % 2, beta)
                    value = orbit[t]
                    switch = int(value[:2] != (a, b))
                    sigma = int(value[3] == -c)
                    decoded = ((switch + sigma) % 2, sigma)
                    require(decoded == t, f"signed-current decoder failed at {t}")

    require(signed_current((1, 2, 3, 5)) == (5, 12, 13, 13), "minimal current hostile failed")
    ba_window = (3, 5, 8, 13)
    require(signed_current(ba_window) == (39, 80, 89, 89), "BA current hostile failed")
    require(signed_current(translate(ba_window, r2)) == (39, 80, 89, -89), "BA correction hostile failed")

    print("THM-3339 VERIFIED-EXACT")
    print("golden_pell_search_n_max=999; exact_converse_proved_by_descent")
    print("fibonacci_indices_checked=2..201")
    print(f"ray_counts={ray_counts}")
    print(f"raw_content_two_count={content_two}")
    print("parameter_identity=G^3=AB; normalized_identity=T*G^3*T^-1=CB")
    print("root_to_child_words=(BA)^r | A(BA)^r | C(BC)^r")
    print("residue_cycle=(b,c,a)->(b,a,c)->(a,b,c)->(a,c,b)->(c,a,b)->(c,b,a)")
    print("cassini_sign=(-1)^k; six_edge_middle_swap=03<->12")
    print("four_vertex_controls=(1,2,3,5),(2,3,5,8): same transitive T4")
    print("root_window_hostile=(1,1,2,3): tie")
    print("s4_edge_images=24; odd_images=0; isolated_middle_swap_in_image=0")
    print("matching_image=6; kernel=4; lifts_per_matching_action=4")
    print(f"owner_fixing_S3_lift_pairs={owner_splittings}; transitive_S4_lift_pairs={transitive_lifts}")
    print(f"lift_profile={dict(sorted(lift_profile.items()))}")
    print("v4_matching_hostile=(1,2,3,5)->(2,1,5,3)")
    print("affine_owner_period=0,p,p,r,q,q; drift=p,0,q,p,0,q")
    print("affine_branch_image_order=8; linear_image=2; translation_kernel=4")
    print("signed_current_stabilizer=1; decoder=t=(s+sigma)p+sigma*q")
    print("BA_current_hostile=(39,80,89,89)->(39,80,89,-89)")
    print("ALL CHECKS PASSED")


if __name__ == "__main__":
    main()
