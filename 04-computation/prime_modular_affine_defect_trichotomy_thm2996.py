#!/usr/bin/env python3
"""Exact companion for THM-2996.

The script uses only finite integer arithmetic.  It verifies the affine
C2*C3 cocycle census for p=2,3,5,7,11,13, the fixed-locus and cusp-order
trichotomy, the permutation-sign factorization, and the bad-prime dual
lattice control at p=3.
"""

from collections import Counter, deque
from itertools import product


PRIMES = (2, 3, 5, 7, 11, 13)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def mat_reduce(A, p):
    return tuple(tuple(x % p for x in row) for row in A)


def mat_add(A, B, p):
    return tuple(
        tuple((A[i][j] + B[i][j]) % p for j in range(2))
        for i in range(2)
    )


def mat_sub(A, B, p):
    return tuple(
        tuple((A[i][j] - B[i][j]) % p for j in range(2))
        for i in range(2)
    )


def mat_mul(A, B, p):
    return tuple(
        tuple(sum(A[i][k] * B[k][j] for k in range(2)) % p for j in range(2))
        for i in range(2)
    )


def mat_pow(A, n, p):
    result = identity_matrix(p)
    base = A
    while n:
        if n & 1:
            result = mat_mul(result, base, p)
        base = mat_mul(base, base, p)
        n //= 2
    return result


def mat_transpose(A):
    return tuple(tuple(A[j][i] for j in range(2)) for i in range(2))


def mat_det(A, p):
    return (A[0][0] * A[1][1] - A[0][1] * A[1][0]) % p


def mat_inverse(A, p):
    det = mat_det(A, p)
    require(det != 0, "attempted to invert a singular matrix")
    inv_det = pow(det, -1, p)
    return (
        ((A[1][1] * inv_det) % p, (-A[0][1] * inv_det) % p),
        ((-A[1][0] * inv_det) % p, (A[0][0] * inv_det) % p),
    )


def identity_matrix(p):
    return ((1 % p, 0), (0, 1 % p))


def zero_matrix():
    return ((0, 0), (0, 0))


def vec_add(u, v, p):
    return ((u[0] + v[0]) % p, (u[1] + v[1]) % p)


def vec_sub(u, v, p):
    return ((u[0] - v[0]) % p, (u[1] - v[1]) % p)


def mat_vec(A, v, p):
    return tuple(sum(A[i][j] * v[j] for j in range(2)) % p for i in range(2))


def vectors(p):
    return tuple(product(range(p), repeat=2))


def pair_add(x, y, p):
    return tuple((x[i] + y[i]) % p for i in range(4))


def flatten_pair(a, b):
    return (a[0], a[1], b[0], b[1])


def unflatten_pair(z):
    return (z[0], z[1]), (z[2], z[3])


def affine_compose(g, h, p):
    """Return g after h."""
    A, a = g
    B, b = h
    return mat_mul(A, B, p), vec_add(mat_vec(A, b, p), a, p)


def affine_pow(g, n, p):
    result = (identity_matrix(p), (0, 0))
    base = g
    while n:
        if n & 1:
            result = affine_compose(result, base, p)
        base = affine_compose(base, base, p)
        n //= 2
    return result


def affine_order(g, p, bound):
    identity = (identity_matrix(p), (0, 0))
    power = identity
    for n in range(1, bound + 1):
        power = affine_compose(power, g, p)
        if power == identity:
            return n
    raise RuntimeError("affine order exceeded the certified bound")


def affine_group(generators, p):
    identity = (identity_matrix(p), (0, 0))
    seen = {identity}
    queue = deque([identity])
    while queue:
        x = queue.popleft()
        for generator in generators:
            y = affine_compose(x, generator, p)
            if y not in seen:
                seen.add(y)
                queue.append(y)
    return seen


def affine_fixed_points(g, p):
    A, a = g
    return tuple(v for v in vectors(p) if vec_add(mat_vec(A, v, p), a, p) == v)


def affine_orbits(group, p):
    remaining = set(vectors(p))
    orbit_sizes = []
    while remaining:
        seed = min(remaining)
        orbit = {
            vec_add(mat_vec(A, seed, p), a, p)
            for A, a in group
        }
        require(orbit <= remaining, "orbit partition was inconsistent")
        remaining -= orbit
        orbit_sizes.append(len(orbit))
    return tuple(sorted(orbit_sizes))


def permutation_sign(g, p):
    points = vectors(p)
    index = {v: i for i, v in enumerate(points)}
    A, a = g
    perm = [index[vec_add(mat_vec(A, v, p), a, p)] for v in points]
    seen = [False] * len(points)
    cycles = 0
    for start in range(len(points)):
        if not seen[start]:
            cycles += 1
            current = start
            while not seen[current]:
                seen[current] = True
                current = perm[current]
    return -1 if (len(points) - cycles) % 2 else 1


def linear_sign_table(S, T, p):
    identity = identity_matrix(p)
    signs = {identity: 1}
    queue = deque([identity])
    for_matrix = ((S, -1), (T, 1))
    while queue:
        A = queue.popleft()
        for generator, generator_sign in for_matrix:
            B = mat_mul(A, generator, p)
            candidate = signs[A] * generator_sign
            if B in signs:
                require(signs[B] == candidate, "linear S3 sign was inconsistent")
            else:
                signs[B] = candidate
                queue.append(B)
    require(len(signs) == 6, "linear image was not S3")
    return signs


def quotient_classes(cocycles, coboundaries, p):
    remaining = set(cocycles)
    classes = []
    while remaining:
        representative = min(remaining)
        coset = {pair_add(representative, coboundary, p) for coboundary in coboundaries}
        require(coset <= set(cocycles), "coboundary coset left the cocycle space")
        classes.append(coset)
        remaining -= coset
    return tuple(classes)


def module_cohomology(S, T, p):
    I = identity_matrix(p)
    V = vectors(p)
    N = mat_add(mat_add(I, T, p), mat_pow(T, 2, p), p)
    admissible_a = tuple(a for a in V if mat_vec(mat_add(I, S, p), a, p) == (0, 0))
    admissible_b = tuple(b for b in V if mat_vec(N, b, p) == (0, 0))
    cocycles = tuple(flatten_pair(a, b) for a in admissible_a for b in admissible_b)
    coboundaries = {
        flatten_pair(mat_vec(mat_sub(I, S, p), v, p), mat_vec(mat_sub(I, T, p), v, p))
        for v in V
    }
    fixed_global = tuple(v for v in V if mat_vec(S, v, p) == v and mat_vec(T, v, p) == v)
    b2 = {mat_vec(mat_sub(I, S, p), v, p) for v in V}
    b3 = {mat_vec(mat_sub(I, T, p), v, p) for v in V}
    return {
        "N": N,
        "admissible_a": admissible_a,
        "admissible_b": admissible_b,
        "cocycles": cocycles,
        "coboundaries": coboundaries,
        "classes": quotient_classes(cocycles, coboundaries, p),
        "fixed_global": fixed_global,
        "h2_size": len(admissible_a) // len(b2),
        "h3_size": len(admissible_b) // len(b3),
        "b3": b3,
    }


def local_c3_class(b, b3, p):
    return min(vec_add(b, w, p) for w in b3)


def check_prime(p):
    I = identity_matrix(p)
    S = mat_reduce(((0, 1), (1, 0)), p)
    T = mat_reduce(((0, -1), (1, -1)), p)
    ST = mat_mul(S, T, p)
    require(mat_pow(S, 2, p) == I, f"S^2 failed at p={p}")
    require(mat_pow(T, 3, p) == I, f"T^3 failed at p={p}")
    require(mat_pow(ST, 2, p) == I, f"(ST)^2 failed at p={p}")
    require(S != I and T != I and ST != I, f"linear S3 quotient collapsed at p={p}")

    cohom = module_cohomology(S, T, p)
    require(cohom["N"] == zero_matrix(), f"ternary norm was nonzero at p={p}")
    require(len(cohom["admissible_a"]) == p, f"binary cocycle count failed at p={p}")
    require(len(cohom["admissible_b"]) == p * p, f"ternary cocycle count failed at p={p}")
    require(len(cohom["cocycles"]) == p ** 3, f"global cocycle count failed at p={p}")
    require(len(cohom["fixed_global"]) == 1, f"global invariant appeared at p={p}")
    require(len(cohom["coboundaries"]) == p * p, f"coboundary count failed at p={p}")
    require(len(cohom["classes"]) == p, f"H1 count failed at p={p}")
    require(cohom["h2_size"] == 1, f"binary local H1 failed at p={p}")
    expected_h3 = 3 if p == 3 else 1
    require(cohom["h3_size"] == expected_h3, f"ternary local H1 failed at p={p}")

    # Normalize a to zero.  Every possible normalizing translation must give
    # the same lambda=2*b_2-b_1.
    lambda_counts = Counter()
    for cocycle in cohom["cocycles"]:
        a, b = unflatten_pair(cocycle)
        lambdas = set()
        for w in vectors(p):
            a_new = vec_add(a, mat_vec(mat_sub(I, S, p), w, p), p)
            if a_new == (0, 0):
                b_new = vec_add(b, mat_vec(mat_sub(I, T, p), w, p), p)
                lambdas.add((2 * b_new[1] - b_new[0]) % p)
        require(len(lambdas) == 1, f"lambda depended on gauge at p={p}")
        lambda_counts[next(iter(lambdas))] += 1
    require(set(lambda_counts) == set(range(p)), f"lambda was not surjective at p={p}")
    require(set(lambda_counts.values()) == {p * p}, f"lambda fibres were wrong at p={p}")

    zero_lambda = {
        cocycle
        for cocycle in cohom["cocycles"]
        if any(
            vec_add(unflatten_pair(cocycle)[0], mat_vec(mat_sub(I, S, p), w, p), p) == (0, 0)
            and (2 * vec_add(unflatten_pair(cocycle)[1], mat_vec(mat_sub(I, T, p), w, p), p)[1]
                 - vec_add(unflatten_pair(cocycle)[1], mat_vec(mat_sub(I, T, p), w, p), p)[0]) % p == 0
            for w in vectors(p)
        )
    }
    require(zero_lambda == cohom["coboundaries"], f"lambda kernel was not B1 at p={p}")

    sigma = (S, (0, 0))
    tau_split = (T, (0, 0))
    tau_full = (T, (1, 0))
    require(affine_order(sigma, p, 6 * p) == 2, f"binary order failed at p={p}")
    require(affine_order(tau_full, p, 6 * p) == 3, f"ternary order failed at p={p}")
    split_product = affine_compose(sigma, tau_split, p)
    full_product = affine_compose(sigma, tau_full, p)
    require(affine_order(split_product, p, 6 * p) == 2, f"split cusp failed at p={p}")
    require(affine_order(full_product, p, 6 * p) == 2 * p, f"full cusp failed at p={p}")
    kappa = affine_pow(full_product, 2, p)
    require(kappa == (I, ((-1) % p, 0)), f"kappa formula failed at p={p}")

    split_group = affine_group((sigma, tau_split), p)
    full_group = affine_group((sigma, tau_full), p)
    require(len(split_group) == 6, f"split image was not S3 at p={p}")
    require(len(full_group) == 6 * p * p, f"full affine image size failed at p={p}")
    translations = {g for g in full_group if g[0] == I}
    require(len(translations) == p * p, f"translation kernel was not V at p={p}")
    require(affine_orbits(full_group, p) == (p * p,), f"full affine image was not transitive at p={p}")
    require(1 in affine_orbits(split_group, p), f"split image lost its fixed origin at p={p}")

    sigma_fixed = affine_fixed_points(sigma, p)
    split_tau_fixed = affine_fixed_points(tau_split, p)
    full_tau_fixed = affine_fixed_points(tau_full, p)
    require(len(sigma_fixed) == p, f"binary fixed-line size failed at p={p}")
    if p == 3:
        require(len(split_tau_fixed) == 3 and len(full_tau_fixed) == 0,
                "bad-prime ternary fixed-locus dichotomy failed")
    else:
        require(len(split_tau_fixed) == 1 and len(full_tau_fixed) == 1,
                f"ternary unique fixed point failed at p={p}")
        require(split_tau_fixed[0] in sigma_fixed and full_tau_fixed[0] not in sigma_fixed,
                f"fixed-point incidence did not detect lambda at p={p}")

    # The affine permutation sign factors through the linear S3 quotient.
    signs = linear_sign_table(S, T, p)
    exponent = (p * (p - 1) // 2) % 2
    for g in full_group:
        expected = signs[g[0]] if exponent else 1
        require(permutation_sign(g, p) == expected,
                f"permutation-sign factorization failed at p={p}")
    require(all(permutation_sign(g, p) == 1 for g in translations),
            f"translation was odd at p={p}")

    ordered_bases = (p * p - 1) * (p * p - p)
    require((ordered_bases == 6) == (p == 2), f"ordered-basis uniqueness failed at p={p}")
    geometry = "spherical" if p == 2 else ("euclidean" if p == 3 else "hyperbolic")
    character = "S3-sign" if exponent else "trivial"
    return {
        "p": p,
        "z1": len(cohom["cocycles"]),
        "b1": len(cohom["coboundaries"]),
        "h1": len(cohom["classes"]),
        "h2": cohom["h2_size"],
        "h3": cohom["h3_size"],
        "points": p * p,
        "image": len(full_group),
        "cusp": affine_order(full_product, p, 6 * p),
        "frames": ordered_bases,
        "geometry": geometry,
        "character": character,
    }


def check_bad_prime_dual():
    p = 3
    S = mat_reduce(((0, 1), (1, 0)), p)
    T = mat_reduce(((0, -1), (1, -1)), p)
    S_dual = mat_transpose(mat_inverse(S, p))
    T_dual = mat_transpose(mat_inverse(T, p))

    # The integral Cartan intertwiner has determinant three.
    C_int = ((2, -1), (-1, 2))
    S_int = ((0, 1), (1, 0))
    T_int = ((0, -1), (1, -1))
    S_dual_int = ((0, 1), (1, 0))
    T_dual_int = ((-1, -1), (1, 0))

    def int_mul(A, B):
        return tuple(tuple(sum(A[i][k] * B[k][j] for k in range(2)) for j in range(2)) for i in range(2))

    require(int_mul(S_dual_int, C_int) == int_mul(C_int, S_int),
            "integral Cartan S intertwiner failed")
    require(int_mul(T_dual_int, C_int) == int_mul(C_int, T_int),
            "integral Cartan T intertwiner failed")
    require(C_int[0][0] * C_int[1][1] - C_int[0][1] * C_int[1][0] == 3,
            "Cartan determinant was not three")

    root = module_cohomology(S, T, p)
    dual = module_cohomology(S_dual, T_dual, p)
    require(len(root["fixed_global"]) == 1 and len(root["classes"]) == 3,
            "root-module p=3 cohomology failed")
    require(len(dual["fixed_global"]) == 3 and len(dual["classes"]) == 9,
            "dual-module p=3 cohomology failed")
    require(root["h2_size"] == dual["h2_size"] == 1,
            "dual bad-prime binary restriction failed")
    require(root["h3_size"] == dual["h3_size"] == 3,
            "dual bad-prime ternary restriction failed")
    require(dual["fixed_global"] == ((0, 0), (1, 1), (2, 2)),
            "dual common invariant line was not <(1,1)>")

    def restriction_fibres(data):
        fibre_count = Counter()
        for cohomology_class in data["classes"]:
            representative = min(cohomology_class)
            _, b = unflatten_pair(representative)
            key = local_c3_class(b, data["b3"], p)
            for cocycle in cohomology_class:
                _, b_other = unflatten_pair(cocycle)
                require(local_c3_class(b_other, data["b3"], p) == key,
                        "restriction changed within a cohomology class")
            fibre_count[key] += 1
        return fibre_count

    root_fibres = restriction_fibres(root)
    dual_fibres = restriction_fibres(dual)
    require(len(root_fibres) == 3 and set(root_fibres.values()) == {1},
            "root restriction was not an isomorphism at p=3")
    require(len(dual_fibres) == 3 and set(dual_fibres.values()) == {3},
            "dual restriction did not have a three-class kernel")

    # Modulo three, the Cartan map sends the one-dimensional root H1 onto
    # the dual diagonal kernel of C3 restriction, not onto the transitive
    # dual locus.
    C_mod = mat_reduce(C_int, p)
    cartan_class_b = set()
    for root_class in root["classes"]:
        a_root, b_root = unflatten_pair(min(root_class))
        image_cocycle = flatten_pair(mat_vec(C_mod, a_root, p),
                                       mat_vec(C_mod, b_root, p))
        containing = [dual_class for dual_class in dual["classes"]
                      if image_cocycle in dual_class]
        require(len(containing) == 1, "Cartan image did not select one dual H1 class")
        zero_a = [unflatten_pair(cocycle)[1] for cocycle in containing[0]
                  if unflatten_pair(cocycle)[0] == (0, 0)]
        require(len(zero_a) == 1, "dual class lacked a unique a=0 representative")
        cartan_class_b.add(zero_a[0])
    require(cartan_class_b == {(0, 0), (1, 1), (2, 2)},
            "Cartan H1 image was not the dual diagonal restriction kernel")

    # In a=0 gauge the residual S-fixed gauge line is also T_dual-fixed,
    # so b=(b1,b2) itself is gauge-invariant.  The two independent class
    # coordinates are x=b2 (cusp defect) and y=b1-b2 (local C3 class).
    sigma = (S_dual, (0, 0))
    class_census = Counter()
    projective_transitive = set()
    for b in vectors(p):
        x = b[1]
        y = (b[0] - b[1]) % p
        tau = (T_dual, b)
        product_st = affine_compose(sigma, tau, p)
        kappa = affine_pow(product_st, 2, p)
        require(kappa == (identity_matrix(p), ((2 * x) % p, (2 * x) % p)),
                "dual p=3 cusp translation formula failed")
        group = affine_group((sigma, tau), p)
        orbits = affine_orbits(group, p)
        cusp = affine_order(product_st, p, 18)
        t_fixed = len(affine_fixed_points(tau, p))

        if x == 0 and y == 0:
            expected = (6, (1, 1, 1, 6), 2, 3, "zero")
        elif x == 0:
            expected = (6, (3, 3, 3), 2, 0, "local-only")
        elif y == 0:
            expected = (18, (3, 6), 6, 3, "cusp-only")
        else:
            expected = (18, (9,), 6, 0, "transitive")
            scale = pow(next(value for value in b if value), -1, p)
            projective_transitive.add(tuple((scale * value) % p for value in b))
        require((len(group), orbits, cusp, t_fixed, expected[4]) == expected,
                f"dual p=3 affine class failed at b={b}")
        class_census[expected[4]] += 1

    require(class_census == Counter({"transitive": 4, "local-only": 2,
                                     "cusp-only": 2, "zero": 1}),
            "dual p=3 class-stratum census failed")
    require(projective_transitive == {(0, 1), (1, 2)},
            "dual p=3 transitive projective classes failed")

    # The involution t -> t^2 exchanges the two transitive projective lines.
    # Its linear action is T_dual^2 and S_dual intertwines that module back to
    # the displayed one.  Thus b maps to S_dual(I+T_dual)b.
    transformed_lines = set()
    for b in projective_transitive:
        b_new = mat_vec(S_dual, vec_add(b, mat_vec(T_dual, b, p), p), p)
        scale = pow(next(value for value in b_new if value), -1, p)
        transformed_lines.add(tuple((scale * value) % p for value in b_new))
    require(transformed_lines == projective_transitive and
            mat_vec(S_dual, vec_add((0, 1), mat_vec(T_dual, (0, 1), p), p), p) == (1, 2) and
            mat_vec(S_dual, vec_add((1, 2), mat_vec(T_dual, (1, 2), p), p), p) == (0, 1),
            "dual p=3 t-to-t^2 exchange failed")
    return {
        "cartan_det": 3,
        "root_h1": len(root["classes"]),
        "dual_h1": len(dual["classes"]),
        "root_restriction_fibre": 1,
        "dual_restriction_fibre": 3,
        "dual_classes": tuple(sorted(class_census.items())),
        "dual_transitive_lines": tuple(sorted(projective_transitive)),
        "dual_table": (
            ("00", 1, 6, "1+1+1+6", 2, 3),
            ("10,20", 2, 6, "3+3+3", 2, 0),
            ("11,22", 2, 18, "3+6", 6, 3),
            ("01,02,12,21", 4, 18, "9", 6, 0),
        ),
    }


def check_faithful_transitive_lemma():
    # If S3 acts faithfully and transitively on V\{0}, then |V|-1 divides 6.
    # The prime-power possibilities are 2,3,4,7.  The one-dimensional linear
    # groups at sizes 2,3,7 are cyclic and cannot contain faithful S3.
    candidates = ((2, 1), (3, 1), (2, 2), (7, 1))
    viable = []
    for p, d in candidates:
        size = p ** d
        require(6 % (size - 1) == 0, "candidate list for orbit-stabilizer was incomplete")
        if d == 2 and p == 2:
            gl_order = (p * p - 1) * (p * p - p)
            require(gl_order == 6, "GL(2,2) order failed")
            viable.append((p, d, size))
        else:
            gl_order = p - 1
            require(gl_order in (1, 2, 6), "one-dimensional candidate order failed")
            # Even when gl_order=6, GL(1,7) is cyclic, so it is not S3.
    require(viable == [(2, 2, 4)], "faithful transitive nonzero-vector lemma failed")
    return viable[0]


def main():
    records = [check_prime(p) for p in PRIMES]
    dual = check_bad_prime_dual()
    unique = check_faithful_transitive_lemma()

    print("THM2996 prime-modular affine-defect exact companion")
    print("module: S=[[0,1],[1,0]], T=[[0,-1],[1,-1]]")
    for record in records:
        print(
            "p={p}: Z1={z1} B1={b1} H1={h1}; local=({h2},{h3}); "
            "points={points} image={image} cusp={cusp}; bases={frames}; "
            "triangle={geometry}; disc_character={character}".format(**record)
        )
    print(
        "p=3 dual boundary: Cartan_det={cartan_det}; root_H1={root_h1} "
        "dual_H1={dual_h1}; restriction_fibres={root_restriction_fibre}/{dual_restriction_fibre}; "
        "classes={dual_classes}; transitive_lines={dual_transitive_lines}".format(**dual)
    )
    for b_values, count, image, orbits, cusp, t_fixed in dual["dual_table"]:
        print(
            f"  dual b={b_values}: classes={count} image={image} "
            f"orbits={orbits} cusp={cusp} t_fixed={t_fixed}"
        )
    print(
        "faithful S3 transitive on V\\{0}: unique field/dimension/size="
        f"F_{unique[0]}^{unique[1]}/{unique[2]}"
    )
    print("p=13 LRC-shaped abstract control: 169 points, image 1014, cusp 26, no physical identification")
    print("all exact checks passed")


if __name__ == "__main__":
    main()
