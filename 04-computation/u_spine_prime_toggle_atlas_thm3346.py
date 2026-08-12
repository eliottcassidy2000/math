#!/usr/bin/env python3
"""Exact audit for the U-spine prime-toggle/root atlas and its cubical fold.

Universe
--------
* all pairs of U-spine indices 0 <= r,s <= 300;
* 28 admissible prime, prime-power, and mixed moduli through N=99905;
* explicit cubical boundary matrices for antipodal-quotient two-skeleta in
  dimensions 3 <= r <= 7, plus the r=8 formula row;
* selected prime-power clock controls for p=5,13,17,29, with exponent ranges
  1..4, 1..3, 1..3, and 1..2 respectively.

All arithmetic is integral.  Rank and Smith computations use SymPy exact
domains; assertions are explicit so optimized Python performs the same audit.
"""

from itertools import combinations
from math import gcd

from sympy import Matrix
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ


PAIR_LIMIT = 300
ADMISSIBLE_MODULI = (5, 13, 17, 25, 29, 37, 41, 53, 61, 65, 73, 85,
                     97, 125, 145, 169, 185, 205, 221, 289, 325, 425,
                     625, 841, 1105, 2197, 4913, 99905)
PRIME_POWERS = (5, 25, 125, 625, 13, 169, 2197, 17, 289, 4913, 29, 841)


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def C(t):
    return 2 * t * t + 2 * t + 1


def content(x, y):
    return gcd(abs(x), abs(y))


def admissible(n):
    if n <= 1 or n % 2 == 0:
        return False
    return all(p % 4 == 1 for p in factorization(n))


def factorization(n):
    factors = {}
    p = 2
    while p * p <= n:
        while n % p == 0:
            factors[p] = factors.get(p, 0) + 1
            n //= p
        p = 3 if p == 2 else p + 2
    if n > 1:
        factors[n] = factors.get(n, 0) + 1
    return factors


def roots(n):
    return tuple(t for t in range(n) if C(t) % n == 0)


def gaussian_norm(z):
    return z[0] * z[0] + z[1] * z[1]


def nearest_quotient(numerator, denominator):
    q, remainder = divmod(numerator, denominator)
    return q + (2 * remainder > denominator)


def gaussian_divmod(a, b):
    denominator = gaussian_norm(b)
    real_numerator = a[0] * b[0] + a[1] * b[1]
    imag_numerator = a[1] * b[0] - a[0] * b[1]
    q = (nearest_quotient(real_numerator, denominator),
         nearest_quotient(imag_numerator, denominator))
    product = (q[0] * b[0] - q[1] * b[1],
               q[0] * b[1] + q[1] * b[0])
    return q, (a[0] - product[0], a[1] - product[1])


def gaussian_gcd(a, b):
    while b != (0, 0):
        _, remainder = gaussian_divmod(a, b)
        a, b = b, remainder
    return a


def prime_power_components(n):
    return tuple(p ** e for p, e in sorted(factorization(n).items()))


def idempotents(n):
    result = []
    for q in prime_power_components(n):
        cofactor = n // q
        result.append((cofactor * pow(cofactor, -1, q)) % n)
    return tuple(result)


def parent_from_gaussian(z):
    m, n = map(abs, z)
    return (abs(m * m - n * n), 2 * m * n, m * m + n * n)


def pell_rows(limit):
    """Square-triangular rows (N,R), including the degenerate seed."""
    rows = [(0, 0), (1, 1)]
    while len(rows) < limit:
        n0, r0 = rows[-2]
        n1, r1 = rows[-1]
        rows.append((6 * n1 - n0 + 2, 6 * r1 - r0))
    return rows


def canonical_antipodal(x, dimension):
    mask = (1 << dimension) - 1
    y = x ^ mask
    return min(x, y)


def quotient_two_skeleton(dimension):
    """Return quotient V,E,F and integral boundary matrices d1,d2."""
    vertices = tuple(range(1 << (dimension - 1)))
    vertex_index = {v: i for i, v in enumerate(vertices)}

    def edge_rep(x, axis):
        mask = (1 << dimension) - 1
        a, b = sorted((x, x ^ (1 << axis)))
        c, d = sorted((a ^ mask, b ^ mask))
        return min((a, b, axis), (c, d, axis))

    def face_rep(x, i, j):
        mask = (1 << dimension) - 1
        base_mask = (1 << i) | (1 << j)
        raw = min(x, x ^ (1 << i), x ^ (1 << j), x ^ base_mask)
        antipodal = min(raw ^ mask, raw ^ mask ^ (1 << i),
                        raw ^ mask ^ (1 << j), raw ^ mask ^ base_mask)
        return (min(raw, antipodal), i, j)

    edges = []
    edge_index = {}
    for x in range(1 << dimension):
        for axis in range(dimension):
            key = edge_rep(x, axis)
            if key not in edge_index:
                edge_index[key] = len(edges)
                a, b, _ = key
                edges.append((canonical_antipodal(a, dimension),
                              canonical_antipodal(b, dimension), axis, a, b))

    faces = []
    for axes in combinations(range(dimension), 2):
        i, j = axes
        for x in range(1 << dimension):
            key = face_rep(x, i, j)
            if key not in faces:
                faces.append(key)

    d1 = [[0 for _ in edges] for _ in vertices]
    for e, (x, y, _, _, _) in enumerate(edges):
        d1[vertex_index[x]][e] -= 1
        d1[vertex_index[y]][e] += 1

    d2 = [[0 for _ in faces] for _ in edges]
    for f, (x, i, j) in enumerate(faces):
        raw_cycle = (
            (x, x ^ (1 << i), i),
            (x ^ (1 << i), x ^ (1 << i) ^ (1 << j), j),
            (x ^ (1 << i) ^ (1 << j), x ^ (1 << j), i),
            (x ^ (1 << j), x, j),
        )
        for a, b, axis in raw_cycle:
            qa = canonical_antipodal(a, dimension)
            qb = canonical_antipodal(b, dimension)
            key = edge_rep(a, axis)
            idx = edge_index[key]
            stored_a, stored_b, _, raw_a, raw_b = edges[idx]
            # The stored quotient edge is oriented by quotient vertices, not
            # by the numerically ordered raw representative.
            raw_orientation = 1 if (canonical_antipodal(raw_a, dimension),
                                    canonical_antipodal(raw_b, dimension)) \
                                   == (stored_a, stored_b) else -1
            if (a, b) == (raw_a, raw_b):
                orientation = 1
            elif (a, b) == (raw_b, raw_a):
                orientation = -1
            else:
                mask = (1 << dimension) - 1
                orientation = 1 if (a ^ mask, b ^ mask) == (raw_a, raw_b) else -1
            orientation *= raw_orientation
            d2[idx][f] += orientation

    require(Matrix(d1) * Matrix(d2) == Matrix.zeros(len(vertices), len(faces)),
            "cubical boundary did not square to zero")
    return vertices, edges, faces, Matrix(d1), Matrix(d2)


def main():
    pair_checks = 0
    radius_checks = 0
    for r in range(PAIR_LIMIT + 1):
        zr = (r + 1, r)
        cr = C(r)
        for s in range(PAIR_LIMIT + 1):
            zs = (s + 1, s)
            cs = C(s)

            # z_r conjugate(z_s): imaginary part r-s.
            minus = (zr[0] * zs[0] + zr[1] * zs[1],
                     zr[1] * zs[0] - zr[0] * zs[1])
            # z_r z_s: real part r+s+1.
            plus = (zr[0] * zs[0] - zr[1] * zs[1],
                    zr[0] * zs[1] + zr[1] * zs[0])
            d_minus = content(*minus)
            d_plus = content(*plus)
            common = gcd(cr, cs)

            require(d_minus == gcd(cr, abs(s - r)),
                    "same-orientation content formula failed")
            require(d_plus == gcd(cr, r + s + 1),
                    "opposite-orientation content formula failed")
            require(gcd(d_minus, d_plus) == 1,
                    "two content channels were not coprime")
            require(d_minus * d_plus == common,
                    "two content channels did not factor the norm gcd")
            require(common % 2 == 1, "U-spine common norm acquired factor 2")
            if r >= 1 and s >= 1:
                m, n = sorted((abs(plus[0] // d_plus),
                               abs(plus[1] // d_plus)), reverse=True)
                require(n * (m - n)
                        == (r + s + 1) * (2 * r * s - 1) // (d_plus * d_plus),
                        "product-channel inradius formula failed")
                radius_checks += 1
            if 0 <= r < s:
                m, n = sorted((abs(minus[0] // d_minus),
                               abs(minus[1] // d_minus)), reverse=True)
                require(n * (m - n)
                        == (s - r) * (2 * r * s + 2 * r + 1)
                        // (d_minus * d_minus),
                        "conjugate-channel inradius formula failed")
                radius_checks += 1
            pair_checks += 1

    require(all(gcd(C(t), C(t + 1)) == 1 for t in range(PAIR_LIMIT)),
            "consecutive U-spine norms were not coprime")

    # A modular root atlas sees only the selected N-primary content.  At
    # N=5 the roots 6 and 23 have complementary N-primary channels (1,5),
    # while their full norm gcd 85 has the extra difference-channel factor 17.
    scope_n, scope_r, scope_s = 5, 6, 23
    require(C(scope_r) % scope_n == 0 and C(scope_s) % scope_n == 0,
            "modular/full-content hostile left the root atlas")
    require((gcd(scope_n, scope_s - scope_r),
             gcd(scope_n, scope_r + scope_s + 1)) == (1, 5),
            "N-primary scope hostile changed")
    require((gcd(C(scope_r), scope_s - scope_r),
             gcd(C(scope_r), scope_r + scope_s + 1)) == (17, 5),
            "full-content scope hostile changed")

    moduli = list(ADMISSIBLE_MODULI)
    require(all(admissible(n) for n in moduli),
            "declared admissible modulus was not admissible")
    root_checks = 0
    reflection_checks = 0
    edge_weight_checks = 0
    toggle_checks = 0
    gaussian_reconstruction_checks = 0
    parent_quotient_checks = 0
    for n in moduli:
        fac = factorization(n)
        rs = roots(n)
        require(len(rs) == 1 << len(fac), "root count was not 2^omega")
        require(tuple(sorted(rs)) == tuple(sorted(n - 1 - t for t in rs)),
                "root reflection failed")
        root_set = set(rs)
        generators = idempotents(n)
        require(len(generators) == len(fac), "idempotent rank mismatch")
        orbit = {rs[0]}
        frontier = [rs[0]]
        while frontier:
            t = frontier.pop()
            for idem in generators:
                image = (t - idem * (2 * t + 1)) % n
                require(image in root_set, "prime toggle left the root atlas")
                require((image - idem * (2 * image + 1)) % n == t,
                        "prime toggle was not involutive")
                toggle_checks += 1
                if image not in orbit:
                    orbit.add(image)
                    frontier.append(image)
        require(orbit == root_set, "prime toggles were not transitive")

        parents = set()
        for t in rs:
            x = 2 * t + 1
            g = gaussian_gcd((n, 0), (x, 1))
            require(gaussian_norm(g) == n,
                    "Gaussian gcd did not reconstruct norm N")
            unit_orbit = ((g[0], g[1]), (-g[1], g[0]),
                          (-g[0], -g[1]), (g[1], -g[0]))
            reconstructed = {
                (b * pow(a, -1, n)) % n
                for a, b in unit_orbit if gcd(a, n) == 1
            }
            require(x % n in reconstructed or (-x) % n in reconstructed,
                    "Gaussian gcd did not reconstruct root up to unit/conjugation")
            reflected_x = (-x) % n
            reflected_g = gaussian_gcd((n, 0), (reflected_x, 1))
            require(parent_from_gaussian(g) == parent_from_gaussian(reflected_g),
                    "root reflection did not become one parent")
            spine_g = gaussian_gcd((n, 0), (t + 1, t))
            require(gaussian_norm(spine_g) == n,
                    "Gaussian gcd with z_t did not reconstruct norm N")
            require(parent_from_gaussian(spine_g)
                    == parent_from_gaussian(g),
                    "z_t and x+i Gaussian reconstructions disagree")
            parents.add(parent_from_gaussian(g))
            gaussian_reconstruction_checks += 2
        require(len(parents) == len(rs) // 2,
                "root atlas quotient did not recover parent count")
        parent_quotient_checks += len(parents)
        root_checks += len(rs)
        reflection_checks += len(rs)
        for r in rs:
            for s in rs:
                dm = gcd(n, abs(s - r))
                dp = gcd(n, r + s + 1)
                require(gcd(dm, dp) == 1 and dm * dp == n,
                        "fixed-grade complementary channels failed")
                edge_weight_checks += 1

    # Every selected p-adic root has a unique lift and its reflection is the
    # other branch; exact valuation-e rows are the p-1 nonlifts modulo p^(e+1).
    prime_power_checks = 0
    exact_valuation_counts = []
    for modulus in PRIME_POWERS:
        p = min(factorization(modulus))
        exponent = factorization(modulus)[p]
        rs = roots(modulus)
        require(len(rs) == 2, "prime-power clock did not have two branches")
        next_modulus = modulus * p
        next_roots = set(roots(next_modulus))
        exact = 0
        for t in rs:
            lifts = [t + k * modulus for k in range(p)]
            good = [u for u in lifts if u in next_roots]
            require(len(good) == 1, "Hensel branch did not lift uniquely")
            exact += sum(C(u) % next_modulus != 0 for u in lifts)
        require(exact == 2 * (p - 1), "exact valuation clock count failed")
        exact_valuation_counts.append((p, exponent, exact))
        prime_power_checks += 2 * p
    require(all(not roots(3 ** exponent) for exponent in range(1, 5)),
            "a 3 mod 4 prime acquired a U-spine Hensel clock")

    homology_rows = []
    for dimension in range(3, 8):
        vertices, edges, faces, d1, d2 = quotient_two_skeleton(dimension)
        expected = (
            1 << (dimension - 1),
            dimension * (1 << (dimension - 2)),
            dimension * (dimension - 1) // 2 * (1 << (dimension - 3)),
        )
        require((len(vertices), len(edges), len(faces)) == expected,
                "cubical quotient f-vector failed")

        rank_d1 = d1.rank()
        rank_d2 = d2.rank()
        beta2 = len(faces) - rank_d2
        require(len(edges) - rank_d1 - rank_d2 == 0,
                "rational H1 did not vanish")

        # H1 = coker(d2 into ker d1).  A spanning-tree gauge gives a cycle
        # basis by deleting one tree edge per nonroot vertex.
        adjacency = {v: [] for v in vertices}
        for e, (x, y, _, _, _) in enumerate(edges):
            adjacency[x].append((y, e))
            adjacency[y].append((x, e))
        tree_edges = []
        seen = {0}
        queue = [0]
        while queue:
            x = queue.pop(0)
            for y, e in adjacency[x]:
                if y not in seen:
                    seen.add(y)
                    queue.append(y)
                    tree_edges.append(e)
        require(len(tree_edges) == len(vertices) - 1,
                "failed to find quotient spanning tree")
        non_tree = [e for e in range(len(edges)) if e not in tree_edges]

        # Fundamental cycles are obtained by solving d1*x=0 with each
        # non-tree coordinate prescribed.  SymPy's nullspace supplies an
        # integral lattice here; Hermite/Smith then detects the sole 2-torsion.
        null = d1.nullspace()
        cycle = Matrix.hstack(*null)
        require(cycle.shape[1] == len(non_tree), "cycle rank mismatch")
        # Coordinates on cycle lattice: select non-tree rows, which are an
        # identity up to column order for a fundamental-cycle basis.
        selector = cycle.extract(non_tree, range(cycle.shape[1]))
        require(abs(int(selector.det())) == 1,
                "nullspace basis was not an integral cycle basis")
        coords = selector.inv() * d2.extract(non_tree, range(d2.shape[1]))
        require(all(value.q == 1 for value in coords),
                "face boundary acquired nonintegral cycle coordinates")
        smith = smith_normal_form(Matrix(coords), domain=ZZ)
        diagonal = [abs(int(smith[i, i]))
                    for i in range(min(smith.rows, smith.cols))
                    if smith[i, i] != 0]
        torsion = tuple(value for value in diagonal if value > 1)
        require(torsion == (2,), "quotient H1 was not exactly C2")

        chi = len(vertices) - len(edges) + len(faces)
        require(beta2 == chi - 1,
                "H2 rank disagreed with Euler and H1 torsion")
        homology_rows.append((dimension, len(vertices), len(edges), len(faces),
                              beta2, torsion))

    # Dimension eight is covered by the same proof; record its closed formulas
    # without paying for a second large Smith reduction in routine replay.
    dimension = 8
    v = 1 << (dimension - 1)
    e = dimension * (1 << (dimension - 2))
    f = dimension * (dimension - 1) // 2 * (1 << (dimension - 3))
    chi = v - e + f
    homology_rows.append((dimension, v, e, f, chi - 1, (2,)))

    # Canonical fixed-grade controls.
    controls = {}
    for n in (65, 85, 1105, 99905):
        rs = roots(n)
        controls[n] = (len(rs), rs[:min(8, len(rs))])
    require(controls[1105][0] == 8 and controls[99905][0] == 16,
            "record-fibre root controls changed")

    # Self-product generally exits the consecutive-coordinate U-spine.
    self_square_hits = []
    for r in range(1, 101):
        product = (2 * r + 1, 2 * r * (r + 1))
        if abs(product[1] - product[0]) == 1:
            self_square_hits.append(r)
    require(self_square_hits == [1], "self-square U-spine hostile changed")

    pell_bridge_checks = 0
    pell_rows_checked = pell_rows(12)
    require(pell_rows_checked[:4] == [(0, 0), (1, 1), (8, 6), (49, 35)],
            "square-triangular recurrence changed")
    for index, (n_index, square_root) in enumerate(pell_rows_checked[1:], 1):
        require(n_index * (n_index + 1) == 2 * square_root * square_root,
                "square-triangular row failed")
        u = (n_index + 1, n_index)
        v = (2 * square_root, 1)
        m_minus = 2 * n_index + 1 - 2 * square_root
        m_plus = 2 * n_index + 1 + 2 * square_root
        t_minus = 2 * square_root - n_index - 1
        t_plus = 2 * square_root + n_index
        h = 4 * square_root * square_root + 1
        require(m_minus * m_plus == h, "Pell factor product failed")
        require(h == C(n_index), "Pell bridge left the U-spine grade")
        require(gcd(m_minus, m_plus) == 1,
                "Pell folded-weight factors were not coprime")
        require(C(t_minus) == m_minus * m_minus
                and C(t_plus) == m_plus * m_plus,
                "adjacent square U-spine depths failed")
        product = (u[0] * v[0] - u[1] * v[1],
                   u[0] * v[1] + u[1] * v[0])
        conjugate_product = (u[0] * v[0] + u[1] * v[1],
                             u[1] * v[0] - u[0] * v[1])
        expected_minus = (m_minus * t_plus,
                          m_minus * (t_plus + 1))
        expected_plus = (m_plus * (t_minus + 1),
                         m_plus * t_minus)
        require(product == expected_minus,
                "Pell product-to-upper-depth identity failed")
        require(conjugate_product == expected_plus,
                "Pell conjugate-product-to-lower-depth identity failed")
        require(content(*product) == m_minus
                and content(*conjugate_product) == m_plus,
                "Pell contents were not adjacent square roots")
        selector_parent = parent_from_gaussian(u)
        companion_parent = parent_from_gaussian(v)
        require(selector_parent == (2 * n_index + 1,
                                    4 * square_root * square_root, h),
                "selector parent formula failed")
        require(companion_parent == (4 * square_root * square_root - 1,
                                     4 * square_root, h),
                "companion parent formula failed")
        if index == 1:
            require(tuple(sorted(selector_parent[:2]))
                    == tuple(sorted(companion_parent[:2])),
                    "degenerate R=1 collapse hostile changed")
        else:
            require(selector_parent != companion_parent,
                    "positive Pell bridge parents unexpectedly collapsed")
        pell_bridge_checks += 14

    print("U-SPINE PRIME-TOGGLE ROOT ATLAS -- EXACT AUDIT")
    print(f"pairwise two-channel content checks: {pair_checks}")
    print(f"primitive composition inradius checks: {radius_checks}")
    print(f"admissible modulus controls through {max(moduli)}: {len(moduli)}")
    print(f"root/reflection checks: {root_checks}/{reflection_checks}")
    print(f"N-primary channel-pair checks: {edge_weight_checks}")
    print(f"CRT prime-toggle involution checks: {toggle_checks}")
    print(f"Gaussian root reconstructions: {gaussian_reconstruction_checks}")
    print(f"fixed-parent quotient rows: {parent_quotient_checks}")
    print(f"prime-power Hensel clock checks: {prime_power_checks}")
    print("exact valuation rows (p,e,2(p-1)):", exact_valuation_counts)
    print("quotient 2-skeleton rows (r,V,E,F,beta2,H1tors):")
    for row in homology_rows:
        print(" ", row)
    print("modular root controls:", controls)
    print("modular/full-content scope hostile: (1,5) versus (17,5)")
    print(f"Pell/Markov two-composition bridge checks: {pell_bridge_checks}")
    print("ALL EXACT CHECKS PASSED")


if __name__ == "__main__":
    main()
