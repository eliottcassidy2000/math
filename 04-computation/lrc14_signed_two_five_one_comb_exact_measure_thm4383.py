from fractions import Fraction as F
from itertools import combinations, permutations, product
from math import gcd

R = F(3, 14)


def pos(x):
    return x if x > 0 else F(0)


def E(t):
    """Period-three quadrature error, evaluated independently by a finite sum."""
    return sum(pos(t - k) for k in range(1, int(t) + 1) if k % 3) - t * t / 3


def component(p, q, k):
    """Endpoint intersection length for determinant K=5k."""
    p, q, k = F(p), F(q), abs(k)
    A = R * abs(q - p) / 5
    B = R * (q + p) / 5
    return F(5, 1) / (p * q) * (pos(B - k) - pos(A - k))


def formula_sum(p, q):
    A = R * abs(q - p) / 5
    B = R * (q + p) / 5
    return 2 * sum(component(p, q, k) for k in range(1, int(B) + 1) if k % 3)


def formula_quadrature(p, q):
    p, q = F(p), F(q)
    A = R * abs(q - p) / 5
    B = R * (q + p) / 5
    return F(6, 245) + F(10, 1) / (p * q) * (E(B) - E(A))


def mod1(x):
    return x - (x.numerator // x.denominator)


def centered(x):
    z = mod1(x)
    return z if z <= F(1, 2) else z - 1


def eligible_j(w, x):
    return [j for j in range(3) if abs(centered(w * (x + F(j, 3)))) < F(1, 14)]


def nint_error(w, y):
    # Called only away from half-integers, on exact open cells.
    z = w * y
    n = (z + F(1, 2)).numerator // (z + F(1, 2)).denominator
    return n, z - n


def physical_walls(speeds):
    walls = {F(0), F(1)}
    for w in speeds:
        for j in range(3):
            # This deliberately over-enumerates m; mod-one deduplication is exact.
            for m in range(-w - 2, 2 * w + 3):
                for sign in (-1, 1):
                    walls.add(mod1((F(m) + sign * F(1, 14)) / w - F(j, 3)))
    return sorted(walls)


def direct_physical_measure(speeds, audit=False):
    walls = physical_walls(speeds)
    total = F(0)
    checked = 0
    for left, right in zip(walls, walls[1:]):
        x = (left + right) / 2
        js = [eligible_j(w, x) for w in speeds]
        fail = all(len(v) == 1 for v in js) and len({v[0] for v in js}) == 3
        if fail:
            total += right - left
        if audit:
            y = mod1(3 * x)
            ell = (3 * x).numerator // (3 * x).denominator
            ns, es, owners = [], [], []
            for w in speeds:
                n, e = nint_error(w, y)
                ns.append(n)
                es.append(e)
                owners.append((-pow(w, -1, 3) * n) % 3)
                assert (abs(e) < R) == (len(eligible_j(w, x)) == 1)
                if abs(e) < R:
                    assert (owners[-1] - ell) % 3 == js[len(owners) - 1][0]
            owner_fail = all(abs(e) < R for e in es) and set(owners) == {0, 1, 2}
            assert fail == owner_fail
            if owner_fail:
                np, nb, nq = ns
                delta = 5 * nb - 2 * np - nq
                k = 11 * np - nb
                assert delta == 0
                assert (53 * np - nq) == 5 * k
                assert k % 3
            checked += 1
    return total, checked


def signed_branches(p, q):
    """All b for which some signed (2,5,1) relation holds."""
    out = set()
    for v in (2 * p + q, abs(2 * p - q)):
        if v and v % 5 == 0:
            out.add(v // 5)
    return sorted(out)


def admissible(p, b, q):
    return (
        len({p, b, q}) == 3
        and all(w > 0 and w % 2 and w % 3 for w in (p, b, q))
        and gcd(gcd(p, b), q) == 1
    )


def finite_sharp_scan(product_limit):
    rows = []
    # p*q < product_limit makes this exhaustive.
    for p in range(1, product_limit, 2):
        if p % 3 == 0:
            continue
        for q in range(1, (product_limit - 1) // p + 1, 2):
            if q % 3 == 0 or p * q >= product_limit:
                continue
            for b in signed_branches(p, q):
                if admissible(p, b, q):
                    rows.append((formula_sum(p, q), p, b, q))
    return sorted(set(rows), reverse=True)


def coefficient_vectors(coefficients):
    """Signed coefficient vectors modulo overall sign."""
    out = set()
    for coeffs in set(permutations(coefficients)):
        for signs in product((-1, 1), repeat=3):
            v = tuple(c * s for c, s in zip(coeffs, signs))
            first = next(x for x in v if x)
            if first < 0:
                v = tuple(-x for x in v)
            out.add(v)
    return out


def positive_intersection_rays(left, right=None):
    """Primitive positive rays cut out by two distinct signed relations."""
    rays = set()
    pairs = combinations(left, 2) if right is None else product(left, right)
    for u, v in pairs:
        if u == v:
            continue
        c = (
            u[1] * v[2] - u[2] * v[1],
            u[2] * v[0] - u[0] * v[2],
            u[0] * v[1] - u[1] * v[0],
        )
        if not all(c) or not (all(x > 0 for x in c) or all(x < 0 for x in c)):
            continue
        c = tuple(abs(x) for x in c)
        d = gcd(gcd(c[0], c[1]), c[2])
        rays.add(tuple(x // d for x in c))
    return rays


def valid_unit_ray(ray):
    return len(set(ray)) == 3 and all(x % 2 and x % 3 for x in ray)


def main():
    target = F(12, 371)
    assert formula_sum(1, 53) == target
    assert formula_quadrature(1, 53) == target
    assert [component(1, 53, k) for k in (1, 2, 3)] == [F(3, 371), F(3, 371), 0]

    direct, checked = direct_physical_measure((1, 11, 53), audit=True)
    assert direct == target

    # Hostile to "the signed speed relation and triple eligibility alone force
    # delta=0": at the first positive q-tooth centre all three tails are
    # eligible, but delta=-1 and the first two physical owners collide.
    hostile_y = F(1, 53)
    hostile = [nint_error(w, hostile_y) for w in (1, 11, 53)]
    hostile_ns = tuple(v[0] for v in hostile)
    hostile_owners = tuple((-pow(w, -1, 3) * n) % 3 for w, n in zip((1, 11, 53), hostile_ns))
    assert hostile_ns == (0, 0, 1)
    assert all(abs(v[1]) < R for v in hostile)
    assert 5 * hostile_ns[1] - 2 * hostile_ns[0] - hostile_ns[2] == -1
    assert hostile_owners == (0, 0, 1)

    # The universal quadrature upper bound is already strict above this product.
    threshold = None
    for N in range(1, 1000):
        if F(6, 245) + F(10, 3 * N) < target:
            threshold = N
            break
    assert threshold == 425
    rows = finite_sharp_scan(threshold)
    assert all(mu == formula_quadrature(p, q) for mu, p, b, q in rows)
    assert rows[0] == (target, 1, 11, 53)
    assert sum(1 for row in rows if row[0] == target) == 1

    # Formula/direct controls on every admissible primitive signed triple
    # with all speeds <= 45 (both plus and difference branches occur).
    direct_controls = []
    for mu, p, b, q in reversed(rows):
        if max(p, b, q) <= 45:
            d, _ = direct_physical_measure((p, b, q))
            assert d == mu == formula_quadrature(p, q)
            direct_controls.append((p, b, q))

    # Explicit controls for the exceptional gcd(p,q)=5 stratum.
    g5_controls = []
    for p, b, q in ((55, 23, 5), (5, 13, 55)):
        assert gcd(p, q) == 5 and gcd(p, b) == 1
        d, _ = direct_physical_measure((p, b, q))
        assert d == formula_sum(p, q) == F(12, 385)
        g5_controls.append((p, b, q))
    # Smallest lexicographic primitive example showing gcd(p,q)=5 really occurs.
    assert admissible(5, 7, 25) and 5 * 7 == 2 * 5 + 25
    assert gcd(5, 25) == 5 and gcd(5, 7) == 1
    assert direct_physical_measure((5, 7, 25))[0] == formula_sum(5, 25) == F(4, 175)
    lex_g5 = sorted((p, b, q) for _, p, b, q in rows if gcd(p, q) == 5)
    assert lex_g5[0] == (5, 7, 25)

    # Every signed branch occurs in the direct-control bank.
    branch_kinds = set()
    for p, b, q in direct_controls + g5_controls:
        if 5 * b == 2 * p + q:
            branch_kinds.add("sum")
        if 5 * b == 2 * p - q:
            branch_kinds.add("2p-q")
        if 5 * b == q - 2 * p:
            branch_kinds.add("q-2p")
    assert branch_kinds == {"sum", "2p-q", "q-2p"}

    # Exact finite coefficient-pair audit: a double presentation would put the
    # speed ray in one of these cross products.  None is distinct, odd, and a
    # three-unit.  This also proves disjointness from the two earlier sectors.
    v251 = coefficient_vectors((5, 2, 1))
    duplicate_rays = positive_intersection_rays(v251)
    cross121 = positive_intersection_rays(v251, coefficient_vectors((2, 1, 1)))
    cross141 = positive_intersection_rays(v251, coefficient_vectors((4, 1, 1)))
    assert len(v251) == 24
    assert len(duplicate_rays) == 45 and not any(map(valid_unit_ray, duplicate_rays))
    assert len(cross121) == 39 and not any(map(valid_unit_ray, cross121))
    assert len(cross141) == 51 and not any(map(valid_unit_ray, cross141))

    print("status=PASS")
    print(f"tuple=(1,11,53) measure={target} physical_cells={checked}")
    print("lattice=(n1,n11,n53)=t*(1,11,53)+k*(0,-1,-5), 3_does_not_divide_k")
    print("nonzero_components=k=-2,-1,1,2; each_length=3/371")
    print("delta_hostile=y=1/53 n=(0,0,1) owners=(0,0,1) delta=-1")
    print("gcd5_hostile=(p,b,q)=(5,7,25) measure=4/175")
    print("formula=6/245+10/(p*q)*(E(B)-E(A))")
    print(f"strict_bulk_threshold=p*q>={threshold}")
    print(f"finite_rows={len(rows)} direct_controls={len(direct_controls)} g5_controls={len(g5_controls)}")
    print(f"branch_kinds={','.join(sorted(branch_kinds))}")
    print(f"collision_rays=251x251:{len(duplicate_rays)},251x121:{len(cross121)},251x141:{len(cross141)}; valid=0")
    print("top10:")
    for mu, p, b, q in rows[:10]:
        print(f"  (p,b,q)=({p},{b},{q}) mu={mu}")


if __name__ == "__main__":
    main()
