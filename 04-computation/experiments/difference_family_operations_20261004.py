"""Exact difference-family arithmetic and typed stabilization controls.

All checks use explicit exceptions and survive python -O. No numerical
limit, prime conjecture, or interpretation of unparenthesized f is assumed.
"""
from itertools import product
from math import gcd


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def gap(pair):
    return pair[1] - pair[0]


def pair_add(p, q):
    return p[0] + q[0], p[1] + q[1]


def pair_mul(p, q):
    """Crossed product; gap(p*q)=gap(p)*gap(q)."""
    a, b = p
    c, d = q
    return a * d + b * c, a * c + b * d


def component_product(p, q):
    return p[0] * q[0], p[1] * q[1]


def coordinates(pair):
    """Nonnegative pair -> (signed gap, diagonal translation)."""
    a, b = pair
    if not isinstance(a, int) or not isinstance(b, int) or min(a, b) < 0:
        raise ValueError("expected two nonnegative integers")
    return b - a, min(a, b)


def from_coordinates(d, t):
    if not isinstance(d, int) or not isinstance(t, int) or t < 0:
        raise ValueError("gap must be integral and diagonal translation nonnegative")
    return max(-d, 0) + t, max(d, 0) + t


def decorated_add(p, q):
    d, t = p
    e, u = q
    carry = (abs(d) + abs(e) - abs(d + e)) // 2
    return d + e, t + u + carry


def decorated_mul(p, q):
    d, t = p
    e, u = q
    return d * e, abs(d) * u + abs(e) * t + 2 * t * u


def triangle(n):
    return n * (n + 1) // 2


def odd_step(a):
    if not isinstance(a, int) or a % 2 != 1:
        raise ValueError("signed odd source required")
    z = 3 * a + 1
    # z cannot vanish at an integer source.
    magnitude = abs(z)
    valuation = (magnitude & -magnitude).bit_length() - 1
    return z // 2 ** valuation


def drop(a):
    return (a - odd_step(a)) // 2


def drop_fibre(d):
    """Complete signed odd fibre of the inherited Syracuse drop K=d."""
    if not isinstance(d, int):
        raise ValueError("integral drop required")
    result = {8 * d + 1, -4 * d - 1}
    n = 6 * d + 1
    v = 3
    while 2 ** v - 3 <= abs(n):
        divisor = 2 ** v - 3
        if n % divisor == 0:
            result.add(2 * d + n // divisor)
        v += 1
    return sorted(result)


def main():
    print("difference_family_operations_20261004: exact arithmetic, explicit universes")
    print("No definition is assigned to the user's unparenthesized f expression.")

    pairs = list(product(range(4), repeat=2))
    quotient_checks = 0
    for p, q in product(pairs, repeat=2):
        need(gap(pair_add(p, q)) == gap(p) + gap(q), "quotient addition")
        need(gap(pair_mul(p, q)) == gap(p) * gap(q), "quotient product")
        need(sum(pair_mul(p, q)) == sum(p) * sum(q), "sum sidecar product")
        need(pair_mul(p, q) == pair_mul(q, p), "commutative product")
        quotient_checks += 1
    algebra_checks = 0
    for p, q, r in product(pairs, repeat=3):
        need(pair_mul(pair_mul(p, q), r) == pair_mul(p, pair_mul(q, r)), "associativity")
        need(pair_mul(p, pair_add(q, r)) == pair_add(pair_mul(p, q), pair_mul(p, r)),
             "distributivity")
        algebra_checks += 1
    for p in pairs:
        need(pair_mul(p, (0, 1)) == p, "multiplicative identity")
        need(pair_add(p, (0, 0)) == p, "additive identity")
        need(pair_mul(p, (1, 0)) == p[::-1], "negative unit")
        need(from_coordinates(*coordinates(p)) == p, "coordinate inverse")
    print(f"Raw pair semiring: {quotient_checks} pair and {algebra_checks} triple controls on {{0..3}}^2")

    decorated = list(product(range(-8, 9), range(5)))
    decorated_checks = 0
    for p, q in product(decorated, repeat=2):
        a, b = from_coordinates(*p), from_coordinates(*q)
        need(from_coordinates(*decorated_add(p, q)) == pair_add(a, b), "decorated sum")
        need(from_coordinates(*decorated_mul(p, q)) == pair_mul(a, b), "decorated product")
        decorated_checks += 2
    print(f"Lossless (gap,translation) laws: {decorated_checks} sum/product controls, gap -8..8 and translation0..4")

    p, shifted, q = (7, 11), (8, 12), (1, 2)
    need(gap(p) == gap(shifted) == 4, "owner's common difference")
    need(gap(component_product(p, q)) == 15 and gap(component_product(shifted, q)) == 16,
         "component product does not descend")
    need(gap(pair_mul(p, q)) == gap(pair_mul(shifted, q)) == 4,
         "crossed repair does descend")
    need(coordinates(p) == (4, 7) and coordinates(shifted) == (4, 8), "translation retained")
    for a, c in product(range(12), repeat=2):
        need(pair_mul((a, a), (c, c)) == (2 * a * c, 2 * a * c), "raw image parity")
    print("Hostile: (7,11)*(1,2) has component-product gap15; (8,12)*(1,2) gap16.")
    print("Crossed products both have gap4; raw zero-fibre product only contains even diagonal pairs.")

    affine_checks = 0
    for a, b, h, k in product(range(-3, 4), repeat=4):
        # Fixed nontrivial affine control; theorem in note covers all coefficients.
        F = lambda x, y: 5 * x - 2 * y + 7
        need(F(a + h, b + k) - F(a, b) == 5 * h - 2 * k, "affine quotient law")
        affine_checks += 1
    print(f"Translation-independent affine controls: {affine_checks}; nonlinear hostile above")

    triangular_checks = local_transfers = 0
    for N in range(3, 201):
        count = 0
        for A in range(N - 1):
            B = N - 2 - A
            if not A < B:
                continue
            d = B - A
            rectangle = (A + 1) * (B + 1)
            pieces = triangle(A) + triangle(B)
            need(d > 0 and d <= N - 2 and d % 2 == N % 2, "gap chart domain")
            need((N - 2 - d) // 2 == A and (N - 2 + d) // 2 == B, "gap chart inverse")
            need(4 * rectangle == N * N - d * d, "rectangle in gap coordinates")
            need(4 * pieces == N * N - 2 * N + d * d, "triangle pieces in gap coordinates")
            need(pieces + rectangle == triangle(N - 1), "owner triangular identity")
            if A + 1 < B - 1:
                new_rectangle = (A + 2) * B
                new_pieces = triangle(A + 1) + triangle(B - 1)
                need(new_rectangle - rectangle == d - 1, "rectangle transfer")
                need(new_pieces - pieces == -(d - 1), "triangle transfer")
                local_transfers += 1
            count += 1
            triangular_checks += 1
        need(count == (N - 1) // 2, "strict nonnegative split count")
    print(f"Triangular gap chart N=3..200: {triangular_checks} splits; {local_transfers} internal transfers")

    scalar_checks = 0
    for A, B, x in product(range(1, 10), range(1, 10), range(-20, 21)):
        need((x == A * B * x) == (x == 0 or A * B == 1), "scalar g classification")
        scalar_checks += 1
    modular_checks = 0
    for m in range(1, 65):
        for c in range(17):
            actual = [x for x in range(m) if (c - 1) * x % m == 0]
            d = gcd(m, c - 1)
            predicted = list(range(0, m, m // d))
            need(actual == predicted and len(actual) == d, "modular fixed subgroup")
            modular_checks += 1
    ray_checks = 0
    for c, root, L in product(range(2, 10), range(1, 7), range(9)):
        ray = {root * c ** k for k in range(L + 1)}
        moved = {c * n for n in ray}
        need(ray - moved == {root}, "root boundary")
        need(moved - ray == {root * c ** (L + 1)}, "finite truncation boundary")
        ray_checks += 1
    print(f"g scalar controls: {scalar_checks}; modular fixed subgroups: {modular_checks}; finite ray boundaries: {ray_checks}")
    print("g boundary: cS=S has no nonempty positive-integer solution for c>1; a forward ray loses its root.")
    print("Universal quotient: c=1 in Z/(c-1); modulo c itself only the zero fixed residue survives.")

    # Independent complete brute list: all fibres |d|<=D lie in |A|<=8D+1.
    D = 1000
    brute = {d: [] for d in range(-D, D + 1)}
    for A in range(-8 * D - 1, 8 * D + 2, 2):
        d = drop(A)
        if -D <= d <= D:
            brute[d].append(A)
    for d, actual in brute.items():
        expected = drop_fibre(d)
        need(actual == expected, "complete signed drop fibre")
        for A in expected:
            need(drop(A) == d, "fibre replay")
    seen = {}
    for A in range(1, 20002, 2):
        signature = drop(A), drop(odd_step(A))
        need(signature not in seen, "inherited positive two-drop injectivity")
        seen[signature] = A
        candidates = [n for n in drop_fibre(signature[0])
                      if n > 0 and drop(odd_step(n)) == signature[1]]
        need(candidates == [A], "two-drop decoder")
    need(drop(13) == drop(33) == 4, "single-drop collision")
    need(drop(39) == -10 and drop(99) == -25, "single-drop product failure")
    need((drop(1), drop(odd_step(1))) == (drop(-1), drop(odd_step(-1))) == (0, 0),
         "signed two-drop collision")
    need(odd_step(7) == 11, "owner's legal rise")
    # Every positive odd rise of gap 2r is the unique branch-one edge.
    for r in range(1, 1001):
        need([a for a in drop_fibre(-r) if a > 0] == [4 * r - 1], "unique positive rise selector")
        need(odd_step(4 * r - 1) == 6 * r - 1, "rise replay")
    print(f"Inherited drop fibres: {2*D+1} complete signed fibres, |drop|<=1000")
    print(f"Inherited two-drop decoding: {len(seen)} positive odd sources1..20001")
    print("Drop hostiles: K(13)=K(33)=4 but K(3*13)=-10, K(3*33)=-25; +/-1 both have signature(0,0).")
    print("D4 contains the unique positive odd Syracuse rise (7,11); translate (8,12) is not an odd-map edge.")
    print("PASS: quotient arithmetic proved; exact orbit guards remain separate; no f/g naming imposed.")


if __name__ == "__main__":
    main()
