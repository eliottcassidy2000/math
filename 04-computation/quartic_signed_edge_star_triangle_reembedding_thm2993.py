"""Exact companion for THM-2993.

The companion verifies the simultaneous three-matching signed-edge cube of a
depressed quartic, its star/triangle covariants, every degeneration factor,
the derivative-square quotient-algebra re-embedding, and the modular
``C2*C3`` compatibility detector.  Every truth-bearing gate is an explicit
``require`` call, so optimized execution remains a genuine verification run.
"""

from itertools import permutations, product

import sympy as sp


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def compose(g, h):
    """First h, then g."""
    return tuple(g[h[i]] for i in range(len(g)))


def inverse(g):
    ans = [None] * len(g)
    for i, image in enumerate(g):
        ans[image] = i
    return tuple(ans)


def permutation_order(g):
    power = tuple(range(len(g)))
    for order in range(1, 100):
        power = compose(g, power)
        if power == tuple(range(len(g))):
            return order
    raise RuntimeError("permutation-order bound exhausted")


def generated_group(generators):
    identity = tuple(range(len(generators[0])))
    seen = {identity}
    stack = [identity]
    while stack:
        h = stack.pop()
        for g in generators:
            gh = compose(g, h)
            if gh not in seen:
                seen.add(gh)
                stack.append(gh)
    return seen


def from_cycles(n, *cycles):
    ans = list(range(n))
    for cycle in cycles:
        for i, x in enumerate(cycle):
            ans[x] = cycle[(i + 1) % len(cycle)]
    return tuple(ans)


# ---------------------------------------------------------------------------
# Cubic coefficient covariants and the two quartic factors.
# ---------------------------------------------------------------------------

T, U, Z = sp.symbols("T U Z")
a, b, c = sp.symbols("a b c")
p, q, r = sp.symbols("p q r")
W, L = sp.symbols("W L")

S = U**3 + a * U**2 + b * U + c
D = sp.expand(sp.discriminant(S, U))
E = sp.expand(
    4 * a**3 * c
    - 3 * a**2 * b**2
    + 18 * a * b * c
    + 4 * b**3
    - 135 * c**2
)
K = sp.expand(b**3 - a**3 * c)

A_star = sp.expand(
    Z**4
    + 4 * (a * b + 7 * c) * Z**3
    - 2 * E * Z**2
    + 4 * (a * b - 9 * c) * D * Z
    + D**2
)
A_triangle = sp.expand(
    Z**4
    + 4 * (a * b - 9 * c) * Z**3
    - 2 * E * Z**2
    + 4 * (a * b + 7 * c) * D * Z
    + D**2
)

require(
    sp.expand(A_star - A_triangle - 64 * c * Z * (Z**2 - D)) == 0,
    "star/triangle difference changed",
)
require(
    sp.cancel(Z**4 * A_star.subs(Z, D / Z) / D**2 - A_triangle) == 0,
    "triangle quartic is not the discriminant reciprocal of the star quartic",
)
require(
    sp.expand(sp.resultant(S, a * U + b, U) - K) == 0,
    "K is not Res(S,aU+b) in the companion convention",
)
require(
    sp.expand(sp.resultant(a * U + b, S, U) - K) == 0,
    "linear-first normalized resultant for K changed",
)

disc_A = sp.factor(sp.discriminant(A_star, Z))
require(
    sp.expand(disc_A - 2**24 * c**2 * K**2 * D**3) == 0,
    "star-quartic discriminant factorization changed",
)
require(
    sp.expand(sp.discriminant(A_triangle, Z) - disc_A) == 0,
    "triangle-quartic discriminant changed",
)

H_cross = sp.expand(
    a**6 * c**2
    - 2 * a**4 * b * c**2
    - 2 * a**3 * b**3 * c
    - 26 * a**3 * c**3
    + 29 * a**2 * b**2 * c**2
    - 2 * a * b**4 * c
    - 18 * a * b * c**3
    + b**6
    - 26 * b**3 * c**2
    + 189 * c**4
)
cross_resultant = sp.factor(sp.resultant(A_star, A_triangle, Z))
require(
    sp.expand(cross_resultant - 2**32 * c**4 * D**4 * H_cross) == 0,
    "star/triangle cross-resultant wall changed",
)
require(
    sp.expand(sp.resultant(A_star, Z**2 - D, Z) - 2**8 * D**2 * H_cross) == 0,
    "H is not the own-complement collision factor",
)
require(
    sp.expand(H_cross.subs(a, 0) - ((b**3 - 13 * c**2) ** 2 + 20 * c**4)) == 0,
    "depressed-cubic real H boundary changed",
)

grade3_disc = sp.expand((2**24 * c**2 * K**2 * D**3).subs(D, -4 * W**2 * L))
require(
    sp.expand(grade3_disc + 2**30 * (c * K * W**3) ** 2 * L**3) == 0,
    "grade-three discriminant corollary changed",
)

L0, B0, c0, Q0 = sp.symbols("L0 B0 c0 Q0", nonzero=True)
pullback_basis_det = sp.expand(
    2**12 * (2 * c0 / L0) * (B0**3 / L0**3) * (-4 * Q0**2 / L0**3)
)
pullback_disc_A = sp.expand(
    2**24 * (4 * c0**2 / L0**2) * (B0**6 / L0**6) * (-4 * Q0**2 / L0**3) ** 3
)
require(
    sp.expand(pullback_basis_det + 2**15 * c0 * B0**3 * Q0**2 / L0**7) == 0,
    "THM2473 chart basis determinant pullback changed",
)
require(
    sp.expand(pullback_disc_A + 2**32 * c0**2 * B0**6 * Q0**6 / L0**17) == 0,
    "THM2473 chart A_star discriminant pullback changed",
)
thm2473_k_hostile = {a: 0, b: 0, c: -sp.Rational(27, 8)}
require(D.subs(thm2473_k_hostile) != 0, "THM2473 K-wall hostile cubic is inseparable")
require(K.subs(thm2473_k_hostile) == 0, "THM2473 K-wall hostile left K=0")
require(
    sp.factor(A_star.subs(thm2473_k_hostile))
    == (8 * Z - 243) ** 3 * (8 * Z - 27) / 4096,
    "THM2473 K-wall hostile factorization changed",
)


# ---------------------------------------------------------------------------
# The depressed quartic realizes A_star as Res(f,Z-f'^2).
# ---------------------------------------------------------------------------

f = T**4 + p * T**2 + q * T + r
f_prime = sp.diff(f, T)
quartic_D = sp.expand(sp.discriminant(f, T))
matching_S = U**3 + 2 * p * U**2 + (p**2 - 4 * r) * U - q**2
require(
    sp.expand(sp.discriminant(matching_S, U) - quartic_D) == 0,
    "quartic and matching-cubic discriminants differ",
)

abc_substitution = {a: 2 * p, b: p**2 - 4 * r, c: -q**2}
A_star_pqr = sp.expand(A_star.subs(abc_substitution))
A_triangle_pqr = sp.expand(A_triangle.subs(abc_substitution))
K_pqr = sp.expand(K.subs(abc_substitution))
H_pqr = sp.expand(H_cross.subs(abc_substitution))

direct_A = sp.expand(sp.resultant(f, Z - f_prime**2, T))
require(
    sp.expand(direct_A - A_star_pqr) == 0,
    "derivative-square resultant does not equal A_star",
)
require(
    sp.expand(sp.Poly(A_star_pqr, Z).coeff_monomial(Z) - quartic_D * (8 * p**3 - 32 * p * r + 36 * q**2))
    == 0,
    "compressed linear coefficient of A_star changed",
)


# ---------------------------------------------------------------------------
# Quotient-algebra basis transition and an explicit inverse coordinate.
# ---------------------------------------------------------------------------


def reduce_mod_f(expr):
    return sp.expand(sp.rem(sp.expand(expr), f, T))


z = reduce_mod_f(f_prime**2)
expected_z = -8 * q * T**3 + 4 * (p**2 - 4 * r) * T**2 + 4 * p * q * T + q**2
require(sp.expand(z - expected_z) == 0, "derivative-square reduction changed")

columns = []
for exponent in range(4):
    reduced = reduce_mod_f(z**exponent)
    columns.append([sp.Poly(reduced, T).coeff_monomial(T**i) for i in range(4)])
basis_matrix = sp.Matrix(4, 4, lambda i, j: columns[j][i])
basis_det = sp.factor(basis_matrix.det())
require(
    sp.expand(basis_det - 2**12 * q**2 * K_pqr * quartic_D) == 0,
    "derivative-square basis determinant changed",
)

P0 = sp.expand(
    (4 * a**3 * c - a**2 * b**2 - 18 * a * b * c + 4 * b**3 + 27 * c**2)
    * (
        8 * a**3 * b * c
        - a**2 * b**3
        - 12 * a**2 * c**2
        + 2 * a * b**2 * c
        - 4 * b**4
        - 9 * b * c**2
    )
)
P1 = sp.expand(
    40 * a**5 * c**2
    - 14 * a**4 * b**2 * c
    + 3 * a**3 * b**4
    + 64 * a**3 * b * c**2
    - 87 * a**2 * b**3 * c
    - 378 * a**2 * c**3
    + 4 * a * b**5
    + 117 * a * b**2 * c**2
    + 100 * b**4 * c
    - 297 * b * c**3
)
P2 = sp.expand(
    -(
        4 * a**3 * b * c
        - 3 * a**2 * b**3
        + 56 * a**2 * c**2
        + 2 * a * b**2 * c
        - 4 * b**4
        + 57 * b * c**2
    )
)
P3 = sp.expand(-(2 * a**2 * c - a * b**2 + 3 * b * c))

inverse_numerator = sum(
    polynomial.subs(abc_substitution) * z**exponent
    for exponent, polynomial in enumerate((P0, P1, P2, P3))
)
require(
    reduce_mod_f(inverse_numerator - 32 * q * K_pqr * quartic_D * T) == 0,
    "explicit inverse coordinate for T changed",
)


# ---------------------------------------------------------------------------
# The signed edge cube and its physical S4 action.
# ---------------------------------------------------------------------------

S4 = tuple(permutations(range(4)))
matchings = (
    ((0, 1), (2, 3)),
    ((0, 2), (1, 3)),
    ((0, 3), (1, 2)),
)


def canonical_edge(edge):
    return tuple(sorted(edge))


def canonical_matching(blocks):
    return tuple(sorted(canonical_edge(block) for block in blocks))


matching_lookup = {canonical_matching(blocks): i for i, blocks in enumerate(matchings)}
edge_bit = []
for blocks in matchings:
    edge_bit.append({canonical_edge(blocks[0]): 0, canonical_edge(blocks[1]): 1})


def quotient_action(g):
    sigma = []
    for blocks in matchings:
        image = canonical_matching(tuple(tuple(g[x] for x in block) for block in blocks))
        require(image in matching_lookup, "matching quotient left its carrier")
        sigma.append(matching_lookup[image])
    return tuple(sigma)


def signed_pair_data(g):
    sigma = quotient_action(g)
    delta = [None, None, None]
    for source_index, blocks in enumerate(matchings):
        target_index = sigma[source_index]
        image_base = canonical_edge(tuple(g[x] for x in blocks[0]))
        delta[target_index] = edge_bit[target_index][image_base]
    return tuple(delta), sigma


def cube_action(g, bits):
    delta, sigma = signed_pair_data(g)
    image = [None, None, None]
    for source_index in range(3):
        target_index = sigma[source_index]
        image[target_index] = bits[source_index] ^ delta[target_index]
    return tuple(image)


signed_images = {signed_pair_data(g) for g in S4}
expected_h0 = {
    (delta, sigma)
    for delta in product((0, 1), repeat=3)
    if sum(delta) % 2 == 0
    for sigma in permutations(range(3))
}
require(signed_images == expected_h0, "physical signed-pair image is not H0")

stars = {
    0: (0, 0, 0),
    1: (0, 1, 1),
    2: (1, 0, 1),
    3: (1, 1, 0),
}
star_set = set(stars.values())
triangle_set = {tuple(bit ^ 1 for bit in bits) for bits in star_set}
require(star_set.isdisjoint(triangle_set), "star and triangle halves overlap")
require(len(star_set | triangle_set) == 8, "signed cube census changed")

for vertex, bits in stars.items():
    for g in S4:
        require(
            cube_action(g, bits) == stars[g[vertex]],
            "S4 action on coherent stars changed",
        )

for bits in star_set:
    require(
        {cube_action(g, bits) for g in S4} == star_set,
        "physical S4 does not act transitively on the star half",
    )
for bits in triangle_set:
    require(
        {cube_action(g, bits) for g in S4} == triangle_set,
        "physical S4 does not act transitively on the triangle half",
    )
for bits in star_set:
    require(
        tuple(bit ^ 1 for bit in bits) in triangle_set,
        "central 111 translation does not exchange parity halves",
    )

point_stabilizers = []
for vertex, bits in stars.items():
    stabilizer = {g for g in S4 if g[vertex] == vertex}
    require(len(stabilizer) == 6, "point stabilizer is not S3-sized")
    require(
        all(cube_action(g, bits) == bits for g in stabilizer),
        "point stabilizer does not fix its star origin",
    )
    require(
        len({quotient_action(g) for g in stabilizer}) == 6,
        "point stabilizer does not split the matching quotient",
    )
    point_stabilizers.append(frozenset(stabilizer))
require(len(set(point_stabilizers)) == 4, "four affine-origin splittings changed")


# ---------------------------------------------------------------------------
# Split quartic control: all eight products and the product sign constraint.
# ---------------------------------------------------------------------------

root_values = (-6, 1, 2, 3)
split_poly = sp.Poly(sp.prod(T - value for value in root_values), T)
split_p = split_poly.coeff_monomial(T**2)
split_q = split_poly.coeff_monomial(T)
split_r = split_poly.coeff_monomial(1)
require((split_p, split_q, split_r) == (-25, 60, -36), "split quartic changed")

base_sums = tuple(root_values[0] + root_values[j] for j in (1, 2, 3))
require(sp.prod(base_sums) == -split_q, "coherent signed-edge product is not -q")
split_D = sp.discriminant(split_poly.as_expr(), T)
derivative_values = tuple(sp.diff(split_poly.as_expr(), T).subs(T, x) ** 2 for x in root_values)
star_products = []
triangle_products = []

for bits in product((0, 1), repeat=3):
    signed_sums = tuple(((-1) ** bits[i]) * base_sums[i] for i in range(3))
    require(
        sp.prod(signed_sums) == (-split_q if sum(bits) % 2 == 0 else split_q),
        "signed-edge cube product parity changed",
    )
    selected_edges = [matchings[i][bits[i]] for i in range(3)]
    block_product = sp.prod(
        (root_values[edge[0]] - root_values[edge[1]]) ** 2 for edge in selected_edges
    )
    if bits in star_set:
        vertex = next(v for v, star_bits in stars.items() if star_bits == bits)
        require(block_product == derivative_values[vertex], "star product is not f'(x)^2")
        star_products.append(block_product)
    else:
        complementary_star = tuple(bit ^ 1 for bit in bits)
        vertex = next(v for v, star_bits in stars.items() if star_bits == complementary_star)
        require(
            block_product == split_D / derivative_values[vertex],
            "triangle product is not Delta/f'(x)^2",
        )
        triangle_products.append(block_product)

require(sorted(star_products) == [64, 196, 324, 254016], "split star-product census changed")
require(sorted(triangle_products) == [4, 3136, 5184, 15876], "split triangle-product census changed")
require(
    sp.Poly(A_star_pqr.subs({p: split_p, q: split_q, r: split_r}), Z)
    == sp.Poly(sp.prod(Z - value for value in star_products), Z),
    "A_star misses the split star products",
)
require(
    sp.Poly(A_triangle_pqr.subs({p: split_p, q: split_q, r: split_r}), Z)
    == sp.Poly(sp.prod(Z - value for value in triangle_products), Z),
    "A_triangle misses the split triangle products",
)


# ---------------------------------------------------------------------------
# The C2*C3 compatibility bit and the cusp-word order detector.
# ---------------------------------------------------------------------------

involution_lifts = [
    g for g in S4 if permutation_order(g) == 2 and permutation_order(quotient_action(g)) == 2
]
ternary_lifts = [
    g for g in S4 if permutation_order(g) == 3 and permutation_order(quotient_action(g)) == 3
]
require(len(involution_lifts) == 6, "order-two affine-lift census changed")
require(len(ternary_lifts) == 8, "order-three affine-lift census changed")

split_pairs = 0
full_pairs = 0
for s_lift in involution_lifts:
    fixed_line = {i for i in range(4) if s_lift[i] == i}
    require(len(fixed_line) == 2, "C2 lift does not fix an affine line")
    for t_lift in ternary_lifts:
        fixed_star = {i for i in range(4) if t_lift[i] == i}
        require(len(fixed_star) == 1, "C3 lift does not fix one star")
        compatible = bool(fixed_line & fixed_star)
        group_order = len(generated_group((s_lift, t_lift)))
        cusp_order = permutation_order(compose(s_lift, t_lift))
        if compatible:
            require(group_order == 6, "compatible free-factor gauges do not give split S3")
            require(cusp_order == 2, "compatible cusp word does not have order two")
            split_pairs += 1
        else:
            require(group_order == 24, "transverse free-factor gauges do not give full S4")
            require(cusp_order == 4, "transverse cusp word does not have order four")
            full_pairs += 1

require((split_pairs, full_pairs) == (24, 24), "C2*C3 compatibility census changed")

s_example = from_cycles(4, (0, 1))
t_split = from_cycles(4, (0, 1, 3))
t_full = from_cycles(4, (1, 2, 3))
require(quotient_action(t_split) == quotient_action(t_full), "example C3 quotients differ")
require(permutation_order(compose(s_example, t_split)) == 2, "split example cusp order changed")
require(permutation_order(compose(s_example, t_full)) == 4, "full example cusp order changed")


# ---------------------------------------------------------------------------
# Sharp degeneration and loss controls.
# ---------------------------------------------------------------------------


def specialize_abc(p_value, q_value, r_value):
    return {
        a: 2 * p_value,
        b: p_value**2 - 4 * r_value,
        c: -(q_value**2),
    }


# K=0 while q and Delta remain nonzero: A_star has a repeated root.
k_hostile = specialize_abc(sp.Integer(1), sp.Integer(1), sp.Rational(3, 4))
require(D.subs(k_hostile) == 125, "K-wall hostile discriminant changed")
require(K.subs(k_hostile) == 0, "K-wall hostile left K=0")
require(
    sp.factor(A_star.subs(k_hostile)) == (Z - 25) ** 2 * (Z**2 + 6 * Z + 25),
    "K-wall hostile factorization changed",
)

# H=0 with q*K*Delta nonzero: the two separable halves share Z=49.
h_hostile = specialize_abc(sp.Integer(-7), sp.Integer(7), sp.Integer(0))
require(D.subs(h_hostile) == 2401, "cross-wall hostile discriminant changed")
require(K.subs(h_hostile) == -16807, "cross-wall hostile K changed")
require(H_cross.subs(h_hostile) == 0, "cross-wall hostile left H=0")
require(
    sp.gcd(sp.Poly(A_star.subs(h_hostile), Z), sp.Poly(A_triangle.subs(h_hostile), Z))
    == sp.Poly(Z - 49, Z),
    "cross-wall hostile common factor changed",
)

# q=0: parity halves coincide and the derivative-square basis loses rank.
q0_hostile = specialize_abc(sp.Integer(-5), sp.Integer(0), sp.Integer(4))
require(D.subs(q0_hostile) != 0 and K.subs(q0_hostile) != 0, "q=0 hostile is not otherwise generic")
require(
    sp.expand(A_star.subs(q0_hostile) - A_triangle.subs(q0_hostile)) == 0,
    "q=0 does not collapse star and triangle covariants",
)
require(
    basis_det.subs({p: -5, q: 0, r: 4}) == 0,
    "q=0 does not collapse derivative-square primitivity",
)

# q-sign reversal preserves the unpointed cover but reverses T.
require(
    sp.expand(A_star_pqr.subs(q, -q) - A_star_pqr) == 0,
    "A_star remembers the sign of q",
)
require(
    sp.expand(A_triangle_pqr.subs(q, -q) - A_triangle_pqr) == 0,
    "A_triangle remembers the sign of q",
)
require(
    sp.expand(z.subs({T: -T, q: -q}) - z) == 0,
    "derivative-square coordinate does not survive T,q sign reversal",
)


print("THM2993 exact quartic star/triangle cube companion")
print("matching_disc_equals_quartic_disc=PASS")
print("signed_cube=8 star=4 triangle=4 physical_H0=24 origins=4")
print("star_products=64,196,324,254016")
print("triangle_products=4,3136,5184,15876")
print("Astar_minus_Atriangle=64*c*Z*(Z^2-Delta)")
print("disc_A=2^24*c^2*K^2*Delta^3")
print("cross_resultant=2^32*c^4*Delta^4*H")
print("basis_det=2^12*q^2*K*Delta")
print("explicit_inverse_degree=3")
print("free_factor_pairs=48 split=24 full=24 cusp_orders=2/4")
print("K_hostile=(p,q,r)=(1,1,3/4),Delta=125")
print("H_hostile=(p,q,r)=(-7,7,0),common=Z-49")
print("q0_hostile=(p,q,r)=(-5,0,4),halves_equal=1")
print("grade3_disc=-2^30*(c*K*W^3)^2*L^3")
print("all_exact_checks=PASS")
