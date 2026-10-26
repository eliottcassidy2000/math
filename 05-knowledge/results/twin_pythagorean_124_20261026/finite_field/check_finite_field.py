#!/usr/bin/env python3
"""Exact finite universes for the conic / paired-design bridge. Stdlib only.

No repository access, network, assert statements, or required output argument.
Use: python3 -I -B check_finite_field.py > receipt.json
"""
from collections import Counter, deque
from itertools import combinations, permutations, product
import json


GATES = 0


def need(condition, description):
    global GATES
    GATES += 1
    if not condition:
        raise RuntimeError(description)


def mul(x, y, degree=3, modulus=11):
    answer = 0
    while y:
        if y & 1:
            answer ^= x
        x <<= 1
        if x & (1 << degree):
            x ^= modulus
        y >>= 1
    return answer


def power(x, exponent, degree=3, modulus=11):
    answer = 1
    while exponent:
        if exponent & 1:
            answer = mul(answer, x, degree, modulus)
        x = mul(x, x, degree, modulus)
        exponent >>= 1
    return answer


def xor(values):
    value = 0
    for item in values:
        value ^= item
    return value


def compose(left, right):
    return tuple(left[x] for x in right)


def mobius(matrix, t, prime=7):
    a, b, c, d = matrix
    numerator, denominator = (a, c) if t == prime else ((a*t+b) % prime, (c*t+d) % prime)
    return prime if denominator == 0 else numerator * pow(denominator, -1, prime) % prime


def permutation(matrix):
    return tuple(mobius(matrix, t) for t in range(8))


def image(design, perm):
    return frozenset(tuple(sorted(perm[t] for t in block)) for block in design)


def generated(generators, size):
    identity = tuple(range(size))
    seen = {identity}
    todo = deque([identity])
    while todo:
        old = todo.popleft()
        for gen in generators:
            new = compose(gen, old)
            if new not in seen:
                seen.add(new)
                todo.append(new)
    return seen


def affine_defects(permutation_on_field):
    f = permutation_on_field
    return [(x, y, f[x] ^ f[y] ^ f[x ^ y] ^ f[0])
            for x in range(len(f)) for y in range(len(f))
            if f[x] ^ f[y] ^ f[x ^ y] ^ f[0]]


def main():
    # Literal eight-point conic and its independently enumerated norm circle.
    infinity = 7
    conic = {}
    for t in range(7):
        z = (1 + t*t) % 7
        need(z != 0, "nonsplit F7 circle chart")
        conic[t] = (((t*t-1)*pow(z, -1, 7)) % 7, 2*t*pow(z, -1, 7) % 7)
    conic[infinity] = (1, 0)
    circle = {(x, y) for x in range(7) for y in range(7) if (x*x+y*y) % 7 == 1}
    need(len(circle) == 8 and set(conic.values()) == circle, "all conic points")
    def circle_product(v, w):
        x, y = v
        a, b = w
        return ((x*a-y*b) % 7, (x*b+y*a) % 7)
    circle_powers = [(1, 0)]
    for unused in range(8):
        circle_powers.append(circle_product(circle_powers[-1], (2, 2)))
    need(len(set(circle_powers[:-1])) == 8 and circle_powers[-1] == (1, 0), "circle is cyclic of order8")
    for t, u in product(range(7), repeat=2):
        x, y = conic[t]
        a, b = conic[u]
        determinant = (x*b-y*a) % 7
        formula = 2*(t-u)*(1+t*u)*pow((1+t*t)*(1+u*u) % 7, -1, 7) % 7
        need(determinant == formula, "circle chord determinant identity")
    need(conic[1] == (0, 1) and conic[6] == (0, 6), "antipodal hostile")
    antipodes = [list(pair) for pair in combinations(range(7), 2)
                 if all((a+b) % 7 == 0 for a, b in zip(conic[pair[0]], conic[pair[1]]))]
    need(len(antipodes) == 3, "punctured circle retains three antipodal ties")

    E = [power(2, t) for t in range(7)] + [0]
    need(E == [1, 2, 4, 3, 6, 7, 5, 0], "fixed exponential marking")
    need(set(E) == set(range(8)), "marking bijective")
    Dplus = frozenset(B for B in combinations(range(8), 4) if xor(E[t] for t in B) == 0)
    Dminus = frozenset(B for B in combinations(range(8), 4) if xor(power(E[t], 3) for t in B) == 0)
    need(len(Dplus) == len(Dminus) == 14 and not Dplus & Dminus, "two disjoint designs")
    for design in (Dplus, Dminus):
        triples = Counter(triple for block in design for triple in combinations(block, 3))
        need(len(triples) == 56 and set(triples.values()) == {1}, "S(3,4,8) exact coverage")
        need(all(tuple(sorted(set(range(8))-set(B))) in design for B in design), "complement closure")
    for x, y, z in combinations(range(8), 3):
        w = x ^ y ^ z
        need(xor(power(v, 3) for v in (x, y, z, w)) == mul(mul(x ^ y, x ^ z), y ^ z), "cube identity")
        need(mul(mul(x ^ y, x ^ z), y ^ z) != 0, "cube identity disjointness")

    pgl = {}
    for matrix in product(range(7), repeat=4):
        a, b, c, d = matrix
        determinant = (a*d-b*c) % 7
        if determinant == 0:
            continue
        perm = permutation(matrix)
        square = determinant in {1, 2, 4}
        if perm in pgl:
            need(pgl[perm] == square, "determinant character independent of scalar representative")
        pgl[perm] = square
    psl = {g for g in pgl if pgl[g]}
    need(len(pgl) == 336 and len(psl) == 168, "projective group orders")
    for g, square in pgl.items():
        need(image(Dplus, g) == (Dplus if square else Dminus), "PGL action on positive design")
        need(image(Dminus, g) == (Dminus if square else Dplus), "PGL action on negative design")
    T = permutation((1, 1, 0, 1))
    S = permutation((0, -1, 1, 0))
    J = permutation((3, 0, 0, 1))
    need(generated([T, S], 8) == psl, "T and S generate PSL")
    def L(x):
        return (x & 1) | ((x & 2) << 1) | ((x & 4) >> 1)
    for t in range(8):
        need(E[T[t]] == mul(2, E[t]), "T is field-linear")
        need(E[S[t]] == L(E[t]) ^ 1, "S is field-affine")
        need(E[J[t]] == power(E[t], 3), "J is field-cube map")
    need(L(mul(2, 2)) != mul(L(2), L(2)), "L is not a field automorphism")

    # Independent full symmetric-group oracle, compared as sets, not just orders.
    single_aut, ordered_aut, pair_swappers = set(), set(), set()
    for perm in permutations(range(8)):
        first = image(Dplus, perm)
        if first == Dplus:
            single_aut.add(perm)
            if image(Dminus, perm) == Dminus:
                ordered_aut.add(perm)
        elif first == Dminus and image(Dminus, perm) == Dplus:
            pair_swappers.add(perm)
    need(len(single_aut) == 1344, "full single-design automorphism count")
    need(ordered_aut == psl, "full ordered-pair automorphism group exactly PSL")
    need(ordered_aut | pair_swappers == set(pgl), "full unordered-pair automorphism group exactly PGL")

    branches = {"U": (2, -1, 1, 0), "A": (2, 1, 1, 0), "D": (1, 2, 0, 1)}
    branch_perms = {name: permutation(M) for name, M in branches.items()}
    reciprocal = permutation((0, 1, 1, 0))
    common_r_matrices = {"U": (0, 1, -1, 2), "A": (0, 1, 1, 2), "D": (1, 0, 2, 1)}
    for name, matrix in common_r_matrices.items():
        need(permutation(matrix) == compose(reciprocal, compose(branch_perms[name], reciprocal)), "common r=n/m gauge conjugacy")
    need(reciprocal[2] == 4 and reciprocal[6] == 6, "common r root and collision gauges")
    need(generated([branch_perms["U"], branch_perms["D"]], 8) == psl, "outer branches generate PSL")
    need(generated(list(branch_perms.values()), 8) == set(pgl), "all branches generate PGL")
    for name, perm in branch_perms.items():
        for t, x in enumerate(E):
            expected = mul(4, x) if name == "D" else mul(4, L(x) ^ 1) if name == "U" else mul(4, power(L(x) ^ 1, 6))
            need(E[perm[t]] == expected, "branch polynomial intertwining")
        need(pgl[perm] == (name != "A"), "branch determinant character")
    def integer_word(word):
        m, n = 2, 1
        for name in word:
            a, b, c, d = branches[name]
            m, n = a*m+b*n, c*m+d*n
        return m, n
    def phi(pair):
        m, n = pair
        return (m*m-n*n, 2*m*n, m*m+n*n)
    shallow = [integer_word(word) for word in ("", "U", "A", "D")]
    need(len({m*pow(n, -1, 7) % 7 for m, n in shallow}) == 4, "no depth-one endpoint collision")
    UU, A = integer_word("UU"), integer_word("A")
    need(UU == (4, 3) and A == (5, 2), "first collision Euclid pairs")
    need(phi(UU) == (7, 24, 25) and phi(A) == (21, 20, 29), "first collision integer triples")
    need(UU[0]*pow(UU[1], -1, 7) % 7 == A[0]*pow(A[1], -1, 7) % 7 == 6, "same endpoint opposite chirality")
    lifted_states = {(2, 0)}
    todo = deque(lifted_states)
    while todo:
        t, bit = todo.popleft()
        for name, perm in branch_perms.items():
            state = (perm[t], bit ^ (name == "A"))
            if state not in lifted_states:
                lifted_states.add(state)
                todo.append(state)
    need(lifted_states == set(product(range(8), range(2))), "all sixteen point-and-chirality states reachable")

    # Puncture at infinity: literal Paley relation and two derived Fano systems.
    lines_plus = frozenset(tuple(t for t in B if t != infinity) for B in Dplus if infinity in B)
    lines_minus = frozenset(tuple(t for t in B if t != infinity) for B in Dminus if infinity in B)
    def edge(u, v):
        return (v-u) % 7 in {1, 2, 4}
    cyclic = frozenset(B for B in combinations(range(7), 3) if all(sum(edge(v, w) for w in B if w != v) == 1 for v in B))
    need(lines_plus | lines_minus == cyclic and len(cyclic) == 14, "derived Fanos equal all Paley cyclic triples")
    for lines in (lines_plus, lines_minus):
        pairs = Counter(pair for B in lines for pair in combinations(B, 2))
        need(len(pairs) == 21 and set(pairs.values()) == {1}, "derived Fano plane")
    paley_aut = {p for p in permutations(range(7)) if all(edge(u, v) == edge(p[u], p[v]) for u in range(7) for v in range(7) if u != v)}
    psl_pole = {g[:7] for g in psl if g[7] == 7}
    need(paley_aut == psl_pole and len(paley_aut) == 21, "full Paley group equals PSL pole stabilizer")
    for g, square in pgl.items():
        if g[7] != 7:
            continue
        for u, v in combinations(range(7), 2):
            need(edge(g[u], g[v]) == (edge(u, v) if square else not edge(u, v)), "pole affine maps preserve or reverse Paley")
    path_count = sum(all(edge(p[i], p[i+1]) for i in range(6)) for p in permutations(range(7)))
    need(path_count == 189 and path_count // len(paley_aut) == 9, "Paley H and normalized count")
    need(edge(1, 6) != edge(6, 1), "Paley distinguishes antipodal hostile")

    # Characteristic two: equal point counts do not identify the equations.
    projective8 = set()
    for v in product(range(8), repeat=3):
        if v == (0, 0, 0):
            continue
        lead = next(x for x in v if x)
        inverse = power(lead, 6)
        projective8.add(tuple(mul(inverse, x) for x in v))
    naive, smooth = [], []
    for x, y, z in sorted(projective8):
        naive_value = mul(x, x) ^ mul(y, y) ^ mul(z, z)
        need(naive_value == mul(x ^ y ^ z, x ^ y ^ z), "characteristic-two double line identity")
        if naive_value == 0:
            naive.append((x, y, z))
        if mul(x, z) ^ mul(y, y) == 0:
            smooth.append((x, y, z))
            need(x != 0 or z != 0, "smooth conic has nonzero gradient")
    need(len(projective8) == 73 and len(naive) == len(smooth) == 9 and naive != smooth, "double line versus smooth conic point counts")

    # Exhaust all primitive-generator exponential markings at 3,7,31 only.
    marking_tests = {}
    for prime, degree, modulus in ((3, 2, 7), (7, 3, 11), (31, 5, 37)):
        base = [power(2, t, degree, modulus) for t in range(prime)]
        need(set(base) == set(range(1, 1 << degree)), "specified polynomial gives primitive generator")
        affine_k, witnesses = [], []
        for k in range(1, prime):
            mark = [base[k*t % prime] for t in range(prime)] + [0]
            unmark = {x: t for t, x in enumerate(mark)}
            transformed = []
            for x in range(1 << degree):
                t = unmark[x]
                u = 0 if t == prime else prime if t == 0 else -pow(t, -1, prime) % prime
                transformed.append(mark[u])
            defects = affine_defects(transformed)
            need(len(set(transformed)) == 1 << degree, "transported inversion is bijective")
            if not defects:
                affine_k.append(k)
            else:
                witnesses.append({"k": k, "first_defect": defects[0]})
            if prime == 31 and k == 1:
                need(defects[0] == (1, 2, 29), "canonical field32 affine defect")
                need([transformed[x] for x in (0, 1, 2, 3)] == [1, 0, 18, 14], "canonical field32 plane image")
        need(len(affine_k) == (prime-1 if prime in (3, 7) else 0), "finite marking boundary")
        marking_tests[str(prime)] = {"degree": degree, "modulus_bits": modulus, "primitive_markings": prime-1, "affine_markings": affine_k, "nonaffine_witnesses": witnesses}

    receipt = {
        "status": "FINITE-EXACT; self-contained elementary proof in README.md",
        "gates": GATES,
        "group_orders": {"PSL2_7": len(psl), "PGL2_7": len(pgl), "single_design_aut": len(single_aut), "ordered_pair_aut": len(ordered_aut), "unordered_pair_aut": len(ordered_aut | pair_swappers), "Paley7_aut": len(paley_aut)},
        "marking_E": E,
        "design_positive": sorted(Dplus), "design_negative": sorted(Dminus),
        "derived_Fano_positive": sorted(lines_plus), "derived_Fano_negative": sorted(lines_minus),
        "circle_order8_powers": circle_powers, "punctured_antipodal_pairs": antipodes,
        "minimal_endpoint_collision": {"UU_pair": UU, "A_pair": A, "UU_triple": phi(UU), "A_triple": phi(A), "r_n_over_m": 6, "t_m_over_n": 6, "bits": [0, 1]},
        "reachable_point_bit_states": len(lifted_states),
        "Paley_H": path_count, "Paley_tc": path_count // len(paley_aut),
        "field32_plane_hostile": {"source": [0, 1, 2, 3], "image": [1, 0, 18, 14], "image_xor": 29},
        "primitive_marking_tests": marking_tests,
        "scope": "No classification beyond fields4,8,32; no twin-prime, abc, or Hamilton-gap theorem."
    }
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
