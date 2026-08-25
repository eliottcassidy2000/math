#!/usr/bin/env python3
"""Characteristic-zero exact affine retained-packet audit over Q(alpha)."""

from hashlib import sha256
from math import comb

from flint import fmpq, fmpq_mat
import sympy as sp


Q0 = fmpq(0)
Q1 = fmpq(1)
REL = (
    fmpq(1276420, 72783360),
    fmpq(-7849770, 72783360),
    fmpq(28419741, 72783360),
    fmpq(77822208, 72783360),
)
K0 = (Q0, Q0, Q0, Q0)
K1 = (Q1, Q0, Q0, Q0)
A = (Q0, Q1, Q0, Q0)


def q(a, b=1):
    return fmpq(a, b)


def kc(a, b=1):
    return (q(a, b), Q0, Q0, Q0)


def ka(a, b):
    return tuple(x + y for x, y in zip(a, b))


def kn(a):
    return tuple(-x for x in a)


def ks(a, b):
    return tuple(x - y for x, y in zip(a, b))


def km(a, b):
    raw = [Q0] * 7
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if y:
                    raw[i + j] += x * y
    for degree in range(6, 3, -1):
        coefficient = raw[degree]
        if coefficient:
            for offset, relation in enumerate(REL):
                raw[degree - 4 + offset] += coefficient * relation
    return tuple(raw[:4])


def kq(a, scalar):
    return tuple(x * scalar for x in a)


def kinv(a):
    basis = (K1, A, km(A, A), km(km(A, A), A))
    columns = [km(a, value) for value in basis]
    matrix = fmpq_mat(4, 4, [columns[c][r] for r in range(4) for c in range(4)])
    solution = matrix.solve(fmpq_mat(4, 1, [Q1, Q0, Q0, Q0]))
    value = tuple(solution[i, 0] for i in range(4))
    if km(a, value) != K1:
        raise RuntimeError("bad inverse")
    return value


def kp(a, n):
    value = K1
    while n:
        if n & 1:
            value = km(value, a)
        a = km(a, a)
        n //= 2
    return value


def padd(left, right):
    result = dict(left)
    for key, value in right.items():
        result[key] = ka(result.get(key, K0), value)
        if result[key] == K0:
            del result[key]
    return result


def pscale(value, scalar):
    if not isinstance(scalar, tuple):
        scalar = kc(scalar)
    return {key: product for key, coefficient in value.items()
            if (product := km(coefficient, scalar)) != K0}


def pmul(left, right, cutoff):
    result = {}
    for (i, j), a in left.items():
        for (u, v), b in right.items():
            if i + j + u + v <= cutoff:
                key = (i + u, j + v)
                result[key] = ka(result.get(key, K0), km(a, b))
    return {key: value for key, value in result.items() if value != K0}


def ppow(value, n, cutoff):
    result = {(0, 0): K1}
    for _ in range(n):
        result = pmul(result, value, cutoff)
    return result


def pdiff(value, axis):
    result = {}
    for (i, j), coefficient in value.items():
        exponent = i if axis == 0 else j
        if exponent:
            key = (i - (axis == 0), j - (axis == 1))
            result[key] = kq(coefficient, exponent)
    return result


def uadd(left, right):
    result = dict(left)
    for degree, value in right.items():
        result[degree] = ka(result.get(degree, K0), value)
        if result[degree] == K0:
            del result[degree]
    return result


def uscale(value, scalar):
    if not isinstance(scalar, tuple):
        scalar = kc(scalar)
    return {degree: product for degree, coefficient in value.items()
            if (product := km(coefficient, scalar)) != K0}


def ushift(value, n):
    return {degree + n: coefficient for degree, coefficient in value.items()}


def around(value, point, cutoff):
    result = {}
    for degree, coefficient in value.items():
        for h in range(min(degree, cutoff) + 1):
            term = kq(coefficient, comb(degree, h) * point ** (degree - h))
            result[(h, 0)] = ka(result.get((h, 0), K0), term)
    return {key: coefficient for key, coefficient in result.items() if coefficient != K0}


Qbase = {0: kc(-3, 4), 1: K1, 2: kc(-27, 4), 3: kc(-2), 4: kc(9, 2), 5: K1}
Pbase = {2: K1, 4: kc(-2), 6: K1}
Q6 = uadd(Qbase, uscale(Pbase, kc(-259, 36)))
R1 = uadd(Pbase, uscale(ushift(Pbase, 2), kc(-1)))
R2 = uadd(uscale(Pbase, kc(4)), uscale(ushift(Pbase, 1), kc(-9)))
parabola = ka(ka(kq(km(A, A), q(520, 9)), kq(A, q(-1688, 81))), kc(-5717, 729))
Qalpha = uadd(Q6, uadd(uscale(R1, parabola), uscale(R2, A)))

POINTS = (-1, 0, 1)
CUT = 5
JET = CUT + 1
MIXED_QUADRATIC = 0
MIXED_POWER = 2
ROWS = tuple((stable, point, source)
             for stable in range(CUT + 1)
             for source in range(CUT + 1 - stable)
             for point in POINTS)
TERMINAL = tuple(i for i, (stable, _, _) in enumerate(ROWS) if stable == CUT)


def build():
    packets = {}
    for point in POINTS:
        vertical = {(0, 1): K1}
        if MIXED_QUADRATIC:
            vertical[(0, MIXED_POWER)] = kc(MIXED_QUADRATIC)
        qjet = padd(around(Qalpha, point, JET), vertical)
        xjet = {(0, 0): kc(point), (1, 0): K1}
        D = padd({(0, 0): K1}, pmul(pmul(xjet, xjet, JET), qjet, JET))
        y = pscale(pmul(pmul(xjet, D, JET), padd(D, {(0, 0): kc(2)}), JET), kc(1, 3))
        z = padd(pmul(qjet, padd(D, {(0, 0): kc(3)}), JET), {(0, 0): kc(3)})
        if y.get((0, 0), K0) != K0 or z.get((0, 0), K0) != K0:
            raise RuntimeError(("retained origin", point))
        yh, yt, zh, zt = pdiff(y, 0), pdiff(y, 1), pdiff(z, 0), pdiff(z, 1)
        area = padd(pmul(yh, zt, CUT), pscale(pmul(yt, zh, CUT), kc(-1)))
        packets[point] = ([ppow(y, n, CUT) for n in range(CUT + 1)],
                          [ppow(z, n, CUT) for n in range(CUT + 1)],
                          (area, yh, zh))
    columns = []
    labels = []
    for kind in range(3):
        for a in range(CUT + 1):
            for b in range(CUT + 1 - a):
                for c in range(CUT + 1 - a - b):
                    values = {}
                    for point in POINTS:
                        yp, zp, bases = packets[point]
                        value = pmul(yp[a], zp[b], CUT)
                        value = pmul(value, {(0, c): K1}, CUT)
                        values[point] = pmul(value, bases[kind], CUT)
                    columns.append([values[point].get((source, stable), K0)
                                    for stable, point, source in ROWS])
                    labels.append((kind, a, b, c))
    return columns, labels


def pseries_inverse(value, cutoff):
    constant = value.get((0, 0), K0)
    if constant == K0:
        raise RuntimeError("nonunit series")
    inverse_constant = kinv(constant)
    remainder = padd(pscale(value, inverse_constant), {(0, 0): kc(-1)})
    answer = {}
    term = {(0, 0): K1}
    for degree in range(cutoff + 1):
        answer = padd(answer, pscale(term, kc(-1 if degree % 2 else 1)))
        term = pmul(term, remainder, cutoff)
    return pscale(answer, inverse_constant)


def build_fixed_a_columns():
    packets = {}
    local_cut = JET + 1
    for point in POINTS:
        vertical = {(0, 1): K1}
        if MIXED_QUADRATIC:
            vertical[(0, MIXED_POWER)] = kc(MIXED_QUADRATIC)
        qjet = padd(around(Qalpha, point, local_cut), vertical)
        xjet = {(0, 0): kc(point), (1, 0): K1}
        D = padd({(0, 0): K1}, pmul(pmul(xjet, xjet, local_cut), qjet, local_cut))
        Dinv = pseries_inverse(D, local_cut)
        ajet = pmul(qjet, pmul(Dinv, Dinv, local_cut), local_cut)
        Alocal = padd(ajet, {(0, 0): kc(3, 4)})
        cjet = pmul(pmul(xjet, D, local_cut), padd(D, {(0, 0): kc(2)}), local_cut)
        packets[point] = (
            [ppow(Alocal, n, CUT + 1) for n in range(CUT + 2)],
            [ppow(cjet, n, CUT + 1) for n in range(CUT + 2)],
            pdiff(ajet, 0),
        )
    hprime = {(0, 0): K1}
    if MIXED_QUADRATIC:
        hprime[(0, MIXED_POWER - 1)] = kc(MIXED_POWER * MIXED_QUADRATIC)
    result = []
    for a in range(CUT + 2):
        for b in range(CUT + 2 - a):
            for c in range(CUT + 2 - a - b):
                if a + b + c == 0:
                    continue
                values = {}
                for point in POINTS:
                    ap, cp, ax = packets[point]
                    density = {}
                    if b:
                        gc = pmul(ap[a], cp[b - 1], CUT)
                        gc = pmul(gc, {(0, c): K1}, CUT)
                        density = padd(density, pscale(pmul(hprime, gc, CUT), kc(-3 * b)))
                    if c:
                        gw = pmul(ap[a], cp[b], CUT)
                        gw = pmul(gw, {(0, c - 1): K1}, CUT)
                        density = padd(density, pscale(pmul(ax, gw, CUT), kc(c)))
                    values[point] = density
                result.append([values[point].get((source, stable), K0)
                               for stable, point, source in ROWS])
    return result


def build_linearized_pair_columns():
    """Tangent image of (a+eps f,-4c+eps g) on the affine H=t source."""
    packets = {}
    local_cut = JET + 1
    for point in POINTS:
        qjet = padd(around(Qalpha, point, local_cut), {(0, 1): K1})
        xjet = {(0, 0): kc(point), (1, 0): K1}
        D = padd({(0, 0): K1}, pmul(pmul(xjet, xjet, local_cut), qjet, local_cut))
        Dinv = pseries_inverse(D, local_cut)
        ajet = pmul(qjet, pmul(Dinv, Dinv, local_cut), local_cut)
        Alocal = padd(ajet, {(0, 0): kc(3, 4)})
        cjet = pmul(pmul(xjet, D, local_cut), padd(D, {(0, 0): kc(2)}), local_cut)
        packets[point] = (
            [ppow(Alocal, n, CUT + 1) for n in range(CUT + 2)],
            [ppow(cjet, n, CUT + 1) for n in range(CUT + 2)],
            pdiff(ajet, 0), pdiff(cjet, 0),
        )
    result = []
    for coordinate in ("f", "g"):
        for a in range(CUT + 2):
            for b in range(CUT + 2 - a):
                for c in range(CUT + 2 - a - b):
                    if a + b + c == 0:
                        continue
                    values = {}
                    for point in POINTS:
                        ap, cp, ax, cx = packets[point]
                        density = {}
                        if coordinate == "f":
                            if a:
                                fa = pmul(ap[a - 1], cp[b], CUT)
                                fa = pmul(fa, {(0, c): K1}, CUT)
                                density = padd(density, pscale(fa, kc(12 * a)))
                            if c:
                                fw = pmul(ap[a], cp[b], CUT)
                                fw = pmul(fw, {(0, c - 1): K1}, CUT)
                                density = padd(density, pscale(pmul(cx, fw, CUT), kc(4 * c)))
                        else:
                            if b:
                                gc = pmul(ap[a], cp[b - 1], CUT)
                                gc = pmul(gc, {(0, c): K1}, CUT)
                                density = padd(density, pscale(gc, kc(-3 * b)))
                            if c:
                                gw = pmul(ap[a], cp[b], CUT)
                                gw = pmul(gw, {(0, c - 1): K1}, CUT)
                                density = padd(density, pscale(pmul(ax, gw, CUT), kc(c)))
                        values[point] = density
                    result.append([values[point].get((source, stable), K0)
                                   for stable, point, source in ROWS])
    return result


def rref_nullspace(matrix):
    rows = [list(row) for row in matrix]
    pivots = []
    row = 0
    for column in range(len(rows[0])):
        selected = next((i for i in range(row, len(rows)) if rows[i][column] != K0), None)
        if selected is None:
            continue
        rows[row], rows[selected] = rows[selected], rows[row]
        inverse = kinv(rows[row][column])
        rows[row] = [km(value, inverse) for value in rows[row]]
        for i in range(len(rows)):
            if i != row and rows[i][column] != K0:
                coefficient = rows[i][column]
                rows[i] = [ks(value, km(coefficient, pivot))
                           for value, pivot in zip(rows[i], rows[row])]
        pivots.append(column)
        row += 1
        if row == len(rows):
            break
    free = [column for column in range(len(rows[0])) if column not in pivots]
    basis = []
    for column in free:
        vector = [K0] * len(rows[0])
        vector[column] = K1
        for i, pivot in enumerate(pivots):
            vector[pivot] = kn(rows[i][column])
        basis.append(vector)
    return pivots, basis


def dot(left, right):
    value = K0
    for a, b in zip(left, right):
        if a != K0 and b != K0:
            value = ka(value, km(a, b))
    return value


def build_closed_columns(columns, labels):
    lookup = {label: columns[index] for index, label in enumerate(labels)}
    result = []

    def combine(terms):
        vector = [K0] * len(ROWS)
        for scalar, label in terms:
            vector = [ka(value, kq(entry, scalar))
                      for value, entry in zip(vector, lookup[label])]
        if any(value != K0 for value in vector):
            result.append(vector)

    for kind in range(3):
        for a in range(CUT + 2):
            for b in range(CUT + 2 - a):
                for c in range(CUT + 2 - a - b):
                    terms = []
                    if kind == 0:
                        if b:
                            terms.append((-b, (0, a, b - 1, c)))
                        if c:
                            terms.append((-c, (1, a, b, c - 1)))
                    elif kind == 1:
                        if a:
                            terms.append((a, (0, a - 1, b, c)))
                        if c:
                            terms.append((-c, (2, a, b, c - 1)))
                    else:
                        if a:
                            terms.append((a, (1, a - 1, b, c)))
                        if b:
                            terms.append((b, (2, a, b - 1, c)))
                    if terms:
                        combine(terms)
    return result


def configure(cut):
    global CUT, JET, ROWS, TERMINAL
    CUT = cut
    JET = cut + 1
    ROWS = tuple((stable, point, source)
                 for stable in range(cut + 1)
                 for source in range(cut + 1 - stable)
                 for point in POINTS)
    TERMINAL = tuple(index for index, (stable, _, _) in enumerate(ROWS)
                     if stable == cut)


def constant_vector():
    return [K1 if stable == 0 and source == 0 else K0
            for stable, _, source in ROWS]


def relation_serialization(relations):
    return "\n".join(
        ";".join(
            f"{stable},{point},{source}:"
            + ",".join(str(coordinate) for coordinate in relation[index])
            for index, (stable, point, source) in enumerate(ROWS)
            if relation[index] != K0
        )
        for relation in relations
    )


def field_norm(value):
    basis = (K1, A, km(A, A), km(km(A, A), A))
    columns = [km(value, element) for element in basis]
    matrix = fmpq_mat(4, 4,
                      [columns[column][row]
                       for row in range(4) for column in range(4)])
    return matrix.det()


def uevaluate(value, point, derivative=False):
    answer = K0
    for degree, coefficient in value.items():
        if derivative:
            if degree:
                answer = ka(answer, kq(coefficient, degree * point ** (degree - 1)))
        else:
            answer = ka(answer, kq(coefficient, point ** degree))
    return answer


def run_form_packet(mixed):
    global MIXED_QUADRATIC
    MIXED_QUADRATIC = mixed
    columns, labels = build()
    pivots, relations = rref_nullspace(columns)
    closed_columns = build_closed_columns(columns, labels)
    closed_pivots, closed_relations = rref_nullspace(closed_columns)
    constant = constant_vector()
    expected_rank = len(ROWS) - 4
    if len(columns) != 168 or len(closed_columns) != 231:
        raise RuntimeError(("column census", len(columns), len(closed_columns)))
    if len(pivots) != expected_rank or len(relations) != 4:
        raise RuntimeError(("rank/nullity", len(pivots), len(relations)))
    if len(closed_pivots) != expected_rank or len(closed_relations) != 4:
        raise RuntimeError(("exact rank/nullity", len(closed_pivots),
                            len(closed_relations)))
    for bank, bank_relations in ((columns, relations),
                                 (closed_columns, closed_relations)):
        if not all(dot(column, relation) == K0
                   for column in bank for relation in bank_relations):
            raise RuntimeError("relation residual")
        if not all(dot(constant, relation) == K0
                   for relation in bank_relations):
            raise RuntimeError("constant outside image")
        if not all(all(relation[index] == K0 for index in TERMINAL)
                   for relation in bank_relations):
            raise RuntimeError("terminal projection nonzero")
    serialization = relation_serialization(relations)
    return {
        "columns": columns,
        "rank": len(pivots),
        "relations": relations,
        "closed_columns": closed_columns,
        "closed_rank": len(closed_pivots),
        "hash": sha256(serialization.encode()).hexdigest(),
        "counts": tuple(sum(value != K0 for value in relation)
                        for relation in relations),
    }


configure(5)
q_values = tuple(uevaluate(Qalpha, point) for point in POINTS)
q_derivatives = tuple(uevaluate(Qalpha, point, True) for point in POINTS)
if q_values != (kc(-3), kc(-3, 4), kc(-3)):
    raise RuntimeError(("exceptional values", q_values))
if q_derivatives != (kc(-9, 2), K1, kc(9, 2)):
    raise RuntimeError(("exceptional derivatives", q_derivatives))

x_symbol, q_symbol = sp.symbols("x q")
d_symbol = 1 + x_symbol**2 * q_symbol
a_symbol = q_symbol / d_symbol**2
b_symbol = (d_symbol - 1) * (d_symbol + 2)**2
c_symbol = x_symbol * d_symbol * (d_symbol + 2)
e_symbol = q_symbol * (d_symbol + 3)
jacobian_ac = sp.factor(sp.diff(a_symbol, x_symbol) * sp.diff(c_symbol, q_symbol)
                        - sp.diff(a_symbol, q_symbol) * sp.diff(c_symbol, x_symbol))
if sp.factor(e_symbol / (b_symbol + 4) - a_symbol) != 0 or jacobian_ac != -3:
    raise RuntimeError(("local Darboux identities", jacobian_ac))

affine = run_form_packet(0)
mixed = run_form_packet(1)
if affine["hash"] != "bb0a2ba7f46b6bf721c875a1a7a19000e1b3b8845ba11bc4707bd299b9740852":
    raise RuntimeError(("affine relation hash", affine["hash"]))
if mixed["hash"] != "bbaf81804c69200a9214711f4f3865acb39831f22af02d623436e788dd6b633e":
    raise RuntimeError(("mixed relation hash", mixed["hash"]))
if affine["counts"] != (3, 7, 14, 23) or mixed["counts"] != (3, 7, 14, 23):
    raise RuntimeError(("relation support counts", affine["counts"], mixed["counts"]))
combined_pivots, _ = rref_nullspace(affine["relations"] + mixed["relations"])
if len(combined_pivots) != 6:
    raise RuntimeError(("combined relation rank", len(combined_pivots)))

MIXED_QUADRATIC = 1
fixed_a_columns = build_fixed_a_columns()
fixed_a_pivots, fixed_a_relations = rref_nullspace(fixed_a_columns)
constant = constant_vector()
fixed_responses = [dot(constant, relation) for relation in fixed_a_relations]
fixed_nonzero = [value for value in fixed_responses if value != K0]
expected_response = (
    q(-2073506706944, 1678822119),
    q(372679949312, 62178597),
    q(-184159683584, 6908733),
    q(-73442787328, 2302911),
)
expected_norm = q(28278059768285603108255604733711150481408,
                  392425657272606224710564875)
if (len(fixed_a_columns), len(fixed_a_pivots), len(fixed_a_relations)) != (83, 57, 6):
    raise RuntimeError("fixed-a census")
if not all(dot(column, relation) == K0
           for column in fixed_a_columns for relation in mixed["relations"]):
    raise RuntimeError("fixed-a bank escaped the mixed two-form image")
if fixed_nonzero != [expected_response]:
    raise RuntimeError(("fixed-a response", fixed_nonzero))
if field_norm(fixed_nonzero[0]) != expected_norm:
    raise RuntimeError("fixed-a response norm")

linearized_columns = build_linearized_pair_columns()
linearized_pivots, linearized_relations = rref_nullspace(linearized_columns)
deformation = [kc(-24) if stable == 1 and source == 0 else K0
               for stable, _, source in ROWS]
linearized_nonzero = [dot(deformation, relation)
                      for relation in linearized_relations
                      if dot(deformation, relation) != K0]
if (len(linearized_columns), len(linearized_pivots),
        len(linearized_relations), len(linearized_nonzero)) != (166, 59, 4, 0):
    raise RuntimeError("linearized-pair gate")
if not all(dot(column, relation) == K0
           for column in linearized_columns for relation in affine["relations"]):
    raise RuntimeError("linearized pair image differs from affine form image")

configure(4)
MIXED_QUADRATIC = 1
fixed_a_cut4 = build_fixed_a_columns()
fixed_a_cut4_pivots, fixed_a_cut4_relations = rref_nullspace(fixed_a_cut4)
fixed_a_cut4_nonzero = [dot(constant_vector(), relation)
                        for relation in fixed_a_cut4_relations
                        if dot(constant_vector(), relation) != K0]
if (len(fixed_a_cut4), len(fixed_a_cut4_pivots),
        len(fixed_a_cut4_relations)) != (55, 40, 5):
    raise RuntimeError("fixed-a cutoff-four census")
if fixed_a_cut4_nonzero:
    raise RuntimeError(("fixed-a cutoff-four response", fixed_a_cut4_nonzero))
configure(5)

print("EXCEPTIONAL AFFINE/SIMPLE-ZERO PACKET -- EXACT OVER Q(alpha)")
print("Q_values=(-3,-3/4,-3)")
print("Q_derivatives=(-9/2,1,9/2);Q_prime_0=1;Q_is_not_even=True")
print("local_identity=a=e/(b+4)=q/D^2;Jac_xq(a,c)=-3")
print("affine_local_pair=(a,-4c);source_Jacobian_for_H=t=12")
print("cut5=rows63;packet=J_n through x-order 5-n for n=0..5")
print(f"H=t:arbitrary_columns=168 rank={affine['rank']} relation_dim=4")
print(f"H=t:polynomial_exact_columns=231 rank={affine['closed_rank']}")
print("H=t:exact_image_equals_arbitrary=True;constant_in_image=True;terminal_J5_projection=0")
print("H=t:relation_sha256=" + affine["hash"])
print(f"H=t+t^2:arbitrary_columns=168 rank={mixed['rank']} relation_dim=4")
print(f"H=t+t^2:polynomial_exact_columns=231 rank={mixed['closed_rank']}")
print("H=t+t^2:exact_image_equals_arbitrary=True;constant_in_image=True;terminal_J5_projection=0")
print("H=t+t^2:relation_sha256=" + mixed["hash"])
print("affine_vs_mixed:combined_relation_rank=6;images_equal=False")
print("affine_mixed-pair_tangent_image_equals_affine_arbitrary_image=True")
print("relation_nonzero_counts_for_both=" + str(affine["counts"]))
print("fixed_F=a:H=t+t^2:cut4_columns=55 rank=40 relation_dim=5;constant_in_span=True")
print("fixed_F=a:H=t+t^2:cut5_columns=83 rank=57 relation_dim=6;constant_in_span=False")
print("fixed_F=a:nonzero_constant_response=" + str(fixed_nonzero[0]))
print("fixed_F=a:response_norm=" + str(expected_norm))
print("mixed_first_order:columns=166 rank=59 relation_dim=4;required_-24t_in_span=True")
print("scope=retained_three-multigerm_cutoff_five_and_dual_numbers_only")
print("RESULT=PASS")
