#!/usr/bin/env python3
"""Exact bounded probe for the THM-3784 codifferent/observable lane.

Universe:
  W1(m) = span_Q{R_0,...,R_m,U,V,E}, m=2,3;
  W2    = span_Q of generator-degree 1 and 2 monomials in W1(2).

The script deliberately separates a constant in the linear image of the
exterior bracket map from a decomposable (rank-two) preimage.
"""

from __future__ import annotations

import platform
from itertools import combinations, combinations_with_replacement

import sympy as sp
from sympy.polys.matrices import DomainMatrix


x, y = sp.symbols("x y")
CHECKS = 0


def check(condition: bool, label: str) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(f"FAILED: {label}")


def poly_dict(f):
    return sp.Poly(f, x, y, domain=sp.QQ).as_dict()


def jac(f, g):
    return sp.Poly(
        sp.diff(f, x) * sp.diff(g, y) - sp.diff(f, y) * sp.diff(g, x),
        x,
        y,
        domain=sp.QQ,
    )


def dm_rank(rows):
    if not rows:
        return 0
    if not rows[0]:
        return 0
    return DomainMatrix.from_list_sympy(len(rows), len(rows[0]), rows).to_field().rank()


def source_basis(polys, names):
    dicts = [poly_dict(f) for f in polys]
    mons = sorted(set().union(*(d.keys() for d in dicts)))
    mat = DomainMatrix.from_list_sympy(
        len(mons), len(polys), [[d.get(q, 0) for d in dicts] for q in mons]
    ).to_field()
    _, pivots = mat.rref()
    return [polys[i] for i in pivots], [names[i] for i in pivots], tuple(pivots)


def bracket_rank_packet(left, right=None, wedge=False):
    cols = []
    labels = []
    monomials = set()
    if wedge:
        pairs = combinations(range(len(left)), 2)
        for i, j in pairs:
            d = jac(left[i], left[j]).as_dict()
            cols.append(d)
            labels.append((i, j))
            monomials.update(d)
    else:
        for i, f in enumerate(left):
            for j, g in enumerate(right):
                d = jac(f, g).as_dict()
                if d:
                    cols.append(d)
                    labels.append((i, j))
                    monomials.update(d)
    nonconstant = sorted(q for q in monomials if q != (0, 0))
    rows = [[d.get(q, 0) for d in cols] for q in nonconstant]
    constant_row = [[d.get((0, 0), 0) for d in cols]]
    return {
        "columns": cols,
        "labels": labels,
        "monomials": monomials,
        "rank_nonconstant": dm_rank(rows),
        "rank_with_constant": dm_rank(rows + constant_row),
    }


def tower(m):
    A = 1 + x * y
    t = x**2 * A
    B = 1 + x * t**m
    U = sp.expand(x * A * B)
    P = 1 / x + sp.Rational(m + 1, m) * t**m
    V = sp.expand(sp.cancel(U * P))
    E = sp.expand(sp.cancel(P * (V - 1)))
    check(sp.denom(V) == 1 and sp.denom(E) == 1, f"m={m} source polynomiality")
    rungs = [sp.expand(x * t**r) for r in range(m + 1)]
    return t, U, P, V, E, rungs


def axis_packet(m, names, funcs):
    z = sp.symbols("z")
    x_on_b = -z ** (-m)
    # t=z and A=t/x^2=z^(2m+1), hence y=(A-1)/x.
    y_on_b = sp.expand((z ** (2 * m + 1) - 1) / x_on_b)
    rows = []
    for name, f in zip(names, funcs):
        d = poly_dict(f)
        ord_x = min(q[0] for q in d)
        on_x = sp.factor(f.subs(x, 0))
        on_a = sp.factor(sp.cancel(f.subs(y, -1 / x)))
        on_b = sp.factor(sp.cancel(f.subs({x: x_on_b, y: y_on_b})))
        rows.append((name, ord_x, on_x, on_a, on_b))
    return rows


print(f"PYTHON={platform.python_version()}")
print(f"SYMPY={sp.__version__}")
print("FIELD=QQ (component statements interpreted over characteristic zero)")


# -------------------------------------------------------------------------
# W1 for m=2 and m=3.
# -------------------------------------------------------------------------

linear_packets = {}
for m in (2, 3):
    t, U, P, V, E, rungs = tower(m)
    funcs = rungs + [U, V, E]
    names = [f"R{r}" for r in range(m + 1)] + ["U", "V", "E"]
    packet = bracket_rank_packet(funcs, wedge=True)
    linear_packets[m] = (t, U, P, V, E, rungs, funcs, names, packet)
    possible = packet["rank_with_constant"] > packet["rank_nonconstant"]
    print(
        f"W1 m={m}: dim={len(funcs)} wedges={len(packet['columns'])} "
        f"monomials={len(packet['monomials'])} "
        f"rank_nonconstant={packet['rank_nonconstant']} "
        f"rank_with_constant={packet['rank_with_constant']} "
        f"pure_constant_possible={possible}"
    )
    check(not possible, f"m={m} W1 has no pure nonzero constant exterior image")
    print(f"W1 m={m} generator sidecars:")
    for name, ord_x, on_x, on_a, on_b in axis_packet(m, names, funcs):
        if name.startswith("R"):
            r = int(name[1:])
            trace = 0 if r < m else -m
            in_target = "NO"
        else:
            trace = f"{m + 1}*{name}"
            in_target = "YES"
        print(
            f"  {name}: target={in_target} trace={trace} ord_x={ord_x} "
            f"components[X,A,B]=[{on_x},{on_a},{on_b}]"
        )

# The only m=3 W1 bracket-kernel relation.
t3, U3, P3, V3, E3, r3, funcs3, names3, packet3 = linear_packets[3]
relation3 = jac(r3[1], r3[2]).as_expr() - sp.Rational(1, 3) * jac(r3[0], r3[3]).as_expr()
check(sp.expand(relation3) == 0, "m=3 rung bracket relation")
non3 = sorted(q for q in packet3["monomials"] if q != (0, 0))
rows3 = [[d.get(q, 0) for d in packet3["columns"]] for q in non3]
check(len(packet3["columns"]) - dm_rank(rows3) == 1, "m=3 W1 kernel is one-dimensional")
print("W1 m=3 sole kernel: J(R1,R2)-(1/3)J(R0,R3)=0")


# -------------------------------------------------------------------------
# W2 for m=2.
# -------------------------------------------------------------------------

t, U, P, V, E, rungs, G, Gnames, _ = linear_packets[2]
raw = list(G)
raw_names = list(Gnames)
for i, j in combinations_with_replacement(range(len(G)), 2):
    raw.append(sp.expand(G[i] * G[j]))
    raw_names.append(f"{Gnames[i]}*{Gnames[j]}")

W, Wnames, Wpivots = source_basis(raw, raw_names)
check(len(raw) == 27 and len(W) == 24, "m=2 W2 source dimension")
omitted = [raw_names[i] for i in range(len(raw)) if i not in Wpivots]
print(f"W2 raw=27 source_dim=24 pivot_omissions={omitted}")

# The first exact component-equalizing mixed identity.  It also proves 1 in W2.
R0, R1, R2 = rungs
constant_identity = sp.expand(V - sp.Rational(3, 2) * R2 - R0 * E + sp.Rational(3, 2) * R2 * V)
check(constant_identity == 1, "mixed W2 constant identity")
print("W2 constant identity: 1=V-(3/2)R2-R0*E+(3/2)R2*V")
print("W2 constant identity component values [X,A,B]=[1,1,1]")

# y is not in W2, whereas 1 and x are.
rank_w = len(W)
rank_w_y = len(source_basis(W + [y], Wnames + ["y"])[0])
rank_w_one = len(source_basis(W + [1], Wnames + ["1"])[0])
check(rank_w_one == rank_w and rank_w_y == rank_w + 1, "constant/y membership")
print(f"W2 membership ranks: W={rank_w}, W+<1>={rank_w_one}, W+<y>={rank_w_y}")

# Quotient constants before applying the Pluecker rank convention.
Wder_raw = [sp.expand(f - f.subs({x: 0, y: 0})) for f in W]
Q, Qnames, _ = source_basis(Wder_raw, Wnames)
check(len(Q) == 23, "derivative quotient dimension")
print(f"PLUCKER_SPACE=W2/QQ dim={len(Q)} wedges={len(Q) * (len(Q) - 1) // 2}")

# Full component/axis table for the chosen source W2 basis.
print("W2 SOURCE BASIS SIDECARS:")
component_rows = axis_packet(2, Wnames, W)


# -------------------------------------------------------------------------
# Target-field membership and field trace on W2.
# -------------------------------------------------------------------------

u0, p0, T = sp.symbols("u0 p0 T")
cover = T**3 - 2 * p0 * T + 2 * u0
g0 = p0 - sp.Rational(3, 2) * T**2
K0 = sp.QQ.frac_field(u0, p0)
xpoly = sp.invert(sp.Poly(g0, T, domain=K0), sp.Poly(cover, T, domain=K0)).as_expr()


def reduce_cover(z):
    return sp.rem(
        sp.Poly(sp.cancel(z), T, domain=K0), sp.Poly(cover, T, domain=K0)
    ).as_expr()


field_generators = [reduce_cover(xpoly * T**r) for r in range(3)] + [
    u0,
    u0 * p0,
    p0 * (u0 * p0 - 1),
]
field_raw = list(field_generators)
for i, j in combinations_with_replacement(range(len(field_generators)), 2):
    field_raw.append(reduce_cover(field_generators[i] * field_generators[j]))
field_W = [reduce_cover(field_raw[i]) for i in Wpivots]

membership_rows = []
for degree in (1, 2):
    coeffs = [sp.Poly(z, T, domain=K0).nth(degree) for z in field_W]
    common_den = sp.lcm([sp.denom(sp.cancel(c)) for c in coeffs])
    nums = [
        sp.Poly(sp.cancel(c * common_den), u0, p0, domain=sp.QQ).as_dict()
        for c in coeffs
    ]
    mons = sorted(set().union(*(d.keys() for d in nums)))
    membership_rows.extend([[d.get(q, 0) for d in nums] for q in mons])
membership_rank = dm_rank(membership_rows)
target_intersection_dim = len(W) - membership_rank
check(target_intersection_dim == 9, "W2 target-field intersection dimension")
print("W2_INTERSECT_k(U,P)=span{1,U,V,E,U^2,UV,UE,VE,E^2} dim=9")


def field_trace(z):
    power_basis = [1, T, T**2]
    cols = []
    for b0 in power_basis:
        q = sp.Poly(reduce_cover(z * b0), T, domain=K0)
        cols.append([q.nth(i) for i in range(3)])
    multiplication = sp.Matrix(3, 3, lambda i, j: cols[j][i])
    return sp.factor(sp.trace(multiplication))


for (name, ord_x, on_x, on_a, on_b), fw in zip(component_rows, field_W):
    q = sp.Poly(fw, T, domain=K0)
    in_target = q.nth(1) == 0 and q.nth(2) == 0
    tr = field_trace(fw)
    print(
        f"  {name}: target={in_target} trace={tr} ord_x={ord_x} "
        f"components[X,A,B]=[{on_x},{on_a},{on_b}]"
    )


# -------------------------------------------------------------------------
# Exterior rank tests and the Pluecker distinction.
# -------------------------------------------------------------------------

cross = bracket_rank_packet(G, Q, wedge=False)
check(
    cross["rank_nonconstant"] == cross["rank_with_constant"] == 56,
    "W1 tensor W2 has no pure constant",
)
print(
    "W1x(W2/QQ): "
    f"rank_nonconstant={cross['rank_nonconstant']} "
    f"rank_with_constant={cross['rank_with_constant']} "
    "pure_constant_possible=False"
)

full = bracket_rank_packet(Q, wedge=True)
check(
    full["rank_nonconstant"] == 111 and full["rank_with_constant"] == 112,
    "full W2 exterior ranks",
)
check(
    len(full["monomials"]) == 621
    and max(sum(q) for q in full["monomials"]) == 91,
    "full W2 coefficient universe",
)
fibre_dim = len(full["columns"]) - full["rank_with_constant"]
print(
    "Lambda2(W2/QQ): "
    f"rank_nonconstant={full['rank_nonconstant']} "
    f"rank_with_constant={full['rank_with_constant']} "
    f"coefficient_monomials={len(full['monomials'])} "
    f"max_total_degree={max(sum(q) for q in full['monomials'])} "
    f"constant_fibre_affine_dim={fibre_dim}"
)

# A sparse exact constant preimage.  Its skew rank is six, so it is not a
# decomposable bivector and therefore is not a Darboux pair.
certificate_terms = [
    (sp.Rational(1), "R0", "E"),
    (sp.Rational(9, 2), "R0", "R1*U"),
    (sp.Rational(18), "R0", "R2*E"),
    (sp.Rational(-30), "R1", "R1*E"),
    (sp.Rational(45, 8), "R1", "U*U"),
    (sp.Rational(9, 2), "R0*R1", "U*V"),
]
certificate_value = 0
skew = sp.zeros(len(Q))
for coefficient, left_name, right_name in certificate_terms:
    i = Qnames.index(left_name)
    j = Qnames.index(right_name)
    certificate_value += coefficient * jac(Q[i], Q[j]).as_expr()
    skew[i, j] += coefficient
    skew[j, i] -= coefficient
check(sp.expand(certificate_value) == 1, "constant exterior certificate")
check(skew.rank() == 6, "constant certificate skew rank")
i0, i1, i2, i3 = [Qnames.index(n) for n in ("R0", "R1", "E", "R1*E")]
pfaffian = sp.expand(
    skew[i0, i1] * skew[i2, i3]
    - skew[i0, i2] * skew[i1, i3]
    + skew[i0, i3] * skew[i1, i2]
)
check(pfaffian == 30, "hostile nonzero four-Pfaffian")
print("EXTERIOR_CERTIFICATE:")
for coefficient, left_name, right_name in certificate_terms:
    print(f"  {coefficient}*J({left_name},{right_name})")
print("  =1")
print("EXTERIOR_CERTIFICATE_SKEW_RANK=6")
print("HOSTILE_PFAFFIAN indices=[R0,R1,E,R1*E] value=30")

# THM-3779's closest target-field pair is only axis-local.
check(sp.expand(jac(U, E).as_expr() - (2 * V - 1)) == 0, "J(U,E)=2V-1")
print("AXIS_LOCAL_PAIR=(U,E): J=2V-1; component values [X,A,B]=[1,-1,-1]")
print("AXIS_LOCAL_PAIR target=[True,True] trace=[3U,3E] ord_x=[1,0]")

# -------------------------------------------------------------------------
# Cheapest normalized decomposability slice.
# -------------------------------------------------------------------------

# The linear-jet map on W2/QQ has rank two.  Its kernel H has dimension 21.
# Every genuine constant-bracket pair has a unique SL2-normalized basis
# F=x+h, G=E+k with h,k in H.  The full chart therefore has 42 variables.
jet = sp.Matrix(
    [
        [sp.Poly(f, x, y).coeff_monomial(x) for f in Q],
        [sp.Poly(f, x, y).coeff_monomial(y) for f in Q],
    ]
)
H_vectors = jet.nullspace()
H = [sp.expand(sum(v[i] * Q[i] for i in range(len(Q)))) for v in H_vectors]
check(jet.rank() == 2 and len(H) == 21, "normalized jet chart dimensions")

# Exhaust all 21^2 one-correction slices
#     F=x+a*h_i,  G=E+b*h_j.
# Linearize with c=ab, solve the coefficient equations exactly, then impose
# the Segre equation c=ab.  Over an algebraically closed field an affine
# polynomial in the remaining free parameters fails to have a zero only
# when it is a nonzero constant.
q00 = jac(x, E).as_dict()
q00[(0, 0)] = q00.get((0, 0), 0) - 1
one_correction_survivors = 0
for hi in H:
    q10 = jac(hi, E).as_dict()
    for hj in H:
        four = [q00, q10, jac(x, hj).as_dict(), jac(hi, hj).as_dict()]
        mons = sorted(set().union(*(d.keys() for d in four)))
        mat = sp.Matrix([[d.get(q, 0) for d in four] for q in mons])
        rhs = -mat[:, 0]
        coeff = mat[:, 1:]
        if coeff.rank() != coeff.row_join(rhs).rank():
            continue
        solution_set = sp.linsolve((coeff, rhs))
        if solution_set is sp.EmptySet:
            continue
        a_value, b_value, c_value = next(iter(solution_set))
        segre = sp.expand(c_value - a_value * b_value)
        free = set().union(
            a_value.free_symbols, b_value.free_symbols, c_value.free_symbols
        )
        # If there are parameters and segre is nonconstant, it has a complex
        # zero.  If it is identically zero, every point works.  Only a forced
        # nonzero constant is impossible.
        if segre == 0 or (free and segre.free_symbols):
            one_correction_survivors += 1
check(one_correction_survivors == 0, "all one-correction Pluecker slices fail")
print(
    "NORMALIZED_PLUCKER_CHART: variables=42 "
    f"four_pfaffians={sp.binomial(len(Q), 4)}"
)
print("ONE_CORRECTION_SLICES=441 survivors=0")

# Positive control: adjoining the missing source coordinate y adds one
# dimension and gives the decomposable rank-two bivector x wedge y.
Qplus, Qplus_names, _ = source_basis(Q + [y], Qnames + ["y"])
check(len(Qplus) == 24, "positive-control enlargement")
positive_skew = sp.zeros(len(Qplus))
ix = Qplus_names.index("R0")
iy = Qplus_names.index("y")
positive_skew[ix, iy] = 1
positive_skew[iy, ix] = -1
check(jac(Qplus[ix], Qplus[iy]).as_expr() == 1, "positive Darboux control")
check(positive_skew.rank() == 2, "positive control rank two")
print("POSITIVE_CONTROL=adjoin_y: J(R0,y)=1 skew_rank=2")

print("FULL_W2_RANK2_INTERSECTION=OPEN")
print(
    "FIRST_UNRESOLVED_CASE=both legs require genuinely quadratic content; "
    "the exact affine constant fibre has dimension 141, while the normalized "
    "rank-two chart has 42 variables."
)
print(f"CHECKS={CHECKS}")
