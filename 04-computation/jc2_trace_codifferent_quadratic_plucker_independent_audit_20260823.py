"""Exact bounded Plucker diagnostics for the m=2 THM-3784 W2 space.

Everything is over QQ except the explicitly labelled F_1009 quadratic-span test.
"""

import itertools
import sympy as s

x, y, t = s.symbols("x y t")
R = [x * t**r for r in range(3)]
U = t / x + t**3
P = 1 / x + s.Rational(3, 2) * t**2
V = s.expand(U * P)
E = s.expand(P * (V - 1))
gens = R + [U, V, E]
gen_names = ["R0", "R1", "R2", "U", "V", "E"]

raw, raw_names = [], []
for name, z in zip(gen_names, gens):
    raw.append(s.expand(z))
    raw_names.append(name)
for i in range(6):
    for j in range(i, 6):
        raw.append(s.expand(gens[i] * gens[j]))
        raw_names.append(f"{gen_names[i]}*{gen_names[j]}")

# Independent source-polynomial basis found by exact QQ elimination.  In this
# basis the additional relation modulo constants has V coefficient one, so
# dropping V gives a basis of W2/QQ.
pivots = list(range(12)) + [13, 14, 15, 16, 17, 19, 20, 21, 22, 23, 25, 26]
basis24 = [raw[i] for i in pivots]
names24 = [raw_names[i] for i in pivots]
keep = [i for i, name in enumerate(names24) if name != "V"]
basis = [basis24[i] for i in keep]
names = [names24[i] for i in keep]


def source(z):
    return s.expand(z.subs(t, x**2 * (1 + x * y)))


def laurent_dict(z):
    ans = {}
    for term in s.Add.make_args(s.expand(z)):
        powers = term.as_powers_dict()
        i, j = int(powers.get(x, 0)), int(powers.get(t, 0))
        ans[(i, j)] = s.cancel(ans.get((i, j), 0) + term / x**i / t**j)
    return {q: c for q, c in ans.items() if c}


def jac(f, g):
    # J_xy = x^3 J_xt because t=x^2(1+xy).
    return laurent_dict(x**3 * (s.diff(f, x) * s.diff(g, t) - s.diff(f, t) * s.diff(g, x)))


def rank_source(polys):
    ps = [s.Poly(source(z), x, y) for z in polys]
    mons = sorted(set().union(*(p.monoms() for p in ps)))
    return s.Matrix([[p.coeff_monomial(q) for p in ps] for q in mons]).rank()


print("raw_count", len(raw), "source_rank", rank_source(raw))
constant_relation = V - gens[0] * E + s.Rational(3, 2) * (gens[2] * V - gens[2])
print("constant_relation", s.expand(source(constant_relation)))
print("quotient_basis_size", len(basis), "names", names)

pairs = [(i, j) for i in range(23) for j in range(i + 1, 23)]
brackets, monomials = [], set()
for i, j in pairs:
    d = jac(basis[i], basis[j])
    brackets.append(d)
    monomials.update(d)
monomials = sorted(monomials)
B = s.polys.matrices.DomainMatrix.from_list_sympy(
    len(monomials), len(pairs),
    [[brackets[c].get(q, 0) for c in range(len(pairs))] for q in monomials],
)
nonconstant = [r for r, q in enumerate(monomials) if q != (0, 0)]
B_non = B.extract(nonconstant, list(range(len(pairs))))
print("exterior_shape", B.shape, "rank_nonconstant", B_non.rank(), "rank_full", B.rank())
print("affine_fibre_dimension", len(pairs) - B.rank())

Bm = B.to_Matrix()
unit = s.Matrix([1 if q == (0, 0) else 0 for q in monomials])
solution = next(iter(s.linsolve((Bm, unit))))
free = set().union(*(z.free_symbols for z in solution))
c0m = s.Matrix([z.subs({u: 0 for u in free}) for z in solution])
C0 = s.zeros(23)
for coefficient, (i, j) in zip(c0m, pairs):
    C0[i, j], C0[j, i] = coefficient, -coefficient
support = [(f"{names[i]}^{names[j]}", z) for z, (i, j) in zip(c0m, pairs) if z]
print("particular_rank", C0.rank(), "particular_support", support)

# Hostile fixed-x slice: J(x,G)=1 is inconsistent in W2.
dy = [s.Poly(s.diff(source(z), y), x, y) for z in basis]
dy_mons = sorted(set().union(*(p.monoms() for p in dy), {(0, 0)}))
Dy = s.Matrix([[p.coeff_monomial(q) for p in dy] for q in dy_mons])
dy_unit = s.Matrix([1 if q == (0, 0) else 0 for q in dy_mons])
print("fixed_x_rank", Dy.rank(), "fixed_x_augmented_rank", Dy.row_join(dy_unit).rank())

# Positive control: after adjoining y, x^y is rank two and maps to 1.
print("positive_control_Jxy", s.diff(x, x) * s.diff(y, y) - s.diff(x, y) * s.diff(y, x))

# Test the lowest possible affine-Pfaffian certificate.  If a linear functional
# ell on Lambda^4 made ell(c^c) a fixed nonzero constant on c0+K, it would
# separate c0^c0 from span(c0^K, K^K).  We test this exact finite-field
# relaxation over F_1009.
Kmat = B.nullspace().to_Matrix()


def two_form(row=None):
    v = Kmat.row(row) if row is not None else c0m.T
    return {pairs[j]: v[j] for j in range(len(pairs)) if v[j]}


def sign4(seq):
    inversions = sum(seq[i] > seq[j] for i in range(4) for j in range(i + 1, 4))
    return -1 if inversions % 2 else 1


def wedge(a, b):
    out = {}
    for (i, j), u in a.items():
        for (k, l), v in b.items():
            if len({i, j, k, l}) < 4:
                continue
            q = tuple(sorted((i, j, k, l)))
            out[q] = s.cancel(out.get(q, 0) + sign4((i, j, k, l)) * u * v)
    return {q: z for q, z in out.items() if z}


prime = 1009


def modq(z):
    z = s.Rational(z)
    return int(z.p) % prime * pow(int(z.q) % prime, -1, prime) % prime


def modvec(v):
    return {q: modq(z) for q, z in v.items() if modq(z)}


def insert(row, piv):
    while row:
        lead = max(row)
        value = row[lead]
        if lead not in piv:
            inv = pow(value, -1, prime)
            piv[lead] = {q: z * inv % prime for q, z in row.items()}
            return True
        old = piv[lead]
        for q, z in old.items():
            new = (row.get(q, 0) - value * z) % prime
            if new:
                row[q] = new
            elif q in row:
                del row[q]
    return False


kernel = [two_form(i) for i in range(Kmat.rows)]
c0 = two_form()
span = {}
row_count = 0
for k in kernel:
    insert(modvec(wedge(c0, k)), span)
    row_count += 1
for i in range(len(kernel)):
    for j in range(i, len(kernel)):
        insert(modvec(wedge(kernel[i], kernel[j])), span)
        row_count += 1
target = modvec(wedge(c0, c0))
independent = insert(target, span)
print(
    "quadratic_relaxation_prime", prime,
    "rows", row_count,
    "span_rank_before_target", len(span) - int(independent),
    "ambient", s.binomial(23, 4),
    "target_independent", independent,
    "target_remainder_terms", len(target),
)
