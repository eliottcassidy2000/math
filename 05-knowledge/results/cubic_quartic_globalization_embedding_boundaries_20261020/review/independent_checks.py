#!/usr/bin/env python3
"""Independent finite-exact audit controls; imports no producer code."""
import hashlib
import itertools
import json
import math
from pathlib import Path

import sympy as s

x, y, u, v, a, b = s.symbols("x y u v a b")
checks = {}


def check(name, value):
    if not bool(value):
        raise RuntimeError(name)
    checks[name] = True


def zero(value):
    return s.cancel(s.expand(value)) == 0


def vector(poly, n):
    pp = s.Poly(poly, y)
    return s.Matrix([pp.nth(j) for j in range(n)])


def reduce(poly, modulus):
    return s.rem(poly, modulus, y)


def quotient_conductor_test(f, g, c, label):
    """Find A/(c), then its largest B/(c)-ideal, independently."""
    c = s.Poly(c, y).monic().as_expr()
    n = s.degree(c, y)
    columns = s.Matrix.hstack(vector(s.Integer(1), n))
    active = [s.Integer(1)]
    while active:
        pp = active.pop()
        for gg in (f, g):
            rr = reduce(pp * gg, c)
            vv = vector(rr, n)
            extended = columns.row_join(vv)
            if extended.rank() > columns.cols:
                columns = extended
                active.append(rr)
    y_matrix = s.Matrix.hstack(*[
        vector(reduce(y**(j+1), c), n) for j in range(n)
    ])
    equations = s.Matrix.vstack(*[q.T for q in columns.T.nullspace()])
    constraints = s.Matrix.vstack(*[
        equations * y_matrix**j for j in range(n)
    ])
    check(label + ":no_larger_conductor_mod_c", constraints.rank() == n)
    _, pivots = columns.T[:, ::-1].rref()
    attained = {n - 1 - i for i in pivots}
    return {"delta": int(n - columns.cols),
            "gaps": [j for j in range(n) if j not in attained],
            "conductor_kernel_dimension": int(n - constraints.rank())}


f = y**4
g = x*y**3 + y**2 + y
c = (x*y**2 + 1)*(x**2*y**4 + 2*x*y**3 + 2*y**2 + 2*y + 1)
H = (v**4 - (4*x+2)*u*v**2 - 4*u*v - 4*x**2*u**2*v
     - x**4*u**3 + (2*x**2-4*x+1)*u**2-u)

# A different elimination order: eliminate y directly, rather than form the
# producer's multiplication matrix.
check("direct_resultant_relation", zero(s.resultant(y**4-u, g-v, y)-H))
check("conductor_derivative_identity",
      zero(s.diff(H, v).subs({u:f, v:g})-s.diff(f, y)*c))

# The monogenic dual Gram matrix is polynomial over Q[x,u], not merely
# invertible over Q(x,u).  Its entry is the leading-remainder functional.
gram = s.Matrix(4, 4, lambda i, j:
                s.Poly(s.rem(v**(i+j), H, v), v).nth(3))
check("two_dimensional_base_trace_gram_determinant", gram.det() == 1)
check("trace_gram_zero_below_antidiagonal",
      all(gram[i,j] == 0 for i in range(4) for j in range(4) if i+j < 3))
check("trace_gram_unit_antidiagonal",
      all(gram[i,3-i] == 1 for i in range(4)))

# Verify conductor containment over Q[x,u] with an independent linear solve.
relation = y**4-u
powers = s.Matrix.hstack(*[
    vector(reduce(g**j, relation), 4) for j in range(4)
])
conductor_coefficients = []
for j in range(4):
    target = vector(reduce(c*y**j, relation), 4)
    sol = powers.inv()*target
    sol = sol.applyfunc(s.cancel)
    check(f"relative_conductor_containment_{j}",
          all(s.denom(t).free_symbols == set() for t in sol))
    check(f"relative_conductor_identity_{j}",
          all(zero(t) for t in powers*sol-target))
    conductor_coefficients.append([str(t) for t in sol])

# Explicit collisions cover every nonzero parameter, not only sampled fibers.
check("all_nonzero_fiber_collision_P", zero(f-f.subs(y,-y)))
check("all_nonzero_fiber_collision_Q",
      zero(s.rem(g-g.subs(y,-y), x*y**2+1, y)))
check("nonzero_collision_image_P",
      zero(s.rem(f-1/x**2, x*y**2+1, y)))
check("nonzero_collision_image_Q",
      zero(s.rem(g+1/x, x*y**2+1, y)))

fibers = {}
for aa in (0,1,-1,2,-2,s.Rational(1,2),s.Rational(-1,2)):
    record = quotient_conductor_test(f, g.subs(x,aa), c.subs(x,aa),
                                    "fiber_"+str(aa))
    check("fiber_"+str(aa)+":codimension",
          record["delta"] == (1 if aa == 0 else 3))
    check("fiber_"+str(aa)+":semigroup_gaps",
          record["gaps"] == ([1] if aa == 0 else [1,2,5]))
    fibers[str(aa)] = record

h, k = y**3, y**4+y
z = 8*h**2+12*h+3
F = v**3-u*(u+1)**3
check("generation_inverse", zero((z*z-64*k**3-9)/8-h))
G = s.resultant(y**4+y-u, 8*y**6+12*y**3+3-v, y)
check("second_kernel_primitive", s.Poly(G,u,v).content() == 1)
check("second_kernel_degree", s.degree(G,v) == 4)
check("second_equation_on_image", zero(G.subs({u:k,v:z})))
check("standard_triple_point_tangent_cone",
      zero(s.expand(F.subs(u,u-1))-(v**3+u**3-u**4)))

# Compute primitives rather than copying the producer's displayed expressions.
L1 = s.integrate(h*s.diff(k,y),y)
L2 = s.integrate(k*s.diff(z,y),y)
ell1, ell2 = reduce(L1,y**3+1), reduce(L2,y**3+1)
check("period_denominator_nonzero", ell1 == s.Rational(9,28)*y)
check("period_ratio", zero(ell2/ell1+s.Rational(28,5)))

# Test every possible abstract automorphism simultaneously.
root_identity = s.Poly((a*y+b)**3+1-a**3*(y**3+1),y)
check("normalization_translation_forced_zero",
      root_identity.nth(2) == 3*a*a*b)
check("normalization_scaling_condition",
      root_identity.nth(0).subs(b,0) == 1-a**3)
local_det = s.det(s.Matrix([[0,a],[16*u+12,0]])).subs(u,-1)
check("all_automorphisms_local_determinant", local_det == 4*a)
check("all_automorphisms_period_obstruction",
      zero(a*ell2-local_det*ell1+s.Rational(108,35)*a*y))
for first, second, label in ((h,k,"standard"),(k,z,"alternate")):
    data = quotient_conductor_test(first,second,(y**3+1)**2,label)
    check(label+":delta_three", data["delta"] == 3)
    check(label+":same_gaps", data["gaps"] == [1,2,5])

# Exhaust by recursive gcd descent, rather than the producer's Cartesian
# sequence-length loop. The initial and later bounds are the cited bounds.
allowed = [3,4,6,7]
found = []
def representable(value, generators):
    reachable = {0}
    for q in range(value+1):
        if q in reachable:
            reachable.update(q+j for j in generators if q+j <= value)
    return value in reachable

def extend(ds):
    d = math.gcd(*ds)
    if d == 1:
        if not {3,4}.issubset(ds):
            return
        previous = ds[0]
        cc = 1-ds[0]
        for j in ds[1:]:
            nextgcd = math.gcd(previous,j)
            cc += (previous//nextgcd-1)*j
            previous = nextgcd
        if cc == 6:
            found.append(ds)
        return
    ni = math.gcd(*ds[:-1])//d
    for q in (3,4):
        nextgcd = math.gcd(d,q)
        if nextgcd == d or q >= ds[-1]*ni:
            continue
        if representable(q*(d//nextgcd),ds):
            extend(ds+[q])

for first, second in itertools.combinations(allowed,2):
    if second % first:
        extend([second,first])
check("source_axioms_recursive_characteristic_patterns",
      sorted(found) == [[4,3],[6,4,3]])

root = Path(__file__).resolve().parents[1]
hashes = {p.name: hashlib.sha256(p.read_bytes()).hexdigest()
          for p in (root/"CUBIC_QUARTIC_GLOBALIZATION_AND_EMBEDDING_BOUNDARIES.md",
                    root/"verify_cubic_quartic_boundaries.py",
                    root/"PRIMARY_SOURCES.md")}
print(json.dumps({"status":"PASS", "check_count":len(checks), "checks":checks,
                  "fibers":fibers, "characteristic_patterns":sorted(found),
                  "trace_dual_gram":str(gram),
                  "relative_conductor_coefficients":conductor_coefficients,
                  "producer_hashes":hashes,
                  "sympy_version":s.__version__},
                 indent=2,sort_keys=True))
