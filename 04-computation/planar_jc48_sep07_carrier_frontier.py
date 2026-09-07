#!/usr/bin/env python3
"""Exact controls for every polynomial-coefficient y-linear carrier.

No research producer is imported. The proof is analytic; this finite universe
checks its discriminant, Newton polygons, exceptional factors, and flow jet.
"""
import hashlib
import json
from math import gcd
import sympy as S

p, y, c, aa, bb, s, tau = S.symbols('p y c aa bb s tau')
gates = 0


def check(v, label):
    global gates
    gates += 1
    if not v:
        raise RuntimeError(label)


def equal(a, b, label):
    check(S.cancel(a - b) == 0, label)


def lower_hull(points):
    hull = []
    for pt in sorted(points):
        while len(hull) >= 2:
            a, b = hull[-2:]
            cross = (b[0]-a[0])*(pt[1]-b[1])-(b[1]-a[1])*(pt[0]-b[0])
            if cross > 0:
                break
            hull.pop()
        hull.append(pt)
    return hull


def newton_ram(points):
    h = lower_hull(points)
    return sum(b[0]-a[0]-gcd(b[0]-a[0], abs(b[1]-a[1]))
               for a, b in zip(h, h[1:]))


def local(a, b, zero):
    # a=None means A=0 identically. b is finite because B is nonzero.
    shift = 2 if zero else 0
    left = None if a is None else 3*a + shift
    right = 2*b
    nu = right if left is None else min(left, right)
    if left is None or left > right:
        ram = 0 if (b+shift) % 3 == 0 else 2
    elif left < right:
        ram = a % 2
    else:
        ram = 0
    pts = [(0, 0), (1, b+5 if zero else b), (3, b+shift)]
    if a is not None:
        pts.append((2, a+shift))
    return nu, ram, newton_ram(pts)


def order(f, q):
    if f == 0:
        return None
    n = 0
    while S.rem(f, q, p) == 0:
        f = S.quo(f, q, p)
        n += 1
    return n


# Literal Sylvester determinant for the cubic and its derivative.
J = p**2*(p**3-y**2)*(aa+bb*y)
f = S.Poly(J-c, y)
fc = f.all_coeffs()
gc = S.Poly(S.diff(f.as_expr(), y), y).all_coeffs()
rows = [[0]*i+fc+[0]*(1-i) for i in range(2)]
rows += [[0]*i+gc+[0]*(2-i) for i in range(3)]
disc = S.cancel(-S.det(S.Matrix(rows))/fc[0])
N = (4*p**7*(bb**2*p**3-aa**2)**2
     +4*aa*c*p**2*(9*bb**2*p**3-aa**2)-27*bb**2*c**2)
equal(disc, p**4*N, 'unnormalized cubic discriminant')
equal(S.diff(J, y).subs(y, -aa/(3*bb)),
      p**2*(aa**2+3*bb**2*p**3)/(3*bb), 'triple-root condition')
q, e, v0, v1, v2 = S.symbols('q e v0 v1 v2')
df = S.discriminant(y**2*(y+q)+e*(v0+v1*y+v2*y**2), y)
equal(S.diff(df, e).subs(e, 0), -4*q**3*v0, 'double-root variation')

# Independent Newton hulls for all local valuation pairs in this finite box.
local_rows = []
for zero in (False, True):
    for a in list(range(13))+[None]:
        for b in range(13):
            nu, ram, nr = local(a, b, zero)
            check(ram == nr, ('Newton hull', a, b, zero))
            check(nu <= 2*b, ('divisor budget', a, b, zero))
            local_rows.append([zero, a, b, nu, ram])

# Named coefficient universe. It includes coprime and repeated common
# factors, both balanced polygons, A=0, roots of B away from zero, and
# high-degree A. B is never zero.
pairs = [
    (0, 1), (1, 1), (1+p**4, 1), (0, p), (0, p**2),
    (1, p), (1+p, p**2), (1, p**4), (p**2, p**4),
    (p**2, p**3), (p**3, p**2), (p**4, p**7),
    ((p-1)**2, (p-1)**3), ((p-1)**3, (p-1)**4),
    ((p-1)**4, (p-1)**7), (0, (p-1)**3),
    (1+p, (p-1)**2*(p+2)),
    (p**2*(p-1)**2, p**4*(p-1)**3),
    (p**3*(p+1)**2, p**2*(p+1)**3),
    ((p**2+1)**2, (p**2+1)**3),
    (p**7+p+1, p**2+p+1), (p**2, p**5+1),
    (0, p**3*(p**2+1)**2),
    (p**2*(p**2+1), p**4*(p**2+1)**2),
]
records = []
for idx, (A, B) in enumerate(pairs):
    A, B = S.expand(A), S.expand(B)
    da = -1000 if A == 0 else int(S.degree(A, p))
    db = int(S.degree(B, p))
    n = max(4*db+13, 4*da+7)
    NN = S.Poly(N.subs({aa:A, bb:B}), p, c).as_expr()
    check(S.degree(NN, p) == n, ('discriminant degree', idx))
    factors = [z[0] for z in S.factor_list(B, p)[1]]
    if p not in factors:
        factors.append(p)
    F = S.Integer(1)
    nu_sum = ram_sum = 0
    local_manifest = []
    for fac in factors:
        a, b = order(A, fac), order(B, fac)
        zero = fac == p
        nu, ram, _ = local(a, b, zero)
        actual_nu = order(NN, fac)
        check(actual_nu == nu, ('fixed discriminant multiplicity', idx, fac))
        degree = int(S.degree(fac, p))
        F *= fac**nu
        nu_sum += degree*nu
        ram_sum += degree*ram
        local_manifest.append([str(fac), a, b, nu, ram])
    reduced = S.cancel(NN/F)
    check(S.denom(reduced) == 1, ('polynomial reduced divisor', idx))
    # A squarefree good specialization certifies nonzero generic resultant
    # for this named control. It is not used for the all-parameter theorem.
    good = None
    for cv in range(1, 13):
        spec = S.Poly(reduced.subs(c, cv), p)
        if spec.degree() == n-nu_sum and S.gcd(spec, spec.diff()).degree() == 0:
            if all(spec.rem(S.Poly(fac, p)).as_expr() != 0 for fac in factors):
                good = cv
                break
    check(good is not None, ('simple residual branch specialization', idx))
    # The complete infinity Newton hull always has total index defect one.
    inf_pts = [(0, -da-5 if A != 0 else 0), (1, -db-5), (3, -db-2)]
    if A != 0:
        inf_pts.append((2, -da-2))
    check(newton_ram(inf_pts) == 1, ('infinity ramification', idx))
    twice_g = n-nu_sum+ram_sum-3
    check(twice_g % 2 == 0 and twice_g >= 12, ('genus and parity', idx))
    if db == 0:
        check(twice_g//2 == max(6, 2*da+3), ('incoming constant slope', idx))
    W = s**8*(A.subs(p, s**2)+s**3*B.subs(p, s**2))
    check(W != 0, ('nonzero logarithmic coefficient', idx))
    II = (p**2*(p**3-y**2)*(A+B*y)).subs({p:s**2+tau, y:s*(s**2+tau)})
    equal(S.diff(II, tau).subs(tau, 0), W, ('actual invariant jet', idx))
    # tau-Poisson response to f(I)=I^k; only the first coefficient is needed.
    for k in (1, 2, 3):
        lead_delta_p = 2*s*k*W**k
        rhs = (tau*(2*s*S.diff(II, tau)-S.diff(II, s)))
        # delta(p) for I^k equals k I^(k-1) delta_I(p).
        lead = S.diff(rhs, tau).subs(tau, 0)*k*W**(k-1)
        equal(lead, lead_delta_p, ('outer response', idx, k))
    records.append(dict(A=str(A), B=str(B), degree=n, removed=nu_sum,
                        local_ram=ram_sum, genus=twice_g//2,
                        good_c=good, local=local_manifest))

# Sharp minimum and corrected-near-miss controls.
check(records[0]['genus'] == 6, 'constant-slope equality')
check(records[3]['genus'] == 6, 'variable-slope equality')
check(records[4]['genus'] == 8, 'old formula hostile A=0 B=p2')
check(records[5]['genus'] == 6, 'coprime B-zero equality')
check(local(2, 4, True)[:2] == (8, 0), 'balanced zero boundary')
check(local(2, 3, False)[:2] == (6, 0), 'balanced nonzero boundary')
blob = json.dumps(dict(local_rows=local_rows, controls=records), sort_keys=True,
                  separators=(',', ':')).encode()
print('Polynomial-coefficient y-linear carrier: exact controls PASS')
print('Always-active gates:', gates)
print('Local Newton cases:', len(local_rows))
print('Named coefficient pairs:', len(records))
print('Genera:', [r['genus'] for r in records])
print('Semantic SHA256:', hashlib.sha256(blob).hexdigest())
print('Scope: analytic all-A,B theorem; finite controls are not its proof.')
