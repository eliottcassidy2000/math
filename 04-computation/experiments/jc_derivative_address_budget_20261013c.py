"""Exact source-polynomial derivative frontiers, not point-jet truncations.

Finite declared controls; a failed search budget is inconclusive.  No JC2 claim.
"""
from dataclasses import dataclass, replace
from itertools import product
from math import comb, factorial
import sympy as s
import jc_bilinear_reconstruction_20261013a as prior

CHECKS = 0


def need(ok, message):
    if not ok:
        raise ValueError(message)


def check(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(message)


def natural(value, minimum=0):
    need(type(value) is int and value >= minimum, "exact integer outside domain")


def address(alpha, n):
    need(type(alpha) is tuple and len(alpha) == n, "complete address tuple")
    for a in alpha:
        natural(a)


def predecessors(alpha):
    for j, a in enumerate(alpha):
        if a:
            yield alpha[:j] + (a-1,) + alpha[j+1:]


def order_ideal(rows, n):
    need(type(rows) is tuple and rows, "nonempty finite address tuple")
    for alpha in rows:
        address(alpha, n)
    need(rows == tuple(sorted(set(rows))), "canonical distinct addresses")
    points = set(rows)
    need((0,)*n in points, "retain zero address")
    need(all(beta in points for a in rows for beta in predecessors(a)),
         "address set is not downward closed")
    return points


def boundary(rows, n):
    points = order_ideal(rows, n)
    candidates = {a[:j]+(a[j]+1,)+a[j+1:]
                  for a in rows for j in range(n)} - points
    return tuple(sorted(a for a in candidates
                        if all(b in points for b in predecessors(a))))


def derivative_reader(v, f):
    """Per-source memoization; public callers validate before entering the cache."""
    f, dual = prior.carrier(v, f)
    cache = {(0,)*len(v): v}
    def read(alpha):
        address(alpha, len(v))
        if alpha not in cache:
            j = next(j for j, a in enumerate(alpha) if a)
            beta = alpha[:j]+(alpha[j]-1,)+alpha[j+1:]
            cache[alpha] = tuple(prior.derive(v, dual, j, h) for h in read(beta))
        return cache[alpha]
    return f, dual, read


@dataclass(frozen=True)
class Frontier:
    variables: tuple
    source_map: tuple
    addresses: tuple
    anchor: tuple


def audit(receipt):
    need(type(receipt) is Frontier, "typed exact frontier receipt")
    v = prior.variables(receipt.variables)
    rows = receipt.addresses
    order_ideal(rows, len(v))
    need(type(receipt.anchor) is tuple and len(receipt.anchor) == len(v),
         "complete retained source point")
    anchor = tuple(prior.rational(a) for a in receipt.anchor)
    f, dual, read = derivative_reader(v, receipt.source_map)
    for alpha in boundary(rows, len(v)):
        need(all(h == 0 for h in read(alpha)), "nonzero source-polynomial boundary")
    return v, f, anchor, read


def reconstruct(receipt):
    v, f, a, read = audit(receipt)
    at = dict(zip(v, a))
    base = tuple(h.subs(at) for h in f)
    out = [s.Integer(0)]*len(v)
    for alpha in receipt.addresses:
        factor = s.prod((v[j]-base[j])**alpha[j]/factorial(alpha[j])
                        for j in range(len(v)))
        for i, h in enumerate(read(alpha)):
            out[i] += h.subs(at)*factor
    out = tuple(s.expand(h) for h in out)
    need(tuple(s.expand(h.subs(dict(zip(v, out)), simultaneous=True)) for h in f) == v,
         "right inverse identity")
    need(tuple(s.expand(h.subs(dict(zip(v, f)), simultaneous=True)) for h in out) == v,
         "left inverse identity")
    return out


def compile_frontier(v, f, anchor=None, max_order=12):
    """Search the literal derivative support; None means only cap exhaustion."""
    v = prior.variables(v)
    natural(max_order)
    if anchor is None:
        anchor = (0,)*len(v)
    need(type(anchor) is tuple and len(anchor) == len(v), "complete source anchor")
    tuple(prior.rational(a) for a in anchor)
    f, dual, read = derivative_reader(v, f)
    waiting, seen, live = {(0,)*len(v)}, set(), set()
    while waiting:
        alpha = min(waiting, key=lambda a: (sum(a), a))
        waiting.remove(alpha)
        seen.add(alpha)
        if all(h == 0 for h in read(alpha)):
            continue
        if sum(alpha) > max_order:
            return None
        live.add(alpha)
        for j in range(len(v)):
            beta = alpha[:j]+(alpha[j]+1,)+alpha[j+1:]
            if beta not in seen:
                waiting.add(beta)
    receipt = Frontier(v, f, tuple(sorted(live)), anchor)
    audit(receipt)
    return receipt


def weighted_monomials(weights, bound):
    need(type(weights) is tuple and weights, "positive source weight tuple")
    for w in weights:
        natural(w, 1)
    natural(bound)
    return tuple(a for a in product(*(range(bound//w+1) for w in weights))
                 if sum(w*b for w,b in zip(weights,a)) <= bound)


def weight_budget(v, f, weights):
    """Sufficient support inequalities. Failure says nothing about invertibility."""
    v = prior.variables(v)
    need(type(weights) is tuple and len(weights) == len(v), "complete source weights")
    for w in weights:
        natural(w, 1)
    f, dual = prior.carrier(v, f)
    for j in range(len(v)):
        for i in range(len(v)):
            for alpha, coefficient in s.Poly(dual[j,i], *v).terms():
                if coefficient:
                    need(sum(a*w for a,w in zip(alpha,weights)) <= weights[i],
                         "dual increases the retained weighted filtration")
    size = len(weighted_monomials(weights, max(weights)))
    return size, size-1


def degree_floor(n, r):
    """Least degree that can hold r+1 independent n-variable polynomials."""
    natural(n, 1)
    natural(r)
    d = 0
    while comb(n+d, n) < r+1:
        d += 1
    return d


def strand(v, f, j, h, depth):
    v = prior.variables(v)
    natural(j)
    need(j < len(v), "dual index outside source")
    natural(depth)
    h = prior.polynomial(h, v)
    f, dual = prior.carrier(v, f)
    values = [h]
    for _ in range(depth):
        values.append(prior.derive(v, dual, j, values[-1]))
    return tuple(values)


def main():
    x,y,z = s.symbols('x y z')
    v = (x,y)
    examples = []
    # A three-shear family, with an independently written inverse, not the old nine-map census.
    for a,b,c in ((2,2,2),(2,3,1),(3,1,2),(1,3,3)):
        u = x+y**a
        w = y+u**b
        f = (s.expand(u+w**c), s.expand(w))
        back_u = x-y**c
        back_y = y-back_u**b
        expected = (s.expand(back_u-back_y**a), s.expand(back_y))
        for anchor in ((0,0),(1,-1)):
            receipt = compile_frontier(v,f,anchor,max_order=a*b*c+1)
            check(receipt is not None, "three-shear frontier closes")
            inv = reconstruct(receipt)
            check(inv == expected, "independent ordered inverse")
            max_order = max(map(sum,receipt.addresses))
            check(max_order == max(s.Poly(h,*v).total_degree() for h in inv),
                  "exact joint termination order equals inverse degree")
            for point in product(range(-1,2),repeat=2):
                image = tuple(h.subs(dict(zip(v,point))) for h in f)
                check(tuple(h.subs(dict(zip(v,image)),simultaneous=True) for h in inv) == point,
                      "marked source round trip")
            if anchor == (0,0):
                examples.append((a,b,c,len(receipt.addresses),len(boundary(receipt.addresses,2)),max_order))
    # Minimal excluded addresses are essential: a rectangle is not a total-degree simplex.
    f = (x-y*y,y)
    triangle = compile_frontier(v,f,max_order=2)
    check(triangle.addresses == ((0,0),(0,1),(0,2),(1,0)), "anisotropic staircase")
    check(boundary(triangle.addresses,2) == ((0,3),(1,1),(2,0)), "minimal boundary")
    oversized = replace(triangle,addresses=tuple(product(range(2),range(3))))
    check(reconstruct(oversized) == (x+y*y,y), "nonminimal rectangle still sufficient")
    check(compile_frontier(v,f,max_order=1) is None, "failed cap inconclusive")

    budgets=[]
    for m in (2,3,5,12,20):
        f=(x-y**m,y)
        size,cut=weight_budget(v,f,(m-1,1))
        check((size,cut)==(m+1,m), "sharp triangular weighted capacity")
        values=strand(v,f,1,x,m+1)
        check(values[m]==factorial(m) and values[m+1]==0, "sharp nilpotent strand")
        check(prior.coefficient_matrix(v,values[:m+1]).rank()==m+1,
              "every nonzero strand entry independent")
        for r in range(m+1):
            d=max(s.Poly(h,*v).total_degree() for h in values[:r+1])
            check(d>=degree_floor(2,r), "degree packing lower bound")
        budgets.append((m,size,cut))
    # Arbitrarily late scalar jets, checked at a declared large degree without a long derivative run.
    m=1024
    large=(x-y**m,y)
    check(weight_budget(v,large,(m-1,1))==(1025,1024), "large symbolic weighted invoice")
    _,dual=prior.carrier(v,large)
    check(dual[1,0]==1024*y**1023, "first source derivative already exposes hidden term")
    for r in (2,3,7,20):
        h=s.prod(range(m-r+1,m+1))*y**(m-r)
        check(h != 0 and h.subs(y,0)==0, "point-jet zero is not a polynomial zero")
    false_point_frontier=Frontier(v,large,((0,0),(0,1),(1,0)),(0,0))
    _,_,read=derivative_reader(v,large)
    check(all(h.subs({x:0,y:0})==0
              for alpha in boundary(false_point_frontier.addresses,2) for h in read(alpha)),
          "every proposed boundary value vanishes at the point, but the frontier is false")

    # Nonincreasing weights are sufficient but not necessary even for an explicit automorphism.
    hard=(x+y*y, y+(x+y*y)**2)
    _,dual=prior.carrier(v,hard)
    check(s.Poly(dual[0,0],*v).coeff_monomial(x*y)!=0,
          "support inequality would force a positive y weight <=0")
    check(compile_frontier(v,hard,max_order=4) is not None,
          "frontier succeeds beyond weighted-support gate")

    # Actual inherited Keller hostile: field-dependent recurrence coefficients hide infinite growth.
    w=(x,y,z)
    f3=prior.thm1300(w)
    values=strand(w,f3,2,x,4)
    check(all(h!=0 for h in values), "inherited nonnilpotent strand finite control")
    check(prior.coefficient_matrix(w,values).rank()==5, "constant coefficient rank retained")
    check(s.cancel(values[1]/values[0]).free_symbols, "rational-function recurrence coefficient is not constant")
    check(s.cancel(values[1]-(values[1]/values[0])*values[0])==0,
          "rank-one over source field is the wrong coefficient ring")
    # Outside Keller: Euler has bounded polynomial degree and a nonzero eigenvalue.
    check([s.expand(x*s.diff(x,x)) for _ in range(3)]==[x,x,x],
          "Keller no-eigenvalue premise cannot be omitted")

    malformed=[lambda:audit(replace(triangle,anchor=(True,0))),
        lambda:audit(replace(triangle,addresses=((0,0),(0,True),(0,2),(1,0)))),
        lambda:audit(replace(triangle,addresses=((0,0),(0,2),(1,0)))),
        lambda:audit(replace(triangle,addresses=((0,0),(1,0),(0,1)))),
        lambda:audit(replace(triangle,addresses=((0,0),(0,1),(1,0)))),
        lambda:audit(replace(triangle,source_map=(x-y**3,y))),
        lambda:audit(false_point_frontier),
        lambda:compile_frontier(v,v,max_order=True),
        lambda:compile_frontier(v,(x+1.0,y)),
        lambda:weight_budget(v,f,(False,1)),
        lambda:weight_budget(v,hard,(100,1)),
        lambda:degree_floor(2,1.0),lambda:strand(v,v,True,x,2)]
    for test in malformed:
        try:test()
        except ValueError:check(True,"malformed or false certificate rejected")
        else:check(False,"hostile accepted")
    print('Exact derivative addresses retain characteristic-zero source polynomials')
    print('Three-shear controls (a,b,c,address_count,boundary_count,inverse_degree):',examples)
    print('Triangular weights (degree,capacity,max_joint_order):',budgets)
    print('Degree1024 source: weighted capacity1025; scalar jets may vanish long before polynomial derivatives')
    print('THM1300 strand0..4: independent overQ; degrees:',[s.Poly(h,*w).total_degree() for h in values])
    print('No universal frontier/weight budget; failed caps and failed weight inequalities are inconclusive')
    print('checks:',CHECKS)


if __name__=='__main__':
    main()
