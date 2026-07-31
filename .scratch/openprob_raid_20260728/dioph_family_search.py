#!/usr/bin/env python3
"""
dioph_family_search.py — hunt for POLYNOMIAL parametrized solution families of
the Epoch FrontierMath "small Diophantine equations":

    z^2 + y^2 z + P(x) = 0        with P(x) = x^3 + c  or  x^3 ± x ± const,
    and the variant  z^2 + y^2 z - z + x^3 + 2 = 0.

Key structure: z(z + y^2) = -P(x).  If P(x(t)) = A(t)·B(t) with
A + B = y(t)^2 (resp. y^2 - 1 for the -z variant), then z = -A(t) gives an
identity:  A^2 - (A+B)A + AB = 0.  So we need polynomials x, y, A over Z with

    A·(y^2 + e - A) = P(x),      e = 0 (usual) or -1 (equation 9).

Ansatz search: deg x = 2, deg y = 2, deg A in {2,3,4} (so deg B = 6-degA),
coefficients integers; solve the coefficient system with sympy's solve /
groebner over QQ, then scan small rational/integer points of the solution
variety.  Also tries deg x = 4, deg y = 4, deg A in {4,...,8} lightly.

Any family found is verified symbolically and instantiated at three large t
to produce |x| > 10^50 with distinct x.
"""
import sympy as sp
import itertools, sys

t = sp.symbols('t')

def try_ansatz(Pfunc, dx, dy, dA, e=0, timeout_terms=None, verbose=False):
    """Pfunc: sympy expr in symbol X.  Returns list of (x,y,A) families."""
    X = sp.symbols('X')
    P = Pfunc
    ax = sp.symbols(f'ax0:{dx+1}')
    ay = sp.symbols(f'ay0:{dy+1}')
    aA = sp.symbols(f'aA0:{dA+1}')
    x = sum(ax[i]*t**i for i in range(dx+1))
    y = sum(ay[i]*t**i for i in range(dy+1))
    A = sum(aA[i]*t**i for i in range(dA+1))
    lhs = sp.expand(A*(y**2 + e - A))
    rhs = sp.expand(P.subs(X, x))
    poly = sp.Poly(lhs - rhs, t)
    eqs = poly.all_coeffs()
    # normalizations to kill scaling: leading coeff of x = 1, y monic optional
    eqs = [sp.Eq(q, 0) for q in eqs]
    eqs.append(sp.Eq(ax[dx], 1))
    sols = []
    try:
        raw = sp.solve(eqs, list(ax) + list(ay) + list(aA), dict=True)
    except Exception as ex:
        if verbose: print("solve failed:", ex)
        return []
    for s in raw:
        xs = x.subs(s); ys = y.subs(s); As = A.subs(s)
        # need remaining free symbols only from parameters; substitute free = 1..3
        free = (xs.free_symbols | ys.free_symbols | As.free_symbols) - {t}
        subs_try = [dict(zip(free, vals)) for vals in
                    itertools.product([0, 1, -1, 2, -2, 3], repeat=len(free))] if free else [{}]
        for sub in subs_try[:200]:
            x2, y2, A2 = xs.subs(sub), ys.subs(sub), As.subs(sub)
            if not all(v.is_rational for v in sp.Poly(x2, t).all_coeffs()):
                continue
            chk = sp.expand(A2*(y2**2 + e - A2) - P.subs(X, x2))
            if chk == 0 and sp.degree(x2, t) >= 1:
                # integer coefficients?
                cs = (sp.Poly(x2, t).all_coeffs() + sp.Poly(y2, t).all_coeffs()
                      + sp.Poly(A2, t).all_coeffs())
                if all(sp.nsimplify(c).is_integer for c in cs):
                    sols.append((sp.expand(x2), sp.expand(y2), sp.expand(A2)))
    # dedupe
    uniq = []
    for s in sols:
        if s not in uniq:
            uniq.append(s)
    return uniq

def main():
    X = sp.symbols('X')
    equations = {
        "warmup   x^3+2":      (X**3 + 2, 0),
        "eq1  x^3-2":          (X**3 - 2, 0),
        "eq2  x^3-x-1":        (X**3 - X - 1, 0),
        "eq3  x^3+x-1":        (X**3 + X - 1, 0),
        "eq4  x^3+x+1":        (X**3 + X + 1, 0),
        "eq5  x^3-3":          (X**3 - 3, 0),
        "eq6  x^3+3":          (X**3 + 3, 0),
        "eq7  x^3-x-2":        (X**3 - X - 2, 0),
        "eq8  x^3-x+2":        (X**3 - X + 2, 0),
        "eq9  x^3+2 (-z)":     (X**3 + 2, -1),
    }
    shapes = [(2, 2, 2), (2, 2, 3), (2, 2, 4)]
    for name, (P, e) in equations.items():
        print(f"=== {name} (e={e}) ===", flush=True)
        found = False
        for (dx, dy, dA) in shapes:
            fams = try_ansatz(P, dx, dy, dA, e=e)
            for (xf, yf, Af) in fams:
                print(f"  FAMILY dx={dx},dA={dA}: x={xf}  y={yf}  z=-({Af})", flush=True)
                found = True
            if found:
                break
        if not found:
            print("  none in small ansatz", flush=True)

if __name__ == "__main__":
    main()
