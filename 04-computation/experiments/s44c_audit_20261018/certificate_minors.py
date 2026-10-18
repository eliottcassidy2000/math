#!/usr/bin/env python3
"""Uniqueness certificates for THM-4710's recursion for ALL N >= 9 (complements witt_symbolic.py).

For L in {3,5,7,9} take the condition matrix K_L(N) of witt_symbolic.closure_symbolic (unknowns: the coefficients of
g in the partition basis and of D_m g in the bases of M_(L-m), 2 <= m <= L; rows: D_1 g = 0 and
D_m g - sum b M_(L-m) = 0).  Generic nullity is 1.  Choose a maximal nonsingular minor at a random large integer
point, clear denominators row by row (polynomial entries in N), bound the degree of its determinant by the sum of
the row degrees, evaluate the determinant EXACTLY (flint) at that many + 1 integer points, and interpolate.  The
interpolated polynomial is the determinant (a polynomial of degree <= D is determined by D+1 values).
At every integer N0 >= 9 that is not a root of this determinant nor of the row denominators, the nullity of K_L(N0)
is <= 1, hence dim S_L(N0) <= 1 (and = 1, since the generic solution specializes).  The integer roots are listed for
direct exact checking (witt_closure_exact.py).
"""
import json, os, time
import sympy as sp
import flint
from witt_symbolic import closure_symbolic, K, Nsym

HERE = os.path.dirname(os.path.abspath(__file__))


def main():
    t0 = time.time()
    S, M, cond = closure_symbolic()
    out = {}
    for L in (3, 5, 7, 9):
        rows, n, nb = cond[L]
        nvar = n + nb
        symrows = [[sp.together(K.to_sympy(x)) for x in r] for r in rows]
        N0 = 1000003
        num = flint.fmpq_mat([[flint.fmpq(int(sp.numer(e).subs(Nsym, N0)), int(sp.denom(e).subs(Nsym, N0))) for e in r]
                              for r in symrows])
        R, rk = num.rref()
        assert rk == nvar - 1
        # pivot columns of rows-space and pivot rows (via transpose)
        Rl = [[R[i, j] for j in range(nvar)] for i in range(rk)]
        pc = [next(j for j in range(nvar) if Rl[i][j] != 0) for i in range(rk)]
        RT, rkT = num.transpose().rref()
        RTl = [[RT[i, j] for j in range(len(rows))] for i in range(rkT)]
        pr = [next(j for j in range(len(rows)) if RTl[i][j] != 0) for i in range(rkT)]
        sub = [[symrows[i][j] for j in pc] for i in pr]
        dens = [sp.lcm([sp.denom(e) for e in r]) for r in sub]
        prow = [[sp.expand(sp.cancel(e * d)) for e in r] for r, d in zip(sub, dens)]
        degs = [max(sp.Poly(e, Nsym).degree() if e != 0 else 0 for e in r) for r in prow]
        D = sum(degs)
        pts = []
        N1 = 10
        while len(pts) < D + 1:
            if all(sp.Poly(d, Nsym).eval(N1) != 0 for d in dens):
                mat = flint.fmpz_mat([[int(sp.Poly(e, Nsym).eval(N1)) if e != 0 else 0 for e in r] for r in prow])
                pts.append((N1, int(mat.det())))
            N1 += 1
        poly = sp.interpolate(pts, Nsym)
        poly = sp.factor(poly)
        roots = sorted(set(int(r) for r in sp.Poly(sp.expand(poly), Nsym).ground_roots() if r.is_integer))
        denroots = sorted(set(int(r) for d in dens for r in sp.Poly(d, Nsym).ground_roots() if r.is_integer)) \
            if any(sp.Poly(d, Nsym).degree() > 0 for d in dens) else []
        out[L] = {'nvar': nvar, 'minor_size': len(pc), 'degree_bound': D, 'det': str(poly),
                  'integer_roots_det': roots, 'integer_roots_row_denominators': denroots}
        print(f'L={L}: minor {len(pc)}x{len(pc)}, degree bound {D}; det = {poly}; integer roots {roots}; '
              f'row-denominator integer roots {denroots}  ({time.time()-t0:.1f}s)', flush=True)
    json.dump(out, open(os.path.join(HERE, 'certificate_minors.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
