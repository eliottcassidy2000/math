#!/usr/bin/env python3
"""
paley_starstar_catalytic_fit_monad.py
monad-explorer-2026-06-07 (8th session)

HUNT for the CATALYTIC (Tutte/BMJ) functional equation of the cycle-rank
generating function  U(x,y) = sum_{k>=1} x^k sum_{m=1}^k t(k,m) y^m,
where  t(k,m) = sum_{even-series sigma, rank m} prod_v(|B_v|-1)!  and
U(x,-1) = F(x)-1,  F = (-1+sqrt(1+4x))/(2x)  solves the loop equation x F^2 + F = 1.

[x^k]U is a polynomial in y of degree k (catalytic variable y marks cycle rank m).
We posit a quadratic-in-U equation, single power of x (root removal lowers k by 1):

   U = x * [ A(y) + B(y) U + C(y) U^2 + D(y) U1 + E(y) (U - U1)/(y-1) ]

with U=U(x,y), U1=U(x,1), and A,B,C,D,E unknown polynomials in y of degree <= DEG.
Matching order x^k (k=1..KMAX) gives a LINEAR system in the unknown coefficients
(U, U^2, U1, DD are known series).  We solve; a consistent low-degree solution that
also PREDICTS k beyond the fitting range is a conjectured catalytic equation.
"""
import sys
import sympy as sp

y = sp.symbols('y')

# triangle t[k][m]  (k=1..6, validated; k=6 from fast enumerator cross-checked vs SVD)
TRI = {
    1: [1],
    2: [1, 3],
    3: [1, 9, 13],
    4: [1, 18, 72, 69],
    5: [1, 30, 230, 580, 421],
    6: [1, 45, 560, 2626, 4845, 2867],
}
KMAX = max(TRI)

# U as list of y-polynomials: Uc[k] = sum_m t(k,m) y^m  (k=1..KMAX), Uc[0]=0
Uc = [sp.Integer(0)] * (KMAX + 1)
for k in TRI:
    Uc[k] = sum(TRI[k][m - 1] * y**m for m in range(1, k + 1))

# U1c[k] = [x^k]U(x,1)
U1c = [sp.expand(Uc[k].subs(y, 1)) for k in range(KMAX + 1)]

# U^2 series coefficients: (sum_k Uc[k] x^k)^2
U2c = [sp.Integer(0)] * (KMAX + 1)
for k in range(KMAX + 1):
    s = sp.Integer(0)
    for i in range(k + 1):
        s += Uc[i] * Uc[k - i]
    U2c[k] = sp.expand(s)

# DD[k] = [x^k] (U - U1)/(y-1)
DDc = [sp.cancel((Uc[k] - U1c[k]) / (y - 1)) for k in range(KMAX + 1)]
DDc = [sp.expand(d) for d in DDc]
# y-derivative series:  Uy[k] = [x^k] d/dy U ;  U1y[k] = its value at y=1
Uyc = [sp.expand(sp.diff(Uc[k], y)) for k in range(KMAX + 1)]
U1yc = [sp.expand(Uyc[k].subs(y, 1)) for k in range(KMAX + 1)]

# basis of known [x^{k-1}] series we may use on the RHS (all multiplied by an x):
#   const-source (only k=1), U, U^2, U1, DD, Uy, U1y
BASIS = "ABCDEFG"   # A=source, B=U, C=U^2, D=U1, E=DD, F=Uy, G=U1y
SER = {"B": Uc, "C": U2c, "D": U1c, "E": DDc, "F": Uyc, "G": U1yc}


def fit(DEG, kfit, use):
    """Fit polynomials (deg<=DEG in y) for the terms in `use` (subset of BASIS)
       using orders x^1..x^kfit. Returns (sol, polys, coeffs)."""
    coeffs, polys = {}, {}
    for name in use:
        cs = sp.symbols(f'{name}0:{DEG+1}')
        coeffs[name] = cs
        polys[name] = sum(cs[j] * y**j for j in range(DEG + 1))
    eqs = []
    for k in range(1, kfit + 1):
        rhs = sp.Integer(0)
        for name in use:
            if name == "A":
                if k == 1:
                    rhs += polys["A"]
            else:
                rhs += polys[name] * SER[name][k - 1]
        expr = sp.expand(Uc[k] - rhs)
        p = sp.Poly(expr, y) if expr != 0 else None
        if p is not None:
            for c in p.all_coeffs():
                eqs.append(sp.expand(c))
    allc = [c for name in use for c in coeffs[name]]
    sol = sp.solve(eqs, allc, dict=True)
    if not sol:
        return None, polys, coeffs
    return sol[0], polys, coeffs


def verify(use, polys, sol, free_zero=True):
    coeffs_used = {name: sp.expand(polys[name].subs(sol)) for name in use}
    if free_zero:
        frees = set()
        for v in coeffs_used.values():
            frees |= v.free_symbols
        frees.discard(y)
        z = {f: 0 for f in frees}
        coeffs_used = {n: sp.expand(v.subs(z)) for n, v in coeffs_used.items()}
    okk = []
    for kk in range(1, KMAX + 1):
        rhs = sp.Integer(0)
        for name in use:
            if name == "A":
                if kk == 1:
                    rhs += coeffs_used["A"]
            else:
                rhs += coeffs_used[name] * SER[name][kk - 1]
        okk.append(sp.expand(Uc[kk] - rhs) == 0)
    return coeffs_used, okk


def main():
    # candidate term-sets to try (A=source must be present)
    trials = [
        ("AF", "source + Uy only (pure ODE in y)"),
        ("ACF", "source + U^2 + Uy"),
        ("ABCF", "source + U + U^2 + Uy"),
        ("ABCDF", "+ U1"),
        ("ABCDEF", "+ divided-diff"),
        ("ABCDEFG", "all incl U1y"),
        ("ABCEF", "U,U^2,DD,Uy"),
        ("ACEF", "U^2,DD,Uy"),
    ]
    for use, desc in trials:
        use = tuple(use)
        for DEG in range(0, 5):
            kfit = KMAX - 1
            sol, polys, coeffs = fit(DEG, kfit, use)
            if sol is None:
                continue
            cu, okk = verify(use, polys, sol)
            holds_all = all(okk)
            predicts6 = okk[KMAX - 1]
            if holds_all:
                print(f"[{use}] {desc}  DEG={DEG}: HOLDS all k<=6 "
                      f"(fit k<=5, predicted k=6 = {predicts6})")
                for name in use:
                    print(f"      {name}: {cu[name]}")
                print()
                break
        else:
            print(f"[{''.join(use)}] {desc}: no consistent solution up to DEG=4")


if __name__ == "__main__":
    main()


if __name__ == "__main__":
    main()
