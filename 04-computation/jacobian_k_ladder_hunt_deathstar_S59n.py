#!/usr/bin/env python3
"""
death-star-2026-07-19-S59n (HYP-8080) -- the k-ladder existence hunt, honest version.

Variables: B (deg<=k+1), C (deg<=k+1), E0 (deg<=k), delta   [A = (1+t)^{k+1},
D = delta (1+t)^2 from the c2-relation; E1 = 1 normalized].
Equations: all t-coefficients of c1 (must vanish) and of c0 - c0(0) (constancy).
Complex damped Newton, many random starts.  VALIDATION: k=2 must rediscover the
known witness (B = 4+7t+3t^2, C = 1+12t+9t^2 scaled, E0 = 2-3t, delta = 3).
"""
import random, cmath

def build_system(k, degB, degC, degE0):
    """returns eval(vars)->list of complex residuals; vars = B..,C..,E0..,delta.
    Implemented by symbolic univariate polys with coefficient slots."""
    from fractions import Fraction as Fr
    NB, NC, NE = degB+1, degC+1, degE0+1
    NUV = NB + NC + NE + 1
    # representation: dict t-deg -> dict of monomials in vars (tuple of (idx,pow)) -> coeff
    # simpler: build c1, c0 as functions numerically at eval time using poly ops on lists
    def polymul(a, b):
        out = [0j]*(len(a)+len(b)-1)
        for i, x in enumerate(a):
            if x == 0: continue
            for j, y in enumerate(b):
                if y == 0: continue
                out[i+j] += x*y
        return out
    def polyadd(*ps):
        n = max(len(p) for p in ps)
        out = [0j]*n
        for p in ps:
            for i, x in enumerate(p): out[i] += x
        return out
    def polyscale(p, s): return [x*s for x in p]
    def polyd(p): return [p[i]*i for i in range(1, len(p))]
    def shift(p, m): return [0j]*m + list(p)
    # fixed
    A = [1.0]
    for _ in range(k+1): A = polymul(A, [1.0, 1.0])
    def system(vars_):
        B = list(vars_[:NB]); C = list(vars_[NB:NB+NC]); E0 = list(vars_[NB+NC:NB+NC+NE])
        delta = vars_[NB+NC+NE]
        D = polyscale(polymul([1.0,1.0],[1.0,1.0]), delta)
        Ap, Bp, Cp, Dp, E0p = polyd(A), polyd(B), polyd(C), polyd(D), polyd(E0)
        M1 = polyadd(polymul(Ap, D), polyscale(polymul(A, Dp), -1))
        c1 = polyadd(
            polymul(E0, M1),
            polyscale(polymul(E0p, polymul(A, D)), k-1),
            polyscale(shift(polymul(Bp, D), k), -2),
            polyscale(shift(polymul(B, D), k-1), -2*k),
            polyscale(shift(polymul(B, Dp), k), k),
            polyscale(polymul(A, C), k+1),
            polyscale(shift(polymul(A, Cp), 1), k+1),
            polyscale(shift(polymul(Ap, C), 1), -1))
        K0 = polyadd(E0, shift(E0p, 1))
        M0 = polyadd(shift(polymul(Bp, D), k), polyscale(shift(polymul(B, D), k-1), k),
                     polyscale(polymul(A, C), -1), polyscale(shift(polymul(A, Cp), 1), -1))
        c0 = polyadd(
            polymul(K0, M0),
            polyscale(polymul(E0p, polyadd(shift(polymul(Bp, D), k+1),
                                           polyscale(shift(polymul(A, Cp), 2), -1))), -1),
            polyscale(shift(polyadd(polymul(Bp, C), polyscale(polymul(B, Cp), -k)), k+1), -1))
        res = [x for x in c1 if True] + [x for x in c0[1:]]
        # c1 all coeffs; c0 all except t^0
        return [x for x in c1] + list(c0[1:])
    return system, NUV

def newton_hunt(k, degB, degC, degE0, ntrials=600, seed=7):
    system, NUV = build_system(k, degB, degC, degE0)
    random.seed(seed)
    hits = []
    for trial in range(ntrials):
        z = [complex(random.uniform(-4,4), random.uniform(-2,2)) for _ in range(NUV)]
        damp = 1e-2
        for it in range(150):
            F = system(z)
            m = len(F)
            # numeric Jacobian (forward diff, complex step not valid for non-analytic? system IS polynomial=analytic)
            h = 1e-7
            J = []
            for j in range(NUV):
                z2 = list(z); z2[j] = z2[j] + h
                F2 = system(z2)
                J.append([(F2[i]-F[i])/h for i in range(m)])
            # solve (J^H J + damp I) dx = -J^H F
            Ac = [[sum(J[i][r].conjugate()*J[jj][r] for r in range(m)) + (damp if i==jj else 0)
                   for jj in range(NUV)] for i in range(NUV)]
            bc = [-sum(J[i][r].conjugate()*F[r] for r in range(m)) for i in range(NUV)]
            M = [row[:] + [bc[i]] for i, row in enumerate(Ac)]
            ok = True
            for c in range(NUV):
                piv = max(range(c, NUV), key=lambda i: abs(M[i][c]))
                if abs(M[piv][c]) < 1e-16: ok = False; break
                M[c], M[piv] = M[piv], M[c]
                pv = M[c][c]
                M[c] = [xx/pv for xx in M[c]]
                for i in range(NUV):
                    if i != c and M[i][c] != 0:
                        f = M[i][c]
                        M[i] = [aa - f*bb for aa, bb in zip(M[i], M[c])]
            if not ok: break
            dx = [M[i][NUV] for i in range(NUV)]
            z = [zi + di for zi, di in zip(z, dx)]
            damp = max(damp*0.7, 1e-12)
        F = system(z)
        res = sum(abs(x)**2 for x in F)
        # reject the degenerate collapse B=C=E0=0 and delta=0
        norm = sum(abs(x)**2 for x in z[:-1])
        if res < 1e-16 and norm > 1e-6 and abs(z[-1]) > 1e-6:
            # also demand C not identically ~0 (F2 must be nondegenerate)
            hits.append((res, [complex(round(v.real,8), round(v.imag,8)) for v in z]))
    return hits

print("=== VALIDATION: k=2 hunt (must find the known witness) ===")
hits2 = newton_hunt(2, 2, 2, 1, ntrials=200)
print(f"k=2 hits: {len(hits2)}/200")
uniq = []
for res, z in sorted(hits2)[:400]:
    key = tuple((round(v.real,4), round(v.imag,4)) for v in z)
    if any(max(abs(complex(*a)-complex(*b)) for a,b in zip(key,s)) < 1e-3 for s in uniq): continue
    uniq.append(key)
print(f"k=2 distinct solutions found: {len(uniq)}")
for s in uniq[:6]:
    print("   ", [f"{a:+.4f}{b:+.4f}i" for a, b in s])

print("\n=== k=3 hunt ===")
hits3 = newton_hunt(3, 4, 4, 3, ntrials=600)
print(f"k=3 hits: {len(hits3)}/600")
uniq3 = []
for res, z in sorted(hits3)[:400]:
    key = tuple((round(v.real,4), round(v.imag,4)) for v in z)
    if any(max(abs(complex(*a)-complex(*b)) for a,b in zip(key,s)) < 1e-3 for s in uniq3): continue
    uniq3.append(key)
print(f"k=3 distinct solutions: {len(uniq3)}")
for s in uniq3[:8]:
    print("   ", [f"{a:+.4f}{b:+.4f}i" for a, b in s])
