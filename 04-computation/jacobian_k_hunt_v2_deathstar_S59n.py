#!/usr/bin/env python3
"""
death-star-2026-07-19-S59n (HYP-8080) -- k-ladder hunt v2: Keller-forced.
Augmented system: vars (B, C, E0, delta, nu); equations c1 (all), c0[1:] (all),
and c0[0]*nu - 1 (forces det = const != 0, deleting the degenerate basin).
VALIDATION: k=2 random starts must now find Keller points.
"""
import random

src = open("04-computation/jacobian_k_ladder_hunt_deathstar_S59n.py").read()
ns = {}
exec(src.split('def newton_hunt')[0], ns)
build_system = ns['build_system']

def hunt(k, degB, degC, degE0, ntrials, seed=11):
    base, NUV0 = build_system(k, degB, degC, degE0)
    NUV = NUV0 + 1
    NB, NC, NE = degB+1, degC+1, degE0+1
    # need c0[0]: rebuild small: base returns c1 coeffs + c0[1:]; we need c0[0] too.
    # trick: c0[0] = det const = evaluate det at t=0? cheaper: recompute via the
    # same polys -- patch: use base() on vars and ALSO compute c0[0] by finite
    # reconstruction: base was built from c1 and c0[1:]; instead rebuild here.
    def system(vars_):
        z = vars_[:NUV0]
        nu = vars_[NUV0]
        res = base(z)
        # c0[0]: reconstruct via polynomial evaluation of det at s=0? c0(t=0):
        # instead compute c0 fully again -- inline minimal recompute:
        B = list(z[:NB]); C = list(z[NB:NB+NC]); E0 = list(z[NB+NC:NB+NC+NE]); delta = z[NB+NC+NE]
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
        A = [1.0]
        for _ in range(k+1): A = polymul(A, [1.0, 1.0])
        D = polyscale(polymul([1.0,1.0],[1.0,1.0]), delta)
        Bp, Cp, E0p = polyd(B), polyd(C), polyd(E0)
        K0 = polyadd(E0, shift(E0p, 1))
        M0 = polyadd(shift(polymul(Bp, D), k), polyscale(shift(polymul(B, D), k-1), k),
                     polyscale(polymul(A, C), -1), polyscale(shift(polymul(A, Cp), 1), -1))
        c0 = polyadd(
            polymul(K0, M0),
            polyscale(polymul(E0p, polyadd(shift(polymul(Bp, D), k+1),
                                           polyscale(shift(polymul(A, Cp), 2), -1))), -1),
            polyscale(shift(polyadd(polymul(Bp, C), polyscale(polymul(B, Cp), -k)), k+1), -1))
        res = list(res) + [c0[0]*nu - 1.0]
        return res
    random.seed(seed)
    hits = []
    for trial in range(ntrials):
        z = [complex(random.uniform(-4,4), random.uniform(-2,2)) for _ in range(NUV)]
        damp = 1e-2
        for it in range(120):
            F = system(z)
            m = len(F)
            h = 1e-7
            J = []
            for j in range(NUV):
                z2 = list(z); z2[j] = z2[j] + h
                F2 = system(z2)
                J.append([(F2[i]-F[i])/h for i in range(m)])
            Ac = [[sum(J[i][r].conjugate()*J[jj][r] for r in range(m)) + (damp if i==jj else 0)
                   for jj in range(NUV)] for i in range(NUV)]
            bc = [-sum(J[i][r].conjugate()*F[r] for r in range(m)) for i in range(NUV)]
            M = [row[:] + [bc[i]] for i, row in enumerate(Ac)]
            ok = True
            for c in range(NUV):
                piv = max(range(c, NUV), key=lambda i: abs(M[i][c]))
                if abs(M[piv][c]) < 1e-18: ok = False; break
                M[c], M[piv] = M[piv], M[c]
                pv = M[c][c]
                M[c] = [x/pv for x in M[c]]
                for i in range(NUV):
                    if i != c and M[i][c] != 0:
                        f = M[i][c]
                        M[i] = [a - f*b for a, b in zip(M[i], M[c])]
            if not ok: break
            dx = [M[i][NUV] for i in range(NUV)]
            z = [zi + di for zi, di in zip(z, dx)]
            damp = max(damp*0.7, 1e-12)
        F = system(z)
        res = sum(abs(x)**2 for x in F)
        if res < 1e-16:
            hits.append((res, z[:NUV0]))
    return hits

print("=== VALIDATION: k=2 Keller-forced hunt ===")
h2 = hunt(2, 2, 2, 1, 60)
print(f"k=2 hits: {len(h2)}/60")
for res, z in h2[:3]:
    print("  sol:", [f"{v.real:+.4f}{v.imag:+.4f}i" for v in z])

print("\n=== k=3 Keller-forced hunt ===")
h3 = hunt(3, 4, 4, 3, 250)
print(f"k=3 hits: {len(h3)}/250")
for res, z in h3[:5]:
    print("  sol:", [f"{v.real:+.4f}{v.imag:+.4f}i" for v in z])
