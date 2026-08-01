"""Boundary test in ONE variable, where FC is known TRUE -- and cheap.
N = d+1 monomials (deg<=d, constant included).  Prediction: the window
m=1..M is solvable iff M < N-1."""
import numpy as np
from math import factorial
from scipy.optimize import least_squares

def Lm(c, m):
    # f = sum c_j s^j ; L(s^k)=k!
    p = np.zeros(1, dtype=complex); p[0] = 1
    for _ in range(m):
        p = np.convolve(p, c)
    return sum(p[k]*factorial(k) for k in range(len(p)))

for d in (2,3,4):
    N = d+1
    print(f"n=1, deg<={d}: N={N} monomials, 2N={2*N} params, predicted solvable iff M < {N-1}")
    for M in range(1, N+2):
        def resid(x):
            c = x[:N] + 1j*x[N:]
            out=[]
            for m in range(1,M+1):
                v = Lm(c,m); out += [v.real, v.imag]
            out.append(float(np.sum(np.abs(c)**2)) - 1.0)
            return out
        rng = np.random.default_rng(7)
        best = 1e9
        for _ in range(60):
            x0 = rng.normal(size=2*N); x0 /= np.linalg.norm(x0)
            r = least_squares(resid, x0, xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=5000)
            best = min(best, float(np.max(np.abs(resid(r.x)))))
            if best < 1e-11: break
        got = "SOLVED" if best < 1e-9 else "NOT FOUND"
        pred = "solvable" if M < N-1 else "NOT solvable"
        flag = "" if (got=="SOLVED")==(M<N-1) else "   <-- PREDICTION BROKEN"
        print(f"   M={M}: {2*M+1:2d} eqs vs {2*N} params  resid {best:.2e}  {got:10s} (pred {pred}){flag}")
    print()
