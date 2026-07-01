#!/usr/bin/env python3
"""
Is the GLOBAL congestion potential enough to force a lonely-measure floor?  mac-mini-S91, HYP-3822.

The potential bound |L| >= 1 - A + Phi/C_max is VALID but treats congestion uniformly.
The SHARP question: given ONLY the global moments A=int C, Phi=int C(C-1), and support {0..n-1},
what is the MINIMUM possible lonely measure m_0=|{C=0}| over all coverage distributions?
  This is a MOMENT LP:  min m_0  s.t.  sum m_k=1, sum k m_k=A, sum k(k-1)m_k=Phi, m_k>=0.
If the LP minimum is 0, then the global potential/PoA argument PROVABLY cannot prove inf meas>0 --
the floor is a LOCAL / arithmetic fact (which t has which coverage), not a global-moment fact.
This is the PoA-lens sharpening of HYP-3817: a MOMENT works only if it is the arithmetic/local
(Gamma_0(N) congruence, HYP-3571) moment; the global empirical moment is blind to the atom, like a
transform.  Also: the potential bound IS positive in the SUB-critical regime A<1 -- find the crossover.
"""
import numpy as np
from scipy.optimize import linprog

def coverage_moments(S, r, G=1_000_000):
    t = (np.arange(G) + 0.5) / G
    C = np.zeros(G, dtype=np.int32)
    for v in S:
        x = (v * t) % 1.0
        C += (np.minimum(x, 1.0 - x) < r).astype(np.int32)
    A = C.mean(); Phi = (C*(C-1)).mean(); Cmax=int(C.max()); L=(C==0).mean()
    hist = np.bincount(C, minlength=Cmax+1) / G
    return A, Phi, Cmax, L, hist

def min_m0_given_moments(A, Phi, N):
    """min m_0 over distributions on {0..N} with mean A and 2nd-factorial-moment Phi."""
    K = np.arange(N+1)
    # variables m_0..m_N ; minimize m_0  => c = e_0
    c = np.zeros(N+1); c[0] = 1.0
    A_eq = np.vstack([np.ones(N+1), K.astype(float), (K*(K-1)).astype(float)])
    b_eq = np.array([1.0, A, Phi])
    res = linprog(c, A_eq=A_eq, b_eq=b_eq, bounds=[(0,1)]*(N+1), method='highs')
    return res

if __name__ == "__main__":
    n = 14; N = n-1
    S_c2 = list(range(1,13)) + [182]
    print("#"*72)
    print("# MOMENT LP: can the GLOBAL potential (A,Phi) force a lonely-measure floor?  (n=14)")
    print("#"*72)
    for r, tag in [(1/14, "r=1/14 (LRC threshold, A~1.86)"), (14/183, "r=14/183 covering-min (A~1.99)")]:
        A, Phi, Cmax, L, hist = coverage_moments(S_c2, r)
        res = min_m0_given_moments(A, Phi, N)
        print(f"\n=== {tag} ===   S={S_c2}")
        print(f"  A={A:.4f}  Phi={Phi:.4f}  C_max={Cmax}  ACTUAL |L|={L:.5f}")
        print(f"  coverage histogram m_k (k=0..):  {np.round(hist,4)}")
        print(f"  --- moment LP: min |L| given (A,Phi,support) ---")
        if res.success:
            m0min = res.x[0]
            witness = {k:round(v,4) for k,v in enumerate(res.x) if v>1e-6}
            print(f"  min m_0 (LP) = {m0min:.6f}   witness dist = {witness}")
            print(f"  => GLOBAL potential {'CANNOT' if m0min<1e-5 else 'CAN'} force a floor "
                  f"(actual |L|={L:.4f} but a same-moment dist with |L|={m0min:.4f} exists)")
        else:
            print(f"  LP failed: {res.message}")

    # crossover: where does the potential bound 1-A+Phi/Cmax go positive? (sub-critical regime)
    print("\n" + "#"*72)
    print("# The potential bound 1-A+Phi/C_max IS positive in the SUB-CRITICAL regime (A<1).")
    print("# Sweep r; report union bound (1-A) vs potential bound (1-A+Phi/Cmax) vs actual |L|.")
    print("#"*72)
    print(f"  {'r':>9} {'A':>7} {'Phi':>8} {'Cmax':>5} {'1-A(union)':>11} {'pot.bound':>10} {'actual|L|':>10}")
    for r in [0.01,0.02,0.03,0.0357,0.04,0.05,0.06,1/14]:
        A, Phi, Cmax, L, hist = coverage_moments(S_c2, r, G=500_000)
        ub = 1-A; pb = 1-A+Phi/Cmax
        print(f"  {r:9.4f} {A:7.4f} {Phi:8.4f} {Cmax:5d} {ub:11.4f} {pb:10.4f} {L:10.5f}")
    print("\n  (A=1 at r=1/(2*13)=0.03846; below it the union/potential bounds are live, above they die.)")
