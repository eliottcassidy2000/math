#!/usr/bin/env python3
"""
verblunsky_lonely_measure_klein.py  --  klein-2026-07-01-S69

VERBLUNSKY COEFFICIENTS of the LRC lonely measure -- the RECURSIVE encoding of a measure on the unit
circle (OPUC / Szego recursion), applied to L_C. Pushes the "runners on a loop" picture to the unit-circle
theory of orthogonal polynomials: the lonely measure mu = (1/L) 1_{L_C} dt is a probability measure on the
circle, so (Verblunsky's theorem) it is EQUIVALENT to a sequence of Verblunsky coefficients alpha_n in the
disk, built one at a time by the Szego recursion. This is a genuine RECURSIVE metaphor for L_C.

KEY EXPECTATIONS (OPUC theory, to verify):
  * alpha_0 = -c_1 = -moment_1 = -(hat1(1)/L). For a measure concentrated at t*, 1-t* (the binding pair),
    c_1 = cos(2 pi t*), so alpha_0 ~ -cos(2 pi t*) -- the S66 TWO-ATOM law, now the FIRST Verblunsky coeff.
  * A measure with EXACTLY N atoms has |alpha_{N-1}| = 1 (termination). The construction's extremal lonely
    measure is 2 atoms {t*, 1-t*} (THM-515 / mac-mini HYP-3789), so as r -> M_C we expect |alpha_1| -> 1.
  * For r < M_C (L_C = union of intervals, a.c.), alpha_n is an infinite DECAYING sequence (Szego).
  * Any Phi6 / prime periodicity in alpha_n mirrors the phase structure (S68, p(w)=nw mod Phi6).

Moments from FFT of 1_{L_C}; Verblunsky via the (real) Levinson-Durbin / Szego recursion.
"""
import numpy as np

def norms(v, t):
    x = (v*t) % 1.0
    return np.minimum(x, 1-x)

def moments_of_LC(C, r, N=600000, K=60):
    t = np.arange(N)/N
    G = np.full(N, 1.0)
    for v in C: G = np.minimum(G, norms(v, t))
    f = (G > r).astype(np.float64)
    L = f.mean()
    Fhat = (np.fft.rfft(f)/N).real          # hat1(k), real (1_{L_C} even)
    c = Fhat[:K+1]/L                          # normalized moments c_k = hat1(k)/L ; c_0 = 1
    return c, L

def verblunsky(c):
    """Real Levinson-Durbin: returns reflection/Verblunsky coeffs a_n (|a_n|=|alpha_n|) and errors E_n."""
    K = len(c)-1
    a = [1.0]           # prediction poly coeffs a[0..n]
    E = c[0]            # = 1
    alphas = []; Es = [E]
    for n in range(1, K+1):
        if E <= 1e-13:              # E->0 => measure is (numerically) atomic: TERMINATION reached
            break
        acc = c[n] + sum(a[j]*c[n-j] for j in range(1, n))
        k = -acc/E
        k = max(-1.0, min(1.0, k))  # clamp tiny numerical overshoot past |k|=1
        newa = [1.0] + [a[j] + k*a[n-j] for j in range(1, n)] + [k]
        a = newa
        E = E*(1-k*k)
        alphas.append(k)            # reflection coeff; Verblunsky alpha_{n-1} = -k (|.| identical)
        Es.append(E)
    return np.array(alphas), np.array(Es)

if __name__ == "__main__":
    n=14; C=list(range(1,n-1))+[n*(n-1)]; Phi6=n*n-n+1; tstar=n/Phi6; Mc=n/Phi6
    print(f"n={n} construction C={C} t*={tstar:.5f}=14/183 M_C={Mc:.5f}  -cos(2π t*)={-np.cos(2*np.pi*tstar):+.4f}")
    print("Verblunsky (reflection) coeffs alpha_n of the lonely measure mu=(1/L)1_{L_C}dt at several levels r:\n")
    for r in [0.05, 0.06, 0.07, 0.074, 0.0762]:
        c, L = moments_of_LC(C, r, K=40)
        al, E = verblunsky(c)
        # alpha_0 vs two-atom law
        a0 = al[0]; two_atom = -np.cos(2*np.pi*tstar)
        print(f" r={r:.4f} (meas L={L:.4f}): alpha_0={a0:+.4f} vs -cos(2πt*)={two_atom:+.4f} (c_1/L two-atom law)")
        # magnitudes of first several, and how fast |alpha_n| grows toward 1 (termination = atomic)
        mags = np.abs(al[:12])
        print(f"    |alpha_n| n=0..11: " + " ".join(f"{m:.3f}" for m in mags))
        print(f"    max|alpha_n| (n<=39) = {np.nanmax(np.abs(al)):.4f}  (->1 signals atomic/termination); E_40={E[-1]:.2e}")
    print()
    # r -> M_C : approach the 2-atom termination |alpha_1| -> 1
    print("APPROACH TO THE BINDING (r -> M_C): does |alpha_1| -> 1 (2-atom termination at t*,1-t*)?")
    def g(al, i): return abs(al[i]) if i < len(al) else float('nan')
    for r in [0.070, 0.073, 0.075, 0.0762, 0.0764]:
        c, L = moments_of_LC(C, r, K=12)
        al,_ = verblunsky(c)
        term = len(al)   # termination index (E hit ~0 here => ~term atoms)
        print(f"    r={r:.4f} L={L:.5f}: alpha_0={al[0]:+.4f} |alpha_1|={g(al,1):.4f} |alpha_2|={g(al,2):.4f} |alpha_3|={g(al,3):.4f}  (terminated at n={term})")
    print()
    # closed-form 2-atom check: mu=(delta_{t*}+delta_{1-t*})/2 => c_k=cos(2πk t*); expect alpha_0=-cos, |alpha_1|=1
    ck = np.array([np.cos(2*np.pi*k*tstar) for k in range(13)])
    al2,_ = verblunsky(ck)
    print(f"CLOSED-FORM 2-atom {{t*,1-t*}}: alpha_0={al2[0]:+.4f} (=-cos(2πt*)={-np.cos(2*np.pi*tstar):+.4f}), |alpha_1|={abs(al2[1]):.4f} (=1 => exactly 2 atoms)")
    print(f"  => the extremal lonely measure (2 atoms at the Phi6-denominator points) TERMINATES at n=1: |alpha_1|=1.")
    print(f"  => the Verblunsky sequence RECURSIVELY refines S66's two-atom law (alpha_0) into the full L_C.")
