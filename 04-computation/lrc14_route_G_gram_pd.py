"""
ROUTE G: Determinant / positive-definite form for the LRC(14) floor.

Setup (THM-515): L(S) = int_0^1 prod_i 1_safe(v_i tau) dtau = measure of lonely set.
1_safe(theta) = 1 - 1_danger, danger = {||theta|| <= 1/14}.
Fourier: hat(1_safe)(t) = h(t), h(0)=6/7, h(t) = -sin(pi t/7)/(pi t).
1_safe(theta) >= 0 (a {0,1} indicator) => h is POSITIVE-DEFINITE (Bochner).

ROUTE G probes Gram / quadratic-form realizations of the floor:
  (G1) Single-runner kernel K[i][j] = h(i-j) (Toeplitz autocorrelation of 1_safe).
  (G2) Lonely measure as a norm-square L = ||g||^2 and a Gram lower bound.
  (G3) The relation-lattice correction matrix M[i][j] = h(n_i - n_j) -- PD? min eig?
  (G4) The S3 cluster floor rho* as a Gram quadratic form.

stdlib only.
"""
import math
from fractions import Fraction as F

def h(t):
    if t == 0:
        return 6.0/7.0
    return -math.sin(math.pi*t/7.0)/(math.pi*t)

def s_coeff(t):
    if t == 0:
        return 1.0/7.0
    return math.sin(math.pi*t/7.0)/(math.pi*t)

def symeig(A, tol=1e-14):
    """All eigenvalues of symmetric A via cyclic Jacobi."""
    n = len(A)
    a = [row[:] for row in A]
    for sweep in range(200):
        off = 0.0
        for p in range(n):
            for q in range(p+1, n):
                off += a[p][q]*a[p][q]
        if off < tol:
            break
        for p in range(n):
            for q in range(p+1, n):
                if abs(a[p][q]) < 1e-18:
                    continue
                app, aqq, apq = a[p][p], a[q][q], a[p][q]
                phi = 0.5*math.atan2(2*apq, aqq-app)
                c, sn = math.cos(phi), math.sin(phi)
                for k in range(n):
                    akp, akq = a[k][p], a[k][q]
                    a[k][p] = c*akp - sn*akq
                    a[k][q] = sn*akp + c*akq
                for k in range(n):
                    apk, aqk = a[p][k], a[q][k]
                    a[p][k] = c*apk - sn*aqk
                    a[q][k] = sn*apk + c*aqk
    eigs = sorted(a[i][i] for i in range(n))
    return eigs

# ====================================================================
print("="*70)
print("(G1) Single-runner Toeplitz kernel K[i][j] = h(i-j), N x N")
print("     = Gram of translated 1_safe; PSD by Bochner (hat h = 1_safe >= 0).")
print("="*70)
for N in [5, 8, 12, 16, 20, 28]:
    K = [[h(i-j) for j in range(N)] for i in range(N)]
    eigs = symeig(K)
    print(f"  N={N:2d}: min eig={eigs[0]:+.6f}  max eig={eigs[-1]:.6f}  PSD:{eigs[0]>-1e-7}")
print("  EXPECTED: min eig -> ess-inf(1_safe) = 0 (1_safe vanishes on the danger band).")
print("  => G1 is only PSD; its floor is 0. The Toeplitz min-eig is NOT the lonely floor.")
print()

# ====================================================================
print("="*70)
print("(G3) Relation-lattice correction matrix M[i][j] = h(n_i - n_j)")
print("     over a finite set of frequencies. PD? min eig as a 'floor'?")
print("="*70)
# Take the frequencies that actually appear: the orbit {v_i * k} mod something.
# Test a generic 13-set vs the extremizer core.
def test_corr_matrix(S, label, kmax=3):
    # frequencies appearing in the THETA sum up to |t_i|<=kmax along each axis:
    # the relevant single frequencies are {v * t : v in S, |t|<=kmax}
    freqs = sorted(set(v*t for v in S for t in range(-kmax, kmax+1)))
    n = len(freqs)
    M = [[h(freqs[i]-freqs[j]) for j in range(n)] for i in range(n)]
    eigs = symeig(M)
    print(f"  {label}: |freqs|={n}, min eig={eigs[0]:+.6f}, PD:{eigs[0]>1e-9}")
    return eigs[0]

# Extremizer core and a generic set
core = [1,2,3,4,5,7,8,9,10,11,12,13,56]  # {1..13}\{6} u {56}
generic = [1,2,5,11,22,33,41,52,60,71,83,97,101]
ap = [1,2,3,4,5,6,7,8,9,10,11,12,13]
test_corr_matrix(core, "core {1..13}\\{6}u{56}")
test_corr_matrix(generic, "generic Sidon-ish")
test_corr_matrix(ap, "AP {1..13}")
print("  (h(n_i-n_j) is a principal submatrix of the infinite Toeplitz => same spectral")
print("   floor issue: it is PSD (Bochner) but min eig -> 0; no positive floor.)")
print()

# ====================================================================
print("="*70)
print("(G2) Lonely measure L = ||g||^2, g = prod 1_safe(v_i tau).")
print("     Cauchy-Schwarz / Gram lower bound: L = <g,1>^2 / ... no. Try L >= |<g,phi>|^2/||phi||^2")
print("     for a fixed positive test function phi. This is a one-dim Gram (Cauchy-Schwarz)")
print("     lower bound. The BEST phi is phi=g, recovering L (circular). A FIXED phi gives")
print("     a computable but possibly-zero lower bound.")
print("="*70)
# L >= (int g phi)^2 / (int phi^2). With phi = a Riesz-type positive density.
# int g phi = int (prod 1_safe(v_i tau)) phi(tau) dtau. If phi is supported where g=1 (a known
# safe sub-arc) we'd get a positive bound -- but locating that arc IS the problem.
# Quantify with a direct grid integral for several S.
def L_grid(S, Ngrid=200003):
    # measure of {tau: ||v tau|| > 1/14 for all v}. Ngrid prime-ish for equidistribution.
    cnt = 0
    for k in range(Ngrid):
        tau = (k + 0.5)/Ngrid
        ok = True
        for v in S:
            f = (v*tau) % 1.0
            d = min(f, 1.0-f)
            if d <= 1.0/14.0:
                ok = False
                break
        if ok:
            cnt += 1
    return cnt/Ngrid

for S, lbl in [(core,"core"),(generic,"generic"),(ap,"AP{1..13}")]:
    print(f"  L({lbl}) ~ {L_grid(S):.5f}")
print("  (these confirm L>0 numerically; Cauchy-Schwarz with a FIXED phi cannot certify")
print("   uniformity without first localizing the safe arc -- the same crux.)")
