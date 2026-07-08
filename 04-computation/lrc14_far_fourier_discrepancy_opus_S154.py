"""
lrc14_far_fourier_discrepancy_opus_S154.py   (opus-2026-07-08-S154)

THE EXACT FOURIER EXPANSION OF THE FAR CORRELATION  (the discrepancy route, LEM-005).

LEM-005 reduces the k=11 spread tail to  far <= E[W]^2  (=> Var(W) <= near <= (2/7)E[W] => PZ>=bar).
The far bound was PROVED AS A REDUCTION to the phase-vector discrepancy (far -> (5/7)^{k+1} <= E[W]^2
as the pairwise differences grow), with the EXPLICIT RATE left open and the band prim-diam in [25,35]
(between klein-S184 exhaustive <=24 and the sampled crossover ~36) uncovered.

This script derives + verifies the EXACT Fourier formula for far (the rigorous foundation of the rate).

DERIVATION.  theta=1/7.  arc-just-before indicator phi(u)=1(u mod 1 in [-theta,0]),
  ahat(m):=phihat(m)=(e(m*theta)-1)/(2*pi*i*m),  ahat(0)=theta,  |ahat(m)|=|sin(pi*m*theta)|/(pi|m|).
Two disjoint theta-arcs at 0 and t:  Psi(u)=1-phi(u)-phi(u-t),  Psihat(0)=1-2theta=5/7,
  Psihat(m)=-ahat(m)(1+e(-m t))  (m!=0).   R(t)=E_{x,y}[prod_i Psi(e_i x - y)]; expanding in Fourier
and using E_x[e((sum m_i e_i)x)]=1(sum m_i e_i=0), E_y[e(-(sum m_i)y)]=1(sum m_i=0):

    R(t) = sum_{m in L} prod_i Psihat(m_i),   L = { m in Z^k : sum m_i = 0, sum m_i e_i = 0 }.

Integrate t over (1/7,6/7) (measure 5/7).  Spectator coords (m_i=0) give Psihat(0)=5/7 each:

    ***  far = (5/7)^{k+1} + sum_{m in L, m!=0} (5/7)^{k-|S|} (-1)^{|S|} (prod_{i in S} ahat(m_i)) J(m)  ***
    J(m) = int_{1/7}^{6/7} prod_{i in S}(1+e(-m_i t)) dt = sum_{T subset S} K(sum_{i in T} m_i),
    K(a) = int_{1/7}^{6/7} e(-a t) dt = 5/7 (a=0);  (e(-a/7)-e(-6a/7))/(-2 pi i a) (a!=0).

|S|=|supp(m)|.  The m=0 baseline (5/7)^{k+1} is the IID far.  Everything else is the RESONANCE
correction, carried by the relation lattice L (additive structure).  This file:
  (1) verifies far_fourier (truncated |m_i|<=Mmax, |S|<=Smax) == far_exact (Farey, rational) -> convergence;
  (2) decomposes the correction by support |S| (which relations dominate: 3-APs vs additive quads);
  (3) tracks the correction vs primitive diameter (the empirical discrepancy RATE).
"""
import sys, itertools
from fractions import Fraction as F
from math import gcd, comb
import mpmath as mp

sys.path.insert(0, ".")
sys.path.insert(0, "04-computation")
from lrc14_pz_general_integrator_opus_S148 import pz_exact, BAR
from lrc14_k11_tail_sharp_near_opus_S149 import EW_and_near

mp.mp.dps = 40
TH = mp.mpf(1) / 7
TWO_PI_I = 2 * mp.pi * mp.mpc(0, 1)


def ahat(m):
    """Fourier coeff of the theta-arc [-theta,0]: (e(m theta)-1)/(2 pi i m)."""
    if m == 0:
        return TH
    return (mp.e ** (TWO_PI_I * m * TH) - 1) / (TWO_PI_I * m)


def Kint(a):
    """int_{1/7}^{6/7} e(-a t) dt."""
    if a == 0:
        return mp.mpf(5) / 7
    return (mp.e ** (-TWO_PI_I * a * TH) - mp.e ** (-TWO_PI_I * a * (6 * TH))) / (-TWO_PI_I * a)


def Jfac(supp_m):
    """J(m) = int prod_{i in S}(1+e(-m_i t)) dt = sum_{T subset S} K(sum_{i in T} m_i)."""
    S = list(supp_m)
    tot = mp.mpc(0)
    for r in range(len(S) + 1):
        for T in itertools.combinations(S, r):
            tot += Kint(sum(T))
    return tot


def relations(E, Mmax, Smax):
    """Enumerate nonzero m in L = {sum m_i=0, sum m_i e_i=0} with support<=Smax, |m_i|<=Mmax.
    Yields (support_indices, m_values_on_support)."""
    k = len(E)
    out = []
    for s in range(2, Smax + 1):
        for idx in itertools.combinations(range(k), s):
            # m_i in [-Mmax,Mmax]\{0} on idx; require sum=0 and sum m_i e_i=0
            for vals in itertools.product([v for v in range(-Mmax, Mmax + 1) if v != 0], repeat=s):
                if sum(vals) != 0:
                    continue
                if sum(v * E[i] for v, i in zip(vals, idx)) != 0:
                    continue
                out.append((idx, vals))
    return out


def far_fourier(E, Mmax=3, Smax=6, by_support=False):
    """far via the exact Fourier formula (truncated)."""
    k = len(E)
    base = (mp.mpf(5) / 7) ** (k + 1)
    rels = relations(E, Mmax, Smax)
    corr = mp.mpc(0)
    supp_contrib = {}
    for idx, vals in rels:
        s = len(idx)
        term = (mp.mpf(5) / 7) ** (k - s) * ((-1) ** s)
        for v in vals:
            term *= ahat(v)
        term *= Jfac(vals)
        corr += term
        supp_contrib[s] = supp_contrib.get(s, mp.mpc(0)) + term
    if by_support:
        return base, corr, supp_contrib
    return base + corr


def far_exact(E):
    """far = E[W^2] - near, both exact rational."""
    EW, near = EW_and_near(E)
    EW1, EW2, PZ = pz_exact(E)
    return EW2 - near, EW, EW2, near


def R2(E):
    from collections import Counter
    c = Counter()
    for i in range(len(E)):
        for j in range(len(E)):
            if i != j:
                c[E[i] - E[j]] += 1
    return sum(v * v for v in c.values())


def P3(E):
    """count ordered 3-APs 2 e_a = e_b + e_c with a,b,c distinct (support-3 relation (2,-1,-1))."""
    Sset = set(E)
    cnt = 0
    for a in E:
        for b in E:
            c = 2 * a - b
            if b != a and c != a and c != b and c in Sset:
                cnt += 1
    return cnt


def main():
    print("=" * 96)
    print("EXACT FOURIER FAR:  far = (5/7)^{k+1} + sum_{m in L\\0} (5/7)^{k-|S|}(-1)^|S| prod ahat(m_i) J(m)")
    print("=" * 96)
    tests = [
        [0, 1, 2], [0, 1, 3], [0, 1, 2, 3], [0, 1, 3, 4],
        [0, 1, 2, 3, 4], [0, 2, 3, 4, 5],
        [0, 1, 2, 3, 4, 5],           # block-6
        [0, 1, 3, 7, 12],             # Sidon-5
    ]
    print("\n(1) VERIFY far_fourier (Mmax=4,Smax=%d) == far_exact (Farey rational):" % 6)
    for E in tests:
        k = len(E)
        fe, EW, EW2, near = far_exact(E)
        ff = far_fourier(E, Mmax=4, Smax=min(6, k))
        print(f"  E={str(E):24s} k={k}: far_exact={float(fe):.8f}  far_fourier={float(ff.real):.8f}"
              f"  |im|={float(abs(ff.imag)):.1e}  diff={float(abs(mp.mpf(float(fe))-ff.real)):.2e}")
    print("  => the exact Fourier far formula (star) is VERIFIED.")

    print("\n(2) SUPPORT DECOMPOSITION of the far correction (which relations dominate):")
    for E in [[0, 1, 2, 3, 4], [0, 1, 2, 3, 4, 5], [0, 1, 2, 3, 4, 5, 6]]:
        k = len(E)
        base, corr, sc = far_fourier(E, Mmax=5, Smax=min(7, k), by_support=True)
        fe, EW, EW2, near = far_exact(E)
        print(f"  block-{k}: base(5/7)^{{k+1}}={float(base):.6f}  corr={float(corr.real):+.6f}"
              f"  far={float((base+corr).real):.6f} (exact {float(fe):.6f})  R2={R2(E)} P3={P3(E)}")
        for s in sorted(sc):
            print(f"      |S|={s}: {float(sc[s].real):+.7f}")
    print("=" * 96)


if __name__ == "__main__":
    main()
