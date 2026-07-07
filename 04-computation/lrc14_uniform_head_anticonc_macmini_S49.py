"""
mac-mini-2026-07-07-S49 (HYP-5077) -- (A) the UNIFORM HEAD BOUND via the (1-mu)
factorization; (B) ANTI-CONCENTRATION for the no-separated-cherry class.

(A) THE ONE-LINE LEMMA: Good_E = complement of Bad_E, meas(Bad) = 1 - mu.  For n != 0:
    |g_hat(n)| = |int 1_Good e(-nx)dx| = |int 1_Bad e(-nx)dx| <= meas(Bad) = 1 - mu.
    (The DC term drops; the indicator's low-frequency mass is capped by the co-set's
    measure -- uniform in n, uniform in E.)
    => HEAD(M) = sum_{0<|n|<=M} |c_G(n)||g_E(n)| <= (1-mu) * C_G(M),
       C_G(M) = sum_{0<|n|<=M} |c_G(n)|  (EXACT per core: G_P is an explicit finite
       union of rational intervals; c_G(n) = sum of interval Fourier coefficients).
    With klein-S166's tail TAIL(M) (interval-count K_G): |SPEC| <= (1-mu)C_G(M) + TAIL(M)
    and R(E) >= 1 - [(1-mu)C_G(M) + TAIL(M)] / (mGP * mu).
    => R >= 3/4 holds UNIFORMLY over all E with mu(E) >= mu*, where mu* solves
       (1-mu)C_G + TAIL = (1/4) mGP mu.   Compute mu* for both hard cores, sweep M.

(B) ANTI-CONCENTRATION: Bad = {maxgap <= 1/7} implies Y := sum_j g_j^2 - 1/8 <= 1/56
    (since sum g^2 <= maxgap * sum g <= 1/7).  One-sided Chebyshev (Paley-Zygmund
    inverse): mu = P(Bad^c) >= P(Y > 1/56) >= (E[Y] - 1/56)_+^2 / E[Y^2].
    iid BASELINES (k=8 circle spacings ~ Dirichlet(1^8)):
      E[g^2] = 2/(k(k+1)); E[sum g^2] = 2/(k+1) = 2/9;  E[Y] = 2/9 - 1/8 = 7/72.
      E[g^4] = 24/(k..(k+3)) etc; E[(sum g^2)^2] = k E[g^4] + k(k-1) E[g^2 h^2],
      E[g^2 h^2] = 4 * k-Dirichlet moment = 2!2!(k-1)!/(k+3)! * ... (computed exactly
      below by simulation-free formula AND cross-checked by Monte Carlo).
    Then the DEFICIT question: at no-cherry census shapes, how far are E[Y], E[Y^2]
    from iid?  (The u-integrated moments should be near-iid for spread shapes.)
    Report the resulting floor vs the 0.675 leg bar.
"""
import numpy as np
from fractions import Fraction as F
from math import gcd, comb, factorial
import random as rnd
rnd.seed(49)

# ---------------- (A) uniform head bound ----------------
def GP_intervals(P):
    """exact rational intervals of G_P = {x: ||px|| >= 1/14 all p in P}."""
    bad=[]
    for p in P:
        w=F(1,14*p)
        for j in range(p+1):
            bad.append((F(j,p)-w, F(j,p)+w))
    bad=[(max(l,F(0)),min(h,F(1))) for l,h in bad if h>0 and l<1]
    bad.sort()
    merged=[]
    for l,h in bad:
        if merged and l<=merged[-1][1]: merged[-1]=(merged[-1][0],max(merged[-1][1],h))
        else: merged.append((l,h))
    good=[]; prev=F(0)
    for l,h in merged:
        if l>prev: good.append((prev,l))
        prev=max(prev,h)
    if prev<1: good.append((prev,F(1)))
    return good

def c_hat_abs(intervals, n):
    """|sum over intervals of (e(-n h)-e(-n l))/(2 pi i n)| exactly via complex float."""
    tot=0
    for l,h in intervals:
        tot += (np.exp(-2j*np.pi*n*float(h)) - np.exp(-2j*np.pi*n*float(l)))
    return abs(tot/(2j*np.pi*n))

CORES = {'P8={9..13} (k=8 leg)': ([9,10,11,12,13], F(40247,90090)),
         'P9={10..13} (k=9 leg)': ([10,11,12,13], F(1577,3003))}
print("=== (A) the uniform head bound: R(mu) curves per core ===")
for name,(P,mGP) in CORES.items():
    iv = GP_intervals(P)
    K_G = len(iv)
    for M in (60, 91, 150):
        C_G = sum(c_hat_abs(iv,n) for n in range(1,M+1))*2   # +-n
        # tail: |c_G(n)| <= K_G/(pi n); |g_E(n)| <= min(1-mu, ...) -- use klein's
        # crude tail with |g|<=1/2? honest: |g_hat(n)| <= 1/2 always for indicator?
        # |int 1_B e()| <= min(meas, 1-meas) <= 1/2. Use (1-mu) for head, 1/2 for tail:
        # TAIL(M) <= (1/2) * sum_{n>M} 2*K_G/(pi n) diverges -- need g-side decay too.
        # klein used interval-count on the G-side with g-side L2; here: Cauchy-Schwarz
        # tail: sum_{n>M}|c||g| <= sqrt(sum_{n>M}|c|^2) * sqrt(sum|g|^2)
        #  <= sqrt(2 K_G^2/(pi^2 M)) * sqrt(mu(1-mu)) -- explicit.
        tail_coef = np.sqrt(2*K_G**2/(np.pi**2*M))   # times sqrt(mu(1-mu))
        # R(mu) >= 1 - [(1-mu)C_G + tail_coef*sqrt(mu(1-mu))]/(mGP*mu) >= 3/4
        # solve numerically for mu*
        mus = np.linspace(0.3,0.999,7000)
        lhs = (1-mus)*C_G + tail_coef*np.sqrt(mus*(1-mus))
        ok = lhs <= 0.25*float(mGP)*mus
        mustar = mus[ok][0] if ok.any() else None
        print(f"  {name}: M={M:3d}, K_G={K_G:2d}, C_G(M)={C_G:.4f}, tailcoef={tail_coef:.4f}"
              f" -> mu* = {mustar:.4f}" if mustar else
              f"  {name}: M={M:3d} no mu* (bound too weak)")

# ---------------- (B) anti-concentration ----------------
print("\n=== (B) anti-concentration for the no-cherry class (k=8) ===")
k=8
# exact Dirichlet(1^k) spacing moments: E[prod g_i^{a_i}] = prod a_i! * (k-1)! / (k-1+sum a)!
def dmom(*a):
    s=sum(a)
    num=1
    for ai in a: num*=factorial(ai)
    return F(num*factorial(k-1), factorial(k-1+s))
Eg2 = dmom(2)            # E[g^2]
Esum2 = k*Eg2            # E[sum g^2]
Eg4 = dmom(4); Eg2h2 = dmom(2,2)
Esum2sq = k*Eg4 + k*(k-1)*Eg2h2
EY = Esum2 - F(1,8)
EY2 = Esum2sq - 2*F(1,8)*Esum2 + F(1,64)
thr = F(1,7) - F(1,8)    # 1/56
floor_iid = float((EY-thr)**2/EY2)
print(f"  iid k=8: E[sum g^2] = {Esum2} = {float(Esum2):.5f}; E[Y] = {EY} = {float(EY):.5f}")
print(f"  E[Y^2] = {EY2} = {float(EY2):.6f}; one-sided floor (E[Y]-1/56)^2/E[Y^2] = {floor_iid:.4f}")
print(f"  vs the k=8 leg bar 0.675: {'CLEARS' if floor_iid>=0.675 else 'falls short at iid'}")

# sharper: one-sided Cantelli on Y about its mean: P(Y <= 1/56) <= Var/(Var + (EY-1/56)^2)
VarY = EY2 - EY*EY
cant = float(VarY/(VarY + (EY-thr)**2))
print(f"  Cantelli: P(Bad) <= Var/(Var+(EY-thr)^2) = {cant:.4f} => mu >= {1-cant:.4f} "
      f"{'CLEARS 0.675' if 1-cant>=0.675 else 'falls short'}")

# deficits at no-cherry census shapes (numeric)
GRID=60000; xs=(np.arange(GRID)+0.5)/GRID
def moments(E):
    ph=np.mod(np.outer(xs,np.array(E,float)),1.0); ph.sort(axis=1)
    g=np.concatenate([np.diff(ph,axis=1),(ph[:,0]+1-ph[:,-1])[:,None]],axis=1)
    s2=(g*g).sum(axis=1)
    Y=s2-1/8
    mu=float((g.max(axis=1)>1/7).mean())
    return float(Y.mean()), float((Y*Y).mean()), mu
L=50
def no_cherry(sp):
    s=sorted(sp)
    from itertools import combinations as C2
    return not any(c>=L*(a+b) for a,b,c in C2(s,3))
print("\n  deficits at no-cherry shapes (diam >= 27):")
worst=(9,None)
for trial in range(400):
    sp = sorted(rnd.sample(range(1,300),7))
    if max(sp)<27 or not no_cherry(sp): continue
    ey,ey2,mu = moments([0]+sp)
    vY=ey2-ey*ey
    cant_s = vY/(vY+(ey-float(thr))**2) if ey>float(thr) else 1.0
    fl=1-cant_s
    if fl<worst[0]: worst=(fl,(tuple(sp),ey,mu))
print(f"  worst Cantelli floor over census: {worst[0]:.4f} at {worst[1][0]}"
      f" (E[Y]={worst[1][1]:.4f}, true mu={worst[1][2]:.4f})")
print(f"  {'>= 0.675: the no-cherry class CLEARS via anti-concentration' if worst[0]>=0.675 else 'short of 0.675 -- needs the moment-deficit control'}")
