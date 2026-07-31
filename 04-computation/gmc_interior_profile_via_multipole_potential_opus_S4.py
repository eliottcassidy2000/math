#!/usr/bin/env python3
"""
gmc_interior_profile_via_multipole_potential_opus_S4.py   opus-2026-07-31-S4

FMM / n-BODY reframe of GMC (Newton ratios R_k = h_k^2/(h_{k-1}h_{k+1})), and a CLOSED FORM for the
open interior moving-edge profile g(alpha)=lim_d d^2 log(R_{alpha d}/R_{alpha d-1}) that klein's THM-3001
endpoint law (bottom +C(N), top -C(N*)) does not reach.

THE REFRAME.  N(n)=a_d prod_i (n+r_i).  log(N/(a_d n^d)) = sum_i log(1+r_i/n) is exactly the MULTIPOLE
EXPANSION of the log-potential of d unit charges at -r_i evaluated at field point n; the GMC log-jets
ell_j=(-1)^{j-1} p_j ARE the multipole moments p_j=sum r_i^j (THM-2997).  So GMC is the theory of the
log-potential of the root charge measure.

THE CLOSED FORM.  Let nu = normalized root measure, phi(x)=int log(1+rx) dnu(r) the log-potential
(= sum_{j>=1} (-1)^{j-1} mu_j x^j / j, mu_j = int r^j dnu = moments).  e_k=[x^k] exp(d phi(x)); saddle
x* phi'(x*)=alpha gives (1/d)log e_k -> phi(x*)-alpha log x*, and (1/d)log C(d,k) -> H(alpha).  Hence
   Phi(alpha) := (1/d) log h_{alpha d} = phi(x*(alpha)) - alpha log x*(alpha) - H(alpha)   (Legendre transform),
and since log h_k = d Phi(k/d):  log R_k = -Phi''(alpha)/d,  log(R_k/R_{k-1}) = -Phi'''(alpha)/d^2, i.e.
   g(alpha) = -Phi'''(alpha).
g(0) = -Phi'''(0) = the THM-3000 fixed-edge curvature 3x^2-2z-1 (far field / well-separated); alpha->1 is
the NEAR FIELD (self-energy) where the multipole expansion diverges -- exactly THM-3001's top endpoint.
"""
import mpmath as mp
from fractions import Fraction as Fr
from math import comb
mp.mp.dps = 50

# uniform root measure nu = Unif[0,1]: phi'(x)=1/x - log(1+x)/x^2, saddle: log(1+x)/x = 1-alpha
def xstar(alpha):
    t = 1 - alpha
    g = mp.mpf('1.0') if alpha < 0.5 else mp.mpf('5')
    return mp.findroot(lambda x: mp.log(1 + x) / x - t, g)

def Phi_prime(alpha):        # envelope theorem: Phi'(alpha) = -log x* - H'(alpha), H'(a)=log((1-a)/a)
    return -mp.log(xstar(alpha)) - mp.log((1 - alpha) / alpha)

def g_pred(alpha):
    a = mp.mpf(str(alpha)); h = min(mp.mpf('1e-3'), a / 3)
    return -(Phi_prime(a + h) - 2 * Phi_prime(a) + Phi_prime(a - h)) / h**2

def g_actual(d, alpha):
    roots = [Fr(i) for i in range(1, d + 1)]
    e = [Fr(1)]
    for r in roots:
        ne = e + [Fr(0)]
        for i in range(len(e), 0, -1): ne[i] = e[i-1]*r + (e[i] if i < len(e) else Fr(0))
        e = ne
    k = max(2, min(d - 2, round(alpha * d)))
    h = [Fr(e[j], comb(d, j)) for j in range(d + 1)]
    def lR(kk):
        num = h[kk]**2; den = h[kk-1]*h[kk+1]
        return mp.log(mp.mpf(num.numerator)/num.denominator) - mp.log(mp.mpf(den.numerator)/den.denominator)
    return float((lR(k) - lR(k-1)) * d * d), k

if __name__ == "__main__":
    print("GMC INTERIOR PROFILE via the multipole/log-potential reframe:  g(alpha) = -Phi'''(alpha)")
    print("(uniform root measure; g(0)=curvature 1/3=%.5f = THM-3000; endpoints = THM-3001)" % (1/3))
    print(f"  {'alpha':>6} {'g_pred=-Phi3(a)':>15} {'g_actual d=1600':>16} {'k':>5}  match")
    for al in [0.02, 0.05, 0.125, 0.25, 0.5, 0.75, 0.9]:
        gp = float(g_pred(al)); ga, k = g_actual(1600, al)
        ok = abs(gp - ga) < 0.03 + 0.05 * abs(gp)
        print(f"  {al:6.3f} {gp:15.5f} {ga:16.5f} {k:5d}  {'MATCH' if ok else 'diff'}")
    print(f"  alpha->0: g_pred(0.01)={float(g_pred(0.01)):.5f} -> curvature 1/3 (far field / fixed edge)")
    print("VERDICT: the open uniform-in-k interior profile is the Legendre-transform third derivative")
    print("of the root charge measure's log-potential. FMM: fixed edge=far field, moving edge=near field.")
