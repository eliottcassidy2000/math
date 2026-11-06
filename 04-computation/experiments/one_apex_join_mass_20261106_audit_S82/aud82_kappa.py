# aud82_kappa.py (independent audit of S82): the constant of Theorem B.
# c' = int_0^1 dlam int_0^lam dmu P(lam*mu*(1 + r s) <= 1), r, s iid U(0,1)   (pairs v < u <= sqrt X)
#    = int_0^1 Gi(q) / (2q) dq,  Gi(q) = int_0^q G(p) dp,  G(p) = P(rs <= 1/p - 1) = 1 (p <= 1/2), w(1 - ln w) (w = 1/p - 1).
# kappa = kappa1 - c', with kappa1 from aud82_kappa.c (limit of D(m)); block mean of r1 = (pi^2/12)(k+2) ln 2 + 2 kappa.
import mpmath as mp
mp.mp.dps = 30
def G(p):
    if p <= mp.mpf(1)/2: return mp.mpf(1)
    w = 1/p - 1
    if w <= 0: return mp.mpf(0)
    return w*(1 - mp.log(w))
def Gi(q):
    if q <= mp.mpf(1)/2: return q
    return mp.mpf(1)/2 + mp.quad(G, [mp.mpf(1)/2, q])
cprime = mp.mpf(1)/4 + mp.quad(lambda q: Gi(q)/(2*q), [mp.mpf(1)/2, 1])
# direct 2-D check of the same constant by Monte-Carlo-free quadrature
c2 = mp.quad(lambda lam: mp.quad(lambda mu: G(lam*mu), [0, lam]), [0, 1])
kappa1 = mp.mpf('0.2474888569')   # limit of D(m), aud82_kappa.out (D(31) = 0.247488858178, remaining block deviations ~ -1.3e-9)
kappa = kappa1 - cprime
print("c' =", mp.nstr(cprime, 15), "  (2-D quadrature:", mp.nstr(c2, 15), ")")
print("kappa1 =", kappa1, "  kappa = kappa1 - c' =", mp.nstr(kappa, 12), "  2 kappa =", mp.nstr(2*kappa, 12))
pl2 = mp.pi**2/12*mp.log(2)
for k in (20, 25, 28, 29):
    print(f"predicted block mean of r1 at k={k}: {mp.nstr(pl2*(k+2) + 2*kappa, 10)}")
