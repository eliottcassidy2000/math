"""Convergents of theta = log2 3 and the Eliahou-type exclusion of Syracuse cycles.
A cycle of the Syracuse map U(n) = (3n+1)/2^{v_2(3n+1)} with k odd elements, total exponent K and minimum m satisfies
   2^K/3^k = prod_i (1 + 1/(3 n_i))  =>  0 < K ln 2 - k ln 3 < k ln(1 + 1/(3m)) < k/(3m).
Hence m < k/(3 delta) with delta = K ln 2 - k ln 3 > 0, and 0 < K/k - theta < 1/(3 m ln 2); if k^2 <= 3 m ln 2 / 2 then
K/k (in lowest terms) is a convergent of theta (Legendre: |theta - K/k| < 1/(2k^2)).  For a convergent p/q > theta, every cycle with K/k = p/q
(any common multiple) has m < q/(3 (p ln 2 - q ln 3)) =: m_max(p/q).  Usage: python3 log23_convergents.py [M0_log2]"""
import sys
from mpmath import mp, mpf, log, floor
mp.dps = 120
theta = log(3)/log(2)
M0 = mpf(2)**int(sys.argv[1]) if len(sys.argv)>1 else mpf(2)**71
# continued fraction
x = theta; a = []
for i in range(60):
    ai = int(floor(x)); a.append(ai); x = 1/(x-ai)
p0,q0,p1,q1 = 1,0,a[0],1
conv=[(p1,q1)]
for ai in a[1:]:
    p0,q0,p1,q1 = p1,q1,ai*p1+p0,ai*q1+q0
    conv.append((p1,q1))
print("partial quotients:", a[:40])
print("%-4s %-28s %-28s %-7s %-24s %-24s %s"%("j","p (=K)","q (=k)","sign","delta=p ln2 - q ln3","m_max = q/(3 delta)","status vs M0=2^%s"%(sys.argv[1] if len(sys.argv)>1 else "71")))
first_alive=None
for j,(p,q) in enumerate(conv):
    delta = p*log(2) - q*log(3)
    if delta <= 0:
        print("%-4d %-28d %-28d %-7s %-24s %-24s %s"%(j,p,q,"below","-","-","K/k < theta: no cycle (2^K < 3^k)")); continue
    mmax = q/(3*delta)
    st = "EXCLUDED (m_max < M0)" if mmax < M0 else "alive"
    if st=="alive" and first_alive is None: first_alive=(j,p,q)
    print("%-4d %-28d %-28d %-7s %-24s %-24s %s"%(j,p,q,"above",mp.nstr(delta,12),mp.nstr(mmax,12),st))
    if q > 10**30: break
j,p,q = first_alive
from mpmath import sqrt
kL = sqrt(3*M0*log(2)/2)
print("\nFirst surviving convergent: j=%d, K/k = %d/%d (m_max = %s)."%(j,p,q,mp.nstr(q/(3*(p*log(2)-q*log(3))),6)))
print("Legendre threshold: k <= sqrt(3 M0 ln 2 / 2) = %s; every cycle with k below it has K/k a convergent, hence k >= q_%d = %d > threshold: contradiction."%(mp.nstr(kL,8),j,q))
print("So (Legendre route) every nontrivial Syracuse cycle has more than %s odd steps; Hercher's Corollary 29 (finer analysis) gives k >= %d."%(mp.nstr(kL,8),q))
# irrationality-measure parametrisation: |theta - K/k| > c k^-mu  and  < 1/(3 m ln 2)  =>  k > (3 c m ln 2)^(1/mu)
for mu in (2, 3, 4.117, 5, 8.616):
    print("  if |theta - p/q| > q^(-%s) for all q, a cycle with minimum >= M0 needs k > (3 M0 ln 2)^(1/%s) = %s"%(mu,mu,mp.nstr((3*M0*log(2))**(1/mpf(mu)),6)))
print("  the actual continued fraction (first 40 partial quotients <= 55) behaves like mu = 2 in this range: m_max(j) ~ a_(j+1) q_j^2 / (3 ln 2)")
