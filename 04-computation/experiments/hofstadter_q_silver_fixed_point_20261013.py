"""
hofstadter_q_silver_fixed_point_20261013.py -- the second-order statistics of the Hofstadter Q jitter as an exact fixed
point.  Exact: d(n) + d(n-1)/2 + d(n-2)/2 = s(n) + O(1/n) with s(n) = (p/n) d(p) + (q/n) d(q), p(n) = n - Q(n-1),
q(n) = n - Q(n-2), and q(n+1) = p(n) + 1.  Model: s is white except for the lag-1 correlation c = rho(1)/2 produced by
the shared adjacent index (weights 1/2); the ARMA(2,1) fixed point then gives rho(1)^2 - 2 rho(1) - 1 = 0, i.e.
rho(1) = 1 - sqrt 2, var(d)/var(s) = 1 + 1/sqrt 2, and (with var(s) = var_prev/2 under white sampling) the sampled
variance of the previous generation exceeds var(d) by the factor 4 - 2 sqrt 2; white sampling predicts
rho_pq = corr(d(p(n)), d(q(n))) = 0 (an earlier draft of this docstring wrongly asked for 3 - 2 sqrt 2).
Measured here to five decimals on the octaves 2^20, 2^22, 2^24, 2^26: rho(1), rho_s(1), var(d)/var(s), rho_pq, the
cross-correlations corr(d(p(n)), d(p(n+1))), corr(d(q(n)), d(p(n+1))), and corr(d(q(n+1)), d(p(n)+1)) [= 1 exactly by
the index identity q(n+1) = p(n) + 1]; the increment statistics of Q (fraction |Q(n) - Q(n-1)| <= 2, rms of (Q(n) - Q(n-1))/n); and the
weights p/n, q/n (mean, std).
"""
import numpy as np, time
from numba import njit

N = 1 << 27

@njit(cache=True)
def hofq(N):
    q = np.zeros(N + 1, dtype=np.int32); q[1] = 1; q[2] = 1
    for n in range(3, N + 1):
        i = n - q[n - 1]; j = n - q[n - 2]
        if i < 1 or j < 1:
            q[0] = n; break
        q[n] = q[i] + q[j]
    return q

@njit(cache=True)
def pieces(q, lo, hi):
    m = hi - lo
    d = np.zeros(m); s = np.zeros(m); dp = np.zeros(m); dq = np.zeros(m); dp1 = np.zeros(m); wp = np.zeros(m)
    small = 0; inc2 = 0.0
    for n in range(lo, hi):
        p = n - q[n - 1]; qq = n - q[n - 2]
        d[n - lo] = q[n] / n - 0.5
        dp[n - lo] = q[p] / p - 0.5; dq[n - lo] = q[qq] / qq - 0.5; dp1[n - lo] = q[p + 1] / (p + 1) - 0.5
        wp[n - lo] = p / n
        s[n - lo] = (p / n) * dp[n - lo] + (qq / n) * dq[n - lo]
        inc = q[n] - q[n - 1]
        if abs(inc) <= 2: small += 1
        inc2 += (inc / n) ** 2
    return d, s, dp, dq, dp1, wp, small / m, np.sqrt(inc2 / m)

def corr(a, b):
    a = a - a.mean(); b = b - b.mean()
    return float((a * b).mean() / np.sqrt((a * a).mean() * (b * b).mean()))

t0 = time.time(); q = hofq(N); print(f"Q to 2^27 in {time.time()-t0:.1f}s")
r2 = np.sqrt(2.0)
print(f"silver-ratio predictions: rho(1) = 1 - sqrt2 = {1-r2:.5f}; var(d)/var(s) = 1 + 1/sqrt2 = {1+1/r2:.5f}; rho_pq (white sampling) = 0; rho_s(1) = (1 - sqrt2)/2 = {(1-r2)/2:.5f}; sampled-variance ratio 4 - 2 sqrt2 = {4-2*r2:.5f}")
print(" octave | rho(1)    rho_s(1)   var(d)/var(s)  rho_pq    corr(dp,dp+)  corr(dq,dp+)  corr(dq(n+1),dp(n)+1)=1  frac|dQ|<=2  rms(dQ/n)  mean p/n  std p/n")
for J in (20, 22, 24, 26):
    lo, hi = 1 << J, 1 << (J + 1)
    d, s, dp, dq, dp1, wp, frac, rinc = pieces(q, lo, hi)
    rho1 = corr(d[:-1], d[1:]); rs1 = corr(s[:-1], s[1:]); vr = d.var() / s.var()
    rpq = corr(dp, dq); c1 = corr(dp[:-1], dp[1:]); c2 = corr(dq[:-1], dp[1:]); c3 = corr(dq[1:], dp1[:-1])
    print(f" 2^{J}   | {rho1:+.5f}  {rs1:+.5f}   {vr:.5f}        {rpq:+.5f}  {c1:+.5f}      {c2:+.5f}      {c3:+.5f}            {frac:.4f}      {rinc:.5f}   {wp.mean():.5f}  {wp.std():.5f}")
print("done")
