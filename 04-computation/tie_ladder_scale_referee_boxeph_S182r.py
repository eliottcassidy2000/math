#!/usr/bin/env python3
"""tie_ladder_scale_referee_boxeph_S182r.py — hostile-referee check on THM-1635 §1/§3
premise beta_j(m) = beta_j(1 + c_j/m + ...). Computes beta_eff(m) = I_m[corrected arc]/
I_m[pure arc] as a Gamma(m/2)-weighted expectation (pure-python Simpson, log-space).
FINDING: a half-step Puiseux correction t = C r^{-1/2}(1 + u r^{-1/2}) — generic, since
inverting r(t) = C^2/t^2 + A/t + ... produces it — gives beta_eff(m) = e^{-u sqrt(2m)}
(1 + O(m^{-1/2})), NOT 1 + c/m: log(beta_eff)/sqrt(2m) -> -u (-0.2957 at m=512, u=0.3).
Integer-step control t = C r^{-1/2}(1 + v/r) DOES give a 1/m ladder, but around the
DRESSED constant e^{-2v} beta (0.5488 at v=0.3). Scale set is {e^{a sqrt m} m^{-k/2}}."""
import math
u = 0.3; v = 0.3
def beta_eff(m, K):
    c = m/2.0; sd = math.sqrt(c)
    lo = max(1e-9, c - 25*sd); hi = c + 25*sd + 80.0
    N = 40001; h = (hi-lo)/(N-1); num = den = 0.0
    for i in range(N):
        r = lo + i*h
        w = 1 if i in (0, N-1) else (2 if i % 2 == 0 else 4)
        g = math.exp(-r + (c-1)*math.log(r) - (c-1)*math.log(c) + (c-1))
        num += w*g*K(r, m); den += w*g
    return num/den
KH = lambda r, m: (1+u/math.sqrt(r))**(-(m+1)) * (1+2*u/math.sqrt(r))
KI = lambda r, m: (1+v/r)**(-(m+1)) * (1+3*v/r)
print(" m    betaH(m)       log(betaH)/sqrt(2m)   betaI(m)    m*(betaI-e2v)/e2v")
e2v = math.exp(-2*v)
for m in [8, 16, 32, 64, 128, 256, 512]:
    bH = beta_eff(m, KH); bI = beta_eff(m, KI)
    print(f"{m:4d}  {bH:12.6e}  {math.log(abs(bH))/math.sqrt(2*m):+.6f}           {bI:.6f}    {m*(bI-e2v)/e2v:+.4f}")
print(f"\nverdict: half-step arcs live at scale e^{{-u sqrt(2m)}} (col3 -> {-u}); THM-1635's")
print(f"1/m-Vandermonde premise is model-specific; int-step base is DRESSED e^-2v = {e2v:.4f}")
