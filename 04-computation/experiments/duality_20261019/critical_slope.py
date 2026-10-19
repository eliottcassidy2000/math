"""Sub/critical slopes: Christoffel (balanced) rational 3x+1 cycles along the convergents of an irrational a converge to a
bounded circle iff a < log_3 2 = 0.63093 (then 3 < 2^(1/a)); at larger a the minima blow up (negative clocks)."""
from fractions import Fraction
from math import log, floor
def mech(p, q, t): return [(((k+1)*p + t)//q) - ((k*p + t)//q) for k in range(q)]
def carry(w):
    P, Bc, Q = 1, 0, 1
    for b in w:
        if b: P, Bc, Q = 3*P, 3*Bc + Q, 2*Q
        else: Q = 2*Q
    return Bc
def convergents(x, n):
    a=[]; h0,h1,k0,k1 = 0,1,1,0; y=x; out=[]
    for _ in range(n):
        ai=floor(y); h0,h1 = h1, ai*h1+h0; k0,k1 = k1, ai*k1+k0; out.append((h1,k1))
        if y-ai < 1e-12: break
        y = 1/(y-ai)
    return out
theta = log(2)/log(3)
for name, a in (("1/phi = 0.6180 (subcritical)", (5**0.5-1)/2), ("0.625 (subcritical)", 0.625), ("log_3 2 (critical)", theta), ("0.64 (supercritical)", 0.64), ("2/3 (supercritical)", 2/3)):
    print(f"slope {name}:")
    for (p, q) in convergents(a, 14)[1:12]:
        if q < 3 or q > 400 or p >= q: continue
        D = 2**q - 3**p; el = [carry(mech(p,q,t)) for t in range(q)]
        print(f"   p/q = {p}/{q}: D = 2^{q}-3^{p} = {D if abs(D) < 10**12 else ('%.3e' % D)}; x_min = {min(el)/D:.4f}, x_max = {max(el)/D:.4f}")
