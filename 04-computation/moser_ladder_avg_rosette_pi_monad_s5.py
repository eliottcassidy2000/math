import sympy, math
from math import isqrt
def is_sq(n):
    if n<0: return -1
    r=isqrt(n); return r if r*r==n else -1
def B_ideal(t):
    res=1
    for p,e in sympy.factorint(t).items():
        if p%3==1: res*=(e+1)
        elif p==3: res*=1
        else: res*=1 if e%2==0 else 0
    return res
def degenerate(t):
    D=4*t-1
    if is_sq(D)>=0: return True
    if D%3==0 and is_sq(D//3)>=0: return True
    return False

# Asymptotic: avg over ALL t of #units = 12 + 6*avg(B). 
# sum_{t<=T} r_E(t) = #Eisenstein pts (nonzero) of norm<=T ~ (2pi/sqrt3) T  (hex Gauss circle)
# => avg r_E -> 2pi/sqrt3 ;  avg #units -> 12 + 2pi/sqrt3.
target = 12 + 2*math.pi/math.sqrt(3)
print(f"predicted asymptotic average rosette size = 12 + 2pi/sqrt3 = {target:.6f}")
print(f"{'T':>8}{'avg #units (all t)':>22}{'avg (non-degen t)':>20}")
for T in (1000,10000,100000,1000000):
    s_all=0; n_all=0; s_nd=0; n_nd=0
    for t in range(2,T+1):
        u=12+6*B_ideal(t)
        s_all+=u; n_all+=1
        if not degenerate(t): s_nd+=u; n_nd+=1
    print(f"{T:>8}{s_all/n_all:>22.6f}{s_nd/n_nd:>20.6f}")

# record rungs: t giving record-high #units (highly-Loeschian)
print("\nRecord rungs (new max #unit-vectors as t increases), t<=200000:")
best=-1
rec=[]
for t in range(2,200001):
    if degenerate(t): continue
    u=12+6*B_ideal(t)
    if u>best:
        best=u; rec.append((t,u,sympy.factorint(t)))
for t,u,f in rec:
    print(f"  t={t:>6}  #units={u:>4}  B={(u-12)//6:>3}  t={dict(f)}")
