from fractions import Fraction as F
import random
def fn(x):
    r=x-int(x)
    if r<0: r+=1
    return min(r,1-r)
lam=F(1,14)
print("PROJECTION LAW under pi: t -> u=2t  (the quotient by the free involution s).")
print("Claim (a) v EVEN, v=2w:  D_v = pi^{-1}(comb of speed w at level lam).")
print("Claim (b) v ODD:         pi(D_v) = comb of speed v at level 2*lam.")
bad_a=bad_b=0; N=6000
for _ in range(N):
    t=F(random.randint(0,9999),10000); u=2*t
    w=random.randint(1,7); v=2*w
    if (fn(v*t)<lam) != (fn(w*u)<lam): bad_a+=1
for _ in range(N):
    u=F(random.randint(0,9999),10000); t=u/2
    v=2*random.randint(0,6)+1
    inproj = (fn(v*t)<lam) or (fn(v*(t+F(1,2)))<lam)
    if inproj != (fn(v*u)<2*lam): bad_b+=1
print(f"  (a) mismatches {bad_a}/{N}   (b) mismatches {bad_b}/{N}")
print()
print("MEASURE ACCOUNTING of the reduction (does it gain?)")
V=list(range(1,14)); ev=[v for v in V if v%2==0]; od=[v for v in V if v%2]
before=len(V)*2*lam
after=len(ev)*2*lam+len(od)*4*lam
print(f"  V={{1..13}}: {len(ev)} even, {len(od)} odd")
print(f"  union bound BEFORE projecting: {len(V)}*2lam = {before} = {float(before):.3f}")
print(f"  union bound AFTER  projecting: {len(ev)}*2lam + {len(od)}*4lam = {after} = {float(after):.3f}")
print(f"  -> the reduction LOSES ({float(after):.3f} > {float(before):.3f}): odd combs DOUBLE.")
print("  Recorded as a dead end for the measure route, per the dead-end rule.")
