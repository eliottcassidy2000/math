import sys
from fractions import Fraction
from functools import reduce
from math import gcd
sys.stdout.reconfigure(encoding='utf-8') if hasattr(sys.stdout,'reconfigure') else None

def covered(S, r):
    """True iff M(S) <= r : danger arcs of half-width r/v around m/v cover [0,1]."""
    arcs=[]
    for v in S:
        for m in range(0, v+1):
            lo=Fraction(m,v)-Fraction(r,v); hi=Fraction(m,v)+Fraction(r,v)
            if hi<0 or lo>1: continue
            arcs.append((lo if lo>0 else Fraction(0), hi if hi<1 else Fraction(1)))
    arcs.sort()
    cur=Fraction(0)
    for lo,hi in arcs:
        if lo>cur: return False
        if hi>cur: cur=hi
        if cur>=1: return True
    return cur>=1

def exact_M(S):
    """M(S) exact. M is rational with denominator dividing some v_i+-v_j (<=2maxS).
       Binary search r in [0,1/2]; snap to denominator L=lcm-ish bound."""
    S=tuple(sorted(set(S)))
    maxs=max(S)
    L=1
    # candidate denominators: all v_i+-v_j
    dens=set()
    for i in range(len(S)):
        for j in range(i+1,len(S)):
            dens.add(S[i]+S[j]); dens.add(abs(S[i]-S[j]))
    dens.discard(0)
    # binary search by float to localize, then test exact candidates near it
    lo,hi=0.0,0.5
    for _ in range(60):
        mid=(lo+hi)/2
        if covered(S, Fraction(mid).limit_denominator(4*maxs+4)):
            hi=mid
        else:
            lo=mid
    # M ~ hi. exact value = smallest r over candidates a/d (d in dens, a< d/2) with covered and not covered just below
    best=None
    for d in dens:
        for a in range(1, d//2+1):
            r=Fraction(a,d)
            if abs(float(r)-hi)<0.01:
                if covered(S,r):
                    # check it's the min such (not covered below by tiny)
                    below=r-Fraction(1,(4*maxs**2*len(S)))
                    if below>0 and not covered(S,below):
                        if best is None or r<best: best=r
    return best, hi

print("=== EXACT M for family F(k,mult)={1..k-2,k,mult*(k-1)} at primorial-1 k ===")
print("    (verify M > floor 1/(k+1): LRC must hold; report achieved level & g*k^2)\n")
cases=[(2,5),(3,7),(3,13),(4,31),(5,211)]
for mult,k in cases:
    S=sorted(set(list(range(1,k-1))+[k,mult*(k-1)]))
    if len(S)!=k:
        print(f"  mult={mult} k={k}: size mismatch {len(S)}"); continue
    M,approx=exact_M(S)
    floor=Fraction(1,k+1)
    if M is None:
        print(f"  mult={mult} k={k}: exact snap failed, M~{approx:.6f} floor={float(floor):.6f}")
        continue
    p,q=M.numerator,M.denominator
    e=p*(k+1)-q
    lvl = p if e==1 else f"p/q={p}/{q}, e={e}"
    above_floor = M>floor
    gk2=float(M-floor)*k*k
    print(f"  mult={mult} k={k:4d}: M={M} ({float(M):.6f})  floor=1/{k+1}({float(floor):.6f})  "
          f"M>floor={above_floor}  level a={lvl}  g*k^2={gk2:.4f}")
