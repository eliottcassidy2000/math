#!/usr/bin/env python3
"""alignment_arithmetic_kps_S128c70.py -- kind-pasteur S128 cont.70.
THE ARITHMETIC BEHIND THE ALIGNMENT BIAS.

Decoding the cont.68/69 worst case (core [1,3,5,6,7,8,11,12], killers 371/374/377/379):
its longest surviving component is [1373/5278, 1385/5306].  Since 5278 = 14*377 and
5306 = 14*379, that runs from the RIGHT edge of tooth j=98 of comb 377 to the LEFT edge of
tooth j=99 of comb 379.  So the gap is not a statistical accident -- it is the interval
between tooth j of one comb and tooth j+1 of another.

THE LAW.  For combs a<b and the gap from tooth j of a to tooth j+1 of b,
    raw length = (j+1)/b - j/a = (a - j(b-a)) / (a b),
so with d = b-a the gap is LINEAR IN j, running from ~1/a at j=0 down to 0 at j=a/d.
Subtracting the two tooth radii 1/(14a)+1/(14b) gives the usable length.
CONSEQUENCE: close moduli (small d) keep the gap near a full period 1/a for MANY j -- that
IS the 'alignment bias', and it is arithmetic, not resonance.  Verify and quantify.
PRINT DATA ONLY."""
import sys, itertools, random
from fractions import Fraction as F
sys.stdout.reconfigure(line_buffering=True)
def nd(x):
    fx=x-int(x); fx=fx+1 if fx<0 else fx; return min(fx,1-fx)
def safe_set(P):
    bps={F(0),F(1)}
    for p in P:
        for j in range(p+1):
            for s in (F(1,14*p),-F(1,14*p)):
                v=F(j,p)+s
                if 0<=v<=1: bps.add(v)
    B=sorted(bps); out=[]
    for i in range(len(B)-1):
        a,b=B[i],B[i+1]
        if b<=a: continue
        if all(nd(p*((a+b)/2))>=F(1,14) for p in P): out.append((a,b))
    mg=[]
    for a,b in out:
        if mg and mg[-1][1]==a: mg[-1]=(mg[-1][0],b)
        else: mg.append((a,b))
    return mg
def rm(iv,k):
    out=[]
    for a,b in iv:
        jlo=int(a*k)-1; jhi=int(b*k)+1; cur=a
        for j in range(jlo,jhi+1):
            x=F(j,k)-F(1,14*k); y=F(j,k)+F(1,14*k)
            if y<=a or x>=b: continue
            x=max(x,a); y=min(y,b)
            if x>cur: out.append((cur,x))
            cur=max(cur,y)
        if cur<b: out.append((cur,b))
    return out
print("### (1) verify the endpoint law on the worst case ###")
P=[1,3,5,6,7,8,11,12]; ks=(371,374,377,379)
iv=safe_set(P)
for k in ks: iv=rm(iv,k)
L=max(b-a for a,b in iv); seg=[s for s in iv if s[1]-s[0]==L][0]
print("  component = [%s, %s]"%(seg[0],seg[1]))
for k in ks:
    for endp,nm in [(seg[0],'left'),(seg[1],'right')]:
        for j in range(int(endp*k)-1,int(endp*k)+2):
            if endp==F(j,k)+F(1,14*k): print("     %s endpoint = RIGHT edge of tooth j=%d of comb %d"%(nm,j,k))
            if endp==F(j,k)-F(1,14*k): print("     %s endpoint = LEFT  edge of tooth j=%d of comb %d"%(nm,j,k))
print("  predicted raw gap (j=98, a=377, b=379): (a - j(b-a))/(ab) = %s"%(F(377-98*2,377*379)))
print("  minus radii 1/(14*377)+1/(14*379) -> %s ; actual L = %s"%(
    F(377-98*2,377*379)-F(1,14*377)-F(1,14*379),L))
print()
print("### (2) the linear-in-j law: gap(a,b,j) = (a - j*d)/(ab) - radii ###")
print("  a=371 b=379 (d=8):  j     raw gap        usable      usable*b")
for j in [0,5,10,20,40,46]:
    raw=F(371-j*8,371*379); use=raw-F(1,14*371)-F(1,14*379)
    print("                       %-5d %-14.8f %-11.8f %.4f"%(j,float(raw),float(use),float(use*379)))
print("  a=371 b=372 (d=1):  j     raw gap        usable      usable*b")
for j in [0,50,150,300,370]:
    raw=F(371-j*1,371*372); use=raw-F(1,14*371)-F(1,14*372)
    print("                       %-5d %-14.8f %-11.8f %.4f"%(j,float(raw),float(use),float(use*372)))
print()
print("### (3) so where do real surviving components sit?  small j = big gap ###")
print("  measure j* = index of the tooth bounding the longest component, vs a/d")
random.seed(70)
C8=[sorted(c) for c in itertools.combinations(range(1,13),8)]
rows=[]
for _ in range(160):
    Pp=random.choice(C8); M=max(Pp); lo=13*M+1
    k1=random.randint(lo,lo+700); kk=[k1]
    for _ in range(3): kk.append(kk[-1]+random.randint(1,5))
    ivv=safe_set(Pp)
    for k in kk: ivv=rm(ivv,k)
    if not ivv: continue
    LL=max(b-a for a,b in ivv); sg=[s for s in ivv if s[1]-s[0]==LL][0]
    mid=(sg[0]+sg[1])/2
    jstar=int(mid*kk[-1])
    d=min(kk[i+1]-kk[i] for i in range(3))
    rows.append((float(LL*kk[-1]),jstar,kk[-1],d,float(jstar*d/kk[-1])))
rows.sort()
print("  worst 6 by L*k4:   L*k4    j*     k4    d_min   j*d/k4")
for r in rows[:6]: print("                     %-8.4f %-6d %-5d %-7d %.4f"%r)
print("  best 3 by L*k4:")
for r in rows[-3:]: print("                     %-8.4f %-6d %-5d %-7d %.4f"%r)
import statistics
print("  correlation check: median j*d/k4 among worst half = %.4f ; best half = %.4f"%(
    statistics.median([r[4] for r in rows[:len(rows)//2]]),
    statistics.median([r[4] for r in rows[len(rows)//2:]])))
print("DONE")
