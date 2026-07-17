#!/usr/bin/env python3
"""bring_jerrard_kps_S128c35.py -- Bring-Jerrard reduction of the level-5 wall quintic.
Exact: depressed form; principal form via quadratic Tschirnhaus (resolvent field Q(sqrt(delta))).
Numeric (50 digits): full Bring-Jerrard y^5+py+q, normalized one-parameter form, invariant hunt."""
import sys
from fractions import Fraction as F
from decimal import Decimal, getcontext
getcontext().prec=60
sys.stdout.reconfigure(line_buffering=True)

# f(x) = 4x^5 -170x^4 +4300x^3 -77680x^2 +930376x -5569395, monicize: divide by 4
f=[F(-5569395,4),F(930376,4),F(-77680,4),F(4300,4),F(-170,4),F(1)]
def polv(c,x):
    v=F(0)
    for co in reversed(c): v=v*x+co
    return v
# 1. depressed: x = y + 17/2
sh=F(17,2)
# compute coefficients of f(y+sh) via synthetic shifts
def shift(c,s):
    c=c[:]
    n=len(c)
    for i in range(n-1):
        for j in range(n-2,i-1,-1):
            c[j]+=s*c[j+1]
    return c
d=shift(f,sh)
print("depressed (monic, x=y+17/2): y^5 + %s y^3 + %s y^2 + %s y + %s"%(d[3],d[2],d[1],d[0]))
assert d[4]==0
# power sums of roots of depressed quintic (Newton's identities), e1=0
e=[F(1),F(0),d[3],-d[2],d[1],-d[0]]  # e0..e5 with signs: y^5 + c3 y^3 + c2 y^2 + c1 y + c0 = sum (-1)^k e_k y^{5-k}: e2=c3, e3=-c2, e4=c1, e5=-c0
p=[F(5)]  # p0 = 5
for k in range(1,9):
    # Newton: p_k = e1 p_{k-1} - e2 p_{k-2} + ... + (-1)^{k-1} k e_k
    s=F(0)
    for i in range(1,min(k,5)+1):
        term=e[i]*(p[k-i] if k-i>0 else F(5) if k-i==0 else F(0))
        s+= (-1)**(i-1)*term
    if k<=5: s+=0  # k e_k term folded: careful -- standard: p_k - e1 p_{k-1} + e2 p_{k-2} - ... + (-1)^{k-1} e_{k-1} p_1 + (-1)^k k e_k = 0
    p.append(s)
# redo Newton carefully
p=[F(5)]
for k in range(1,9):
    s=F(0)
    for i in range(1,k):
        if i<=5: s+=(-1)**(i-1)*e[i]*p[k-i]
    if k<=5: s+=(-1)**(k-1)*F(k)*e[k]
    p.append(s)
print("power sums p1..p6:",p[1:7])
# 2. principal form via z = y^2 + a y + b: need sum z_i = 0 and sum z_i^2 = 0
# sum z = p2 + a p1 + 5b = p2 + 5b (p1=0) -> b = -p2/5
b=-p[2]/5
# sum z^2 = sum (y^2+ay+b)^2 = p4 + 2a p3 + (a^2+2b) p2 + 2ab p1 + 5 b^2 = quadratic in a
A2=p[2]
A1=2*p[3]
A0=p[4]+2*b*p[2]+5*b*b
print("principal-form quadratic for a: (%s) a^2 + (%s) a + (%s) = 0"%(A2,A1,A0))
delta=A1*A1-4*A2*A0
print("delta =",delta)
num=delta.numerator; den=delta.denominator
def sqfree(n):
    s=1; d=2; nn=abs(n)
    while d*d<=nn:
        while nn%(d*d)==0: nn//=d*d; s*=d
        d+=1
    return nn,s
sf_num,sq_num=sqfree(num); sf_den,sq_den=sqfree(den)
print("delta squarefree kernel: %d / %d  (i.e. Q(sqrt(delta)) = Q(sqrt(%d)))"%(sf_num,sf_den,sf_num*sf_den))
k14={14:'14',183:'183',91:'91',7:'7',13:'13',2:'2'}
kern=sf_num*sf_den
for q,name in k14.items():
    if kern%q==0: print("   kernel divisible by",name)
# 3. numeric BJ: roots of depressed quintic to 50 digits (Newton from float seeds)
co=[Decimal(str(c.numerator))/Decimal(str(c.denominator)) for c in d]
import cmath
fl=[float(c) for c in d]
seeds=[]
ws=[complex(0.4,0.9)**i*20 for i in range(5)]
for _ in range(2000):
    new=[]
    for i,wi in enumerate(ws):
        num_=sum(fl[j]*wi**j for j in range(6)); den_=1.0
        for j,wj in enumerate(ws):
            if j!=i: den_*=(wi-wj)
        new.append(wi-num_/den_ if den_!=0 else wi)
    ws=new
# refine with complex Decimal Newton (use complex float 'refined' via mpc-free: use python complex with extra iterations - fine at ~1e-15; for 50 digits use Decimal for the real root only)
print("depressed roots (float):",["%.8f%+.8fi"%(w.real,w.imag) for w in ws])
# principal-form roots: z_i = y_i^2 + a y_i + b with a from the quadratic (pick + branch)
import math
aplus=(-float(A1)+math.copysign(1,1)*math.sqrt(abs(float(delta))))/(2*float(A2)) if float(delta)>=0 else None
if aplus is None:
    sq=cmath.sqrt(complex(float(delta)))
    a1=(-float(A1)+sq)/(2*float(A2)); branches=[a1,(-float(A1)-sq)/(2*float(A2))]
else:
    branches=[(-float(A1)+math.sqrt(float(delta)))/(2*float(A2)),(-float(A1)-math.sqrt(float(delta)))/(2*float(A2))]
print("a branches:",branches,"(delta sign %s)"%("+" if float(delta)>0 else "-"))
a=branches[0]; bb=float(b)
zs=[w*w+a*w+bb for w in ws]
# principal quintic coefficients from roots (elementary symmetric)
def esym(rs):
    es=[complex(1)]
    for r in rs:
        new=[complex(1)]*(len(es)+1)
        new[0]=complex(1)
        prev=es
        cur=[complex(0)]*(len(es)+1)
        for i,c in enumerate(prev): cur[i]+=c
        for i,c in enumerate(prev): cur[i+1]+=c*(-r)
        es=cur
    return es  # coeffs of prod (x - r), low->high? build: es holds poly coeffs highest first? let's do simple
def poly_from_roots(rs):
    c=[complex(1)]
    for r in rs:
        c=[0]+c
        c=[c[i]- r*(c[i+1] if i+1<len(c) else 0) for i in range(len(c))]
    return c  # low->high
pc=poly_from_roots(zs)
pc=[x/pc[-1] for x in pc]
print("principal quintic (numeric): z^5 + (%.6g%+.6gi) z^2 + (%.6g%+.6gi) z + (%.6g%+.6gi); z^4,z^3 coeffs: %.2e %.2e"%(
    pc[2].real,pc[2].imag,pc[1].real,pc[1].imag,pc[0].real,pc[0].imag,abs(pc[4]),abs(pc[3])))
# 4. BJ: quartic Tschirnhaus w = z^4 + c3 z^3 + c2 z^2 + c1 z + c0 killing w^4,w^3,w^2 -- solve via power sums of z
P=[complex(5)]
zz=zs
for k in range(1,21):
    P.append(sum(z**k for z in zz))
# want sum w = 0, sum w^2 = 0, sum w^3 = 0 with w_i = z_i^4 + c3 z_i^3 + c2 z_i^2 + c1 z_i + c0
# c0 fixed by sum w = 0 given others; sum w^2 = 0 and sum w^3 = 0: 2 eqs in c1,c2,c3 (one dof free: set c3 = 0)
# use c3 as free = 0; unknowns c1, c2; eqs are polynomial (degree 2 and 3) -- solve numerically by 2d Newton
def wsum(c1,c2,c3,k):
    ws_=[z**4 + c3*z**3 + c2*z**2 + c1*z for z in zz]
    m=sum(ws_)/5
    return sum((w-m)**k for w in ws_)
def eqs(v):
    c1,c2=v
    return [wsum(c1,c2,0,2), wsum(c1,c2,0,3)]
def jac(v,h=1e-6):
    c1,c2=v
    e0=eqs(v)
    J=[[0,0],[0,0]]
    for j,dv in enumerate([(h,0),(0,h)]):
        e1=eqs([c1+dv[0],c2+dv[1]])
        for i in range(2): J[i][j]=(e1[i]-e0[i])/h
    return e0,J
v=[complex(1,0.5),complex(1,-0.5)]
for it in range(200):
    e0,J=jac(v)
    det=J[0][0]*J[1][1]-J[0][1]*J[1][0]
    if abs(det)<1e-30: break
    dv0=( e0[0]*J[1][1]-e0[1]*J[0][1])/det
    dv1=(-e0[0]*J[1][0]+e0[1]*J[0][0])/det
    v=[v[0]-dv0,v[1]-dv1]
    if abs(dv0)+abs(dv1)<1e-13: break
c1,c2=v
wlist=[z**4 + c2*z**2 + c1*z for z in zz]
m=sum(wlist)/5
wlist=[w-m for w in wlist]
bj=poly_from_roots(wlist)
bj=[x/bj[-1] for x in bj]
print("BJ check: w^4..w^2 coeffs: %.2e %.2e %.2e"%(abs(bj[4]),abs(bj[3]),abs(bj[2])))
pco=bj[1]; qco=bj[0]
print("Bring-Jerrard: w^5 + p w + q with p = %.10g%+.10gi , q = %.10g%+.10gi"%(pco.real,pco.imag,qco.real,qco.imag))
# one-parameter normal form: w = s u with s^4 = -p (if we want u^5 - u - t) -> t = -q / s^5
s4=-pco
s=s4**0.25
t=-qco/ s**5
print("normal form u^5 - u - t: t = %.12g%+.12gi ; |t| = %.12g"%(t.real,t.imag,abs(t)))
for name,val in [("14",14),("183",183),("91",91),("183/14",183/14),("91/14",91/14),("1/91",1/91),("1/14",1/14),("35035/371293",35035/371293)]:
    r=abs(t)/val
    print("   |t| / %s = %.8f"%(name,r))
print("DONE")
