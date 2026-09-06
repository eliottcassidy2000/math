#!/usr/bin/env python3
"""Exact bounded hostiles to over-abstracting A2B3 signed transport.

Standard-library Fraction arithmetic only. No imported repository producers,
root finders, floating-point signs, or optimized-away assert statements.
The only universe is two declared hostile carriers and two actual controls.
"""
from fractions import Fraction as Q
from math import comb, lcm
import hashlib
import json

GATES = 0

def need(test, label):
    global GATES
    GATES += 1
    if not test:
        raise RuntimeError(label)

def trim(p):
    p = list(map(Q, p))
    while len(p) > 1 and not p[-1]:
        p.pop()
    return p

def add(a, b):
    c = [Q(0)] * max(len(a), len(b))
    for i, v in enumerate(a): c[i] += v
    for i, v in enumerate(b): c[i] += v
    return trim(c)

def scale(a, s): return trim([s*v for v in a])

def conv(a, b):
    c = [Q(0)] * (len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b): c[i+j] += x*y
    return trim(c)

def at(p, x):
    r = Q(0)
    for c in reversed(p): r = r*x+c
    return r

def shift(p, x):
    c = [Q(0)] * len(p)
    for i, a in enumerate(p):
        for j in range(i+1): c[j] += a*comb(i,j)*x**(i-j)
    return trim(c)

def rem(a, b):
    a = trim(a); b = trim(b)
    while len(a) >= len(b) and any(a):
        d = len(a)-len(b); t = a[-1]/b[-1]
        for j, bj in enumerate(b): a[d+j] -= t*bj
        a = trim(a)
    return a

def integral_form(p):
    den = lcm(*(v.denominator for v in p))
    nums = [int(v*den) for v in p]
    return nums, den

def fmt(p): return [str(v) for v in p]

def bracket(p, a, b, label):
    need(a < b and at(p,a)*at(p,b) < 0, label)

def kernels(roots):
    q = len(roots)
    b = [Q(1)]
    for a in roots: b = conv(b, [-a,1])
    c = [Q(3*(i+1), q+2+2*i)*b[i+1] for i in range(q)]
    d = [Q(2+3*i, q+1+2*i)*c[i] for i in range(q)]
    return b,c,d

def even_kernel(p, s):
    # p(v) has degree d. K_s(u)=s^d p(u^2/s).
    d = len(p)-1
    a = [Q(0)] * (2*d+1)
    for i,v in enumerate(p): a[2*i] = v*s**(d-i)
    return a

def carrier(p, s, g):
    return conv(even_kernel(p,s), [Q(comb(g,j)) for j in range(g+1)])

def coeff(p, j): return p[j] if 0 <= j < len(p) else Q(0)

def positive_g(p):
    d = len(p)-1
    return [p[d-j]*(-1)**j for j in range(d+1)]

def raw_paths(b,c,d,h,x):
    # Polynomial in positive s=-t; every returned response is multiplied by s.
    q=h+x; g=x+3*h+1
    B,C,D = map(positive_g,(b,c,d))
    need(all(v>0 for a in (B,C,D) for v in a), 'full beta coefficients positive')
    BB,CD=conv(B,B),conv(C,D)
    O=[Q(comb(g,2*j+1)) for j in range(g//2+1) if 2*j+1<=g]
    E=[Q(comb(g,2*j)) for j in range(g//2+1)]
    OO,EE=conv(O,O),conv(E,E)
    # Fully retain raw B indices -x..h and C,D indices -x..h-1.
    # alpha_double includes j=-1: binom(2g,2+2j).
    V=[Q(0)]*(2*h+2); W=[Q(0)]*(2*h+2); R=[Q(0)]*(2*h+2)
    for j in range(-1,2*h+1):
        bj=coeff(BB,j+2*x)
        cdj=coeff(CD,j-1+2*x)
        aj=Q(comb(2*g,2+2*j)) if 0<=2+2*j<=2*g else Q(0)
        need(aj == coeff(OO,j)+coeff(EE,j+1), 'complete alpha parity')
        sign = -1 if j % 2 else 1
        V[j+1]=sign*coeff(OO,j)*bj
        W[j+1]=sign*aj*bj
        R[j+1]=sign*aj*2*cdj
    return trim(V),trim(W),trim(R)

# Main hostile: q6,h3,x3,g13 with distinct integer squared-root parameters.
roots=list(map(Q,[1,2,3,5,9,12])); q=6; h=3; x=3; g=13; k=7
b,c,d=kernels(roots)
need(b == list(map(Q,[3240,-7218,5739,-2110,380,-32,1])), 'factored B')
need(c == [Q(-10827,4),Q(17217,5),Q(-3165,2),Q(2280,7),Q(-30),Q(1)], 'C bank')
need(d == [Q(-10827,14),Q(1913),Q(-12660,11),Q(25080,91),Q(-28),Q(1)], 'D bank')
C_intervals=[(Q(19,10),Q(2)),(Q(28,10),Q(29,10)),(Q(47,10),Q(48,10)),(Q(85,10),Q(86,10)),(Q(118,10),Q(119,10))]
D_intervals=[(Q(5,10),Q(6,10)),(Q(31,10),Q(32,10)),(Q(44,10),Q(45,10)),(Q(81,10),Q(82,10)),(Q(116,10),Q(117,10))]
for label,p,intervals in [('C',c,C_intervals),('D',d,D_intervals)]:
    need(len(intervals)==len(p)-1, label+' degree exhaustion')
    for i,(a,z) in enumerate(intervals):
        bracket(p,a,z,label+' root bracket')
        need(a>0 and (i==0 or intervals[i-1][1]<a),label+' positive disjoint')
for i,(a,z) in enumerate(C_intervals):
    need(roots[i]<a and z<=roots[i+1] and not at(c,z)==0,'C interlaces B')
need(sum(3<a<z<5 for a,z in D_intervals)==2,'D has two roots in (3,5)')
need(all(z<=2 or a>=3 for a,z in D_intervals),'D has no root in (2,3)')
need(at(c,1)*at(d,1)<0 and at(c,3)*at(d,3)<0,'opposite residues at B roots')

# Exact all-s identities, checked coefficientwise (no interpolation).
# K_B'=(2/3)u[(q+2)K_C+u K_C']; (q+1)K_D+u K_D'=2K_C+(3/2)u K_C'.
for i in range(q):
    need(2*(i+1)*b[i+1] == Q(2,3)*(q+2+2*i)*c[i], 'first Euler identity')
    need((q+1+2*i)*d[i] == (2+3*i)*c[i], 'second Euler identity')
P=[b[h-j]*comb(g,k-2*(h-j)) for j in range(h+1)]
p=list(map(Q,[-1055,63129,-357291,213840]))
need(P==scale(p,26),'same literal first coefficient')
for a,z in [(Q(0),Q(1,10)),(Q(1,10),Q(1)),(Q(1),Q(3,2))]:
    bracket(p,a,z,'first-root degree exhaustion')

cd,bb=conv(c,d),conv(b,b)
T=[Q(0)]*7
for i in range(7): T[6-i]=cd[i]*comb(26,12-2*i)
sW=[Q(0)]*8
for i in range(8): sW[7-i]=bb[i]*comb(26,14-2*i)
sR=scale(T,-2)
sV,rawW,rawR=raw_paths(b,c,d,h,x)
need(sW==rawW and sR==rawR,'literal carriers versus full Laurent Hadamard path')
sQ=add(sW,sR)
expected_T=([-12900924131028983563588405,762374596098079635101930499,-3799751043113696471769929961],920602755072000)
need(integral_form(rem(T,p))==expected_T,'decisive exact remainder')
remainders={}
for name,poly in [('T',T),('sW',sW),('sQ',sQ),('sV',sV),('s(Q-V)',add(sQ,scale(sV,-1)))]:
    rr=rem(poly,p)
    need(all(v<0 for v in shift(rr,Q(1))),name+' remainder negative on s>=1')
    remainders[name]=integral_form(rr)
need(coeff(sW,0)!=0 and coeff(sR,0)!=0,'load-bearing lower carry retained')
for phase in [Q(1),Q(3,2),Q(2)]:
    HB,HC,HD=[carrier(a,phase,g) for a in (b,c,d)]
    need(coeff(HB,k)==phase**x*at(P,phase),'direct first extraction')
    need(coeff(conv(HC,HD),4*h)==phase**4*at(T,phase),'direct mixed extraction')
    need(coeff(conv(HB,HB),2*k)==phase**(2*x-1)*at(sW,phase),'direct hit extraction')
    KB,KC=[even_kernel(a,phase) for a in (b,c)]
    lowered=conv([Q(comb(g-1,j)) for j in range(g)],add(KB,scale([0]+KC,Q(2,3))))
    need(k*coeff(HB,k)==g*coeff(lowered,k-1),'same-zero coefficient lowering')

# Actual positive control: h1,x1,g5. Full beta F6,F5,F4.
actual_b=[Q(1),Q(-4),Q(1)]
actual_c=[Q(-3),Q(1)]
actual_d=[Q(-2),Q(1)]
actualV,actualW,actualR=raw_paths(actual_b,actual_c,actual_d,1,1)
HB,HC,HD=[carrier(a,Q(2),5) for a in (actual_b,actual_c,actual_d)]
need(coeff(HB,3)==0,'actual (-9,1,6) first root')
need(coeff(conv(HC,HD),4)>0,'actual control mixed positive')
need(at(actualW,2)<0 and at(actualR,2)<0,'actual control both groups negative')
actual_control={'s':2,'sV':str(at(actualV,2)),'sW':str(at(actualW,2)),'sR':str(at(actualR,2)),'mixed':str(coeff(conv(HC,HD),4))}

# Second hostile: strict generic interlacing and the origin condition are insufficient.
Fpoly=[Q(1)]
for r in [Q(-3),Q(-2),Q(1),Q(1,4)]: Fpoly=conv(Fpoly,[1,r])
J=conv([0,1],conv([1,Q(-5,2)],[1,Q(3,4)]))
need(coeff(Fpoly,2)==0,'generic first zero retained')
need(coeff(conv(J,J),4)==Q(-11,16),'strict-interlacer square hostile')
Froots=[Q(-4),Q(-1),Q(1,3),Q(1,2)]
Jroots=[Q(-4,3),Q(0),Q(2,5)]
for i,r in enumerate(Froots):need(at(Fpoly,r)==0,'F exact root')
for i,r in enumerate(Jroots):
    need(at(J,r)==0 and Froots[i]<r<Froots[i+1],'strict J interlacing')
need(coeff(conv(Fpoly,Fpoly),4)==Q(-351,8),'pencil coefficient bank')
# F +/- aJ have roots by four disjoint interlacing sign changes for every nonzero a;
# the analytic proof is in the note. Two rational named controls independently
# verify this root count by sign brackets between J roots and suitable endpoints.
for a in [Q(1),Q(10)]:
    for sign in [-1,1]:
        pencil=add(Fpoly,scale(J,sign*a))
        intervals=[(Q(-100),Jroots[0]),(Jroots[0],Jroots[1]),(Jroots[1],Jroots[2]),(Jroots[2],Q(100))]
        for lo,hi in intervals:bracket(pencil,lo,hi,'pencil exact RR control')

# One actual same-degree beta triple distinguishes the path kernel from the hostile.
def actual_kernel(n):
    q=n//3
    return [Q(comb(n-2*j,j))*(-1)**j for j in reversed(range(q+1))]
AB,AC,AD=[actual_kernel(n) for n in [18,17,16]]
need(AB!=b and AC!=c and AD!=d,'hostile is not actual composition row')
for i in range(6):
    need(2*(i+1)*AB[i+1]==Q(2,3)*(8+2*i)*AC[i],'actual first coupling')
    need((7+2*i)*AD[i]==(2+3*i)*AC[i],'actual second coupling')
# The actual beta midpoint identity, in the squared-root variable, is literal.
AE=actual_kernel(36)
need(AE==add(conv(AB,AB),scale([0]+conv(AC,AD),-2)),'actual doubled-beta identity')

bank={'roots':fmt(roots),'B':fmt(b),'C':fmt(c),'D':fmt(d),'p':fmt(p),'remainders':remainders,'actual_control':actual_control,'F':fmt(Fpoly),'J':fmt(J)}
semantic=hashlib.sha256(json.dumps(bank,sort_keys=True,separators=(',',':')).encode()).hexdigest()
print('STATUS: EXACT REFUTATIONS; actual all-height transport remains OPEN')
print('UNIVERSE: one q6/h3 Euler-coupled PF hostile; one strict quartic interlacer hostile; two actual controls')
print('B,C,D ascending:',fmt(b),fmt(c),fmt(d))
print('first p ascending:',fmt(p),'largest root in (1,3/2)')
print('C root brackets:',[(str(a),str(z)) for a,z in C_intervals])
print('D root brackets:',[(str(a),str(z)) for a,z in D_intervals])
for name,data in remainders.items():print(name,'remainder numerator ascending / denominator:',data)
print('MAIN HOSTILE: mixed<0, skip>0, but formal Q<V<0 at the largest first root')
print('ACTUAL CONTROL:',actual_control)
print('STRICT INTERLACER: F=',fmt(Fpoly),'J=',fmt(J),'[u4]J2=-11/16')
print('PENCIL: [u4](F2-lambda J2)=-351/8+11lambda/16; RR for every lambda>=0')
print('LOWERING: k[u^k](1+u)^g K_B = g[u^(k-1)](1+u)^(g-1)(K_B+(2/3)uK_C), g=q+k')
print('gates:',GATES)
print('semantic_sha256:',semantic)
