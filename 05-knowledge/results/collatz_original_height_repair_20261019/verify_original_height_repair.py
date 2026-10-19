#!/usr/bin/env python3
"""Exact finite certificate for an all-lift original-height repair theorem.
No orbit census: 533 bounded offsets are certified by affine inequalities.
Prints JSON; no file writes, network, or assert statements.
"""
from fractions import Fraction as F
from collections import Counter
import json

checks=0
def check(v,msg):
    global checks
    checks+=1
    if not v: raise RuntimeError(msg)

def branch(w):
    a,b,m=1,0,1
    for bit in w:
        if bit=='1': a,b=3*a,3*b+m
        elif bit!='0': raise ValueError(bit)
        m*=2
    return a,b,m

def cylinder(w):
    a,b,m=branch(w)
    return (-b*pow(a,-1,m))%m

def T(n): return (3*n+1)//2 if n%2 else n//2

def word(n,L):
    out=''
    for _ in range(L):
        out+=str(n%2); n=T(n)
    return out,n

def ff(f): return f'{f.numerator}/{f.denominator}'

K='110110110110110101011011101010'
Ku='001111111111111101010110000100'
R='10'; Ru='01'
H='110101000'; Hu='001101001'
V=K*4+R; Vu=Ku*4+Ru
P,C0,Q=branch(K); D=P-Q
r=cylinder(K)
check((P,C0,Q,r)==(1162261467,4912087879,1073741824,815464443),'K data')
for name,w,u in [('K',K,Ku),('F',R,Ru),('H',H,Hu),('V',V,Vu)]:
    a,b,m=branch(w); au,bu,mu=branch(u)
    check((au,mu,bu)==(a,m,b+m-a),name+' independent upper branch')
    check(cylinder(u)==cylinder(w)+1,name+' independent upper guard')
    check(cylinder(w)<m-1,name+' no residue wrap')
for name,w in [('K',K),('V',V)]:
    for j in range(1,len(w)+1):
        a,b,m=branch(w[:j]); check(a>m,name+' expanding prefix')
for name,w in [('F',R),('H',H)]:
    a,b,m=branch(w)
    check(F(a,m)<1 and F(a+b,m)<=1,name+' reset nonincrease above one')
    for j in range(len(w)):
        aj,bj,mj=branch(w[:j])
        check(F(aj,mj)>F(a,m),name+' proper-prefix slope above endpoint')
        check(F(aj+bj,mj)>=F(a+b,m),name+' proper prefix at one above endpoint')
# Two independent constructions of the forced prefix cylinder.
a4,b4,m4=branch(K*4)
r4=cylinder(K*4)
check(r4==588943162011782273059357097228697595,'forced source')
check(m4==2**120,'forced modulus')
check(b4==C0*(P**4-Q**4)//D,'repetition carry')
check((D*r4+C0)%Q**4==0,'center valuation guard')
check(0<r4<m4 and r4>9000,'source threshold')
# Small expansion cases are already descending at their first reset.
for ell in (1,2,3):
    a,b,m=branch(K*ell+R)
    check(a<m and a*r+b<m*r,'small expansion F descent')
aH,bH,mH=branch(H); aF,bF,mF=branch(R)
check(F(aH,mH)<F(aF,mF) and F(bH,mH)<F(bF,mF),'H below F')
# Bounded carry: inherited bank proof, with constants authenticated here.
aV,bV,mV=branch(V)
center=F(bV,aV-mV)
check(center>F(C0,D),'V is maximal center')
check(F(bH,mH-aH)<1,'reset fixed points bounded by one')
E=center+1
check(533<E<534,'integer displacement bound')
# Discover first contracting prefix on the anchor, then certify every lift
# independently by the full parity-word congruence and an affine inequality.
rows=[]; classes={}; max_threshold=F(0)
for d in range(1,534):
    n=r4+d; a,b,m=1,0,1; w=''
    for j in range(1,121):
        bit=n%2; w+=str(bit)
        if bit: a,b=3*a,3*b+m
        m*=2; n=T(n)
        if a<m: break
    else: raise RuntimeError('uncertified offset '+str(d))
    check(j<=92,'continuation deadline')
    check(branch(w)==(a,b,m),'independent word coefficients')
    check((r4+d-cylinder(w))%m==0,'all-lift source guard')
    check(m4%m==0,'prefix modulus divides inherited modulus')
    threshold=F(a*d+b,m-a)
    check(threshold<9000,'uniform exact affine threshold')
    max_threshold=max(max_threshold,threshold)
    # Actual replay is a separate positive control, never the all-lift proof.
    for t in (0,1,17):
        x=r4+t*m4
        wr,y=word(x+d,j)
        check(wr==w and y==(a*(x+d)+b)//m and y<x,'literal lifted replay')
    key=(j,d%m)
    if key not in classes:
        classes[key]={'time':j,'offset_residue':d%m,'offset_modulus':m,
                      'parity_word':w,'multiplier_numerator':a,'carry':b,
                      'offset_count':0,'min_offset':d,'max_offset':d}
    row=classes[key]; row['offset_count']+=1; row['max_offset']=d
    rows.append((d,j,ff(threshold)))
check(len(rows)==533,'complete positive offset range')
check(max(j for _,j,_ in rows)==92,'sharp certificate depth')
check(max_threshold==F(116927,13),'maximum threshold')
check(len(classes)==54,'compressed class count')
check(sum(z['offset_count'] for z in classes.values())==533,'class coverage')
# Ensure the compressed residue description has exactly the same coverage.
for d in range(1,534):
    hits=[v for v in classes.values() if d%v['offset_modulus']==v['offset_residue']]
    check(len(hits)==1,'disjoint exhaustive offset classes')
    row=hits[0]
    check(row['min_offset']<=d<=row['max_offset'],'class endpoint bounds')
check(word(r4,120)[0]==K*4,'zero offset expanding hostile')
check(word(r4+4,91)[1]>r4 and word(r4+4,92)[1]<r4,'longest repair hostile')
# A known small periodic adjacent return explains why zero displacement
# cannot be declared impossible from the affine-pair identity alone.
check(T(T(1))==1 and T(T(2))==2,'F equality control')
print(json.dumps({'status':'PASS','explicit_checks':checks,
 'scope':'533 offset certificates, each valid on all inherited positive source lifts; no global source coverage',
 'bank_center':ff(center),'carry_upper_bound':ff(E),'max_positive_displacement':533,
 'forced_source':r4,'forced_modulus':m4,'repair_max_steps':92,
 'max_affine_threshold':ff(max_threshold),'compressed_classes':len(classes),
 'long_repairs':[{'offset':d,'steps':j} for d,j,_ in rows if j>=30],
 'offset_certificate':sorted(classes.values(),key=lambda z:(z['time'],z['offset_residue']))},indent=2,sort_keys=True))
