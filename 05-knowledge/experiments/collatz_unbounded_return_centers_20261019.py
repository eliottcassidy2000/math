#!/usr/bin/env python3
"""Read-only exact resonant controls; no floating-point decisions or search."""
from fractions import Fraction as F
import hashlib,json

checks=0

def check(ok,label):
    global checks
    checks+=1
    if not ok:raise RuntimeError(label)

def coeff(bits):
    a,b,m=1,0,1
    for bit in bits:
        if bit=='1':a,b=3*a,3*b+m
        elif bit!='0':raise ValueError('binary word required')
        m*=2
    return a,b,m

def walk(x,bits):
    y=x
    for bit in bits:
        check(y%2==int(bit),'independent literal parity guard')
        y=(3*y+1)//2 if y%2 else y//2
        check(y>x,'no lower prefix descent')
    return y

def digest(n):
    return hashlib.sha256(n.to_bytes((n.bit_length()+7)//8,'big')).hexdigest()

P,Q,C=3**19,2**30,4912087879
D=P-Q;c=F(C,D)
lo='110110110110110101011011101010'
hi='001111111111111101010110000100'
check(coeff(lo)==(P,C,Q),'K lower coefficients')
check(coeff(hi)==(P,C+Q-P,Q),'K upper coefficients')
for t in range(1,31):
    a,b,m=coeff(lo[:t]);check(a>m,'K prefix slope')
rows=[];prev_gap=None;prev_center=None
for L,j in [(1,0),(4,1),(11,3),(40,11),(69,19),(3341,920),(6613,1821)]:
    pl,ql=P**L,Q**L
    cl=C*(pl-ql)//D
    sj,tj=3**j,4**j
    S,M=pl*sj,ql*tj
    H=sj*cl+ql*(tj-sj)
    Hu=sj*(cl+ql-pl)+ql*2*(tj-sj)
    check(3*S<4*M and S>M,'exact floor j: 1<A<4/3')
    A=F(S,M);gap=A-1;e=F(H,M);center=F(H,S-M)
    check(e==c*A+1-(c+1)*F(3,4)**j,'normalized intercept identity')
    check(Hu==H+M-S and Hu>0,'independent upper composition')
    rho=(-H*pow(S,-1,M))%M
    ru=(-Hu*pow(S,-1,M))%M
    check(0<rho<M-1 and ru==rho+1,'joint nonwrapping source cylinder')
    check((S*rho+H)%M==0,'exact lower endpoint')
    check((S*(rho+1)+Hu)//M==(S*rho+H)//M+1,'paired endpoint')
    if prev_gap is not None:
        check(gap<prev_gap,'listed resonance gaps decrease')
        check(center>prev_center,'listed centers increase')
    prev_gap,prev_center=gap,center
    # Every F block-boundary slope is >= final A>1; odd prefixes multiply by3/2.
    check(A>1 and A*F(3,2)>1,'F prefix envelope')
    if L<=69:
        lower=lo*L+'10'*j;upper=hi*L+'01'*j
        check(coeff(lower)==(S,H,M),'independent full lower word')
        check(coeff(upper)==(S,Hu,M),'independent full upper word')
        y=walk(rho,lower)
        check(y==(S*rho+H)//M,'literal lower endpoint')
        # Upper prefixes may descend: only check their parity and final adjacency.
        z=rho+1
        for bit in upper:
            check(z%2==int(bit),'upper literal parity guard')
            z=(3*z+1)//2 if z%2 else z//2
        check(z==y+1,'literal paired adjacent endpoint')
    scale=10**12
    floor_gap=(gap.numerator*scale)//gap.denominator
    floor_center=center.numerator//center.denominator
    check(F(floor_gap,scale)<=gap<F(floor_gap+1,scale),'rational slope-gap enclosure')
    check(floor_center<=center<floor_center+1,'rational center enclosure')
    rows.append(dict(L=L,j=j,steps=30*L+2*j,gap_lower_numerator=floor_gap,
       gap_upper_numerator=floor_gap+1,gap_denominator=scale,
       center_floor=floor_center,source_bit_length=rho.bit_length(),source_sha256=digest(rho)))
print(json.dumps(dict(status='PROVED unbounded compressed-return centers; no source-entry claim',
 checks=checks,controls=rows),indent=2,sort_keys=True))
