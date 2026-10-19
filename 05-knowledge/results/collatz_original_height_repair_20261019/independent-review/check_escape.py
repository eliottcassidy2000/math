from fractions import Fraction as F
from collections import Counter
import json

K='110110110110110101011011101010'
Ku='001111111111111101010110000100'
Fword='10'; Fu='01'; H='110101000'; Hu='001101001'
V=K*4+Fword; Vu=Ku*4+Fu

def T(n): return (3*n+1)//2 if n&1 else n//2

def follow(n,L):
    bits=[]
    for _ in range(L):
        bits.append(str(n&1)); n=T(n)
    return ''.join(bits),n

# Coefficients from the nonrecursive odd-position sum, independently of
# the writer's branch recurrence.
def coeff(w):
    ones=w.count('1')
    b=sum((1<<j)*3**w[j+1:].count('1') for j,v in enumerate(w) if v=='1')
    return 3**ones,b,1<<len(w)

# Recover cylinder by binary lifting and direct parity checks, without
# the writer's affine modular-inverse cylinder formula.
def lift_cylinder(w):
    r=0
    for j,bit in enumerate(w):
        n=r
        for _ in range(j): n=T(n)
        if n%2 != int(bit): r+=1<<j
        if follow(r,j+1)[0] != w[:j+1]: raise ValueError('binary lift')
    return r

r4=lift_cylinder(K*4); M4=1<<120
if r4!=588943162011782273059357097228697595: raise ValueError('r4')
r1=lift_cylinder(K)
centers=[]
for name,w,u in [('K',K,Ku),('V',V,Vu),('F',Fword,Fu),('H',H,Hu)]:
    a,b,m=coeff(w); au,bu,mu=coeff(u)
    if (au,bu,mu)!=(a,b+m-a,m): raise ValueError('companion affine '+name)
    r=lift_cylinder(w); ru=lift_cylinder(u)
    if ru!=r+1: raise ValueError('companion guard '+name)
    if name in ('K','V'):
        centers.append(F(b,a-m))
        for j in range(1,len(w)+1):
            aa,bb,mm=coeff(w[:j])
            if not aa>mm: raise ValueError('expansion prefix')
    else:
        A=F(a,m); B=F(b,m)
        if not 0<A<1 or A+B>1: raise ValueError('reset decreasing')
        for j in range(len(w)):
            aa,bb,mm=coeff(w[:j]); AA=F(aa,mm); BB=F(bb,mm)
            if not AA>A or AA+BB<A+B: raise ValueError('reset prefix endpoint')

C=max(centers)
if not F(533)<C+1<F(534): raise ValueError('carry bound')
for ell in (1,2,3):
    for reset in (Fword,H):
        a,b,m=coeff(K*ell+reset)
        if not (a<m and a*r1+b<m*r1): raise ValueError('small expansion descent')

rows=[]
for d in range(1,534):
    n=r4+d; bits=''; odds=0
    for j in range(1,121):
        bit=n&1; bits+=str(bit); odds+=bit; n=T(n)
        if 3**odds < 1<<j: break
    else: raise ValueError('no contracting prefix')
    a,b,m=coeff(bits)
    if (a*(r4+d)+b)%m: raise ValueError('guard')
    threshold=F(a*d+b,m-a)
    if j>92 or threshold>=9000 or not m<=M4: raise ValueError('escape')
    if not n<r4: raise ValueError('anchor replay')
    rows.append((d,j,str(threshold),bits))

hist=Counter(j for _,j,_,_ in rows)
maxth=max((F(t),d) for d,_,t,_ in rows)
longest=[d for d,j,_,_ in rows if j==92]
print(json.dumps({'status':'PASS','method':'closed-form odd-position carry; binary-lift cylinder reconstruction; exact 533 offset affine certificates', 'forced_source':r4,'forced_modulus':M4,'carry_bound':str(C+1),'max_steps':max(hist),'longest_offsets':longest,'maximum_threshold':str(maxth[0]),'maximum_threshold_offset':maxth[1],'depth_histogram':dict(sorted(hist.items()))},indent=2))
