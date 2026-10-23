#!/usr/bin/env python3
"""Exact finite controls; no floating-point comparisons, no assert statements."""
from fractions import Fraction as F
from math import gcd
from pathlib import Path
from collections import defaultdict
import json

CHECKS=0

def ck(ok,label):
    global CHECKS
    CHECKS+=1
    if not ok: raise RuntimeError(label)

def vbase(k,b):
    j=0
    while k%b==0:k//=b;j+=1
    return j

def rotate(x,k,M):
    k%=M
    return ((x<<k)|(x>>(M-k)))&((1<<M)-1)

def combmask(b,n,M):
    ck(M%(2*b*n)==0,'unaligned comb grid')
    # Cell midpoint (2r+1)/(2M), strict condition by exact integers.
    z=0
    for r in range(M):
        v=((2*r+1)*n)%(2*M)
        if min(v,2*M-v)*b<M:z|=1<<r
    ck(z.bit_count()*b==M,'comb mass')
    return z

def components(x,M):
    if x==0:return 0
    if x==(1<<M)-1:return 1
    return (x & ~rotate(x,1,M)).bit_count()

def rows(b,J,M):
    return [rotate(combmask(b,b**j,M),s,M) for j in range(J) for s in range(M)]

def signature(x,rs):return tuple((x&r).bit_count() for r in rs)

def rank_mod(rs,M,p=1000003):
    a=[[int(bool(row>>j&1)) for j in range(M)] for row in rs]
    r=0
    for col in range(M):
        pivot=next((i for i in range(r,len(a)) if a[i][col]%p),None)
        if pivot is None:continue
        a[r],a[pivot]=a[pivot],a[r]
        inv=pow(a[r][col],-1,p)
        a[r]=[(x*inv)%p for x in a[r]]
        for i in range(r+1,len(a)):
            c=a[i][col]%p
            if c:a[i]=[(x-c*y)%p for x,y in zip(a[i],a[r])]
        r+=1
        if r==M:break
    return r

def matrix_controls():
    out=[]
    for b,J,mul in [(2,1,4),(2,2,2),(2,3,1),(3,1,2),(3,2,1),(4,1,2),(4,2,1),(7,1,2),(7,2,1),(7,2,2)]:
        q=b**J;M=2*q*mul;p=M//q;rs=rows(b,J,M)
        rank=rank_mod(rs,M)
        ck(rank==M-p+1,'rank formula')
        for a in range(p-1):
            u=sum(1<<(a+r*p) for r in range(q))
            v=sum(1<<(p-1+r*p) for r in range(q))
            for row in rs:ck((row&u).bit_count()==(row&v).bit_count(),'explicit kernel basis')
        out.append(dict(base=b,depth=J,cells=M,rows=len(rs),rank_mod_prime=rank,rational_rank=rank,kernel_dimension=p-1))
    return out

def exhaustive_controls():
    out=[]
    for b,J,M in [(2,1,8),(3,1,12),(4,1,8),(7,1,14),(2,2,16),(2,3,16)]:
        q=b**J;p=M//q;rs=rows(b,J,M);seen={};collisions=0;least=None
        for x in range(1<<M):
            sig=signature(x,rs)
            if sig not in seen:seen[sig]=x;continue
            y=seen[sig];collisions+=1
            ck(x.bit_count()==y.bit_count(),'collision mean')
            h=[((x>>i)&1)-((y>>i)&1) for i in range(M)]
            ck(all(h[i]==h[(i+p)%M] for i in range(M)),'collision periodicity')
            cx=components(x,M);cy=components(y,M)
            ck(cx>=q and cy>=q,'component threshold')
            least=min(cx,cy) if least is None else min(least,cx,cy)
        ck(least==q,'threshold attained in exhaustive grid')
        out.append(dict(base=b,depth=J,cells=M,binary_carriers=1<<M,response_classes=len(seen),nonrepresentative_collisions=collisions,min_collision_components=least))
    return out

def sharp_hostiles():
    out=[]
    for b,J in [(2,1),(2,3),(3,2),(4,2),(7,1),(7,2),(7,3)]:
        q=b**J;M=4*q
        K=sum(1<<(4*r) for r in range(q));L=rotate(K,2,M)
        ck(K&L==0,'disjoint sharp hostile')
        ck(components(K,M)==q and components(L,M)==q,'sharp components')
        for j in range(J):
            d=combmask(b,b**j,M)
            for s in range(M):
                a=(K&rotate(d,s,M)).bit_count();z=(L&rotate(d,s,M)).bit_count()
                ck(a==z and F(a,M)==F(1,4*b),'full phase response sharp hostile')
        # The finest original grid need not align the next probe; use 4*b*q.
        M2=4*b*q
        K2=sum(((1<<b)-1)<<(4*b*r) for r in range(q))
        L2=rotate(K2,2*b,M2);d=combmask(b,q,M2)
        k0=F((K2&d).bit_count(),M2);l0=F((L2&d).bit_count(),M2)
        ck(k0==F(1,2*b) and l0==0,'next probe distinguishes')
        out.append(dict(base=b,depth=J,components=q,all_old_phase_values=str(F(1,4*b)),next_probe_values=[str(k0),str(l0)],phase_knots=M))
    return out

def mod1(x):return x-x.numerator//x.denominator

def dist(x):
    x=mod1(x);return min(x,1-x)

def comb_intervals(b,n,s):
    result=[]
    r=F(1,2*b*n)
    for k in range(n):
        c=mod1(s+F(k,n));a=c-r;z=c+r
        if a<0:result.extend([(F(0),z),(a+1,F(1))])
        elif z>1:result.extend([(a,F(1)),(F(0),z-1)])
        else:result.append((a,z))
    return result

def overlap(K,D):return sum((max(F(0),min(z,d)-max(a,c)) for a,z in K for c,d in D),F(0))

def scalar_sample_hostile():
    bank=[(1,F(0)),(7,F(1,31)),(49,F(3,17)),(3,F(1,5)),(2,F(2,11))]
    walls={F(0),F(1),F(1,2),F(5,6)}
    Ds=[]
    for n,s in bank:
        D=comb_intervals(7,n,s);Ds.append(D)
        for a,z in D:walls.update([a,z])
    ws=sorted(walls);a,z=max(zip(ws,ws[1:]),key=lambda az:az[1]-az[0]);w=z-a
    K=[(a+w/8,a+3*w/8)];L=[(a+5*w/8,a+7*w/8)]
    answers=[]
    for D in Ds:
        k=overlap(K,D);l=overlap(L,D)
        ck(k==l and k in (0,w/4),'finite scalar hostile')
        answers.append(str(k))
    ck(K[0][1]<L[0][0],'distinct one-arcs')
    return dict(bank=[[n,str(s)] for n,s in bank],cell=[str(a),str(z)],K=[[str(x) for x in K[0]]],L=[[str(x) for x in L[0]]],responses=answers,common_measure=str(w/4))

def ap_skeleton_controls():
    out=[]
    for m in range(1,21):
        delta=F(1,m+1);walls={F(0),F(1)}
        for v in range(1,m+1):
            for k in range(v):
                walls.add(mod1((F(k)+delta)/v));walls.add(mod1((F(k)-delta)/v))
        ws=sorted(walls)
        safe=lambda t:all(dist(v*t)>=delta for v in range(1,m+1))
        actual=[t for t in ws[:-1] if safe(t)]
        expected=[F(a,m+1) for a in range(1,m+1) if gcd(a,m+1)==1]
        ck(actual==expected,'tight AP weak-safe skeleton')
        for a,z in zip(ws,ws[1:]):ck(not safe((a+z)/2),'AP has no positive interval')
        out.append(dict(speeds=m,safe_points=[str(x) for x in actual],wall_count=len(ws)-1))
    return out

# Integer cyclotomic arithmetic, independent of the convolution signature.
def trim(a):
    while len(a)>1 and a[-1]==0:a.pop()
    return a

def divide(a,b):
    a=a[:];q=[0]*max(1,len(a)-len(b)+1)
    while len(a)>=len(b):
        k=len(a)-len(b);c=a[-1]
        q[k]=c
        for j,z in enumerate(b):a[k+j]-=c*z
        trim(a)
        if a==[0]:break
    return trim(q),trim(a)

def cyclotomic(n,cache={1:[-1,1]}):
    if n not in cache:
        a=[-1]+[0]*(n-1)+[1]
        for d in range(1,n):
            if n%d==0:
                a,r=divide(a,cyclotomic(d));ck(r==[0],'cyclotomic division')
        cache[n]=a
    return cache[n]

def moment_controls():
    M=12;phi=cyclotomic(M);deg=len(phi)-1;out=[]
    rvec=[]
    for k in range(1,4):
        row=[]
        for i in range(M):
            power=(k*i)%M;a=[0]*power+[1]
            _,r=divide(a,phi);row.append(tuple(r+[0]*(deg-len(r))))
        rvec.append(row)
    for C in (1,2,3):
        seen=defaultdict(list)
        for x in range(1<<M):
            sig=[x.bit_count()]
            for k in range(C):sig.extend(sum(rvec[k][i][j] for i in range(M) if x>>i&1) for j in range(deg))
            seen[tuple(sig)].append(x)
        for group in seen.values():
            if len(group)>1:ck(all(components(x,M)>C for x in group),'finite low moment uniqueness')
        low=sum(components(x,M)<=C for x in range(1<<M))
        out.append(dict(cells=M,first_complex_moments=C,carriers=1<<M,classes=len(seen),low_component_carriers=low))
    return out

def valuation_bank_controls():
    out=[]
    for J in (1,2,3):
        q=7**J
        speeds=[1,2,3,7,14,21,49,98,147,343,686,1029,13,91,637]
        admitted=[n for n in speeds if vbase(n,7)<J]
        tested=0
        for n in admitted:
            for ell in range(-200,201):
                if ell==0:continue
                k=ell*q
                detected=(k%n==0 and (k//n)%7!=0)
                ck(not detected,'valuation-bounded bank kills q modes');tested+=1
        out.append(dict(depth=J,period_denominator=q,admitted_speeds=admitted,nonzero_mode_tests=tested))
    return out

def stability_controls():
    b,J,M=2,2,8;q=b**J;rs=rows(b,J,M)
    sigs=[signature(x,rs) for x in range(1<<M)]
    counts=defaultdict(int)
    for K in range(1<<M):
        c=components(K,M)
        if c>=q:continue
        for G in range(1<<M):
            d=(K^G).bit_count()
            z=max(abs(x-y) for x,y in zip(sigs[K],sigs[G]))
            if K in (0,(1<<M)-1):
                ck(d<=b*z,'empty/full mean stability')
            else:
                # epsilon^(2c+1) <= (2c+1)^(2c+1)*(11bc/7)*rho.
                ck(7*d**(2*c+1)<=(2*c+1)**(2*c+1)*11*b*c*z*M**(2*c),'rational uniform response stability')
            counts[str(c)]+=1
    return dict(base=b,depth=J,cells=M,exact_pairs_by_components=dict(counts),rational_multiplier='11*b*c/7',scope='Uniform phase errors are exact: every response is affine between the aligned grid knots.')

def main():
    result={'scope':'Exact finite-grid and rational-interval controls for an abstract observer theorem, not a native LRC census.',
            'matrix_controls':matrix_controls(),'binary_exhaustions':exhaustive_controls(),
            'sharp_periodic_hostiles':sharp_hostiles(),'finite_phase_hostile':scalar_sample_hostile(),
            'ap_boundary_controls':ap_skeleton_controls(),'low_moment_controls':moment_controls(),'valuation_bank_controls':valuation_bank_controls(),'stability_controls':stability_controls()}
    result['always_active_checks']=CHECKS
    txt=json.dumps(result,sort_keys=True,indent=2)+'\n'
    Path(__file__).with_name('receipt.json').write_text(txt)
    print(txt,end='')
if __name__=='__main__':main()
