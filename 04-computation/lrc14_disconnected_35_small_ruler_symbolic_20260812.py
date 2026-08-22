#!/usr/bin/env python3
"""Exact all-dilation 3:5 floor on the 29 small LRC14 rulers.

Universe: the 2,530 ordered upper-median contexts with L<4592 from the
canonical THM-3350 atlas.  For every such (L,j,e,f), prove exactly for every
g>=1 that I(L,j;e,3g;f,5g)>=1/294.  The proof uses affine branch stabilization
on each residue class modulo |5e-3f| (period one if zero), followed by exact
quadratic integer minimization.  Body-safe endpoint inequalities make every
physical tooth full, so the digital tent is literal rather than an unbounded
proxy.
"""
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction as F
from hashlib import sha256
from importlib.util import module_from_spec,spec_from_file_location
from pathlib import Path
import argparse

ROOT=Path(__file__).resolve().parents[1]
TAIL=ROOT/'04-computation/lrc14_connected_low_uniform_high_forest_tail_thm3350.py'
EXPECTED_TAIL='78daaf73966d283c0c0bafa1c0975684e6167d2ef6375a3abeece4e00cdc87f9'
P,Q=3,5

def require(x,d):
    if not x:raise RuntimeError(d)
def load(n,p):
    s=spec_from_file_location(n,p);m=module_from_spec(s);s.loader.exec_module(m);return m
def filehash(p):return sha256(p.read_bytes().replace(b'\r\n',b'\n')).hexdigest()
def ceildiv(a,b):return -((-a)//b)

def floor_moments(n,m,a,b):
    if n==0:return 0,0,0
    s1=n*(n-1)//2;s2=n*(n-1)*(2*n-1)//6
    qa,a0=divmod(a,m);qb,b0=divmod(b,m)
    base0=qa*s1+qb*n;base1=qa*s2+qb*s1;base2=qa*qa*s2+2*qa*qb*s1+qb*qb*n
    if a0==0:return base0,base1,base2
    height=(a0*(n-1)+b0)//m
    if height==0:return base0,base1,base2
    u0,u1,u2=floor_moments(height,a0,m,m-b0+a0-1)
    r0=n*height-u0;r1=height*s1-(u2-u0)//2;r2=n*height*height-2*u1-u0
    return base0+r0,base1+r1,base2+2*qa*r1+2*qb*r0+r2

def residue_prefix(n,m,a,b,threshold,base):
    shifted=floor_moments(n,m,a,b+m-threshold)
    d0=shifted[0]-base[0];d1=shifted[1]-base[1]
    y0d=(shifted[2]-base[2]-d0)//2
    high_sum=a*d1+b*d0-m*y0d
    total=a*n*(n-1)//2+b*n-m*base[0]
    return n-d0,total-high_sum

def triangle_sum(L,n,m,a,b,peak,base,total):
    if peak<=0:return 0
    radius=(peak-1)//L
    require(2*radius<m,('overlapping triangle tails',L,n,m,peak))
    lc,ls=residue_prefix(n,m,a,b,radius+1,base)
    bc,bs=residue_prefix(n,m,a,b,m-radius,base)
    hc,hs=n-bc,total-bs
    return peak*lc-L*ls+(peak-L*m)*hc+L*hs

def fast_numden(L,cell,e,p,f,q):
    require(p<=q,(p,q));z,w=L*p-e,L*q-f
    r,s=e*cell%L,f*cell%L;det=r*w-s*z
    require(det%L==0,('nonintegral phase',L,cell,e,p,f,q,det))
    b=(det//L)%z;a=w%z;base=floor_moments(p,z,a,b)
    total=a*p*(p-1)//2+b*p-z*base[0];unit=L//14
    num=triangle_sum(L,p,z,a,b,unit*(z+w),base,total)
    num-=triangle_sum(L,p,z,a,b,unit*(w-z),base,total)
    return num,z*w

def term_candidates(L,cell,e,f,residue_index,shift_index,kernel,g):
    z,w=L*g*P-e,L*g*Q-f;r0,s0=e*cell%L,f*cell%L
    base=(r0*w-s0*z)//L;cross=P*w-Q*z
    lo=max(0,ceildiv(-shift_index,Q));hi=min(g-1,(g*Q-1-shift_index)//Q)
    affine=base+residue_index*w-shift_index*z
    peak=(L//14)*((z+w) if kernel==0 else (w-z));radius=(peak-1)//L
    if cross<0:affine,cross=-affine,-cross
    if cross:
        plo=ceildiv(-radius-affine,cross);phi=(radius-affine)//cross
        neg=ceildiv(-affine,cross)-1
        return (lo,plo,hi,phi,neg,neg+1),cross
    return (lo,hi,affine,peak-L*abs(affine)),0

def stable_term(L,cell,e,f,ri,si,kernel,residue,period,h):
    here,cross=term_candidates(L,cell,e,f,ri,si,kernel,residue+period*h)
    nxt,cross2=term_candidates(L,cell,e,f,ri,si,kernel,residue+period*(h+1))
    require(cross==cross2,('cross instability',L,cell,e,f,residue,period,h))
    slopes=tuple(y-x for x,y in zip(here,nxt))
    nxt2,_=term_candidates(L,cell,e,f,ri,si,kernel,residue+period*(h+2))
    require(tuple(y-x for x,y in zip(nxt,nxt2))==slopes,('endpoint nonaffinity',L,cell,e,f,residue,period,h))
    if cross:
        lo,plo,hi,phi,neg,_=here;slo,splo,shi,sphi,sneg,_=slopes
        left,sleft=(plo,splo) if plo>lo or (plo==lo and splo>slo) else (lo,slo)
        if not(sleft>=slo and sleft>=splo):return False
        right,sright=(phi,sphi) if phi<hi or (phi==hi and sphi<shi) else (hi,shi)
        if not(sright<=shi and sright<=sphi):return False
        if left>right:return sleft>=sright
        if neg<left:return sneg<=sleft
        if neg>=right:return sneg>=sright
        return sleft<=sneg<=sright
    lo,hi,affine,value=here;slo,shi,saffine,svalue=slopes
    if hi<lo or shi<slo:return False
    if affine>0 and saffine<0:return False
    if affine<0 and saffine>0:return False
    if affine==0 and saffine!=0:return False
    if value>0:return svalue>=0
    if value<0:return svalue<=0
    return svalue==0

def ray_is_stable(L,cell,e,f,residue,period,h):
    return all(stable_term(L,cell,e,f,ri,si,k,residue,period,h)
               for ri in range(P) for si in range(-2,Q+3) for k in (0,1))

def margin(L,cell,e,f,g):
    num,den=fast_numden(L,cell,e,g*P,f,g*Q)
    return 294*num-den

def quadratic_minimum(v0,v1,v2):
    second=v2-2*v1+v0;require(second%2==0,('nonintegral quadratic',v0,v1,v2))
    a=second//2;b=v1-v0-a;require(a>=0,('concave tail',v0,v1,v2))
    candidates={0}
    if a==0:require(b>=0,('decreasing linear tail',v0,v1,v2))
    elif b<0:
        c=(-b)//(2*a);candidates.update((max(0,c-1),c,c+1,c+2))
    return min(a*n*n+b*n+v0 for n in candidates),(a,b,v0)

def prove_context(task):
    L,cell,e,f=task;unit=L//14;r=e*cell%L;s=f*cell%L
    require(unit<=r and r+e<=L-unit,(task,'unsafe e endpoint',r))
    require(unit<=s and s+f<=L-unit,(task,'unsafe f endpoint',s))
    period=abs(Q*e-P*f) or 1;heads=classes=0;maxg=0;minclear=None;sem=sha256()
    for residue in range(1,period+1):
        h=1
        while not ray_is_stable(L,cell,e,f,residue,period,h):
            h*=2;require(h<=1<<20,('stabilization ceiling',task,residue,period))
        maxg=max(maxg,residue+period*h)
        for hh in range(h):
            v=margin(L,cell,e,f,residue+period*hh)
            require(v>=0,('negative head',task,residue,period,hh,v))
            minclear=v if minclear is None else min(minclear,v);heads+=1
        vals=tuple(margin(L,cell,e,f,residue+period*(h+k)) for k in range(4))
        minimum,poly=quadratic_minimum(*vals[:3]);a,b,c=poly
        require(vals[3]==9*a+3*b+c,('fourth point',task,residue,h,vals,poly))
        require(minimum>=0,('negative tail',task,residue,h,poly,minimum))
        minclear=min(minclear,minimum);classes+=1;sem.update(repr((residue,h,poly,minimum)).encode())
    return task,period,heads,classes,maxg,minclear,sem.hexdigest()

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--workers',type=int,default=10);a=ap.parse_args()
    require(filehash(TAIL)==EXPECTED_TAIL,('tail hash',filehash(TAIL)))
    T=load('small35tail',TAIL)
    contexts=set()
    for body,L in T.SEL.MS.body_universe():
        cell,*_=T.SEL.body_geometry(body,L)
        if L<4592:
            for e in body:
                for f in body:
                    if e!=f:contexts.add((L,cell,e,f))
    contexts=tuple(sorted(contexts));require(len(contexts)==2530,('contexts',len(contexts)))
    with ProcessPoolExecutor(max_workers=a.workers) as pool:rows=tuple(pool.map(prove_context,contexts,chunksize=1))
    require(all(x[5]>=0 for x in rows),('negative summary',tuple(x for x in rows if x[5]<0)[:10]))
    weakest=min((x[5],x[0],x[1],x[4]) for x in rows)
    print('LRC14 DISCONNECTED 3:5 SMALL-RULER SYMBOLIC FLOOR')
    print('contexts',len(contexts),'rulers',len(set(x[0] for x in contexts)),'cells',len(set(x[:2] for x in contexts)))
    print('finite_heads',sum(x[2] for x in rows),'affine_residue_classes',sum(x[3] for x in rows),'maximum_stable_g',max(x[4] for x in rows))
    print('negative',sum(x[5]<0 for x in rows),'zero',sum(x[5]==0 for x in rows),'weakest_cleared',weakest)
    print('semantic_sha256',sha256(repr(rows).encode()).hexdigest())
    print('conclusion=for every context and every g>=1, I(3g,5g)>=1/294')

if __name__=='__main__':main()
