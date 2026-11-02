#!/usr/bin/env python3
"""Exact controls for scoped geometric/arithmetic connections."""
from fractions import Fraction as F
from itertools import combinations
from math import isqrt
import json,subprocess
from audit import decode,counts,need,det,hp,cycles,matrix_from_pairs

def point(word,b):
    k=len(word)
    return sum((F((b-1)*int(e),b**(i+1)) for i,e in enumerate(word)),F(1,2*b**k)),F(1,b**k)
def cosh(p,q):
    x,y=p;X,Y=q
    return ((x-X)**2+y*y+Y*Y)/(2*y*Y)
def common(w,v):
    h=0
    while h<min(len(w),len(v)) and w[h]==v[h]:h+=1
    return h
def trees():
    words=['']+[format(j,'0'+str(k)+'b') for k in range(1,7) for j in range(1<<k)]
    rows=[]
    for b in [2,3,4,5]:
        pts={w:point(w,b) for w in words};L=1+F(5*(b-1)**2,8*b);edges=0;checks=0
        for w,v in combinations(words,2):
            t=len(w)+len(v)-2*common(w,v);d=cosh(pts[w],pts[v]);edge=t==1
            need(d>=L and ((d==L)==edge),'sharp distance '+str((b,w,v)))
            edges+=edge;checks+=1
            if b>=3:
                need(d<=F(3,2)*b**t,'upper metric bound')
                if common(w,v)<min(len(w),len(v)):
                    need(d>=F((b-2)**2,2*b*b)*b**t,'separated branch bound')
        rows.append({'base':b,'vertices':len(words),'pairs':checks,'edges_at_minimum':edges,'cosh_minimum':str(L),'depth':6})
    hostile=[]
    for k in range(1,13):
        a='0'+'1'*(k-1);b='1'+'0'*(k-1)
        need(cosh(point(a,2),point(b,2))==F(3,2),'dyadic inverse-distance hostile')
        hostile.append({'level':k,'tree_distance':2*k,'cosh_distance':'3/2'})
    return {'finite_controls':rows,'dyadic_hostiles':hostile,'infinite_theorem':'analytic proof in GEOMETRY.md'}
def ball(neighbors,root,r):
    visited={root};level={root};sizes=[1]
    for i in range(r):
        level={v for u in level for v in neighbors(u)}-visited;visited|=level;sizes.append(len(visited))
    return sizes
def capacity():
    def H(p):
        a,b,s=p
        return [(a,b,1),(a-1,b,1),(a,b-1,1)] if s==0 else [(a,b,0),(a+1,b,0),(a,b+1,0)]
    def T(p):
        a,b=p
        return [(a+c,b+d) for c,d in [(1,0),(-1,0),(0,1),(0,-1),(1,-1),(-1,1)]]
    def C(n):return n//2 if n%2==0 else 3*n+1
    levels=[{5}];seen={5};sizes=[1]
    for r in range(1,27):
        cur=set()
        for m in levels[-1]:
            pred=[2*m]+([(m-1)//3] if m%6==4 else [])
            for n in pred:need(C(n)==m and n not in seen,'native source or clock collision');cur.add(n)
        need(len(cur)==sum(1+(m%6==4) for m in levels[-1]),'shared immediate predecessor')
        levels.append(cur);seen|=cur;sizes.append(len(seen))
    h=ball(H,(0,0,0),26);t=ball(T,(0,0),26)
    for r in range(27):need(h[r]==1+3*r*(r+1)//2 and t[r]==1+3*r*(r+1),'lattice ball')
    need(sizes[21]==704 and h[21]==694 and sizes[26]==2263 and t[26]==2107,'capacity control')
    need(sizes[21]-h[21]==10 and sizes[26]-t[26]==156,'literal deficit subtraction; twenty is false')
    return {'first_honeycomb_violation':next(r for r in range(27) if sizes[r]>h[r]),'first_triangular_violation':next(r for r in range(27) if sizes[r]>t[r]),'honeycomb_deficit_at21':sizes[21]-h[21],'triangular_deficit_at26':sizes[26]-t[26],'root':5,'clock':'ordinary Collatz','inverse_sizes':sizes,'honeycomb_sizes':h,'triangular_sizes':t}
def restricted():
    best=None;which=[];total=0
    for s in subprocess.check_output(['geng','-q','8'],text=True).splitlines():
        n,e=decode(s);E={frozenset(x) for x in e}
        if any(all(frozenset(pair) in E for pair in combinations(v,2)) for v in combinations(range(n),4)):continue
        total+=1;m,w,p=counts(n,e);value=512+33*w-32*m[2]+7*m[3]
        if best is None or value<best:best=value;which=[s]
        elif value==best:which.append(s)
    need(total==6431 and best==240,'K4-free universe')
    return {'universe':total,'minimum':best,'witnesses':which,'remaining_triangle_reserve':str(F(best,512)-F(7,40))}
def restricted_consequences():
    # Coordinates a+b*sqrt(3), kept exactly; two diamonds translated by (0,1).
    def sub(a,b):return (a[0]-b[0],a[1]-b[1])
    def sq(a):return (a[0]**2+3*a[1]**2,2*a[0]*a[1])
    diamond=[((F(-1,2),0),(0,0)),((F(1,2),0),(0,0)),
             ((0,0),(0,F(1,2))),((0,0),(0,F(-1,2)))]
    points=diamond+[(x,(y[0]+1,y[1])) for x,y in diamond]
    need(len(set(points))==8,'diamond source collision')
    D={(0,1),(0,2),(0,3),(1,2),(1,3)}
    E=D|{(a+4,b+4) for a,b in D}|{(2,6),(3,7)}
    for i,j in E:
        x=sq(sub(points[i][0],points[j][0]));y=sq(sub(points[i][1],points[j][1]))
        need((x[0]+y[0],x[1]+y[1])==(1,0),'diamond unit edge')
    m,w,p=counts(8,sorted(E))
    need(512+33*w-32*m[2]+7*m[3]==240,'geometric bracket equality')
    cases=[('unit_diamond_support',E),('fair_pair_equality',{(2*i,2*i+1) for i in range(5)}),
           ('unit_triangle_support',{(0,1),(0,2),(1,2)})]
    rows=[]
    for name,edges in cases:
        need(not any(all(tuple(pair) in edges for pair in combinations(v,2)) for v in combinations(range(11),4)),'K4-free source guard')
        values=[F(1,2) if e in edges else F(1) for e in combinations(range(11),2)]
        a=matrix_from_pairs(values,11);C=cycles(a);H=hp(a)
        u=det([[a[i][j]+int(i==j) for j in range(11)] for i in range(11)]);v=-det(a)
        V=sum(x*(1-x) for x in values)
        rhs=F(781,1280)*V+F(47,160)*C[3]+C[4]/64+C[5]+C[6]/16+C[7]+C[8]/4+C[9]+C[10]
        slack=H-u-abs(v)-rhs;need(slack>=0,'restricted matrix consequence '+name)
        if name=='fair_pair_equality':need(slack==0,'restricted sharp variance')
        rows.append({'name':name,'H':str(H),'u':str(u),'v':str(v),'V':str(V),'tau':str(C[3]),'slack':str(slack)})
    return {'unit_support_edges':sorted(E),'auxiliary_bracket':'15/32','coordinates':'D={(+-1/2,0),(0,+-sqrt(3)/2)} and D+(0,1); selected edges only','matrices':rows}
def arithmetic():
    def divs(n):return [d for d in range(1,n) if n%d==0]
    vals=[4,9,2,8,49,121,27,169,16,10,15,14,21,35,22,33,12,18,20,28,24]
    rows=[]
    for n in vals:
        d=divs(n);s=sum(d);w=[(a,n//a) for a in range(1,isqrt(n)+1) if n%a==0 and s==1+a+n//a]
        rows.append({'n':n,'aliquot':s,'proper_divisors':d,'compressed_factor_pairs':w})
    need([r['compressed_factor_pairs'] for r in rows if r['n']==20]==[[(1,20)]],'twenty witness')
    need([r['compressed_factor_pairs'] for r in rows if r['n']==16]==[[]],'sixteen boundary')
    tail=[20]
    while tail[-1]:tail.append(sum(divs(tail[-1])))
    fiber=[n for n in range(1,42) if n+sum(divs(n))==42]
    need(fiber==[20,26,41],'sigma42 source fibre')
    return {'listed_sources':rows,'twenty_tail':tail,'sigma42_fibre':fiber}
if __name__=='__main__':
    print(json.dumps({'status':'passed','trees':trees(),'capacity':capacity(),'restricted_stability':restricted(),'restricted_consequences':restricted_consequences(),'arithmetic':arithmetic()},indent=2,sort_keys=True))
