#!/usr/bin/env python3
"""Independent norm-five profile proof and compact global sharp finite bases.

Exact Bernstein sign certificates verify all affine-envelope pieces.
Finite projections come from interval intersections, and physical masses
also use the literal three-sheet breakpoint predicate. No producer imports.
"""
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations,permutations,product
from math import gcd,lcm,comb
from pathlib import Path
import csv
import argparse
import sys
import sympy as sp

sys.stdout.reconfigure(newline='\n')
t,z,u=sp.symbols('t z u')
R=sp.Rational
GATES=0
REPO=Path(__file__).resolve().parents[1]


def need(ok,label):
    global GATES
    GATES+=1
    if not ok:
        raise RuntimeError(label)


def bernstein_nonnegative(poly,left,right,depth=0):
    P=sp.Poly(sp.expand(poly.subs(t,left+(right-left)*u)),u)
    degree=P.degree()
    if P.is_zero:
        return
    aa=[P.nth(j) for j in range(degree+1)]
    bb=[sum(aa[j]*R(comb(k,j),comb(degree,j)) for j in range(k+1)) for k in range(degree+1)]
    if all(x>=0 for x in bb):
        need(True,'exact Bernstein coefficients certify numerator sign')
        return
    need(depth<8,('Bernstein subdivision depth',poly,left,right,bb))
    middle=(left+right)/2
    bernstein_nonnegative(poly,left,middle,depth+1)
    bernstein_nonnegative(poly,middle,right,depth+1)


def nonnegative(expr,left,right,label):
    numerator,denominator=sp.cancel(expr).as_numer_denom()
    if numerator==0:
        need(True,label)
        return
    denpoly=sp.Poly(denominator,t)
    roots=denpoly.all_roots() if denpoly.degree()>0 else []
    need(all(x.is_real in (True,False) for x in roots),('exact denominator-root classification',label))
    need(not any(x.is_real and left<x<right for x in roots),('denominator nonzero in parameter sector',label))
    sign=sp.sign(denominator.subs(t,(left+right)/2))
    need(sign in (-1,1),('denominator sign',label))
    bernstein_nonnegative(sign*numerator,left,right)


def integral(line,left,right):
    P=sp.Poly(line,z)
    need(P.degree()<=1,'affine profile integration')
    return sp.cancel(P.nth(0)*(right-left)+P.nth(1)*(right**2-left**2)/2)


def symbolic_profiles():
    sectors=(('I',R(1,2)-t,(2,2,1),R(0),R(1,4),R(1,2),R(29,32)),
             ('II',1-t/2,(1,2,2),R(0),R(2,3),(2+t)/4,R(121,128)),
             ('III',2-2*t,(2,1,2),R(1,2),R(2,3),(2-t)/2,R(121,128)),
             ('IV',t+R(1,2),(2,2,1),R(0),R(1,2),(1+t)/2,R(15,16)))
    selector_domains={'I':[(R(0),R(1,8),0),(R(1,8),R(1,4),2)],
                      'II':[(R(0),R(1,8),0),(R(1,8),R(2,3),2)],
                      'III':[(R(1,2),R(9,16),0),(R(9,16),R(2,3),2)],
                      'IV':[(R(0),R(1,4),0),(R(1,4),R(1,2),0)]}
    physical_domains={
        'I':[(R(0),R(1,4),[(0,R(0),t/2+R(1,4)),(1,t/2+R(1,4),(1-t)/2),(3,(1-t)/2,R(1,2))])],
        'II':[(R(0),R(1,2),[(0,R(0),t/2),(1,t/2,(1-t)/2),(3,(1-t)/2,(2+t)/4)]),
              (R(1,2),R(2,3),[(0,R(0),(1-t)/2),(2,(1-t)/2,t/2),(3,t/2,(2+t)/4)])],
        'III':[(R(1,2),R(2,3),[(0,R(0),t-R(1,2)),(1,t-R(1,2),1-t),(3,1-t,1-t/2)])],
        'IV':[(R(0),R(1,2),[(0,R(0),R(1,4)-t/2),(1,R(1,4)-t/2,R(1,2)),(2,R(1,2),(1+t)/2)])]
    }
    records=[]
    for name,b,ds,left,right,cutoff,maximum in sectors:
        speeds=(t,b,R(1))
        As=[sum(speeds)-x for x in speeds]
        Bs=[speeds[(i+1)%3]*speeds[(i+2)%3] for i in range(3)]
        lines=[R(2)]+[(As[i]-ds[i]*z)/Bs[i] for i in range(3)]
        plateau=[sp.cancel((As[i]-2*Bs[i])/ds[i]) for i in range(3)]
        for i in range(3):
            nonnegative(As[i]/ds[i]-cutoff,left,right,(name,'complete cutoff',i))
            nonnegative(plateau[i],left,right,(name,'initial cap plateau',i))
        need(any(sp.cancel(As[i]/ds[i]-cutoff)==0 for i in range(3)),(name,'some roof attains complete cutoff'))
        selected=[]
        for lo,hi,winner in selector_domains[name]:
            Js=[]
            for i in range(3):
                if (plateau[i]-cutoff).subs(t,(lo+hi)/2)>=0:
                    nonnegative(plateau[i]-cutoff,lo,hi,(name,'cap throughout',i))
                    value=2*cutoff
                else:
                    nonnegative(cutoff-plateau[i],lo,hi,(name,'roof after cap',i))
                    value=2*plateau[i]+integral(lines[i+1],plateau[i],cutoff)
                Js.append(sp.cancel(value))
            for j in range(3):
                nonnegative(Js[j]-Js[winner],lo,hi,(name,'selected integral winner',winner,j))
            nonnegative(maximum-Js[winner],lo,hi,(name,'global continuum maximum'))
            selected.append((lo,hi,sp.factor(Js[winner])))
        for lo,hi,pieces in physical_domains[name]:
            need(pieces[0][1]==0 and sp.cancel(pieces[-1][2]-cutoff)==0,'full physical z-domain')
            need(all(sp.cancel(a[2]-b[1])==0 for a,b in zip(pieces,pieces[1:])),
                 'physical envelope has no missing interval')
            total=R(0)
            for index,a,c in pieces:
                nonnegative(c-a,lo,hi,(name,'ordered physical breakpoints'))
                for j in range(4):
                    nonnegative((lines[j]-lines[index]).subs(z,a),lo,hi,(name,'active line left endpoint',index,j))
                    nonnegative((lines[j]-lines[index]).subs(z,c),lo,hi,(name,'active line right endpoint',index,j))
                total+=integral(lines[index],a,c)
            need(sp.cancel(total-R(7,8))==0,(name,'exact complete physical integral'))
        records.append((name,selected,maximum))
        print('CERTIFIED_SELECTOR',name,selected,'maximum_or_supremum',maximum)
    need(R(3,49)*R(121,128)==R(363,6272),'normalized selected bulk')
    need(R(3,49)*R(7,8)==R(3,56),'normalized physical bulk')
    need(2*R(3,14)**2*R(7,8)==R(9,112),'full raw physical profile integral')
    for c,target,bulk in ((28,R(11,140),R(363,6272)),(51,R(46,665),R(363,6272)),
                          (46,R(51,770),R(3,56))):
        need(bulk+R(4,7*c)<target,('strict layer-cake tail cutoff',c,target))
    print('CERTIFIED_PHYSICAL_J all_four_sectors=7/8; raw_integral=9/112; bulk=3/56')
    print('STRICT_CUTOFFS selected11/140:28 selected46/665:51 physical51/770:46')
    return records


def cross(a,b):
    return (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])


def egcd(a,b):
    x,y,X,Y=1,0,0,1
    while b:
        q=a//b
        a,b=b,a-q*b
        x,X=X,x-q*X
        y,Y=Y,y-q*Y
    return a,x,y


def ray_interval_values(w,v):
    g,s,h=egcd(w[0],w[1])
    one,m,n=egcd(g,w[2])
    bezout=(m*s,m*h,n)
    need(one==1 and sum(x*y for x,y in zip(bezout,w))==1,'primitive full-lattice lift')
    lift=tuple(-x for x in cross(bezout,v))
    need(cross(w,lift)==v,'integral carrier lift')
    E=[Q(0)]*3
    mass=Q(0)
    live=[]
    for k in range(1,w[2]+1):
        if k%3==0:
            continue
        intervals=[(Q(14*k*ni-3,14*wi),Q(14*k*ni+3,14*wi)) for ni,wi in zip(lift,w)]
        joint=min(b for a,b in intervals)-max(a for a,b in intervals)
        if joint<=0:
            continue
        live.append(k)
        cap=min(b-a for a,b in intervals)
        for omitted in range(3):
            others=[intervals[j] for j in range(3) if j!=omitted]
            E[omitted]+=2*min(cap,min(b for a,b in others)-max(a for a,b in others))
        mass+=2*joint
    return tuple(E),mass,live


def literal_tail(w):
    D=42*lcm(*w)
    end=D//3
    points={0,end}
    for speed in w:
        unit=D//(14*speed)
        for k in range(speed):
            for sign in (-1,1):
                for j in range(3):
                    point=((14*k+sign)*unit-j*end)%D
                    if point<=end:
                        points.add(point)
    ordered=sorted(points)
    length=0
    for a,b in zip(ordered,ordered[1:]):
        x=a+b
        bad=True
        for j in range(3):
            if not any(14*min((speed*(x+2*j*end))%(2*D),(-speed*(x+2*j*end))%(2*D))<2*D for speed in w):
                bad=False
                break
        if bad:
            length+=b-a
    return Q(3*length,D)


def finite_bases():
    tsv=REPO/'05-knowledge/results/overnight4_20260906_lrc_parityfree_probe_head63.tsv'
    with tsv.open(newline='',encoding='utf-8') as f:
        table={tuple(int(row[x]) for x in ('a','b','c')):row for row in csv.DictReader(f,delimiter='\t')}
    vectors=sorted({tuple(a*b for a,b in zip(p,signs)) for p in set(permutations((1,2,2)))
                    for signs in product((-1,1),repeat=3) if signs[0]==1})
    expected={}
    for w in combinations([x for x in range(1,51) if x%3],3):
        if gcd(*w)!=1:
            continue
        relations=[v for v in vectors if sum(a*b for a,b in zip(v,w))==0]
        if not relations:
            continue
        a,b,c=w
        sectors=(c==2*(a+b),2*c==a+2*b,2*c==2*a+b,c==2*(b-a))
        need(sum(sectors)==1 and len(relations)==1,'complete disjoint normalized sign sectors')
        E,mass,live=ray_interval_values(w,relations[0])
        need(literal_tail(w)==mass,('literal sheet predicate versus complete ray interval mass',w))
        need(w in table,('frozen full-head membership',w))
        row=table[w]
        D=int(row['denominator'])
        need(E==tuple(Q(int(row['E'+str(i)+'_numerator']),D) for i in range(3))
             and mass==Q(int(row['mass_numerator']),D) and 2*len(live)==int(row['raw_carriers']),
             ('all interval projections and physical mass versus native-audited complete head',w))
        need(min(E)<=Q(46,665) and mass<=Q(51,770),('global sharp finite control',w))
        expected[w]=(E,mass)
    need(len(expected)==174,'entire compact selected H50 base')
    head27={w:values for w,values in expected.items() if w[2]<=27}
    head45={w:values for w,values in expected.items() if w[2]<=45}
    need(len(head27)==48,'entire weaker H27 base')
    need(all(min(E)<Q(11,140) for E,mass in head27.values()),'strict weaker selected H27 base')
    selected=[w for w,(E,mass) in expected.items() if min(E)==Q(46,665)]
    physical=[w for w,(E,mass) in head45.items() if mass==Q(51,770)]
    need(selected==[(2,19,20)] and physical==[(1,11,20)],'unique global equality loci')
    E,mass=expected[(10,11,16)]
    need(E==(Q(17,176),Q(9,140),Q(3,55)) and mass==Q(331,6160) and min(E)-mass==Q(1,1232),
         'fixed-roof hostile retained')
    print('COMPACT_BASES H27',len(head27),'H45',len(head45),'H50',len(expected))
    print('SHARP_SELECTED 46/665 equality',selected,'SHARP_PHYSICAL 51/770 equality',physical)
    print('SWITCHING_ROOF_HOSTILE',(10,11,16),E,mass,'gap',min(E)-mass)
    print('NATIVE_AUDITED_HEAD_SHA256',sha256(tsv.read_bytes()).hexdigest())


def main():
    global REPO
    parser=argparse.ArgumentParser()
    parser.add_argument('--repo',type=Path,default=REPO,help='Repository root containing the frozen H63 table')
    REPO=parser.parse_args().repo
    symbolic_profiles()
    finite_bases()
    print('SCOPE independent concurrent confirmation of incoming THM4441; no new sharp-closure priority')
    print('PASS_OPTIMIZATION_LIVE_GATES',GATES)


if __name__=='__main__':
    main()
