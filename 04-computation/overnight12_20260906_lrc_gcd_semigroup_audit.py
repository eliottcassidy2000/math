"""Independent referee for the twelfth signed-box/actual-entry theorem.

No producer imports: enumerate signed coefficients, every outside coefficient,
and literal phases. The all-height result is audited in the companion prose.
"""
from fractions import Fraction as F
from functools import reduce
from itertools import combinations
from math import gcd, prod
from pathlib import Path
import hashlib
import json
import sys

sys.stdout.reconfigure(newline='\n')
GATES=0
Q=91**6

def need(ok,label):
    global GATES
    GATES+=1
    if not ok:
        raise RuntimeError(label)

def cgcd(v):
    return reduce(gcd,v,0)

def ceildiv(a,b):
    return -((-a)//b)

def interval_member(q,a,b,x):
    r=pow(a,-1,b)*x%b
    s=(x-a*r)//b
    lo=max(ceildiv(-q-r,b),ceildiv(s-q,a))
    hi=min((q-r)//b,(s+q)//a)
    if lo>hi:
        return False,None
    return True,(r+lo*b,s-lo*a)

def box(q,a,b):
    return {a*r+b*s for r in range(-q,q+1) for s in range(-q,q+1)}

def reduced_cross(q,A,B,Y):
    d=gcd(A,B);a=A//d;b=B//d
    delta=gcd(d,Y)
    return d//delta<=q and interval_member(q,a,b,Y//delta)[0]

def norm(y):
    y=y%1
    return min(y,1-y)

def allowed_sum(n):
    if n>356:
        return False
    for p in range(2,n+1):
        if n%p:
            continue
        e=0
        while n%p==0:
            n//=p;e+=1
        if p%3!=2 or e>2:
            return False
        if n==1:
            return True
    return n==1

def atlas_components(row):
    # Union-find, independently of the producer's graph traversal.
    parent=list(range(len(row)))
    def find(i):
        while parent[i]!=i:
            parent[i]=parent[parent[i]];i=parent[i]
        return i
    for i,j in combinations(range(len(row)),2):
        d=gcd(row[i],row[j])
        if allowed_sum((row[i]+row[j])//d):
            parent[find(i)]=find(j)
    return sorted(sorted(i for i in range(len(row)) if find(i)==r)
                  for r in {find(i) for i in range(len(row))})

def complete_grid_phase(V,g):
    # Independently find a tiny-denominator core phase, then select one g-lift
    # of pair phase 1/4 by nearest-integer arithmetic, never enumerate the grid.
    targets=[F(a,d) for d in range(2,101) for a in range(1,d)
             if all(norm(v*F(a,d))>=F(1,12) for v in V)]
    need(bool(targets),'independent lower-dimensional literal phase')
    y=min(targets,key=lambda x:(x.denominator,x.numerator))
    z=F(1,4)
    v=g*y-z
    j=(v+F(1,2)).numerator//(v+F(1,2)).denominator
    x=(z+j)/g
    need(abs(x-y)<=F(1,2*g),'nearest complete g-grid point')
    need(all(norm(v*x)>=F(1,14) for v in V+[g,3*g]),'literal full thirteen-speed weak safety')
    return str(y),str(x)

def main():
    here=Path(__file__).parent
    expected={
      'overnight12_20260906_lrc_gcd_semigroup.py':'15a3fa8f511640c2fee3404cb5cc5b377af4164f8db3dcec7aa1702c2f15e4c3',
      'overnight12_20260906_lrc_gcd_semigroup.out':'b3f47f0f979c9b126dd6a67b2a83c0dc8dc30fe46b9834d389d3425f82f66ab8'}
    for name,h in expected.items():
        source_path=here/name
        if not source_path.exists() and name.endswith('.out'):
            source_path=here.parent/'05-knowledge'/'results'/name
        need(hashlib.sha256(source_path.read_bytes()).hexdigest()==h,'frozen producer input '+name)
    boxes=points=0
    for q in range(2,14):
        for b in range(2,q+1):
            for a in range(1,b):
                if gcd(a,b)!=1:
                    continue
                raw=box(q,a,b);L=q*(a+b);R=L-(a-1)*(b-1)
                first=next(x for x in range(L+2) if x not in raw)
                need(first==R+1,'literal first positive gap')
                need(-first not in raw,'negative first gap')
                need(q*b<=R<=L and 2*R>L,'central interval exceeds half support')
                for x in range(-L-2,L+3):
                    yes,witness=interval_member(q,a,b,x)
                    need(yes==(x in raw),'independent set versus interval membership')
                    if yes:
                        r,s=witness
                        need(abs(r)<=q and abs(s)<=q and a*r+b*s==x,'constructed signed coefficient witness')
                    points+=1
                boxes+=1
    triple_cases=0
    for q in range(2,8):
        for b in range(2,q+1):
            for a in range(1,b):
                if gcd(a,b)!=1:
                    continue
                for d in range(1,10):
                    A=d*a;B=d*b
                    signed={A*r+B*s for r in range(-q,q+1) for s in range(-q,q+1)}
                    for Y in range(1,46):
                        raw=any(C*Y in signed for C in range(1,q+1))
                        need(raw==reduced_cross(q,A,B,Y),'complete outside-coefficient enumeration versus gcd reduction')
                        triple_cases+=1
    need(14 not in box(3,2,3) and 15 in box(3,2,3),'first gap is not a complete support cutoff')
    need(3 not in box(2,1,6) and 6 in box(2,1,6),'b>Q destroys minimal outside-coefficient principle')
    need(any(C*3 in box(2,1,6) for C in (1,2)),'missing minimal point rescued outside hypothesis')
    need(1 in box(2,1,2) and not reduced_cross(2,3,6,1),'minimal outside coefficient gate cannot be dropped')
    discrete=0
    for R in range(1,38):
        for delta in range(1,16):
            for s in range(1,16):
                z=gcd(delta,s);d=delta//z;e=s//z
                B=d*(R//e+1)
                # Literal divisibility/strict inequality at every previous integer.
                first=next(g for g in range(1,B+1) if g*s%delta==0 and g*s>R*delta)
                need(first==B,'exact first scale from strict radius and retained divisibility')
                discrete+=1
    H=Q//(42*177)
    need(Q==567869252041 and H==76388115,'inherited coefficient height and endpoint cutoff')
    need(677*H==51714753855<Q,'native coefficient gate in unresolved gluing branch')
    need((21*355-1)//11==677,'strict negation of the first gluing branch')
    atlas=[(a,b) for b in range(2,356) for a in range(1,b)
           if gcd(a,b)==1 and allowed_sum(a+b)]
    need(len(atlas)==5855 and max(a for a,b in atlas)==177,'literal inherited atlas and small pair coordinate')
    need(max(b for a,b in atlas)==355,'large pair coordinate')
    records=[]
    controls=[([1,4,6,8,10,12,14,15,16,18,22],2**45),
              ([2,3,4,5,6,10,12,14,15,20,30],60*Q+1)]
    for V,g in controls:
        K=max(V);row=V+[g,3*g]
        need(cgcd(V)==1 and cgcd(row)==1 and len(set(row))==13,'primitive and distinct actual-entry control')
        need(sum(row)<=Q*Q and g>2*Q*K,'physical box and uniform all-crossing dominance')
        comps=atlas_components(row)
        need(comps==[list(range(11)),[11,12]],'literal actual decoder graph')
        need(sum(len(c)-1 for c in comps)==11,'connected internal relation span has exact decoder rank')
        # Every crossing has a nonzero outside partial sum, a nonzero integer
        # multiple of g. Its absolute value is >=g, while at most two core
        # labels contribute <=2QK. This covers both support orientations.
        for u,v in combinations(V,2):
            for s in (1,3):
                need(not reduced_cross(Q,u,v,g*s),'all 110 two-core/one-outside selected supports')
        y,x=complete_grid_phase(V,g)
        records.append(dict(core=V,g=g,K=K,sum=sum(row),core_phase=y,full_phase=x))
    V=controls[1][0]
    scores=[]
    for u,v in combinations(V,2):
        d=gcd(u,v);a=u//d;b=v//d
        scores.append((Q*(a+b)-(a-1)*(b-1),u,v))
    need(max(scores)==(29*Q-182,14,15),'exact best interior pair over all 55 choices')
    need(max(x for x in scores if x[2]==30)==(22*Q-84,14,30),'exact best maximum pair')
    hostile=list(range(1,11))+[85]
    need(all(norm(v*F(1,12))>=F(1,12) for v in hostile),'global arc hostile starts at full 1/12 safety')
    need(abs(F(7,85)-F(1,12))==F(1,1020)<F(1,840),'selected-pair radius is too large')
    need(norm(85*F(7,85))==0,'global maximum destroys selected-pair protected arc')
    primes=[37,43,61,67,73,79,97,103,127]
    P=15*prod(primes);V=[2*P//r for r in primes]+[P//3,P//5]
    need(cgcd(V)==1 and min(V)>1 and max(V)==237907127334685115>Q,'primorial is a valid normalization hostile')
    need(max(max(u,v)//gcd(u,v) for u,v in combinations(V,2))==127,'all primitive pair heights stay small')
    need([cgcd(V[:7]),cgcd(V[:8]),cgcd(V[:9]),cgcd(V[:9]+[V[9]])]==[392430,3810,30,5],
         'primorial violates the current 90/30/9/4 subset caps')
    print('STATUS: INDEPENDENT PASS; exact signed boxes, selected three-label iff, and actual-entry closure')
    print('UNIVERSES:',boxes,'full signed boxes;',points,'literal test integers;',triple_cases,'full three-label cases;',discrete,'discrete scale cases')
    print('ACTUAL ENTRY CONTROLS:',json.dumps(records,sort_keys=True))
    print('CUT-OFF:',H,'; all-scale atlas pairs:',len(atlas))
    print('SCOPE: 110 selected two-core/one-outside supports; not the eleven opposite-orientation supports')
    print('SCOPE: source lower-dimensional LRC is CITED; incoming 90/30 hierarchy is INHERITED, not rerun')
    print('SEMANTIC SHA256:',hashlib.sha256(json.dumps(records,sort_keys=True).encode()).hexdigest())
    print('ACTIVE GATES:',GATES)

if __name__=='__main__':
    main()
