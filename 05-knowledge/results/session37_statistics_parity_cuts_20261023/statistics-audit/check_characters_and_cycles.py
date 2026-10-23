#!/usr/bin/env python3
"""Independent, standard-library exact checks. No Session37 worker imports."""
from itertools import permutations, combinations
from collections import Counter
from fractions import Fraction
from math import factorial
import json
from pathlib import Path
OUT=Path(__file__).resolve().parent

def check(b,m):
    if not b: raise AssertionError(m)
def partitions(n,lo=1):
    if n==0:yield ()
    for x in range(lo,n+1):
        for rest in partitions(n-x,x): yield (x,)+rest

def pair_action(mu):
    n=sum(mu);p=[];s=0
    for k in mu:p.extend(list(range(s+1,s+k))+[s]);s+=k
    edges=list(combinations(range(n),2));idx={e:i for i,e in enumerate(edges)}
    image=[];sgn=[]
    for i,j in edges:
        a,b=p[i],p[j];image.append(idx[tuple(sorted((a,b)))]);sgn.append(1 if a<b else -1)
    done=set();orbits=[]
    for e in range(len(edges)):
        if e in done:continue
        x=e;l=0;eps=1
        while x not in done:done.add(x);l+=1;eps*=sgn[x];x=image[x]
        orbits.append((l,eps))
    return orbits

def poly_factors(factors):
    p=[1]
    for k,a in factors:
        q=p+[0]*k
        for i,c in enumerate(p):q[i+k]+=a*c
        p=q
    return p

def chars():
    rows=[];tests=0
    for n in range(13):
        size=n*(n-1)//2;avg_e=[Fraction(0)]*(size+1);avg_w=avg_e.copy()
        for mu in partitions(n):
            orb=pair_action(mu);ex=poly_factors((k,(-1)**(k+1)*eps)for k,eps in orb);wa=poly_factors(orb)
            fixed=0 if any(eps<0 for _,eps in orb) else 2**len(orb)
            check(sum(ex)==fixed==sum(wa),('characters',n,mu))
            if n>=3:check(sum((-1)**k*a for k,a in enumerate(ex))==0,('super',n,mu))
            cm=Counter(mu);z=1
            for k,m in cm.items():z*=k**m*factorial(m)
            for k in range(size+1):avg_e[k]+=Fraction(ex[k],z);avg_w[k]+=Fraction(wa[k],z)
            tests+=1
        check(all(x.denominator==1 and x>=0 for x in avg_e+avg_w),('integrality',n))
        e=list(map(int,avg_e));w=list(map(int,avg_w))
        rows.append({'n':n,'types':sum(e),'exterior_superdimension':sum((-1)**k*a for k,a in enumerate(e)),'exterior_invariant_degrees':e,'walsh_invariant_degrees':w})
    check(rows[3]['exterior_invariant_degrees']==[1,0,0,1],'n3ex')
    check(rows[3]['walsh_invariant_degrees']==[1,0,1,0],'n3wa')
    return {'conjugacy_class_checks':tests,'rows':rows}

def tournament(mask,n):
    A=[[0]*n for _ in range(n)]
    for e,(i,j) in enumerate(combinations(range(n),2)):
        A[i][j]=(mask>>e)&1;A[j][i]=1-A[i][j]
    return A

def nine_fixture():
    A=[[0]*9 for _ in range(9)]
    for i,j in combinations(range(9),2):
        if i//3<j//3 or (i//3==j//3 and (j-i)==1):A[i][j]=1
        else:A[j][i]=1
    A[0][8]=0;A[8][0]=1
    return A

def cycle_cover_census(A):
    n=len(A);d=[0]*(n+1);per=d.copy();pack=[0]*(n//3+1);single=Counter();h=0;valid=0
    for p in permutations(range(n)):
        if all(A[p[i]][p[i+1]] for i in range(n-1)):h+=1
        if not all(i==p[i] or A[i][p[i]] for i in range(n)):continue
        valid+=1;seen=set();lengths=[]
        for i in range(n):
            if i in seen:continue
            j=i;l=0
            while j not in seen:seen.add(j);l+=1;j=p[j]
            lengths.append(l)
        fixed=lengths.count(1);per[fixed]+=1;d[fixed]+=(-1)**(n-len(lengths))
        if all(l%2==1 for l in lengths):
            k=len(lengths)-fixed;pack[k]+=1
            if k==1:single[next(l for l in lengths if l>1)]+=1
    check(h==sum(c*2**k for k,c in enumerate(pack)),('OCF',h,pack))
    return {'hamiltonian_paths':h,'det_xI_A_low_to_high':d,'per_xI_A_low_to_high':per,'odd_cycle_packing':pack,'single_odd_cycles_by_length':dict(single),'valid_cycle_covers':valid,'outmasks':[sum(A[i][j]<<j for j in range(n))for i in range(n)]}

def main():
    c=chars();(OUT/'characters.json').write_text(json.dumps(c,indent=2)+'\n')
    cases=[('real_root_counterexample',nine_fixture()),('det_per_collision_a',tournament(29759988540,9)),('det_per_collision_b',tournament(18875161993,9))]
    r={}
    for name,A in cases:r[name]=cycle_cover_census(A);print(name,r[name],flush=True)
    a=r['det_per_collision_a'];b=r['det_per_collision_b'];q=r['real_root_counterexample']
    check(a['det_xI_A_low_to_high']==b['det_xI_A_low_to_high'],'det collision');check(a['per_xI_A_low_to_high']==b['per_xI_A_low_to_high'],'per collision');check((a['hamiltonian_paths'],b['hamiltonian_paths'])==(2843,2841),'H distinction')
    check(q['odd_cycle_packing']==[1,49,7,1],'root counterexample');check(q['hamiltonian_paths']==135,'H135')
    _,c0,b0,a0=q['odd_cycle_packing'];discriminant=b0*b0*c0*c0-4*a0*c0**3-4*b0**3-27*a0*a0+18*a0*b0*c0
    check(discriminant==-348172,'discriminant');q['cubic_discriminant']=discriminant
    r['status']='PASS';(OUT/'cycles.json').write_text(json.dumps(r,indent=2)+'\n')
    print('PASS',c['conjugacy_class_checks'],'conjugacy class checks')
if __name__=='__main__':main()
