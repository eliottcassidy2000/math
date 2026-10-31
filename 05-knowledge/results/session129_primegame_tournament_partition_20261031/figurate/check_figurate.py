#!/usr/bin/env python3
"""Exact Session129 triangular/Fibonacci/tournament observer bridges.
Standard library only; explicit checks survive python -O; no upstream imports.
"""
from fractions import Fraction
from itertools import combinations, permutations, product
from math import gcd, isqrt, prod
import json

P=(2,3,5,7,11,13,17,19,23,29)
RENAMED=(31,37,41,43,47,53,59,61,67,71)
RULES=((17,91),(78,85),(19,51),(23,38),(29,33),(77,29),(95,23),(77,19),(1,17),(11,13),(13,11),(15,2),(1,7),(55,1))
checks=0

def need(ok, msg):
    global checks
    checks += 1
    if not ok: raise RuntimeError(msg)

def vec(n):
    ans=[]
    for p in P:
        e=0
        while n%p==0:
            n//=p; e+=1
        ans.append(e)
    need(n==1,'number supported on ten primes')
    return tuple(ans)

def integer(v, primes=P): return prod(p**e for p,e in zip(primes,v))
def parity(v): return sum((e%2)<<i for i,e in enumerate(v))
def step(n, rules=RULES):
    for i,(a,b) in enumerate(rules):
        if n%b==0: return i,n//b*a
    raise RuntimeError('unreachable: denominator1')
def tri(n): return n*(n+1)//2
def Q(x,y): return x*x+x*y-y*y
def face(n): return n,n+1,n*(n+1)
def window(n): return 1,n,n+1,2*n+1

def istri(n): return n>=0 and isqrt(8*n+1)**2==8*n+1

def fib_index(n):
    a,b=0,1;k=0
    while a<n: a,b=b,a+b;k+=1
    return k if a==n else None

def orbit_census(n):
    pairs=list(combinations(range(n),2)); positions={p:i for i,p in enumerate(pairs)}
    permlist=list(permutations(range(n))); full=(1<<len(pairs))-1
    def arc(mask,i,j):
        return bool((mask>>positions[min(i,j),max(i,j)])&1)==(i<j)
    def relabel(mask,p):
        return sum(1<<k for k,(i,j) in enumerate(pairs) if arc(mask,p[i],p[j]))
    canonical={}; sizes={}
    for mask in range(full+1):
        if mask in canonical: continue
        orbit={relabel(mask,p) for p in permlist}; c=min(orbit)
        canonical.update((x,c) for x in orbit);sizes[c]=len(orbit)
    merged={m:min(canonical[m],canonical[full^m]) for m in range(full+1)}
    SC=sum(canonical[full^c]==c for c in sizes)
    def hp(c):return sum(all(arc(c,p[i],p[i+1]) for i in range(n-1)) for p in permlist)
    merged_rows=[{'representative':c,'H':hp(c),'labelled_size':sum(merged[m]==c for m in merged)} for c in sorted(set(merged.values()))]
    return {'n':n,'edge_slots':len(pairs),'fixed_path_free_slots':len(pairs)-(n-1),'labelled':full+1,'isoclasses':len(sizes),'self_converse':SC,'merged':len(set(merged.values())),'rows':merged_rows}, merged

census=[]; maps={}
for n in (3,4,5):
    row,m=orbit_census(n);census.append(row);maps[n]=m
need([(r['isoclasses'],r['self_converse'],r['merged']) for r in census]==[(2,2,2),(4,2,3),(12,8,10)],'census and Burnside counts')
for r in census:
    need(r['merged']==(r['isoclasses']+r['self_converse'])//2,'converse orbit count')
    need(r['edge_slots']==tri(r['n']-1),'triangular edge slots')
    need(r['fixed_path_free_slots']==tri(r['n']-2),'triangular free slots')

# Exact state bridge and all factor-support toggle masks on a whole finite box.
alpha=[vec(a) for a,b in RULES];beta=[vec(b) for a,b in RULES]
toggles=[parity(tuple(a-b for a,b in zip(x,y))) for x,y in zip(alpha,beta)]
need(all(max(v)<=1 for v in alpha+beta),'all fraction parts squarefree')
renamed_rules=[(integer(a,RENAMED),integer(b,RENAMED)) for a,b in zip(alpha,beta)]
for exponents in product(range(2),repeat=10):
    n=integer(exponents); i,nxt=step(n); v=vec(nxt)
    need(all(exponents[j]>=beta[i][j] for j in range(10)),'selected divisibility guard')
    need(v==tuple(exponents[j]+alpha[i][j]-beta[i][j] for j in range(10)),'integer/exponent exact update')
    need(parity(v)==parity(exponents)^toggles[i],'K5 toggle cocycle')
    u,w,z=face(n);nextu,nextw,nextz=face(nxt)
    need(w*z==u*(w+z),'harmonic triangle')
    need((gcd(u,w),gcd(u,z),gcd(w,z))==(1,n,n+1),'gcd decoder')
    need(z==2*tri(n),'triangular coordinate')
    a,b=RULES[i]
    need((nextu,nextw,nextz)==(a*n//b,a*n//b+1,(a*n//b)*(a*n//b+1)),'lift conjugacy')
    need(Fraction(tri(nxt))==Fraction(a*a,b*b)*tri(n)+Fraction(a*(b-a),2*b*b)*n,'triangular polarization law')
    ww=window(n)
    need(ww[2]==ww[0]+ww[1] and ww[3]==ww[1]+ww[2],'generalized Fibonacci recurrence window')
    need(ww[0]*ww[3]-ww[1]*ww[2]==Q(1,n),'antipodal Cassini defect')
    j, renamed_nxt=step(integer(exponents,RENAMED),renamed_rules)
    need(j==i and renamed_nxt==integer(v,RENAMED),'prime-renaming conjugacy')

# Minimized parity-observer failure: search by larger integer then smaller one.
def support_integer(n):
    t=n
    for p in P:
        while t%p==0:t//=p
    return t==1
states=[n for n in range(1,101) if support_integer(n)]
parity_hostile=None;merged_hostile=None
for b in states:
    for a in states:
        if a>=b:break
        va,vb=vec(a),vec(b);ia,na=step(a);ib,nb=step(b)
        pa,pb=parity(va),parity(vb);qa,qb=parity(vec(na)),parity(vec(nb))
        if parity_hostile is None and pa==pb and qa!=qb:
            parity_hostile={'states':[a,b],'parity':pa,'rules':[ia+1,ib+1],'next_states':[na,nb],'next_parities':[qa,qb],'next_merged_representatives':[maps[5][qa],maps[5][qb]]}
        if merged_hostile is None and maps[5][pa]==maps[5][pb] and maps[5][qa]!=maps[5][qb]:
            merged_hostile={'states':[a,b],'merged_representative':maps[5][pa],'parities':[pa,pb],'rules':[ia+1,ib+1],'next_states':[na,nb],'next_parities':[qa,qb],'next_merged_representatives':[maps[5][qa],maps[5][qb]]}
need(parity_hostile['states']==[1,4],'smallest same-parity unequal update hostile')
need(parity_hostile['next_merged_representatives'][0]!=parity_hostile['next_merged_representatives'][1],'same-parity hostile survives merged quotient')
need(merged_hostile['states']==[1,4],'smallest merged nonautonomy hostile')
need([abs(Q(1,n)) for n in (1,2,15)]==[1,1,209],'unit Cassini first step hostile')
need(step(2)==(11,15),'seed2 first transition')
need([n for n in range(1,1001) if abs(Q(1,n))==1]==[1,2],'finite unit Cassini audit')

# Positive-ratio compactification a/b -> a/(a+b), with seed (b,a).
fraction_rows=[]
for a,b in RULES:
    x,y=b,a;w=(x,y,x+y,x+2*y);h=(x*y,x*(x+y),y*(x+y))
    need(gcd(a,b)==1,'printed fractions reduced')
    need(h[1]*h[2]==h[0]*(h[1]+h[2]),'fraction compactification harmonic')
    need((gcd(h[0],h[1]),gcd(h[0],h[2]),gcd(h[1],h[2]))==(x,y,x+y),'fraction face gcd decoder')
    need(w[0]*w[3]-w[1]*w[2]==Q(x,y),'fraction Cassini norm')
    need(abs(Q(x,y))!=1,'no literal instruction compactification is unit Cassini')
    fraction_rows.append({'fraction':f'{a}/{b}','seed':[x,y],'window':w,'harmonic_face':h,'Cassini_norm':Q(x,y)})
need(Q(7,1)==55==tri(10),'eraser7 norm is restart numerator and T10')
need(fib_index(55)==10,'55=F10')
need(Q(17,1)==305,'eraser17 norm hostile')
need(Q(11,13)==95 and Q(13,11)==-Q(2,15)==191,'other exact norm hits')
for x in range(-8,9):
    for y in range(-8,9):
        need(Q(y,x+y)==-Q(x,y),'Fibonacci matrix norm reversal')
        need(Q(x+y,-y)==Q(x,y),'quadratic conjugation')
        need(Q(x-y,2*x-y)==-Q(x,y),'projective norm involution')
for q in range(1,101):
    k=Q(q,1)
    need(k==2*tri(q)-1,'eraser norm triangular law')
    need(4*k+5==(2*q+1)**2 and (isqrt(4*k+5)-1)//2==q,'inverse eraser norm square gate')
norm_groups={}
for a,b in RULES:norm_groups.setdefault(abs(Q(b,a)),[]).append(f'{a}/{b}')
need(len(norm_groups)==13,'thirteen distinct absolute norm values')
need([g for g in norm_groups.values() if len(g)>1]==[['11/13','15/2']],'unique equal-norm instruction pair')
need((13-11,2*13-11)==(2,15),'F cubed conjugation realizes unique orbit pair')
need((8,1) not in RULES,'F maps1/7 to8/1 outside rack')
renamed55=RENAMED[2]*RENAMED[4]
need(renamed55==1927 and not istri(renamed55) and fib_index(renamed55) is None,'numeric55 identities fail renaming')
need(Q(RENAMED[3],1)!=renamed55,'norm relation itself not prime-renaming equivariant')

# A choice-free type boundary: a finite observer cannot recover arbitrary states.
need(maps[5][0]==maps[5][1023],'converse quotient identifies global reversal')
# Declared lexicographic gauge identifies register slots with EDGE slots, not nodes.
prime_node_rows=[{'prime':p,'edge':list(e),'singleton_node':maps[5][1<<i]} for i,(p,e) in enumerate(zip(P,combinations(range(5),2)))]
need(len(set(r['singleton_node'] for r in prime_node_rows))<10,'prime singleton observer is not a node bijection')
output={
 'status':'PASS','explicit_checks':checks,
 'fractions':[f'{a}/{b}' for a,b in RULES],
 'prime_edge_gauge':prime_node_rows,'register_singleton_image_size':len(set(r['singleton_node'] for r in prime_node_rows)),
 'census':census,'factor_toggle_masks':toggles,
 'smallest_equal_parity_hostile':parity_hostile,'smallest_equal_merged_hostile':merged_hostile,
 'fraction_compactification':fraction_rows,'absolute_norm_groups':norm_groups,'Fibonacci_unit_conjugation_orbits_met':13,
 'state_lift_samples':[{'N':n,'window':window(n),'face':face(n),'triangular_number':tri(n),'Cassini':Q(1,n)} for n in (1,2,15,55,65)],
 'prime_renaming_numeric_hostile':{'55_image':renamed55,'norm_of_renamed_eraser7':Q(43,1),'is_triangular':istri(renamed55),'Fibonacci_index':fib_index(renamed55)},
 'scope':'Exact proofs in THEORY.md; finite replay supplements them. No intrinsic edge gauge, autonomous finite quotient, unit-Cassini dynamics, or fraction-block node determination asserted.'
}
print(json.dumps(output,sort_keys=True,indent=2))
