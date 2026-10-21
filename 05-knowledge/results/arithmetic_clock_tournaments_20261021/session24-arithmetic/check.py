#!/usr/bin/env python3
"""Exact Session24 arithmetic/signed-section and Fano controls. No repo imports."""
from collections import Counter
from itertools import combinations, permutations, product
from math import gcd, prod
import json
from pathlib import Path

checks=0
def need(condition, message):
    global checks
    checks+=1
    if not condition: raise RuntimeError(message)

def order(a,m):
    need(gcd(a,m)==1,'unit order requested')
    x=a%m;k=1
    while x!=1: x=x*a%m;k+=1
    return k

def subgroup(gens,m):
    H={1};q=[1]
    while q:
        x=q.pop()
        for g in gens:
            y=x*g%m
            if y not in H:H.add(y);q.append(y)
    return H

def orbits(H,X,m):
    left=set(X);out=[]
    while left:
        x=min(left);O={h*x%m for h in H};out.append(O);left-=O
    return out

def signed_section_data(H,X,m):
    oo=orbits(H,X,m);flips=[O for O in oo if {-x%m for x in O}==O]
    return {'orbit_lengths':dict(sorted(Counter(map(len,oo)).items())),
            'flipped_orbits':len(flips),
            'section_count':0 if flips else 2**(len(oo)//2)}

N=341
need(N==11*31==(4**5-1)//3,'341 arithmetic')
need(3*N+1==2**10,'trunk terminal exponent')
need(pow(2,N-1,N)==1,'base2 pseudoprimality')
need(gcd(3**8-3**4-1,3**60-1)==N,'golden common quotient gcd')
need((81*81-81-1)%N==0 and pow(81,15,N)==1,'golden ring map')
orders={str(m):{str(a):order(a,m) for a in (2,3,4,81)} for m in (11,31,341)}
for p in (11,31): need(order(4,p)==5,'primitive divisor at exponent5')
need(pow(2,5,N)*31%N==(-31)%N,'doubling antipodal obstruction')
need(pow(3,15,N)*11%N==(-11)%N,'tripling antipodal obstruction')
units={x for x in range(1,N) if gcd(x,N)==1}
clock_data={str(a):{'full':signed_section_data(subgroup([a],N),range(1,N),N),
                    'units':signed_section_data(subgroup([a],N),units,N)} for a in (2,3,4,81)}
need([clock_data[str(a)]['full']['section_count'] for a in (2,3,4,81)]==[0,0,2**34,2**12], 'clock section counts')
R={p:{x*x%p for x in range(1,p)} for p in (11,31)}
H=subgroup([4,81],N)
need(H=={x for x in units if x%11 in R[11] and x%31 in R[31]},'joint group product')
need(len(H)==75,'joint group order')
OO=orbits(H,range(1,N),N);pairs=[];done=set()
for i,O in enumerate(OO):
    if i in done:continue
    j=next(j for j,P in enumerate(OO) if P=={-x%N for x in O})
    need(i!=j,'no flipped joint orbit');pairs.append((i,j));done|={i,j}
need(len(pairs)==4,'four orientation choices')
full=(1<<N)-1

def rows_for(S):return [sum(1<<((v+d)%N) for d in S) for v in range(N)]
def module_closure(rows,d):
    M=1|(1<<d);q=[d]
    while q:
        v=q.pop();D=(rows[0]^rows[v])&~M&full;M|=D
        while D:
            bit=D&-D;D-=bit;q.append(bit.bit_length()-1)
    return M

def naive_closure(S,d):
    M={0,d}
    while True:
        D={x for x in range(N) if x not in M and any(((x-v)%N in S)!=((x)%N in S) for v in M)}
        if not D:return M
        M|=D

sets=[];records=[]
for bits in product((0,1),repeat=4):
    S=set().union(*(OO[pair[b]] for pair,b in zip(pairs,bits)));sets.append(frozenset(S))
    need(len(S)==170 and not(S&{-x%N for x in S}),'skew connection set')
    need(all({a*x%N for x in S}==S for a in (4,81)),'both clocks act')
    rows=rows_for(S)
    closures=[module_closure(rows,d) for d in range(1,N)]
    # Every pair is translation equivalent to 0,d; all 340 such pairs checked.
    sizes=Counter(M.bit_count() for M in closures)
    for O in OO:
        d=min(O);naive=naive_closure(S,d);fast=closures[d-1]
        need(naive=={x for x in range(N) if fast>>x&1},'independent naive/bitset closure')
    mods={p:all(len({x in S for x in range(r,N,p)})==1 for r in range(1,p)) for p in (11,31)}
    need((not any(mods.values()))==(set(sizes)=={N}),'axis coset constancy agrees with complete pair closures')
    # Full graph invariant: profile of common out-neighbours over directed arcs,
    # divided by N using translation invariance (one representative per difference).
    profile=Counter((rows[0]&rows[d]).bit_count() for d in S)
    records.append({'bits':bits,'prime':set(sizes)=={N},'pair_closure_sizes':dict(sorted(sizes.items())),
                    'CRT_coset_module':mods,'arc_common_out_profile':dict(sorted(profile.items()))})
need(sum(r['prime'] for r in records)==8,'8 prime tournaments')
need(sum(r['CRT_coset_module'][11] for r in records)==4,'4 P11[P31] gauges')
need(sum(r['CRT_coset_module'][31] for r in records)==4,'4 P31[P11] gauges')
# Four normal forms, each with independent coordinate sign gauges.
normal=[]
for p,q in ((11,31),(31,11)):
    for delta in (1,-1):
        gauges=[]
        for e,f in product((1,-1),repeat=2):
            def sign(x):
                a,b=x%p,x%q
                if a==0:return f*(1 if b in R[q] else -1)
                z=e*(1 if a in R[p] else -1)
                return delta*z if b==0 else z
            gauges.append(frozenset(x for x in range(1,N) if sign(x)==1))
        need(set(gauges)<=set(sets) and len(set(gauges))==4,'normal-form gauges')
        S=set(gauges[0]);rows=rows_for(S)
        profile=Counter((rows[0]&rows[d]).bit_count() for d in S)
        top=max(profile)
        need(all(((rows[0]&rows[d]).bit_count()==top)==(d%p==0) for d in range(1,N)),'max common-out relation recovers CRT partition')
        normal.append({'dominant_prime':p,'other_prime':q,'axis_sign':delta,
                       'prime':delta==-1,'arc_common_out_profile':dict(sorted(profile.items()))})
need(len({tuple(x['arc_common_out_profile'].items()) for x in normal})==4,'four nonisomorphic normal forms')
# Multiplier orbits already show at most four isomorphism classes.
seen=set();mclasses=[]
for S in sets:
    if S in seen:continue
    orb={frozenset(u*x%N for x in S) for u in units};need(len(orb)==4,'multiplier orbit length')
    seen|=orb;mclasses.append(orb)
need(len(mclasses)==4 and seen==set(sets),'four multiplier classes, complete')
markers={0,155,187}
need({(x%11,x%31) for x in markers}=={(0,0),(1,0),(0,1)},'three-corner CRT labels')
marker_stabilizer=[(a,b) for a in H for b in range(N) if {(a*x+b)%N for x in markers}==markers]
need(marker_stabilizer==[(1,0)],'three-corner regular subset for full prime-form automorphism group')
signed_H=H|{-a%N for a in H}
signed_marker_stabilizer=[(a,b) for a in signed_H for b in range(N) if {(a*x+b)%N for x in markers}==markers]
need(signed_marker_stabilizer==[(1,0)],'three-corner regular subset for ordinary subset action of full signed group')

# Recovered rank-three squareclass/Fano mechanism, independent of the special
# bracket predicate. Basis coordinates are 2,3,11 and multiplication is XOR.
primes=(2,3,11);divisor={x:prod(p for i,p in enumerate(primes) if x>>i&1) for x in range(1,8)}
lines={a:tuple(x for x in range(1,8) if (a&x).bit_count()%2==0) for a in range(1,8)}
need(all(len(L)==3 and L[0]^L[1]^L[2]==0 for L in lines.values()),'Fano lines')
need(all(sum(set(pair)<=set(L) for L in lines.values())==1 for pair in combinations(range(1,8),2)),'Steiner decomposition')
code=set()
for g in range(128):
    c=sum((sum((g>>(x-1))&1 for x in lines[a])%2)<<(a-1) for a in range(1,8));code.add(c)
need(len(code)==16,'incidence rank4')
need(Counter(c.bit_count() for c in code)=={0:1,3:7,4:7,7:1},'Hamming weight enumerator')
base=96;gauge={base^c for c in code}
def fanorows(m):
    rows=[0]*7
    for a,(x,y,z) in lines.items():
        edges=((x,y),(y,z),(z,x)) if m>>(a-1)&1 else ((y,x),(z,y),(x,z))
        for u,v in edges:rows[u-1]|=1<<(v-1)
    return rows

def hp_dp(rows):
    n=len(rows);dp=[[0]*n for _ in range(1<<n)]
    for v in range(n):dp[1<<v][v]=1
    for M in range(1,1<<n):
        for v in range(n):
            if not dp[M][v]:continue
            D=rows[v]&~M
            while D:
                b=D&-D;D-=b;w=b.bit_length()-1;dp[M|b][w]+=dp[M][v]
    return sum(dp[-1])

fanoH=Counter();syndromes=Counter()
P7=[sum(1<<y for y in range(7) if (y-x)%7 in {1,2,4}) for x in range(7)]
perms=list(permutations(range(7)))
for m in range(128):
    rows=fanorows(m);h=hp_dp(rows)
    h2=sum(all(rows[v[i]]>>v[i+1]&1 for i in range(6)) for v in perms)
    need(h==h2,'Fano DP/permutation HP agreement')
    s=0
    for a in range(1,8):
        if (m^base)>>(a-1)&1:s^=a
    need((m in gauge)==(s==0),'Hamming syndrome test')
    need(h==(189 if s==0 else 171),'recovered Fano H split')
    if s:need((m^(1<<(s-1))) in gauge,'unique line correction')
    fanoH[h]+=1;syndromes[s]+=1
# Paley isomorphism and its automorphism count at the stated convention.
rows=fanorows(base)
iso=[v for v in perms if all(((rows[i]>>j)&1)==((P7[v[i]]>>v[j])&1) for i in range(7) for j in range(7))]
need(len(iso)==21,'base Paley7 identification and 21 automorphisms')
for m in gauge:
    r=fanorows(m)
    need(any(all(((r[i]>>j)&1)==((P7[v[i]]>>v[j])&1) for i in range(7) for j in range(7)) for v in perms),'every gauge is Paley7')
# Exact bracket escape proof reduces to n=1,2; enumeration is a control.
qual=[]
for p in range(2,10000):
    if any(p%d==0 for d in range(2,int(p**0.5)+1)):continue
    n=1
    while p>(2*n+1)**2:n+=1
    if 2*p<=(2*n+1)**2:qual.append(p)
need(qual==[2,3,11],'bracket escape control')
result={'status':'FINITE-EXACT controls; elementary proofs in README', 'checks':checks,
        'orders':orders,'clock_sections':clock_data,'joint_group_order':len(H),
        'joint_orbits':[(min(o),len(o)) for o in OO], 'joint_orbit_pairs':pairs,
        'joint_tournaments':records,'normal_forms':normal,'isomorphism_class_count':4,
        'automorphism_subgroup_order':N*len(H),'three_corner_subset':sorted(markers),'three_corner_stabilizer':marker_stabilizer,'signed_three_corner_stabilizer':signed_marker_stabilizer,'signed_group_order_prime_forms':N*len(signed_H),'Fano_divisor_labels':divisor,
        'Fano_lines':{a:[divisor[x] for x in L] for a,L in lines.items()},
        'Fano_H_histogram':dict(sorted(fanoH.items())), 'Fano_syndrome_histogram':dict(sorted(syndromes.items())),
        'bracket_escape_control_under10000':qual}
print(json.dumps(result,indent=2,sort_keys=True))
