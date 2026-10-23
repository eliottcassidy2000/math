#!/usr/bin/env python3
"""Exact independent controls for square-reset cancellation; no assertions.
All-order argument is in PELL_RESET_ADDENDUM.md. No source/Git writes.
"""
from math import isqrt
from pathlib import Path
import json, hashlib
checks=0

def require(ok,label):
    global checks
    checks+=1
    if not ok:raise RuntimeError(label)

def degree(n):return isqrt(2*n-1)-isqrt(n)

def channels(n):return list(range(isqrt(n)+1,isqrt(2*n-1)+1))

limit=10000
negative=[];m,k=1,1
while k<=limit:
    require(m*m-2*k*k==-1,'negative Pell generator')
    negative.append((m,k));m,k=3*m+4*k,2*m+3*k
positive=[];k,j=1,0
while k<=limit:
    require(k*k-2*j*j==1,'positive Pell generator')
    positive.append((k,j));k,j=3*k+4*j,2*k+3*j
expected={k for m,k in negative if k>=2}|{k for k,j in positive if k>=2}
resets=[]
for k in range(2,limit+1):
    n=k*k
    oldup=isqrt(2*n-3);newup=isqrt(2*n-1)
    up=newup-oldup;down=isqrt(n)-isqrt(n-1)
    require(up in (0,1) and down==1,'square jump bounds')
    require(degree(n)-degree(n-1)==up-down,'difference identity')
    silent=(degree(n)==degree(n-1))
    candidate=any(isqrt(2*n-e)**2==2*n-e for e in (1,2))
    require(silent==candidate,'Pell square-reset equivalence')
    require(silent==(k in expected),'Pell sequence exhausts bounded controls')
    if silent:
        m=newup;rhs=m*m-2*k*k
        require(rhs in (-1,-2),'signed Pell equation')
        require((rhs==-1 and m%2==1) or (rhs==-2 and m%2==0),'correct parity branches')
        old=set(channels(n-1));new=set(channels(n))
        require(old-new=={k} and new-old=={m},'one outgoing and one incoming square channel')
        resets.append({'n':n,'k':k,'incoming_square_root':m,'pell_rhs':rhs,
                       'arrival_degree':degree(n),'outgoing_square':n,'incoming_square':m*m})
# Bounded all-n difference control; the interval proof gives the unbounded claim.
for n in range(2,100001):
    up=isqrt(2*n-1)-isqrt(2*n-3)
    down=isqrt(n)-isqrt(n-1)
    require(up in (0,1),'no multiple upper jumps')
    require(down==int(isqrt(n)**2==n),'only square lower jumps')
    require(degree(n)-degree(n-1)==up-down,'all-n jump identity')
require([r['n'] for r in resets[:4]]==[9,25,289,841],'first four silent resets')
require(channels(24)==[5,6] and channels(25)==[6,7],'n25 square channel exchange')
root=Path(__file__).resolve().parent.parent
source_files=['PELL_ARRIVAL_CLOCKS.md',
              'check_pell_arrivals.py','check_pell_arrivals.json',
              'check_arrival_polynomials.py','check_arrival_polynomials.json']
source_pins={p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in source_files}
producer=json.loads((root/'check_pell_arrivals.json').read_text())
require(producer['silent_square_births']==[r['n'] for r in resets if r['n']<=producer['n_max']],
        'root finite Pell list')
require(producer['Pell_witnesses']==[
    {'n':r['n'],'m':r['incoming_square_root'],'k':r['k'],'norm':r['pell_rhs']}
    for r in resets if r['n']<=producer['n_max']], 'root exact Pell witnesses')
# Certify all provided interval rows by endpoint and maximality checks, without
# repeating the producer's all-n enumeration for each degree.
poly=json.loads((root/'check_arrival_polynomials.json').read_text())
polynomial_rows=[]
for row in poly['rows']:
    j=row['degree'];rmax=j+1+isqrt(2*(j+1)**2)
    supplied={isqrt(L):(L,U) for L,U in row['intervals']}
    require(len(supplied)==len(row['intervals']),'one supplied interval per square band')
    for r in range(1,rmax+1):
        lo,hi=r*r,(r+1)**2-1
        if r not in supplied:
            require(degree(lo)>j or degree(hi)<j,'omitted band is empty')
            continue
        L,U=supplied[r]
        require(lo<=L<=U<=hi,'interval lies in its square band')
        require(degree(L)==j and degree(U)==j,'closed interval endpoints')
        require(L==lo or degree(L-1)!=j,'left endpoint maximal')
        require(U==hi or degree(U+1)!=j,'right endpoint maximal')
    require(all(r<=rmax for r in supplied),'finite cutoff')
    count=sum(U-L+1 for L,U in supplied.values())
    require(count==row['count'],'generating polynomial at one')
    polynomial_rows.append({'degree':j,'count':count,'interval_count':len(supplied)})
require([r['count'] for r in polynomial_rows[:10]]==[3,12,24,34,47,59,68,84,92,105],
        'first ten arrival-level cardinalities')
out={'status':'PASS independent exact Pell and arrival-polynomial controls',
     'explicit_checks':checks,'integer_square_roots_through_k':limit,
     'all_n_difference_control_through':100000,
     'negative_Pell_pairs_m_k':negative,'positive_Pell_pairs_k_j':positive,
     'silent_square_resets':resets,
     'scope':'All-order endpoint and jump proof in PELL_RESET_ADDENDUM.md; bounded controls do not replace it',
     'n25_channels':{'old_squares':[25,36],'new_squares':[36,49],
                      'old_neighbours':[1,12],'new_neighbours':[11,24]}}
out['root_source_sha256']=source_pins
out['arrival_polynomial_endpoint_checks']=polynomial_rows
print(json.dumps(out,indent=2))
