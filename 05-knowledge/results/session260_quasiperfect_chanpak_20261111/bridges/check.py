"""Small hostile probes for the universal clone congruence; no exhaustive claim beyond stated ranges."""
from functools import lru_cache
from math import factorial,comb
import json
from pathlib import Path

def skeleton(A,v):
    outs=[x for x in range(len(A)) if x!=v]; full=(1<<len(outs))-1
    result=[]
    for r in range(1,len(outs)//2+2):
        @lru_cache(None)
        def f(mask,last,k):
            ans=int(mask==full and k==r)
            if k<r and (last<0 or A[last][v]):
                ans+=f(mask,v,k+1)
            for j,x in enumerate(outs):
                if not mask>>j&1 and (last<0 or A[last][x]):
                    ans+=f(mask|1<<j,x,k)
            return ans
        result.append(f(0,-1,0))
    return result

def F(m,r,mod=None):
    z=sum((-1)**(r-j)*comb(r,j)*(pow(j,m,mod) if mod else j**m) for j in range(r+1))
    return z%mod if mod else z

def root_test(r):
    n=2*r-1; A=[[int(i<j) for j in range(n)] for i in range(n)]
    for j in range(1,n): A[0][j]=int(j%2==1); A[j][0]=1-A[0][j]
    return A

def clone(A,v,m):
    old=[i for i in range(len(A)) if i!=v]; names=old+[v]*m
    return [[int(i<j) if a==b else A[a][b] for j,b in enumerate(names)] for i,a in enumerate(names)]

def hp(A):
    n=len(A)
    @lru_cache(None)
    def f(mask,last):
        if mask==(1<<n)-1:return 1
        return sum(f(mask|1<<j,j) for j in range(n) if not mask>>j&1 and A[last][j])
    return sum(f(1<<i,i) for i in range(n))

def factors(m):
    out=[]; p=2
    while p*p<=m:
        if m%p==0:
            e=0
            while m%p==0:m//=p;e+=1
            out.append((p,e))
        p+=1
    if m>1:out.append((m,1))
    return out

def check(q,msg):
    if not q:raise RuntimeError(msg)

rows=[]; literal=0
for r in range(2,7):
    A=root_test(r); a=skeleton(A,0)
    check(len(a)==r and a[-1]==factorial(r-1),('top',r,a))
    rows.append({'r':r,'coefficients':a})
    if r<=4:
        for m in range(1,5):
            pred=sum(x*F(m,j+1) for j,x in enumerate(a))
            check(pred==hp(clone(A,0,m)),('literal',r,m));literal+=1
# Every labelled base through order four, every root, four clone sizes.
all_controls=0
for n in range(1,5):
    pairs=[(i,j) for i in range(n) for j in range(i+1,n)]
    for mask in range(1<<len(pairs)):
        A=[[0]*n for _ in range(n)]
        for b,(i,j) in enumerate(pairs):A[i][j]=(mask>>b)&1;A[j][i]=1-A[i][j]
        for v in range(n):
            a=skeleton(A,v)
            cyc=any(A[v][i] and A[i][j] and A[j][v] for i in range(n) for j in range(n))
            check((len(a)>1 and a[1]>0)==cyc,('cyclic',n,mask,v))
            for m in range(1,5):
                check(sum(x*F(m,j+1) for j,x in enumerate(a))==hp(clone(A,v,m)),('all',n,mask,v,m))
                all_controls+=1
# Finite-difference criterion independently matches Korselt on 2..600.
# r up to largest prime divisor suffices by theorem's triangular proof.
accepted=[]
for m in range(2,601):
    ff=factors(m); R=max(p for p,e in ff)
    got=all(F(m,r,m)==0 for r in range(2,R+1))
    want=all(e==1 and (m-1)%(p-1)==0 for p,e in ff)
    check(got==want,('criterion',m));
    if got:accepted.append(m)
controls=[]
for m in [4,9,15,20,104,341,464,561,650,1105,1729,1952]:
    witness=next((r for r in range(2,7) if sum(a*F(m,j+1,m) for j,a in enumerate(rows[r-2]['coefficients']))%m != rows[r-2]['coefficients'][0]%m),None)
    controls.append({'m':m,'korselt':all(e==1 and (m-1)%(p-1)==0 for p,e in factors(m)),'first_test_r_through_6':witness})
out={'status':'FINITE-EXACT','alternating_root_coefficients':rows,'direct_alternating_controls':literal,'all_labelled_root_clone_controls_n_le_4':all_controls,'finite_difference_criterion_range':[2,600],'accepted_composites':[m for m in accepted if len(factors(m))>1 or factors(m)[0][1]>1],'hostiles_and_positive_controls':controls}
Path(__file__).with_name('receipt.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
