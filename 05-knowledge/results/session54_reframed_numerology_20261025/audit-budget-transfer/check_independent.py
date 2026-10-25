from fractions import Fraction as Q
from itertools import product
from math import comb
import json

def need(ok,label):
    if not ok:raise RuntimeError(label)
def matvec(a,x):return [sum(row[j]*x[j] for j in range(len(x))) for row in a]
def rank(a):
    a=[[Q(v) for v in row] for row in a];i=0
    for j in range(len(a[0])):
        pivot=next((k for k in range(i,len(a)) if a[k][j]),None)
        if pivot is None:continue
        a[i],a[pivot]=a[pivot],a[i];v=a[i][j];a[i]=[x/v for x in a[i]]
        for k in range(len(a)):
            if k!=i:
                v=a[k][j];a[k]=[x-v*y for x,y in zip(a[k],a[i])]
        i+=1
        if i==len(a):break
    return i

def solve(a,rhs):
    a=[[Q(v) for v in row]+[Q(rhs[i])] for i,row in enumerate(a)]
    for j in range(len(a)):
        p=next(k for k in range(j,len(a)) if a[k][j]);a[j],a[p]=a[p],a[j]
        v=a[j][j];a[j]=[x/v for x in a[j]]
        for k in range(len(a)):
            if k!=j:
                v=a[k][j];a[k]=[x-v*y for x,y in zip(a[k],a[j])]
    return [row[-1] for row in a]
# Unit digit strings, not producer modular orbit scan.
exact={}
for r in range(1,7):
    counts=[0]*(r+1)
    for ds in product(range(6),repeat=r):
        if ds[0]%3:counts[sum(d>=3 for d in ds)]+=1
    exact[r]=counts
rows=[]
for b in range(9):
    n=b+1
    T=[[3*int(i==j or i==j+1) for j in range(n)] for i in range(n)]
    init=[2*int(j<=1) for j in range(n)];v=init;totals=[]
    for r in range(1,31):
        want=[2*3**(r-1)*(comb(r,j) if j<=r else 0) for j in range(n)]
        need(v==want,'matrix powers and closed formula')
        if r<=6:need(v==(exact[r]+[0]*n)[:n],'independent unit digit enumeration')
        totals.append(sum(v));v=matvec(T,v)
    # Exact Hankel rank proves finite scalar data genuinely excites all dimensions;
    # the general proof is separately recorded in REPORT.md.
    hankel=[[totals[i+j] for j in range(n)] for i in range(n)]
    need(rank(hankel)==n,'full scalar Hankel rank')
    d=totals[:]
    for _ in range(n):d=[y-3*x for x,y in zip(d,d[1:])]
    need(not any(d),'annihilator (E-3)^(b+1)')
    d=totals[:]
    for _ in range(b):d=[y-3*x for x,y in zip(d,d[1:])]
    need(all(d),'no lower annihilator power')
    # Rational resolvent uses generic exact Gaussian elimination.
    for z in [Q(0),Q(1,11),Q(-2,7),Q(2,5)]:
        matrix=[[int(i==j)-z*T[i][j] for j in range(n)] for i in range(n)]
        observed=sum(solve(matrix,init))
        x=3*z
        formula=2/(1-x)+sum(2*x**(j-1)/(1-x)**(j+1) for j in range(1,n))
        need(observed==formula,'rational GF resolvent')
    N=[[int(i==j+1) for j in range(n)] for i in range(n)]
    excited=init[:]
    for _ in range(b):excited=matvec(N,excited)
    need(sum(excited)==2 and not any(matvec(N,excited)),'top Jordan direction')
    rows.append({'budget':b,'initial':init,'hankel_rank':n,'first_four':totals[:4]})
need(exact[1]==[2,2],'unit initial not invented r0')
need(exact[2]==[6,12,6],'all exact states')
need(sum(matvec([[3,0],[0,3]],[2,2]))!=18,'diagonal hostile')
need(sum(matvec([[3,3],[3,3]],[2,2]))!=18,'wraparound hostile')
need(sum(matvec([[3,0],[3,3]],[2,0]))!=18,'wrong initialization hostile')
need(4-6*Q(1,3)==2,'F1 numerator no cancellation')
print(json.dumps({'status':'PASS','digit_enumeration_depth':6,'matrix_levels':30,'budgets':rows,'GF_rational_sample_count_per_budget':4,'hostile_models_rejected':3},sort_keys=True,indent=2))
