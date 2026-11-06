#!/usr/bin/env python3
"""Independent literal permutation controls, no author imports."""
from itertools import permutations,product
from math import factorial
import json

def check(x,s):
    if not x: raise RuntimeError(s)
def h(a):
    return sum(all(a[p[i]]>>p[i+1]&1 for i in range(len(p)-1)) for p in permutations(range(len(a))))
def tournament(n,code):
    a=[0]*n;k=0
    for i in range(n):
        for j in range(i+1,n):
            if code>>k&1:a[i]|=1<<j
            else:a[j]|=1<<i
            k+=1
    return a
def extend(a,z):
    b=a[:]+[z];n=len(a)
    for i in range(n):
        if not(z>>i&1):b[i]|=1<<n
    return b
parents=cuts=capacities=0
for n in range(1,5):
    for code in range(1<<(n*(n-1)//2)):
        a=tournament(n,code);H=h(a);Q=[[0]*n for _ in range(n)]
        for p in permutations(range(n)):
            failures=[k for k in range(n-1) if not(a[p[k]]>>p[k+1]&1)]
            if len(failures)>1:continue
            exposed=range(n-1) if not failures else failures
            for k in exposed:Q[p[k]][p[k+1]]+=1
        for i in range(n):
            for j in range(i+1,n):
                left=[v for v in range(n) if v not in (i,j)]
                b=[sum((a[v]>>w&1)<<k for k,w in enumerate(left)) for v in left]
                check(Q[i][j]+Q[j][i]>=2*h(b),'capacity floor');capacities+=1
        for z in range(1<<n):
            predicted=H+sum((Q[i][j]+Q[j][i])*((z>>i)&1)*(1-((z>>j)&1)) for i in range(n) for j in range(n) if a[i]>>j&1)
            check(predicted==h(extend(a,z)),'positive cut identity');cuts+=1
        parents+=1
rows=[];cases=0
for m in range(2,6):
    eligible=sum(1 for bins in product(range(3),repeat=m) if bins.count(0)>0 and bins.count(2)>0 and bins.count(1)>=2)
    formula=3**m-3*2**m+3-m*(2**(m-1)-2)
    check(eligible==formula,'three-bin count')
    for pa,pc,ac in product(range(2),repeat=3):
        total=0
        for bits in product(range(2),repeat=m):
            a=[0]*(m+2)
            for i in range(m):
                for j in range(i+1,m):a[i]|=1<<j
            a[m if ac else m+1]|=1<<(m+1 if ac else m)
            for q,polarity in ((m,pa),(m+1,pc)):
                for i,b in enumerate(bits):
                    if b^polarity:a[q]|=1<<i
                    else:a[i]|=1<<q
            total+=h(a);cases+=1
        check(total*16>=formula*2**m,'mean lower bound')
        rows.append({'m':m,'polarities':[pa,pc],'arc_a_to_c':ac,'sum_H':total,'words':2**m,'bin_count':formula})
m=13;N=3**m-3*2**m+3-m*(2**(m-1)-2)
check(N==1516528 and N//16==94783 and N%16==0,'m13 exact')
check(N>16*(6*(3*2**m+1)//2+1),'m13 upper-sign contradiction')
print(json.dumps({'status':'PASS','method':'literal permutations, independent of author subset DP','parents':parents,'cuts':cuts,'capacity_pairs':capacities,'two_apex_graphs':cases,'two_apex_controls':rows,'m13_bin_count':N,'m13_average_floor':N//16,'plus_target_average':73732,'minus_target_average':73730},indent=2,sort_keys=True))
