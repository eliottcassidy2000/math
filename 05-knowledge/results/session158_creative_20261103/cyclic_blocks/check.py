#!/usr/bin/env python3
from math import factorial, comb

def surj(n,k):
    return sum((-1)**(k-j)*comb(k,j)*j**n for j in range(k+1))

def formula(a,b,c):
    p,q,r=[[surj(n,k) for k in range(max(a,b,c)+2)] for n in (a,b,c)]
    return sum(3*p[k]*q[k]*r[k] + p[k+1]*q[k]*r[k]+p[k]*q[k+1]*r[k]+p[k]*q[k]*r[k+1] +p[k+1]*q[k+1]*r[k]+p[k]*q[k+1]*r[k+1]+p[k+1]*q[k]*r[k+1] for k in range(max(a,b,c)+1))

def dp(a,b,c):
    sizes=[a,b,c]; n=sum(sizes); block=[i for i,size in enumerate(sizes) for _ in range(size)]
    out=[0]*n
    for i in range(n):
        for j in range(n):
            if (block[i]==block[j] and i<j) or (block[i]!=block[j] and (block[j]-block[i])%3==1):out[i]|=1<<j
    d=[[0]*n for _ in range(1<<n)]
    for i in range(n): d[1<<i][i]=1
    for s in range(1,1<<n):
        for i in range(n):
            if d[s][i]:
                opts=out[i] & ~s
                while opts:
                    bit=opts & -opts; j=bit.bit_length()-1
                    d[s|bit][j]+=d[s][i]; opts-=bit
    return sum(d[-1])

for a in range(1,5):
 for b in range(1,5):
  for c in range(1,5):
    x=dp(a,b,c); y=formula(a,b,c)
    if x!=y: raise Exception((a,b,c,x,y))
    if c==2 and x!= 2*3**(a+b)-4*(3**a*2**b+2**a*3**b)+13*2**(a+b)+2*(3**a+3**b)-11*(2**a+2**b)+15: raise Exception(("simplified",a,b,c,x))
    if c==1 and x!=2+(2**a-1)*(2**b-1): raise Exception(('singleton',a,b,c,x))
print('Literal subset-endpoint DP agrees for all 64 triples 1..4; also the c=1 and c=2 formulas')
for n in range(1,8): print('equal blocks:',n,formula(n,n,n))
for b in range(2,9): print('(2,b,2):', b, formula(2,b,2))
for e in range(1,7):
    m=3**e; period=2*3**(e-1)
    found={(4*pow(3,b,m)+5*pow(2,b,m)-11)%m for b in range(max(2,e),max(2,e)+period)}
    expected={r for r in range(m) if r%3!=1}
    if found!=expected: raise Exception(('3-adic',e,found^expected))
print('3-adic residue family verified modulo 3^e for 1 <= e <= 6')
