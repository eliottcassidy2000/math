#!/usr/bin/env python3
"""Independent vertex-subset DP for U_m, without root/worker imports."""
from math import comb
from pathlib import Path
import json
OUT=Path(__file__).resolve().parent

def check(ok,msg):
    if not ok:raise AssertionError(msg)

def matrix(m):
    n=3*m;A=[[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1,n):
            if i//3<j//3 or j-i==1:A[i][j]=1
            else:A[j][i]=1
    A[0][-1]=0;A[-1][0]=1
    return A

def formula(m):
    r=m-2;p=[comb(m,k)for k in range(m+1)]
    for k in range(r+1):
        a=comb(r,k)*(9*10**(r-k)-(-2)**(r-k))
        check(a%2==0,'integer formula');p[k+1]+=a//2
    return p

def direct(A):
    n=len(A);N=1<<n;dp=[[0]*n for _ in range(N)];cycles=[0]*N
    for i in range(n):dp[1<<i][i]=1
    for mask in range(1,N):
        root=(mask&-mask).bit_length()-1;row=dp[mask]
        if mask.bit_count()>=3 and mask.bit_count()%2:cycles[mask]=sum(row[j]for j in range(n)if A[j][root])
        available=((N-1)^mask)&~((1<<(root+1))-1)
        for last,ways in enumerate(row):
            if not ways:continue
            rem=available
            while rem:
                bit=rem&-rem;rem^=bit;j=bit.bit_length()-1
                if A[last][j]:dp[mask|bit][j]+=ways
    p=[[1]]+[None]*(N-1)
    for mask in range(1,N):
        bit=mask&-mask;rest=mask^bit;q=list(p[rest]);sub=rest
        while sub:
            cm=sub|bit;w=cycles[cm]
            if w:
                smaller=p[mask^cm]
                if len(q)<=len(smaller):q +=[0]*(len(smaller)+1-len(q))
                for k,a in enumerate(smaller):q[k+1]+=w*a
            sub=(sub-1)&rest
        p[mask]=q
    # Separate ordinary Hamiltonian-path DP, now no minimum-root restriction.
    hp=[[0]*n for _ in range(N)]
    for i in range(n):hp[1<<i][i]=1
    for mask in range(1,N):
        for i,w in enumerate(hp[mask]):
            if not w:continue
            rem=(N-1)^mask
            while rem:
                b=rem&-rem;rem^=b;j=b.bit_length()-1
                if A[i][j]:hp[mask|b][j]+=w
    return p[-1],sum(hp[-1]),sum(cycles)

def main():
    results=[]
    for m in range(2,6):
        p,h,c=direct(matrix(m));f=formula(m)
        check(p==f,('formula',m,p,f));check(h==sum(a*2**k for k,a in enumerate(p)),('OCF',m))
        check(h==(17 if m==2 else 3**m*(1+4**(m-2))),('H',m))
        check(p[-2]==m+4 and p[-3]==comb(m,2)+46*(m-2),'top')
        r={'m':m,'n':3*m,'coefficients':p,'H':h,'odd_cycles':c,'newton_obstruction_numerator':-4*(21*m*m-48*m+4)}
        results.append(r);print(r,flush=True)
    for m in range(3,101):
        p=formula(m)
        check(p[-2]**2*comb(m,2)<m*m*p[-3],'Newton inequality')
    (OUT/'triangle_family.json').write_text(json.dumps({'status':'PASS','direct_DPs':results,'exact_Newton_checks_m3_through100':98},indent=2)+'\n')
if __name__=='__main__':main()
