#!/usr/bin/env python3
from math import factorial,gcd
import json
from itertools import permutations

def check(ok,label):
    if not ok:raise RuntimeError(label)

def profile(a):
    n=len(a);N=1<<n;dp=[[0]*n for _ in range(N)];h=[0]*N;h[0]=1
    for i in range(n):dp[1<<i][i]=1
    for S in range(1,N):
        for v in range(n):
            bits=a[v]&~S
            while bits:
                b=bits&-bits;bits-=b;dp[S|b][b.bit_length()-1]+=dp[S][v]
        h[S]=sum(dp[S])
    pc=[[0]*(n+1) for _ in range(N)];pc[0][0]=1
    for S in range(1,N):
        least=S&-S;A=S
        while A:
            if A&least:
                for c in range(1,n+1):pc[S][c]+=h[A]*pc[S^A][c-1]
            A=(A-1)&S
    return pc[-1][1:]

def cyclic_blocks(r,s):
    n=r+s+1;a=[0]*n
    for i in range(r):a[i]=sum(1<<j for j in range(i+1,r+s))
    for i in range(r,r+s):a[i]=sum(1<<j for j in range(i+1,r+s+1))
    a[-1]=(1<<r)-1
    return a

def join(a,b):
    n=len(a);m=len(b)
    return [x|(((1<<m)-1)<<n) for x in a]+[x<<n for x in b]

def content(p):return gcd(*(factorial(c)*factorial(c-1)*v for c,v in enumerate(p,1)))

def main():
    a,b=cyclic_blocks(2,1),cyclic_blocks(2,2)
    records=[]
    for label,T,want in [('W5',a,1),('W11',b,1),('Z5',join(a,b),5)]:
        p=profile(T);g=content(p);check(g==want,label)
        records.append({'name':label,'order':len(T),'profile':p,'context_content':g})
    rotation=0;small=[]
    for n in range(3,34,2):
        T=[sum(1<<((i+s)%n) for s in range(1,(n+1)//2)) for i in range(n)]
        check(all(((T[i]>>j)&1)==((T[(i+1)%n]>>((j+1)%n))&1) for i in range(n) for j in range(n)),('rotation',n));rotation+=1
        if n<=9:
            p=profile(T);g=content(p);check(g%n==0,('cyclic zero',n));rec={'order':n,'H':p[0],'context_content':g}
            if n==5:
                aut=sum(all(((T[i]>>j)&1)==((T[perm[i]]>>perm[j])&1) for i in range(n) for j in range(n)) for perm in permutations(range(n)))
                check(aut==5 and g==15 and p==[15,45,35,10,1],'strict symmetry hostile')
                rec['automorphism_order']=aut;rec['profile']=p
            small.append(rec)
    n=11
    R=[sum(1<<((i+s)%n) for s in range(1,6)) for i in range(n)]
    rp=profile(R);rg=content(rp);check(rg==33 and rg%55!=0,'R11 content boundary')
    S={1,3,4,5,9};P=[sum(1<<((i+s)%11) for s in S) for i in range(11)]
    affines=[tuple((a*x+b)%11 for x in range(11)) for a in S for b in range(11)]
    check(len(set(affines))==55,'affine order')
    check(all(((P[i]>>j)&1)==((P[g[i]]>>g[j])&1) for g in affines for i in range(11) for j in range(11)),'Paley affine action')
    print(json.dumps({'status':'PASS','R11':{'profile':rp,'context_content':rg},'Paley11_affine_subgroup':55,'primitive_product':records,'rotation_orders':rotation,'small_cyclic_zero_states':small},indent=2,sort_keys=True))
if __name__=='__main__':main()
