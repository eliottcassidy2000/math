from math import gcd
W=[(11,1158,32,220,721),(11,1270,32,472,721),(11,1244,478,32,959),
   (11,1574,988,32,959),(13,4614,956,64,2881),(13,5110,1976,64,2881)]
def need(x):
    if not x:raise RuntimeError("rectangle or calibration mismatch")
def h(w):
    a=0
    while w:
        w,b=divmod(w,2)
        if not b:a+=w+(w&1)
    return 1+a
def dp(w,n):
    row=[0]*(n+1)
    for i in range(n):
        row[i]=((1<<n)-1)^((1<<(i+1))-1)
        if not (w>>(n-1-i))&1:row[i]|=1<<n
        else:row[n]|=1<<i
    A=[[0]*(n+1) for _ in range(1<<(n+1))]
    for i in range(n+1):A[1<<i][i]=1
    for m in range(1,len(A)):
        for i,x in enumerate(A[m]):
            if not x: continue
            b=row[i]&~m
            while b:
                k=b&-b;b^=k
                A[m|k][k.bit_length()-1]+=x
    return sum(A[-1])
def sig(w,A,n):
    f=lambda x:int(f"{x:0{n}b}"[::-1],2)
    return f(w)-f(w^A),(w^A)-w
for n,w,A,B,H in W:
    need((A&B)==0 and [dp(w^m,n) for m in (0,A,B,A^B)]==[H,H+2,H+4,H+6])
lanes=[[(0,2,9,10,0,0,2,32),(0,-2,5,2,1,0,30,96)],
       [(0,2,11,14,0,0,2,48),(1,2,13,10,12,9,92,96)],
       [(0,2,9,10,0,0,2,32),(0,-2,11,14,7,9,66,96)]]
for (n,w,A,B,H),L in zip(W[::2],lanes):
    f=[(*sig(w,A,n),2),(*sig(w,B,n),4)]
    for order in range(2):
        a,b,c=f[order];d,e,f0=f[1-order]
        for k in (-2,2):
            x,y,z=d-k*a,e-k*b,k*c-f0
            g=gcd(x,y)
            l=[r for r in L if r[:2]==(order,k)]
            need(len(l)==int(z%g==0))
            if not l:continue
            _,_,A0,B0,aa,bb,m0,m1=l[0]
            need(x*aa+y*bb==z and (A0,B0)==(abs(y)//g,abs(x)//g))
            need(aa>=0 and bb>=0 and (aa<A0 or bb<B0))
            for t in (0,1,2,100,1<<130):
                u,v=aa+A0*t,bb+B0*t
                P=2+sum(x*y for x,y in zip(sig(w,A,n),(u,v)))
                Q=4+sum(x*y for x,y in zip(sig(w,B,n),(u,v)))
                m=m0+m1*t
                need(sorted([0,P,Q,P+Q])==list(range(min(0,P,Q,P+Q),min(0,P,Q,P+Q)+4*m,m)))
print("six independent Hamilton-DP ledgers and all exterior calibration equations")
