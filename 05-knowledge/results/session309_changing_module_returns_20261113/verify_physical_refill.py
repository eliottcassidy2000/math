from math import isqrt

def require(b):
    if not b: raise ValueError('certificate failed')

def responses(rows):
    n=len(rows); f=(1<<n)-1
    for i in range(n):
        require(not (rows[i]>>i&1))
        for j in range(i):require(((rows[i]>>j&1)+(rows[j]>>i&1))==1)
    d=[[0]*n for _ in range(f+1)]
    for i in range(n):d[1<<i][i]=1
    for s in range(1,f+1):
        for i,v in enumerate(d[s]):
            if v:
                t=rows[i]&~s
                while t:
                    b=t&-t;t^=b;d[s|b][b.bit_length()-1]+=v
    h=[sum(x) for x in d];h[0]=1
    q=[sum(h[s]*h[m^s] for s in range(m+1) if not s&~m) for m in range(f+1)]
    return h[f],q[f],sum(h[s]*q[f^s] for s in range(f+1))

def strong(rows):
    n=len(rows)
    for i in range(n):
        seen=1<<i
        for _ in range(n):
            for j in range(n):
                if seen>>j&1: seen|=rows[j]
        if seen!=(1<<n)-1:return False
    return True

M=[2,12,9,17,39,15];N=[2,4,25,35,75,87,15];R=[14,28,104,176,133,211,155,7]
h,p,w=responses(M);j,q,z=responses(N);r,s,t=responses(R)
require((h,w,j,q,r,s)==(23,1647,53,844,499,7240))
require(strong(N) and strong(R))
def prime(a):return a>=2 and all(a%d for d in range(2,isqrt(a)+1))
count=[499+180*7240*2**(6-k)*3**k for k in range(7)]
require(count[0]==60*w*q+h*j==83405299)
require(all(2*y==3*x-499 for x,y in zip(count,count[1:])))
require([i for i,z in enumerate(count) if prime(z)]==[0,1,3,6])

# Names v0..v7 for R, 20..25 for the movable right tail, 26 for old u.
def make(k):
    tail=list(range(20,26))
    v=list(range(8))+tail[:6-k];u=tail[6-k:]+[26]
    # Quotient positions a,v,c,u,d0,d1,x0,f0,f1.
    order=[30]+v+[31]+u+[32,33,34,35,36]
    arcs={(a,b) for i,a in enumerate(order) for b in order[i+1:]}
    for i in range(8):
        for j in range(i+1,8):
            if R[j]>>i&1:arcs.remove((i,j));arcs.add((j,i))
    for a,b in ((30,32),(30,33),(31,35),(31,36)):
        arcs.remove((a,b));arcs.add((b,a))
    return arcs
for k in range(6):
    a=make(k);b=make(k+1);w=25-k
    require(a-b=={(w,31)} and b-a=={(31,w)})
print('PASS: prime rebase, four module responses, six physical flips')
