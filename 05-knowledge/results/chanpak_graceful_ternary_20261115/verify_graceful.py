"""Independent finite controls for the graceful buffer and single-apex formulas."""
from itertools import product

def need(x):
    if not x: raise RuntimeError('failed identity')

def cap(m,k):
    if m<=0:return 0
    r=(m+k-1)//k
    return r*m-k*r*(r-1)//2

def decode(m,k,x):
    ans=[]
    while x:
        if x<=cap(m-1,k):m-=1
        else:
            ans.append(m)
            x-=m;m-=k
    return ans

def hneedle(w):
    n=len(w)+1
    out=[0]*n
    for i,b in enumerate(w):
        out[i]=((1<<(n-1))-(1<<(i+1)))|((1<<(n-1)) if not b else 0)
        out[n-1]|=b<<i
    dp=[[0]*n for _ in range(1<<n)]
    for i in range(n):dp[1<<i][i]=1
    for s in range(1,1<<n):
        for v,t in enumerate(dp[s]):
            if not t:continue
            candidates=out[v]&~s
            while candidates:
                b=candidates&-candidates
                dp[s|b][b.bit_length()-1]+=t
                candidates-=b
    return sum(dp[-1])

def c(d):return 1 if d==1 else 1<<(d-2)

def h(w):
    return 1+2*sum(c(j-i)*w[i]*(1-w[j]) for i in range(len(w)) for j in range(i+1,len(w)))

def bval(w):
    v=0
    for z in w:v=2*v+z
    return v

def hint(a):
    return h(tuple(map(int,bin(a)[2:])))

for m in range(1,24):
    for k in range(1,m+2):
        for r in range(cap(m,k)+1):
            ans=decode(m,k,r)
            need(sum(ans)==r and all(y-x>=k for y,x in zip(ans,ans[1:])))
for n in range(1,8):
    for w in product((0,1),repeat=n):
        H=h(w)
        need(H==hneedle(w))
        d=bval(tuple(1-b for b in w[::-1]))
        need(h((0,0)+w)==H and h((1,0)+w)==H+2+2*d)
        for i in range(1,n):
            a=bval(w[:i])
            b=bval(tuple(1-b for b in w[i:][::-1]))
            need(H==hint(a)+hint(b)-1+a*b+(a&1)*(b&1))
for b in range(1,512,2):
    d=hint(b+1)-hint(b)
    need(d>=2 and (d==2)==(b in (1,5)))
for a in range(0,64,2):
    for b in range(1,64,2):
        u,v=b+1-2*((a+3)//4),a+hint(b+1)-hint(b)
        need(set((abs(u),abs(v)))!={2,4})
print('graceful bounded-spacing decoders, literal Hamilton DP, exact readers, and additive-rectangle controls passed')
