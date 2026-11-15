import math, random
B,Q,T,S,h = 2855,5760,4034,3914,45

def req(value, *details):
    if not value:
        raise AssertionError(details)

def p(x):
    total=0
    while x:
        u,b=divmod(x,2)
        if not b:
            total+=u+(u&1)
        x=u
    return total

def reverse(x,r):
    z=0
    for _ in range(r):
        z=2*z+(x&1)
        x//=2
    return z

def output(x,y):
    return B+Q*x*y+T*x+S*y+h*(p(x)+p(y))

def place(N):
    req(N>=100_000_000)
    a=math.isqrt((N-B)//Q)
    a+=Q*a*a<N-B
    r=a.bit_length()//4
    q=1<<r
    u,v=divmod(a,q)
    base=B+Q*a*a+(T+S)*a+h*(2*p(v)+2*u*(q-1-reverse(v,r)))
    err=h*(2*u+2)
    def f(i):
        z=q*i
        return base-Q*z*z+(T-S)*z
    def first_below(target):
        req(f(0)>=target and f(u)<target,target,N,a,q)
        lo,hi=0,u
        while hi-lo>1:
            mid=(lo+hi)//2
            if f(mid)>=target: lo=mid
            else: hi=mid
        return hi
    lo=first_below(N-err)
    hi=first_below(N)-1
    counts=[]
    for i in (lo,hi):
        x,y=a+q*i,a-q*i
        value=output(x,y)
        req(0<=value-f(i)<=err,N,i,value,f(i),err)
        counts.append((value,x,y))
    req(counts[0][0]<N<=counts[1][0],counts,N)
    return counts

def lit(x,y):
    # Fixed core, two guards, and canonical address tapes.
    rows=(4,1,130,7,143,159,63,75)
    left=bin(x)[2:] if x else ''
    right=''.join('0' if (y>>i)&1 else '1' for i in range(y.bit_length()))
    marks=left+'1'+''.join(str((35>>i)&1) for i in range(8))+'0'+right
    n=len(marks)+1
    start=len(left)+1
    adj=[0]*n
    for i in range(n-1):
        for j in range(i+1,n-1):
            f=not(start<=i<j<start+8) or ((rows[i-start]>>(j-start))&1)
            a,b=(i,j) if f else (j,i)
            adj[a]|=1<<b
        a,b=(n-1,i) if marks[i]=='1' else (i,n-1)
        adj[a]|=1<<b
    size=1<<n
    dp=[None]*size
    for i in range(n):
        dp[1<<i]=[0]*n
        dp[1<<i][i]=1
    for mask in range(1,size-1):
        if dp[mask] is None: continue
        for i,d in enumerate(dp[mask]):
            if not d: continue
            nxt=adj[i]&~mask
            while nxt:
                bit=nxt&-nxt
                nxt-=bit
                new=mask|bit
                if dp[new] is None: dp[new]=[0]*n
                dp[new][bit.bit_length()-1]+=d
    return sum(dp[-1])

if __name__=='__main__':
    tests=[(0,0),(2,3),(3,2),(5,7),(8,9),(15,15),(17,3),(2,21)]
    for x,y in tests:
        req(lit(x,y)==output(x,y),x,y)
    for n in (100_000_001,100_000_000_001,10**30+39):
        print('target',n,'bracket',place(n))
    rng=random.Random(45)
    for k in range(27,1025):
        for _ in range(5):
            place((1<<k)+rng.getrandbits(k)|1)
    print('literal graph checks',len(tests))
    print('exact brackets over 998 magnitudes, 4990 targets, through 1025 bits')

def complement_check(q,u):
    r=q.bit_length()-1
    # Center the pair at q*u+(q-1)/2.
    base=h*(u*(q-1)+2*p(u))
    values=[]
    for v in range(q//2,q):
        x=q*u+v
        y=q*u+q-1-v
        F=B+Q*x*y+T*x+S*y+base
        value=output(x,y)
        req(value-F==h*(p(v)+p(q-1-v)+(u&1)),q,u,v)
        req(0<=value-F<=h*q,q,u,v)
        values.append(value)
    return values

if __name__=='__main__':
    n_direct=0
    for r in range(1,15):
        q=1<<r
        for u in (0,1,q,q+1,2*q-1):
            outputs=complement_check(q,u)
            req(all(abs(s-t)<=(Q+abs(T-S)+h)*q for s,t in zip(outputs,outputs[1:])),q,u)
            n_direct+=len(outputs)
    print('complemented low digits and step bounds:',n_direct,'physical formula checks')
