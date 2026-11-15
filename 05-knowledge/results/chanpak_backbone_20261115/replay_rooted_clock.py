"""Literal subset-DP checks for both rooted clocks and arbitrary cores."""
from itertools import combinations, product
from fractions import Fraction
from random import Random
from replay_backward_chain import literal, formula, upper

def require(x, *details):
    if not x:
        raise RuntimeError(details)

def subpaths(rows):
    n=len(rows)
    a=[[0]*n for _ in range(1<<n)]
    for v in range(n):
        a[1<<v][v]=1
    for s in range(1,1<<n):
        for u,x in enumerate(a[s]):
            if not x: continue
            nxt=rows[u]&~s
            while nxt:
                b=nxt&-nxt;nxt^=b
                a[s|b][b.bit_length()-1]+=x
    return a

def all_h(rows):
    return [1]+[sum(e) for e in subpaths(rows)[1:]]

def forward(n):
    return [((1<<n)-1)^((1<<(i+1))-1) for i in range(n)]

def reverse(rows,a,b):
    rows[a]^=1<<b; rows[b]^=1<<a

def extension(rows,m,r,u=None):
    n=len(rows); N=n+m+1+(u is not None)
    res=forward(N); mask=(1<<n)-1
    for i,row in enumerate(rows):
        res[i]=(res[i]&~mask)|row
    if u is not None:
        reverse(res,u,N-2)
    reverse(res,r,N-1)
    return res

def chart(gs):
    last=0; pairs=[]
    for m in gs:
        pairs.append((last,last+1+m))
        last+=1+m
    rows=forward(last+1)
    for a,b in pairs: reverse(rows,a,b)
    return rows

def update(m,state):
    f,u,r=state
    a,b=2**m,3**m
    return (2*a-1)*f+2*b*u+(2*b+2)*r,a*f,b*(u+r)

def rooted3(gs):
    state=(1,0,0)
    for m in gs: state=update(m,state)
    return state

def twoport(rows,r,u):
    mask=(1<<len(rows))-1
    h=all_h(rows)
    starts=[]
    for x in (r,u):
        p=[[0]*len(rows) for _ in h]; p[1<<x][x]=1
        for s in range(1,mask+1):
            for y,c in enumerate(p[s]):
                if not c:continue
                nxt=rows[y]&~s
                while nxt:
                    t=nxt&-nxt;nxt^=t
                    p[s|t][t.bit_length()-1]+=c
        starts.append(p)
    one,two=[[sum(v) for v in z] for z in starts]
    qr=sum(one[s]*h[mask^s] for s in range(1,mask+1))
    qu=sum(two[s]*h[mask^s] for s in range(1,mask+1))
    a=b=c=0
    for s in range(1,mask+1):
        if not one[s]:continue
        sg=(-1)**s.bit_count()
        if sg<0:c+=starts[0][s][u]*h[mask^s]
        t=mask^s
        while t:
            q=one[s]*two[t]*h[mask^s^t]
            a+=q;b+=q*(sg+(-1)**t.bit_count())
            t=(t-1)&(mask^s)
    return 2*qr+qu,2*a,h[-1]+b-2*c

def slots(rows,k):
    h=all_h(rows)
    d=[0]*len(h);d[0]=1
    for _ in range(k):
        l=[0]*len(h)
        for m in range(len(h)):
            s=m
            while True:
                l[m]+=h[s]*d[m^s]
                if not s: break
                s=(s-1)&m
        d=l
    return d[-1]

def ordinal(a,b):
    N=len(a)
    return [x|(((1<<len(b))-1)<<N) for x in a]+[x<<N for x in b]

def ceilpow(num,den,p,q,scale=1000000):
    x=scale**q*num**p;y=den**p
    lo,hi=0,scale+1
    while lo+1<hi:
        m=(lo+hi)//2
        if m**q*y>=x:hi=m
        else:lo=m
    require((hi-1)**q*y<x<=hi**q*y,num,den)
    return hi

def cert(p,q,R,M,limit):
    z=Fraction(3,5)
    require(2**p*3**q>5**q)
    s=1000000
    subtotal=0
    state=(1,0,0)
    for r in range(R):
        subtotal+=sum(ceilpow(1,update(m,state)[0],p,q,s)
                      for m in range(1,M))
        state=update(0,state)
    tail=2*z**(R//2)/(1-z)*z/(1-z)+2/(1-z)*z**M/(1-z)
    result=Fraction(subtotal,s)+tail
    require(result<limit,subtotal,tail,limit)
    print(f"rooted word block mass ({p}/{q}): finite upward sum {subtotal}/{s}, tail bound {tail}, total <{limit}")

def main():
    ctr=0
    for length in range(7):
        for gs in product(range(4),repeat=length):
            if sum(gs)+length+2>13:continue
            rows=chart(gs)
            h=literal(gs)
            n=all_h(extension(rows,0,0))[-1]
            state=rooted3(gs)
            require(n==h+state[0]-1,gs,h,n,state)
            ss=state
            for m in reversed(gs):
                f,u,r=ss
                a=(u&-u).bit_length()-1
                require(a==m,gs,a,m)
                f0=u>>a
                r0=(f-((2<<a)-1)*f0-2*r)//2
                u0=r//3**a-r0
                ss=(f0,u0,r0)
            require(ss==(1,0,0),gs,ss)
            ctr+=1
    print(f"literal rooted path DP and backward state decoding: {ctr} words")
    rng=Random(80903)
    ctr=0
    for n in range(2,8):
        for _ in range(20):
            r,u=rng.sample(range(n),2)
            rows=forward(n)
            for a,b in combinations(range(n),2):
                if rng.randrange(2):reverse(rows,a,b)
            A,B,C=twoport(rows,r,u)
            for m in (1,2,3):
                got=all_h(extension(rows,m,r,u))[-1]
                require(got==C+2**m*A+3**m*B,rows,r,u,m,A,B,C,got)
                ctr+=1
    print(f"arbitrary-core two-port literal controls: {ctr}")
    ctr=0
    T=chart((1,))
    for n in range(1,6):
        for _ in range(8):
            rows=forward(n)
            for a,b in combinations(range(n),2):
                if rng.randrange(2):reverse(rows,a,b)
            for u,sg in ((1,1),(2,-1)):
                for m in (1,2):
                    got=all_h(extension(ordinal(T,rows),m,0,u))[-1]
                    need=9*2**m*slots(rows,2)+4*3**m*slots(rows,3)+sg*all_h(rows)[-1]
                    require(got==need,rows,u,sg,m,got,need)
                    ctr+=1
    print(f"ordinal-tail clock literal controls: {ctr}")
    cert(24,25,40,20,Fraction(99,100))

if __name__=='__main__':main()
