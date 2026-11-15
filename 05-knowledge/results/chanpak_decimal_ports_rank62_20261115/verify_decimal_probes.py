from math import factorial, comb
from functools import lru_cache

def build(n, arc):
    return tuple(sum((1<<j) for j in range(n) if i!=j and arc(i,j)) for i in range(n))

def probes():
    yield build(3, lambda i,j:(j-i)%3==1),0
    yield build(4, lambda i,j:i<j if {i,j}!={0,3} else i==3),1
    for r in [3,4,5]:
        t=r-1; n=2*t+1
        def arc(i,j):
            if i==0:return j<=t
            if j==0:return i>t
            if (i<=t)==(j<=t):return i<j
            if i<=t:return j==t+i
            return i!=t+j
        yield build(n,arc),0

def coefficients(out,root):
    n=len(out); v=[i for i in range(n) if i!=root]
    R=(n+1)//2
    FULL=((1<<n)-1)^(1<<root)
    @lru_cache(None)
    def go(mask,last,r):
        h=[0]*5
        if mask==FULL and r:
            h[r-1]=1
        for b in v+[root]:
            if last!=-1 and not (out[last]>>b&1):continue
            if b==root:
                if r>=R:continue
                x=go(mask,b,r+1)
            else:
                if mask>>b&1:continue
                x=go(mask|(1<<b),b,r)
            for k in range(5):h[k]+=x[k]
        return tuple(h)
    return go(0,-1,0)

def compute():
    M=[];c=[]
    for out,p in probes():
        row=coefficients(out,p);c.append(row)
        assert all(row[r]%factorial(r)==0 for r in range(5))
        M.append([row[r]//factorial(r) for r in range(5)])
    # row substitution from top formal degree down; bottom 2x2 inverse has det 1.
    v=[24*comb(5,i+1)//factorial(i) for i in range(5)]
    z=[0]*5
    for r in range(4,1,-1):
        z[r]=v[r]
        for j in range(5):v[j]-=z[r]*M[r][j]
    z[0]=2*v[0]-5*v[1]; z[1]=3*v[1]-v[0]
    assert [sum(z[i]*M[i][j] for i in range(5)) for j in range(5)]==[120,240,120,20,1]
    print('templates',c,'matrix',M,'probe coefficients',z)

def fill(out, root, m):
    # Copy of root replaced by an increasing transitive interval; all other
    # labels retain the probe orientation.
    inputs=[i for i in range(len(out)) if i!=root]+[root]*m
    def arc(i,j):
        u,v=inputs[i],inputs[j]
        return i<j if u==v else (out[u]>>v)&1
    return build(len(inputs),arc)

def count(out):
    n=len(out)
    dp=[[0]*n for _ in range(1<<n)]
    for v in range(n):dp[1<<v][v]=1
    for S in range(1,1<<n):
        for v,x in enumerate(dp[S]):
            if not x:continue
            mask=out[v]&~S
            while mask:
                b=mask&-mask
                mask-=b
                dp[S|b][b.bit_length()-1]+=x
    return sum(dp[-1])

def literal_check():
    probes0=list(probes())
    left=[count(fill(a,r,3)) for a,r in probes0[:2]]
    right=[count(fill(a,r,6)) for a,r in probes0]
    charge=sum(a*b for a,b in zip((652,-566,149,-20,1),right))
    assert left[1]-left[0]==8
    assert charge==24*5**6
    assert (left[1]-left[0])*charge==3*1000**2
    print("Literal ports at decimal depth 2 passed:",left,right)

if __name__=='__main__':
    compute()
    literal_check()
