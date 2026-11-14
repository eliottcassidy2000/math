from math import comb, factorial
from random import Random

def surjections(n,k):
    return sum((-1)**i*comb(k,i)*(k-i)**n for i in range(k+1)) if k else 0

def trans(n):
    return tuple(sum(1<<j for j in range(i+1,n)) for i in range(n))

def join(a,b,c):
    blocks=(a,b,c)
    starts=(0,len(a),len(a)+len(b))
    out=[]
    for i,block in enumerate(blocks):
        j=(i+1)%3
        cross=((1<<len(blocks[j]))-1)<<starts[j]
        out.extend((row<<starts[i])|cross for row in block)
    return tuple(out)

def paths(rows):
    n=len(rows)
    dp=[[0]*n for _ in range(1<<n)]
    for i in range(n): dp[1<<i][i]=1
    for s in range(1,1<<n):
        for i in range(n):
            if not(s>>i&1): continue
            dp[s][i]+=sum(dp[s^(1<<i)][j] for j in range(n) if rows[j]>>i&1)
    return [1]+[sum(v) for v in dp[1:]]

def profile(rows):
    H=paths(rows)
    step=[1]+[0]*(len(H)-1)
    G=[0]
    for _ in rows:
        new=[0]*len(H)
        for s in range(1,len(H)):
            t=s
            while t:
                new[s]+=step[s^t]*H[t]
                t=(t-1)&s
        step=new
        G.append(step[-1])
    return G

def weights(n,L):
    S=[surjections(L,j) for j in range(n+2)]
    return [0]+[S[j+1]**2+2*S[j+1]*S[j]+3*S[j]**2+2*S[j]*S[j-1]+S[j-1]**2 for j in range(1,n+1)]

def digits(vec,w):
    H=sum(a*b for a,b in zip(vec,w))
    out=[0]*len(w)
    for j in range(len(w)-1,0,-1):
        out[j],H=divmod(H,w[j])
    assert H==0
    return out

def coefficient_table(rows):
    n=len(rows)
    counts=[0]*(n+1)
    def dfs(visited,last,j):
        if visited==((1<<n)-2) and j:
            counts[j]+=1
        if j<n and (last<0 or rows[last]&1): dfs(visited,0,j+1)
        for i in range(1,n):
            if not(visited>>i&1) and (last<0 or rows[last]>>i&1):
                dfs(visited|1<<i,i,j)
    dfs(0,-1,0)
    return counts

if __name__=='__main__':
    rg=Random(705)
    physical=0
    for n in range(1,8):
        B=factorial(n)*2**(n-1)
        K=1+(72*n*(n+1)*B-1).bit_length()
        w=weights(n,(n+1)*K)
        assert all(B*sum(w[:j])<w[j] for j in range(1,n+1))
        for _ in range(10):
            rows=[0]*n
            for i in range(n):
                for j in range(i+1,n):
                    a,b=(i,j) if rg.randrange(2) else (j,i)
                    rows[a]|=1<<b
            G=profile(tuple(rows))
            assert all(0<=a<=B for a in G)
            assert digits(G,w)==G
            trial=[0]+[rg.randrange(B+1) for _ in rows]
            assert digits(trial,w)==trial
            if n<=5:
                for L in range(1,4):
                    assert paths(join(rows,trans(L),trans(L)))[-1]==sum(a*b for a,b in zip(G,weights(n,L)))
                    physical+=1
    rows=(0,)
    for d in range(1,5):
        rows=join(rows,(0,),(0,))
        co=coefficient_table(rows)
        for k in range(d+1):
            assert co[k+1]>=factorial(k)*comb(d,k)*3**(d-k)
        assert co[d+1]==factorial(d) and not any(co[d+2:])
    print('PASS:',physical,'literal reader controls, 70 common profiles and four nesting depths')
