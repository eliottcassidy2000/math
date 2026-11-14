from collections import defaultdict
from itertools import permutations
from random import Random

def orient(ap, word):
    k=len(ap)
    n=k+len(word)
    a=[[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1,n):
            if j<k:
                yes=ap[i][j]
            elif i<k:
                yes=word[j-k]>>i&1
            else:
                yes=1
            a[i][j],a[j][i]=yes,1-yes
    return a

def subset_dp(a):
    n=len(a)
    if not n:
        return 1
    dp=[[0]*n for _ in range(1<<n)]
    for i in range(n):
        dp[1<<i][i]=1
    for s in range(1,1<<n):
        rem=((1<<n)-1)^s
        for i,f in enumerate(dp[s]):
            if f==0:
                continue
            v=rem
            while v:
                bit=v&-v
                j=bit.bit_length()-1
                if a[i][j]:
                    dp[s|bit][j]+=f
                v-=bit
    return sum(dp[-1])

def strong(a):
    for i in range(len(a)):
        seen={i}
        queue=[i]
        for j in queue:
            for k,x in enumerate(a[j]):
                if x and k not in seen:
                    seen.add(k);queue.append(k)
        if len(seen)!=len(a):
            return False
    return True

def count(ap, word):
    k=len(ap)
    ans=0
    for p in permutations(range(k)):
        left=(None,*p)
        right=(*p,None)
        empty=tuple(u is None or v is None or ap[u][v]
                    for u,v in zip(left,right))
        states={(0,)*(k+1):1}
        for m in word:
            new=defaultdict(int)
            for i,(u,v) in enumerate(zip(left,right)):
                start=u is None or (m>>u&1)
                close=v is None or not (m>>v&1)
                for s,c in states.items():
                    if s[i] or start:
                        z=2 if close else 1
                        new[s[:i]+(z,)+s[i+1:]]+=c
            states=new
        ans+=sum(c for s,c in states.items()
                 if all(z==2 or (z==0 and ok) for z,ok in zip(s,empty)))
    return ans

def apex_orient(rng,k):
    a=[[0]*k for _ in range(k)]
    for i in range(k):
        for j in range(i+1,k):
            a[i][j]=rng.randrange(2)
            a[j][i]=1-a[i][j]
    return a

def annihilator(k,L):
    p=[1]
    K=k+1
    for j in range(K+1):
        for _ in range(K-j+1):
            q=[0]*(len(p)+1)
            for i,b in enumerate(p):
                q[i]-=b*j**L
                q[i+1]+=b
            p=q
    return p

rng=Random(87421)
sample=0
for k in range(4):
    for _ in range(40):
        a=apex_orient(rng,k)
        w=tuple(rng.randrange(1<<k) for i in range(rng.randrange(8)))
        T=orient(a,w)
        h=count(a,w)
        assert h==subset_dp(T)
        n=len(T)
        if n>=3 and n>=3*k and strong(T):
            assert h>=1+(1<<(n-3*k))
        sample+=1
    for L in (1,2,3):
        coeff=annihilator(k,L)
        seq_len=len(coeff)+4
        for _ in range(4):
            a=apex_orient(rng,k)
            prefix,suffix=[tuple(rng.randrange(1<<k) for _ in range(2))
                           for side in range(2)]
            period=tuple(rng.randrange(1<<k) for _ in range(L))
            vals=[count(a,prefix+period*t+suffix) for t in range(seq_len)]
            assert all(sum(vals[i+j]*x for j,x in enumerate(coeff))==0
                       for i in range(seq_len-len(coeff)+1))
print(f"{sample} exact transfer comparisons and 48 repeat templates passed")
