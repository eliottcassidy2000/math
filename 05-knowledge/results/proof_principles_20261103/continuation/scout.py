from itertools import combinations

def tournament(n, mask):
    rows = [0]*n
    for i, (u,v) in enumerate(combinations(range(n),2)):
        if mask>>i&1: rows[u] |= 1<<v
        else: rows[v] |= 1<<u
    return rows

def h(rows):
    n=len(rows)
    dp=[[0]*n for _ in range(1<<n)]
    for u in range(n): dp[1<<u][u]=1
    for s in range(1,1<<n):
        for u in range(n):
            z=dp[s][u]
            if not z: continue
            poss=rows[u]&~s
            while poss:
                b=poss&-poss; poss-=b
                dp[s|b][b.bit_length()-1]+=z
    return sum(dp[-1])

def endpoints(rows):
    n=len(rows); dp=[[0]*n for _ in range(1<<n)]
    for u in range(n):dp[1<<u][u]=1
    for s in range(1,1<<n):
        for u in range(n):
            z=dp[s][u]
            if not z: continue
            poss=rows[u]&~s
            while poss:
                b=poss&-poss;poss-=b
                dp[s|b][b.bit_length()-1]+=z
    return dp

def root_profiles(rows):
    n=len(rows);e=endpoints(rows);b=endpoints(reverse(rows))
    hs=[sum(z) for z in e];hs[0]=1
    q=[[0]*(n+1) for _ in range(1<<n)];q[0][0]=1
    for s in range(1,1<<n):
        t=s
        while t:
            rem=s^t
            for k in range(1,s.bit_count()+1):q[s][k]+=hs[t]*q[rem][k-1]
            t=(t-1)&s
    bb=[[0]*(n+1) for _ in range(n)];ee=[[0]*(n+1) for _ in range(n)]
    allmask=(1<<n)-1
    for s in range(1,1<<n):
        rem=allmask^s
        for u in range(n):
            for k in range(1,n+1):
                bb[u][k]+=b[s][u]*q[rem][k-1]
                ee[u][k]+=e[s][u]*q[rem][k-1]
    return bb,ee

def add(rows, mask):
    n=len(rows)
    return [r | ((not(mask>>u&1))<<n) for u,r in enumerate(rows)]+[mask]

def inflate(q, module):
    k=len(module); n=k+len(q)-1
    rows=list(module)+[0]*(len(q)-1)
    for a in range(len(q)):
        aa=range(k) if a==0 else [k+a-1]
        for b in range(len(q)):
            if q[a]>>b&1:
                bb=range(k) if b==0 else [k+b-1]
                for u in aa:
                    for v in bb:rows[u]|=1<<v
    return rows

def cyclic(k):
    return [sum(1<<((u+d)%k) for d in range(1,(k+1)//2)) for u in range(k)]

def reverse(rows):
    n=len(rows)
    return [((1<<n)-1)^r^(1<<u) for u,r in enumerate(rows)]

if __name__=='__main__':
    count=0
    for k in [3,5]:
        m=cyclic(k)
        for qn in range(1,4):
            for qi in range(1<<(qn*(qn-1)//2)):
                q=tournament(qn,qi)
                a=inflate(q,m);b=inflate(q,reverse(m))
                for x in range(1<<len(a)):
                    ha=h(add(a,x));hb=h(add(b,x))
                    if ha!=hb: raise RuntimeError((k,qn,qi,x,a,b,ha,hb))
                    count+=1
        print(k,count,flush=True)
