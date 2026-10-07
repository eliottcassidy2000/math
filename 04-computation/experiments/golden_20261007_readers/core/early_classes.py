# list descent survivors mod 2^K and which are merge-certified (n vs n-d at equal time, uniformly) - direct integer check
def T(x): return x//2 if x%2==0 else (3*x+1)//2
def word(r,K):
    w=[];x=r
    for _ in range(K): w.append(x&1); x=T(x)
    return w
def surv(r,K):
    a=0
    for t,b in enumerate(word(r,K),1):
        a+=b
        if 3**a < 2**t: return False
    return True
for K in range(5,11):
    S=[r for r in range(2**K) if surv(r,K)]
    big=2**(K+60)
    cert=[]
    for r in S:
        n=big+r
        for d in range(1,2**K):
            u,v=n-d,n; hit=None
            for t in range(1,K+1):
                u,v=T(u),T(v)
                if u==v: hit=t;break
            if hit:
                # check uniformity: also for n+2^K*s
                ok=True
                for s in (1,7,12345):
                    u,v=n+(s<<K)-d,n+(s<<K)
                    for t in range(hit): u,v=T(u),T(v)
                    ok &= (u==v)
                if ok: cert.append((r,d,hit)); break
    print(K, "survivors",len(S), S if len(S)<=40 else "", "\n   merge-certified (r, d, t):", cert)
