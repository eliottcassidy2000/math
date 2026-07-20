import random
from math import comb, gcd
PR=[(1<<61)-1,(1<<62)-57]
def seq(M,D,rc,K):
    a=[1];Rp=[1]
    for m in range(1,K+1):
        nx=[0]*(len(Rp)+D)
        for i,ci in enumerate(Rp):
            if ci:
                for j,rj in enumerate(rc): nx[i+j]+=ci*rj
        Rp=nx; a.append(Rp[M*m] if M*m<len(Rp) else 0)
    return a
def rk(rows,p):
    rows=[r[:] for r in rows];nc=len(rows[0]);r=0
    for c in range(nc):
        pv=next((i for i in range(r,len(rows)) if rows[i][c]%p),None)
        if pv is None: continue
        rows[r],rows[pv]=rows[pv],rows[r]
        inv=pow(rows[r][c],p-2,p);rows[r]=[(x*inv)%p for x in rows[r]]
        for i in range(len(rows)):
            if i!=r and rows[i][c]%p:
                f=rows[i][c];rows[i]=[(rows[i][k]-f*rows[r][k])%p for k in range(nc)]
        r+=1
        if r==len(rows): break
    return r
def hasrec(a,r,s,p,margin=22):
    nc=(r+1)*(s+1);mm=len(a)-1-r
    if mm+1<nc+margin: return None
    rows=[]
    for m in range(mm+1):
        row=[]
        for i in range(r+1):
            b=a[m+i]%p;mj=1
            for j in range(s+1): row.append(mj*b%p);mj=mj*m%p
        rows.append(row)
    return rk(rows,p)<nc
# discriminating cells: min(M,N) != gcd(M,N)
for (M,N) in [(4,6),(3,6),(4,8)]:
    D=M+N; g=gcd(M,N); mn=min(M,N)
    s_gcd=comb(D,2)-g+1; s_min=comb(D,2)-mn+1
    K=(D+1)*(s_gcd+3)+D+35
    random.seed(99+M*5+N)
    rc=[random.choice([-3,-2,-1,1,2,3])]+[random.randint(-3,3) for _ in range(D-1)]+[random.choice([-3,-2,-1,1,2,3])]
    a=seq(M,D,rc,K)
    at_gcd=all(hasrec(a,D,s_gcd,p) for p in PR)
    below=all(hasrec(a,D,s_gcd-1,p) for p in PR)   # should be False (s_gcd minimal)
    print("(%d,%d) D=%d gcd=%d min=%d | s_gcd_pred=%d s_min_pred=%d | rec@%d:%s rec@%d:%s -> minimal s=%d confirms %s"%(
        M,N,D,g,mn,s_gcd,s_min,s_gcd,at_gcd,s_gcd-1,below,s_gcd,
        "GCD" if (at_gcd and not below) else "??"))
