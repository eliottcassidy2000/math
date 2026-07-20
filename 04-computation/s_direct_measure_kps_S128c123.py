import random
from math import comb
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
def hasrec(a,r,s,p,margin=25):
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
def mins(M,N):
    D=M+N;smax=comb(D,2)+4;K=(D+1)*(smax+2)+D+40
    vals=[]
    random.seed(1234+M*13+N)
    for _ in range(3):
        rc=[random.choice([-3,-2,-1,1,2,3])]+[random.randint(-3,3) for _ in range(D-1)]+[random.choice([-3,-2,-1,1,2,3])]
        a=seq(M,D,rc,K)
        s=None
        for ss in range(smax+1):
            h=[hasrec(a,D,ss,p) for p in PR]
            if None in h: continue
            if all(h): s=ss;break
        vals.append(s)
    return vals
for (M,N) in [(2,2),(2,3),(2,4),(2,5),(3,3),(3,4)]:
    D=M+N;pred=comb(D,2)-min(M,N)+1
    v=mins(M,N)
    mode=max(set(v),key=v.count)
    print("(%d,%d) D=%d  measured s over 3 R: %-14s  mode=%s  pred(C(D,2)-min+1)=%d  match=%s"%(M,N,D,str(v),mode,pred,mode==pred))
