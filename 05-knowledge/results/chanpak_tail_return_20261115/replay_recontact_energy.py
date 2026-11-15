from itertools import product
def paths(r):
    n=len(r); f=(1<<n)-1
    a=[[0]*n for _ in range(f+1)]; h=[1]+[0]*f
    for s in range(1,f+1):
        q=s
        while q:
            b=q&-q;i=b.bit_length()-1; t=s^b
            q^=b
            if not t:a[s][i]=1
            else:
                v=t&r[i]
                while v:
                    b2=v&-v;v^=b2;a[s][i]+=a[t][b2.bit_length()-1]
        h[s]=sum(a[s])
    return a,h
def cap(r,m):
    n=len(r)
    return tuple(v|(0 if m>>i&1 else (1<<n)) for i,v in enumerate(r))+(m,)
def ordinal(a,b):
    n=len(a)
    return tuple(x|(((1<<len(b))-1)<<n) for x in a)+tuple(x<<n for x in b)
def st(r,m):
    ap,h=paths(r); n=len(r);f=len(h)-1
    s=d=0
    for z in range(1,f+1):
        v=sum(ap[z][i] for i in range(n) if m>>i&1)*h[f^z]
        s+=v;d+=v if z.bit_count()%2 else -v
    return h[-1],paths(cap(r,m))[1][-1],s,d
def kernel(r):
    ap,h=paths(r);n=len(r);f=len(h)-1
    aa=[[[0]*n for _ in r] for _ in range(f+1)]
    K=[[0]*n for _ in r]
    for m in range(1,f+1):
        for i in range(n):
            if not m>>i&1:continue
            q=m^(1<<i)
            if not q:aa[m][i][i]=1
            else:
                for j in range(n):
                    if (q&r[i])>>j&1:
                        for k in range(n):aa[m][i][k]+=aa[q][j][k]
            if m.bit_count()%2==0:
                for k in range(n):K[i][k]+=2*aa[m][i][k]*h[f^m]
    return K
def reach(r,i):
    m=1<<i
    for _ in r:
        for j,x in enumerate(r):
            if m>>j&1:m|=x
    return m==(1<<len(r))-1
def main():
    A=(82,0,35,87,6,27,54)
    B=(6,24,18,5,9)
    for r in (A,B):
        K=kernel(r)
        for m in range(1<<len(r)):
            h,n,s,d=st(r,m)
            L=[i for i in range(len(r)) if m>>i&1]
            n1=h+s-d-sum(K[i][j]+K[j][i] for i in L for j in L if i<j)
            if n!=n1:raise AssertionError((r,m,n,n1))
    print('160 literal mask/profile energy controls')
    got=[st(A,m) for m in (48,36,80,68)]
    if got!=[(29,x,162,0) for x in (175,177,169,179)]:raise AssertionError(got)
    print('Same-spine exact offsets and mixed second difference:',[t[1] for t in got],8)
    T=ordinal(A,B)
    for ma,mb,desired in [(48,3,13643),(68,5,13645)]:
        m=ma|(mb<<len(A))
        u=st(T,m); t=ordinal(T,(0,)); w=st(t,m)
        if u!=(435,desired,20658,174) or w!=(435,desired+20832,41316,0) or not all(reach(cap(t,m),i) for i in range(len(t)+1)):
            raise AssertionError((u,w))
        print('Bézout matched literal, then strong guard:',u,w)
if __name__=='__main__':main()
