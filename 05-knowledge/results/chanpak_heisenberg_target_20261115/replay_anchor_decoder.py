from random import Random
def count(word):
    h,u,t=1,0,0
    for b in word:
        h+=(1-b)*(u+t)
        u=2*u+b
        t=b
    return h,u
def word(n,k):
    return tuple((n>>i)&1 for i in range(k-1,-1,-1))
def target(seq,x,C):
    a,u=count(seq)
    return a+1+C*(u+1)-x
def need(ok):
    if not ok:raise AssertionError('decoder replay')
def literal(rows):
    n=len(rows)
    table=[[0]*n for _ in range(1<<n)]
    for i in range(n):table[1<<i][i]=1
    for s in range(1,1<<n):
        for j,m in enumerate(table[s]):
            if not m:continue
            v=rows[j]&~s
            while v:
                b=v&-v
                table[s|b][b.bit_length()-1]+=m
                v^=b
    return sum(table[-1])
def glue(w,m,s,wr=None):
    n=len(w)
    W=range(n+1,n+1+s)
    U=range(n+1+s,n+1+s+m)
    r=[0]*(n+1+m+s)
    for i,b in enumerate(w,1):
        if b:r[0]|=1<<i
        else:r[i]|=1
        for j in range(i+1,n+1):r[i]|=1<<j
        for j in U:r[i]|=1<<j
    for i in range(n+1):
        for j in W:r[i]|=1<<j
    for i in W:
        for j in list(W)+list(U):
            if (j in U or (i<j if wr is None else wr[i-(n+1)]>>(j-(n+1))&1)):r[i]|=1<<j
    for i in U:
        r[i]|=1
        for j in U:
            if i<j:r[i]|=1<<j
    return r
def check():
    rng=Random(4221)
    for j in range(3000):
        a=rng.randrange(1,10000)
        k=rng.randrange(1,15)
        b=rng.randrange(1<<k)
        p,v=word(a,a.bit_length()),word(b,k)
        x=(1<<rng.randrange(1,8))-1
        C=x<<rng.randrange(1,8)
        r=count(v[::-1])[1]
        D=(C+1)*(1<<k)-1-r
        e=p[-1]*(1-v[0])
        A=count(p)[0]
        N=target(p+v,x,C)
        need(N==target(v,x,C)+a*D+A-1+e)
        if a<1<<k:need((N-target(v,x,C))//D==a)
        g=k+(C+1).bit_length()
        N=target((1,)+(0,)*g+v,x,C)
        need(N//(1<<g)==(C+1)*(1<<k)-r)
        need(N%(1<<g)==C*b+count(v)[0]+1+C-x)
    print('3000 exact quotient, seam, and checksum controls')
    for j in range(120):
        n=rng.randrange(7)
        w=word(rng.randrange(1<<n),n)
        m,s=rng.randrange(1,4),rng.randrange(1,4)
        x=(1<<m)-1
        need(literal(glue(w,m,s))==target(w,x,x<<s))
    print('120 literal Hamilton path computations')
    for k,a,b1,b5,c1,c5 in ((3,1,6,7,229,261),(5,5,28,27,951,921)):
        L=32*(1<<k)+5*a
        for b,t,s,wr,offset in ((b1,1,5,None,c1),(b5,5,4,(8,1,3,6),c5)):
            g=1
            seq=(1,)+(0,)*g+word(b,k)
            H=literal(glue(seq,1,s,wr))
            need(H==L*(1<<g)+offset)
    for k in range(5,16):
        g=k+6
        for a,b1,b5,c1,c5 in (((1<<(k-3)),6,7,229,261),(5*(1<<(k-5)),28,27,951,921)):
            L=32*(1<<k)+5*a
            for b,t,offset in ((b1,1,c1),(b5,5,c5)):
                A,u=count((1,)+(0,)*g+word(b,k))
                need((t*A+32*(u+1))==(L<<g)+offset)
    print('Four literal profile-switch controls and 44 symbolic high-window controls')
if __name__=='__main__':check()
