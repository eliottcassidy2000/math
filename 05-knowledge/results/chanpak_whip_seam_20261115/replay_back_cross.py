from random import Random

def table(rows):
    n=len(rows);f=(1<<n)-1
    a=[[0]*n for _ in range(f+1)]; h=[1]+[0]*f
    for s in range(1,f+1):
        for i in range(n):
            if s>>i&1:
                t=s^(1<<i)
                a[s][i]=sum(a[t][j] for j in range(n) if (t&rows[i])>>j&1) if t else 1
        h[s]=sum(a[s])
    return a,h

def split(rows,mask):
    tab,h=table(rows);f=len(h)-1
    tots=[0,0]
    for s,ends in enumerate(tab):
        if s:
            tots[s.bit_count()%2]+=h[f^s]*sum(v for i,v in enumerate(ends) if (mask>>i)&1)
    return tots

def arc(row,i,j):
    return ((row[i]>>j)&1)!=0

def converse(rows):
    return tuple(sum((1<<i) for i in range(len(rows)) if arc(rows,i,j)) for j in range(len(rows)))

def ordinal(a,b):
    n=len(a);q=((1<<len(b))-1)<<n
    return tuple(v|q for v in a)+tuple(v<<n for v in b)

def rev(rows,i,j):
    q=list(rows)
    q[i]^=1<<j;q[j]^=1<<i
    return tuple(q)

def glue(rows,mask):
    n=len(rows)
    return tuple(v|(0 if mask>>i&1 else 1<<n) for i,v in enumerate(rows))+(mask,)

def state(rows,mask):
    h=table(rows)[1][-1]
    H=table(glue(rows,mask))[1][-1]
    E,O=split(rows,mask)
    return h,H,E+O,O-E

def tt(n):
    return tuple((1<<n)-(1<<(i+1)) for i in range(n))

def strong(rows):
    f=(1<<len(rows))-1
    for i in range(len(rows)):
        c=1<<i
        for _ in rows:
            for j,v in enumerate(rows):
                if c>>j&1:c|=v
        if c!=f:return False
    return True

def random_t(n,rng):
    r=[0]*n
    for i in range(n):
        for j in range(i+1,n):
            if rng.randrange(2):r[i]|=1<<j
            else:r[j]|=1<<i
    return tuple(r)

def p(t):
    if t==0:return 0
    u,b=divmod(t,2)
    return p(u)+(0 if b else u+(u%2))

def main():
    g=Random(509)
    for _ in range(200):
        a=random_t(g.randrange(1,6),g);b=random_t(g.randrange(1,6),g)
        x=g.randrange(len(a));y=g.randrange(len(b))
        E,O=split(a,1<<x);e,o=split(converse(b),1<<y)
        old=ordinal(a,b);new=rev(old,x,len(a)+y)
        delta=table(new)[1][-1]-table(old)[1][-1]
        if delta!=2*(E*o+O*e):raise AssertionError((a,b,x,y,delta))
        if (delta==0)!=((a[x]==0) and (converse(b)[y]==0)):raise AssertionError(delta)
    print('200 literal rank-two seam-gain and positivity checks')
    A=(38,16,2,7,13,30)
    for j in range(7):
        base=ordinal(A,tt(j+1));new=rev(base,2,6+j)
        if table(new)[1][-1]-23!=(18 if not j else 27*(1<<j)):raise AssertionError(j)
    print('Seven exact geometric depth controls')
    for a,m,pre,d in ((2,13,221,17),(1,112,219,23)):
        base=rev(ordinal(A,(0,)),a,6)
        if state(base,m)!=(41,pre,295,d):raise AssertionError((a,m))
        mid=ordinal(base,(0,));m1=m|128
        if state(mid,m1)!=(41,pre,631,41):raise AssertionError(a)
        ready=ordinal(mid,(0,))
        if not strong(glue(ready,m1)) or state(ready,m1)!=(41,pre+672,1262,0):raise AssertionError(a)
        for t in (1,2,3,5):
            r=list(ready);mask=m1
            for j,z in enumerate(bin(t)[:1:-1]):
                if z=='0':mask|=1<<(len(ready)+j)
            r=ordinal(r,tt(t.bit_length()))
            H=table(glue(r,mask))[1][-1]
            if H!=pre+672+1262*t+41*p(t):raise AssertionError((a,t,H))
        print('Compensated literal arm:',base,m,'guarded H',pre+672)
if __name__=='__main__':main()
