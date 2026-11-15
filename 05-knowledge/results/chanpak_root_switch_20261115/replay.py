"""Verify the four-state outlet gluing recurrence against literal subset DP."""
from random import Random

def counts(adj):
    n=len(adj); F=(1<<n)-1
    p=[[0]*n for _ in range(F+1)]
    for i in range(n):
        p[1<<i][i]=1
    for m in range(1,F+1):
        for i in range(n):
            if p[m][i]:
                todo=adj[i]&~m
                while todo:
                    x=todo&-todo
                    p[m|x][x.bit_length()-1]+=p[m][i]
                    todo-=x
    h=[sum(row) for row in p]; h[0]=1
    return h,p

def transpose(adj):
    return [sum(1<<i for i in range(len(adj)) if (adj[i]>>j)&1) for j in range(len(adj))]

def slot(h,m):
    s=m
    answer=0
    while True:
        answer+=h[s]*h[m^s]
        if s==0:
            return answer
        s=(s-1)&m

def state(adj,r):
    h,_=counts(adj)
    _,p=counts(transpose(adj)); F=len(h)-1; R=1<<r
    s=sum(row[r]*h[F^m] for m,row in enumerate(p))
    d=sum(row[r]*h[F^m]*(-1 if m.bit_count()%2==0 else 1) for m,row in enumerate(p))
    return h[F],s,d,h[F^R]

def sew(W,u,R):
    n=len(W); m=len(R)
    return [row|(((1<<m)-1)<<n) if i!=u else row for i,row in enumerate(W)] + [(row<<n)|(1<<u) for row in R]

def solve(W,u,R,r):
    h,s,d,g=state(W,u); b=h-s-d+2*g; z=s-g
    hr,_=counts(R); _,p=counts(transpose(R)); F=len(hr)-1; c=1<<max(r,0)
    A=hr[F]; Q=slot(hr,F); H=A*b+Q*z
    if r<0:
        return H,A*g+Q*z,A*g,A*g
    s1=d1=s2=d2=0
    for m,row in enumerate(p):
        if not row[r]:
            continue
        e=(-1 if m.bit_count()%2==0 else 1)
        s1+=row[r]*hr[F^m]; d1+=e*row[r]*hr[F^m]
        s2+=row[r]*slot(hr,F^m); d2+=e*row[r]*slot(hr,F^m)
    return (H, (b+g)*s1+2*z*s2+p[F][r]*(d-g),
            (b-g)*d1+z*d2, b*hr[F^c]+z*slot(hr,F^c)+(d-g if len(R)==1 else 0))

def rand_t(n,rng):
    a=[0]*n
    for i in range(n):
        for j in range(i+1,n):
            if rng.randrange(2):
                a[i]|=1<<j
            else:
                a[j]|=1<<i
    return a

def main():
    rng=Random(731)
    tests=0
    for n in range(1,7):
        for m in range(1,7):
            for _ in range(15):
                W=rand_t(n,rng);R=rand_t(m,rng);u=rng.randrange(n)
                F=sew(W,u,R)
                for r in range(-1,m):
                    actual=state(F,u if r<0 else n+r)
                    predicted=solve(W,u,R,r)
                    if actual!=predicted:
                        raise AssertionError((W,R,u,r,actual,predicted))
                    tests+=1
    pair = [( [12,21,24,18,1], 0, (13,15,3,3)),
            ( [14,20,24,18,1], 4, (15,15,3,3))]
    for W,u,P in pair:
        if state(W,u)!=P:
            raise AssertionError((W,u,state(W,u),P))
        for ell in range(1,9):
            R=[sum(1<<j for j in range(i+1,ell)) for i in range(ell)]
            expected = 12*(1<<ell) + P[0]-12
            actual = state(sew(W,u,R),u)[0]
            if expected != actual:
                raise AssertionError((W,u,ell,expected,actual))
    print('unit repair pair: both roots and 16 literal outlet counts checked')
    print(tests,'literal subset-DP comparisons succeeded')
if __name__=='__main__':
    main()
