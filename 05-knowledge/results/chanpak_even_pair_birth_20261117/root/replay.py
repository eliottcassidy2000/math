from itertools import combinations
from functools import cache

def ck(x,*e):
    if not x:raise ArithmeticError(e)

def table(a):
    n=len(a);D=[[0]*n for _ in range(1<<n)]
    for i in range(n):D[1<<i][i]=1
    for m,row in enumerate(D):
        for i,v in enumerate(row):
            if not v:continue
            w=a[i]&~m
            while w:
                b=w&-w;w-=b;D[m|b][b.bit_length()-1]+=v
    return [1]+[sum(r) for r in D[1:]]

def slots(a):
    h=table(a);full=len(h)-1
    @cache
    def p(q,S):
        if S==0:return 1
        if q==0:return 0
        v=0;s=S
        while True:
            v+=h[s]*p(q-1,S^s)
            if not s:break
            s=(s-1)&S
        return v
    return lambda q:p(q,full)

def tower(a):
    n=len(a)
    return tuple(x|1<<(n+1) for x in a)+((1<<n)-1,1<<n)

def pad(a):
    return tuple(x|1<<len(a) for x in a)+(0,)

def cpoly(q,p):
    return q*(p(q+1)-p(q-1)+q*p(q)) if q else 0

def evendeck(a):
    h=table(a)
    return tuple(tuple(sorted(x for s,x in enumerate(h) if s.bit_count()==i)) for i in range(0,len(a)+1,2))

def main():
    A=(6,12,48,5,11,27); B=(2,12,25,49,35,7)
    ck(evendeck(A)==evendeck(B))
    D1,D2=table(A),table(B);F=len(D1)-1
    ck(D1[-1]==D2[-1]==33)
    ck([sorted(d[F^(1<<i)] for i in range(6)) for d in (D1,D2)]==[[3,11,11,11,15,15],[5,9,9,11,13,15]])
    for a in [(0,), (2,0), (2,4,1), (6,12,8,1), A,B]:
        p=slots(a);r=slots(tower(a))
        s,t=slots(tower(pad(a))),slots(pad(tower(a)))
        for q in range(1,5):
            ck(r(q)==cpoly(q,p),q,a)
            ck(s(q)-t(q)==q*(p(q+1)+p(q-1)),q,a)
        # n parity is transparent modulo 3 in two-vertex growth.
        ck(r(1)%3==((1+(-1)**len(a))*p(1))%3,a)
        d=table(a);M=len(d)-1
        for subset in (tuple(range(len(a))),tuple(range(min(3,len(a)))),()):
            z=0
            for k in range(len(subset)+1):
                for U in combinations(subset,k):z+=(-1)**k*d[M^sum(1<<i for i in U)]
            ck(z>=0,a,subset,z)
    C=(2,4,1)
    ck(len(tower(C))==5 and slots(tower(C))(1)==15)
    a=(3910,2044,3576,3057,4001,4033,3728,3841,3648,3076,2056,2)
    d=table(a);M=len(d)-1
    ck((d[M],d[M^(1<<5)],d[M^(1<<8)],d[M^((1<<5)|(1<<8))])==(97309,23395,24041,5803))
    ck(sum(d[M^(1<<i)]>15745 for i in range(12))==7)
    ck(sum(d[M^((1<<i)|(1<<j))]>3357 for i,j in combinations(range(12),2))==27)
    print('PASS exact sorted even decks, odd-deck separation, cyclic transfer, commutator, order-five compression.')
if __name__=='__main__':main()
