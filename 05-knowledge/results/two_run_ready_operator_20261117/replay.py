from itertools import permutations
from collections import Counter
from functools import cache
E=(18,20,17,7,8)

def hp(t):
    n=len(t);N=1<<n
    z=[[0]*n for _ in range(N)]
    for j in range(n):z[1<<j][j]=1
    for s in range(1,N):
        for j,h in enumerate(z[s]):
            if not h:continue
            r=t[j]&~s
            while r:
                b=r&-r;r-=b;z[s|b][b.bit_length()-1]+=h
    return [1]+[sum(s) for s in z[1:]]

def move(t):
    n=len(t)
    return tuple(x | (1<<n) | (1<<(n+1)) | (1<<(n+3)) | (1<<(n+4)) for x in t)+tuple(((x<<n) | (((1<<n)-1) if q==2 else 0)) for q,x in enumerate(E))

def small():
    n=len(E)
    cy=Counter()
    for s in range(1,1<<n):
        I=[i for i in range(n) if s>>i&1]
        if len(I)<=1 or len(I)%2==0:continue
        for a in permutations(I[1:]):
            a=(I[0],)+a
            if all(E[i]>>j&1 for i,j in zip(a,a[1:]+a[:1])):cy[s]+=1
    @cache
    def P(s):
        if not s:return Counter({0:1})
        b=s&-s;z=Counter({i+1:c for i,c in P(s^b).items()})
        for u,c in cy.items():
            if u&b and u&s==u:
                for i,d in P(s^u).items():z[i+1]+=2*c*d
        return z
    sol=[Counter(),Counter()]
    for s in range(1,1<<n):
        if not s>>2&1 or s.bit_count()<2:continue
        z=[i for i in range(n) if s>>i&1 and i!=2]
        k=sum(all(E[i]>>j&1 for i,j in zip(a,a[1:]+(2,))) for a in permutations(z))
        for j,d in P(((1<<n)-1)^s).items():sol[s.bit_count()%2][j]+=d*k
    return P((1<<n)-1),*sol

if __name__=='__main__':
    print(small())
    for B in ((0,),(2,0),(2,4,1),(2,4,9,3)):
        p=hp(B);u=p[-1];R=sum(h*p[len(p)-1-i] for i,h in enumerate(p))
        T=move(B)
        q=hp(T)
        assert q[-1]==3*u+14*R
        assert q[(len(q)-1)^(1<<(len(B)+3))]==u+2*R
        print('literal',len(B),u,R,q[-1])

def rooted(t,r):
    H=hp(t);n=len(t);N=1<<n
    d=[[0]*n for _ in range(N)];d[1<<r][r]=1
    for s in range(1,N):
        for j,h in enumerate(d[s]):
            if not h:continue
            o=t[j]&~s
            while o:
                b=o&-o;o-=b;d[s|b][b.bit_length()-1]+=h
    eo=[0,0]
    for s in range(1,N):eo[s.bit_count()%2]+=sum(d[s])*H[(N-1)^s]
    return H[-1],H[(N-1)^(1<<r)],*eo

def exhaustive():
    from itertools import combinations
    done=0
    for n in range(1,5):
        C=list(combinations(range(n),2))
        for a in range(1<<len(C)):
            B=[0]*n
            for k,(i,j) in enumerate(C):
                if a>>k&1:B[i]|=1<<j
                else:B[j]|=1<<i
            h=hp(B);u=h[-1];R=sum(t*h[(len(h)-1)^s]for s,t in enumerate(h))
            G=u+2*R;E=2*u+8*R
            assert rooted(move(B),n+3)==(3*u+14*R,G,E,E+G)
            done+=1
    print("full profile checks, all labeled inputs through 4:",done)

if __name__=='__main__':
    assert small()==({1:6,3:8,5:1},{1:6,3:2},{0:3,2:3})
    exhaustive()
