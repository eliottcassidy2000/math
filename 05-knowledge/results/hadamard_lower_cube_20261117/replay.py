from itertools import permutations, product
from functools import cache
def need(x):
    if not x: raise AssertionError

def parity(x): return x.bit_count()&1
def sector(x):return parity(x&6)

def default_frame(n):
    out=[]
    for side in (1,0):
        a=[x for x in range(n) if sector(x)==side and parity(x)]
        b=[x for x in range(n) if sector(x)==side and not parity(x)]
        out.extend(x for y in zip(a,b) for x in y)
    return tuple(out)

def blocks(p):
    n=len(p);half=n//2
    reps=sorted({min(x^y for y in (0,1,6,7)) for x in p[half:]})
    forbidden=(min(p[half]^y for y in (0,1,6,7)),p[half-1]&1)
    return [(r,c) for r in reps for c in (0,1) if (r,c)!=forbidden]

def overlay(p,choices):
    n=len(p)
    M=[[(-1)**(parity(x&y)^parity(y)) for y in p] for x in p]
    for r,c in choices:
        for i,x in enumerate(p):
            if x^r in (0,1,6,7):
                for j,y in enumerate(p):
                    if sector(y)==1 and y&1==c:
                        M[i][j]*=-1
    d=[1]
    for i in range(n-1): d.append(d[-1]*M[i][i+1])
    return [[M[i][j]*d[i]*d[j] for j in range(n)] for i in range(n)]

def halves(M):
    n=len(M);A=[0]*n;B=[0]*n
    for i in range(n):
        for j in range(i+1,n):
            x,y=(i,j) if M[i][j]>0 else (j,i)
            A[x]|=1<<y
            x,y=(i,j) if M[j][i]<0 else (j,i)
            B[x]|=1<<y
    return tuple(A),tuple(B)

@cache
def profiles(rows):
    n=len(rows); full=(1<<n)-1
    @cache
    def d(S,v):
        R=S^(1<<v)
        return 1 if not R else sum(d(R,u) for u,r in enumerate(rows) if R>>u&1 and r>>v&1)
    @cache
    def start(S,v):
        R=S^(1<<v)
        return 1 if not R else sum(start(R,u) for u in range(n) if rows[v]>>u&1 and R>>u&1)
    @cache
    def h(S):
        return sum(d(S,v) for v in range(n) if S>>v&1) if S else 1
    return [(h(full),sum(d(S,v)*h(full^S) for S in range(1,full+1) if S>>v&1),start(full,v))
            if d(full,v) else None for v in range(n)]

def audit(M,M0,p,n):
    need(all(sum(a*b for a,b in zip(M[i],M[j]))==n*int(i==j) for i in range(n) for j in range(n)))
    need(all(M[i][i]==1 for i in range(n)))
    need(all(M[i][i+1]==1 and M[i+1][i]==-1 for i in range(n-1)))
    need(all(M[i][j]==M0[i][j] for i in range(n) for j in range(i-1,n) if j>=0))

def H(A):
    n=len(A); D=[[0]*n for _ in range(1<<n)]
    for i in range(n): D[1<<i][i]=1
    for s in range(1,1<<n):
        for i in range(n):
            x=D[s][i]
            if x:
                R=A[i]&~s
                while R:
                    bit=R&-R; R^=bit
                    D[s|bit][bit.bit_length()-1]+=x
    return sum(D[-1])

if __name__=='__main__':
    for n in (8,16,32):
        p=default_frame(n);B=blocks(p);M0=overlay(p,[])
        need(len(B)==n//4-1)
        for i in range(1<<len(B)):
            M=overlay(p,[b for j,b in enumerate(B) if i>>j&1])
            audit(M,M0,p,n)
            need(sum(M[r][c]!=M0[r][c] for r in range(n) for c in range(n))==n*i.bit_count())
        print(n,'bits',len(B),flush=True)
    p=(2,3,4,5,1,0,7,6)
    A,X=halves(overlay(p,[]));AA,Y=halves(overlay(p,[(0,0)]))
    need(A==AA)
    print('literal',A,X,Y)
    print('profiles',profiles(X),profiles(Y))
    for (_,b,c) in profiles(X):
        for (_,e,f) in profiles(Y):
            corners=(e-b,e-b+16,e-b+16-f+c)
            need(min(corners)>=6 or max(corners)<=-6)
    p=default_frame(16);B=blocks(p)
    menu=[H(halves(overlay(p,[b for j,b in enumerate(B) if i>>j&1]))[1]) for i in range(8)]
    need(menu==[612097097,612081469,612051653,612047197,606270393,606264149,606257381,606246749])
    print('order-16 lower counts in mask order',menu)
