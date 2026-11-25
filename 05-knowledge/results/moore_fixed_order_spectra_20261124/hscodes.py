from collections import Counter
from itertools import combinations, permutations, product
from math import factorial

def word_h(w):
    if not w: return 1
    ans=w[0]+1-w[-1]
    b=0
    for i,x in enumerate(w):
        if not x: ans+=b
        if i and not w[i-1] and x: ans+=1
        b=2*b+x
    return ans

def value(w):
    b=0
    for x in w:b=2*b+x
    return b

def mate_h(w):
    assert w and w[-1]==0
    empty=(0,0,0,1) # clipped run count, top run singleton, pending bad descent, endpoint permitted
    D={(empty,empty,2):1}
    for i,a in enumerate(w):
        out=1-a if i+1<len(w) else 0 # flipped apex outcontact a -> i
        N=Counter()
        for (l,r,last), num in D.items():
            for j,old in enumerate((l,r)):
                c,s,b,ok=old
                # For left lane keep the low run's last vertex, for right the high run's first.
                good=(1-out) if j==0 else out
                z=[]
                if not c:
                    z.append((1,1,0,good))
                else:
                    if j==last:
                        z.append((c,0,0,good if j==0 and c==1 else ok))
                    if not b:
                        z.append((2,1,int(j==last and s),good if j==1 else ok))
                for new in z:
                    N[(new,r,j) if j==0 else (l,new,j)]+=num
        D=N
    return sum(n for (l,r,j),n in D.items() if not l[2] and not r[2] and l[3] and r[3])

def tournament(n, rev=0):
    A=[[0]*n for i in range(n)]
    for k,(i,j) in enumerate(combinations(range(n),2)):
        A[j][i]=(rev>>k)&1
        A[i][j]=1-A[j][i]
    return A

def hp(A):
    n=len(A); M=1<<n
    dp=[[0]*n for _ in range(M)]
    for s in range(1,M):
        for i in range(n):
            if s>>i&1:
                t=s^(1<<i)
                dp[s][i]=sum(dp[t][j]*A[j][i] for j in range(n)) if t else 1
    return sum(dp[-1])

def direct(w,Y=(),contacts=(),K=0,flip=False):
    ny=len(Y); L=len(w); n=L+ny+K+(1 if K else 0)+1
    A=tournament(n)
    if K: word=tuple(w)+tuple(contacts)+(1,)*K+(0,)
    else: word=tuple(w)
    for i in range(n-1):A[-1][i]=word[i];A[i][-1]=1-word[i]
    for i in range(ny):
        for j in range(ny):
            A[L+i][L+j]=Y[i][j]
    if flip:
        for i,j in combinations(range(n),2):
            if j>i+1:A[i][j],A[j][i]=A[j][i],A[i][j]
    return A

def build_graph(d):
    if d in (2,3):
        G=[set() for _ in range(1+d*d)]
        def add(a,b):G[a].add(b);G[b].add(a)
        for i in range(5):
            add(i,(i+1)%5)
            if d==3:add(i,i+5);add(i+5,5+(i+2)%5)
    else:
        assert d==7
        G=[set() for _ in range(50)]
        def add(a,b):G[a].add(b);G[b].add(a)
        for h in range(5):
            for x in range(5):
                add(5*h+x,5*h+(x+1)%5)
                add(25+5*h+x,25+5*h+(x+2)%5)
                for i in range(5):add(5*h+x,25+5*i+(x+h*i)%5)
    return G


from math import factorial
def unrank(items,k):
    items=list(items);out=[]
    while items:
        a,k=divmod(k,factorial(len(items)-1))
        out.append(items.pop(a))
    assert not k
    return out

def phase_chart(a):
    sides=[set() for i in range(50)]
    def edge(i,j):sides[i].add(j);sides[j].add(i)
    for i,x in product(range(5),repeat=2):
        edge(5*i+x,5*i+(x+1)%5)
        edge(25+5*i+x,25+5*i+(x+2)%5)
        for j in range(5):edge(5*i+x,25+5*j+(x+a[i][j])%5)
    return sides

def check(a):
    G=phase_chart(a)
    assert all(len(N)==7 for N in G)
    assert all(len(G[i]&G[j])==int(j not in G[i]) for i,j in combinations(range(50),2))

def phase14(k):
    u,v=divmod(k,24)
    x=[0,1]+unrank([2,3,4],u)
    y=[0]+unrank(range(1,5),v)
    return tuple(tuple(a*b%5 for b in y) for a in x)
def bits(s):return tuple(map(int,s))
def extremes(prefix,width,code):return bits(prefix+format(code,f'0{width}b')+'0')
G={phase14(k) for k in range(144)}; assert len(G)==144
for a in G:check(a)
for n,w in [(14,extremes('1010',8,0)),(16,extremes('1010',10,0)),
            (50,extremes('1010101010',38,47503389336))]:
    h=word_h(w);mate=mate_h(w)
    assert len(w)+1==n
    if n<20:
        assert hp(direct(w))==h
        assert hp(direct(w,flip=True))==mate
    print(n,''.join(map(str,w)),h,mate)
print('14-site 144 decodable classes, number of distinct source/mate heights:',
      len({word_h(extremes('1010',8,k)) for k in range(144)}),
      len({mate_h(extremes('1010',8,k)) for k in range(144)}))
from collections import Counter
pairs=[(word_h(extremes('1010',8,k)),mate_h(extremes('1010',8,k))) for k in range(144)]
print('source min/max and multiplicity histogram',min(p[0] for p in pairs),max(p[0] for p in pairs),dict(Counter(Counter(p[0] for p in pairs).values())))
print('mate min/max',min(p[1] for p in pairs),max(p[1] for p in pairs))
def join(A,a,B,b):
    n=len(A); l=len(B)
    T=tournament(n+l)
    for i in range(n):
        for j in range(n):T[i][j]=A[i][j]
    for i in range(l):
        for j in range(l):T[n+i][n+j]=B[j][i] # second factor is taken in converse
    T[a][n+b]=0;T[n+b][a]=1
    return T
def aw(u):return direct(tuple(map(int,bin(u)[2:])) if u else ())
A=aw(100);B=aw(126)
assert hp(join(A,len(A)-1,B,len(B)-1))==18481
A=aw(149);B=tournament(7,891796)
assert hp(join(A,len(A)-1,B,1))==20281
print('literal F1 and F2 controls checked: 18481, 20281')
