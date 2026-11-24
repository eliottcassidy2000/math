from collections import Counter
from itertools import combinations, permutations, product
from math import factorial
from gadget import profiles

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

def rank(p):
    r=0;a=list(range(len(p)))
    for i in p:
        r=len(a)*r+a.index(i);a.remove(i)
    return r

def unrank(n,r):
    a=list(range(n));res=[]
    for k in range(n-1,-1,-1):
        j,r=divmod(r,factorial(k));res.append(a.pop(j))
    assert r==0
    return res

def encode(G):
    d=len(G[0]); N=sorted(G[0]);branches=[sorted(G[v]-{0}) for v in N]
    def partner(x, j):
        a, = G[x].intersection(branches[j])
        return a
    B=[branches[0]]+[[partner(x,j) for x in branches[0]] for j in range(1,d)]
    C=0;f=factorial(d-1);a=0;pp=[]
    for i,j in combinations(range(1,d),2):
        q=[B[j].index(partner(x,j)) for x in B[i]]
        C=C*f+rank(q);pp.append(q);a+=1
    r=C
    flat=[]
    for _ in range(a):
        r,x=divmod(r,f);flat.append(unrank(d-1,x))
    assert flat[::-1]==pp and r==0
    return C,(f**a-1).bit_length()

if __name__=='__main__':
 for n in range(1,9):
    for b in product((0,1),repeat=n-1):
        w=b+(0,)
        assert word_h(w)==hp(direct(w))
        assert mate_h(w)==hp(direct(w,flip=True))
 print("word and complement transfer checked on 255 small tilings")
 for mask,e in ((177,111),(181,113)):
    Y=tournament(6,mask);contact=[(19>>i)&1 for i in range(6)]
    assert list(profiles(6,mask))[19]==(29,320,162,158,e)
    a=[p for p in permutations(range(6)) if all(Y[i][j] for i,j in zip(p,p[1:]))]
    assert len(a)==29
    assert sum(all(Y[i][j]==Y[p[i]][p[j]] for i,j in combinations(range(6),2)) for p in permutations(range(6)))==1
    for i in range(6):
        S={i}
        for _ in range(5):S|={j for k in S for j in range(6) if Y[k][j]}
        assert len(S)==6
    for n in range(1,4):
      for b in product((0,1),repeat=n-1):
        w=b+(0,);B=value(w);K=2
        assert hp(direct(w,Y,contact,K))==2**K*(320*B+158)+29*(word_h(w)-1)+162*B+e
 print("guarded pair checked by subset DP; gadget common path:",next(p for p in permutations(range(6)) if all(tournament(6,m)[i][j] for m in (177,181) for i,j in zip(p,p[1:]))))
 for d in (2,3,7):
    C,W=encode(build_graph(d))
    w=tuple(map(int,'10'*5+(format(C,f'0{W}b') if W else '')+'0'))
    K=len(w)+7; B=value(w)
    simple=w+(1,)*K+(0,)
    h=word_h(simple);m=mate_h(simple)
    assert h==2**K*(B+1)+word_h(w)
    F=2**K*(320*B+158)+29*(word_h(w)-1)+162*B+111
    print("degree",d,"payload bits",W,"rank",C)
    print(" black tiling order",len(simple)+1,"H",h,"mate",m)
    print(" single arc fork order",len(simple)+7,"H lower",F,"H upper",F+2)
 G=build_graph(7)
 c=0;N=sorted(G[c])
 rho=lambda v: (5*((-(v//5))%5)+(-(v%5))%5 if v<25 else 5*(v//5)+(-(v%5))%5)
 assert all({rho(j) for j in G[i]}==G[rho(i)] for i in range(50))
 u,v=[i for i in N if rho(i)!=i]
 U=sorted(G[u]-{c})
 for a in U: assert rho(a) in G[a] and rho(a) in G[v]
 f=[]
 for i in N:
    if rho(i)!=i:continue
    row=[next(iter(G[x]&(G[i]-{c}))) for x in U]
    f.append(tuple(row.index(rho(z)) for z in row))
 assert len(f)==5
 for i in range(6):assert {g[i] for g in f}==set(range(6))-{i}
 print("Hoffman-Singleton star model: five matching rows on six sites",f)
 d=57;W=(factorial(56)**1540-1).bit_length();L=W+11;K=L+7
 print("degree 57 conditional",W,"payload bits; order",L+K+2,"marked line;",L+K+8,"fork")
