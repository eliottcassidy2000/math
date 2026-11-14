from collections import Counter
from itertools import combinations,permutations,product
from math import comb,factorial,gcd,isqrt

def req(ok,msg):
    if not ok:raise RuntimeError(msg)
def E(f):
    return {tuple(sorted((i,f[i]))) for i in range(1,len(f))}
def chart(f):
    n=len(f)
    residual={e for e in combinations(range(n),2) if e not in E(f)}
    pairs={e:e if e[0] else tuple(sorted((e[1],f[e[1]]))) for e in residual}
    req(set(pairs.values())==set(combinations(range(1,n),2)),"diagram not bijection")
    # Transfer a directed ordered pair, not just a bit in numeric coordinates.
    def lift(a,b):
        if a and b:return (a,b)
        return (b,f[b]) if a==0 else (f[a],a)
    return lift

def valid_functions(n):
    for x in product(*(tuple(j for j in range(n) if j!=i) for i in range(1,n))):
        f=(-1,)+x
        if all(not y or f[y]!=i for i,y in enumerate(f) if i):yield f

def partitions(m,start=1):
    if m==0:
        yield ()
    else:
        for i in range(start,m+1,2):
            for p in partitions(m-i,i):yield (i,)+p

def fix_f(lens):
    M=Counter(lens);Z=1
    for l,t in M.items():
        X=sum(d*c for d,c in M.items() if l%d==0)
        Z*=sum((-1)**k*factorial(t)//(factorial(t-2*k)*factorial(k)*2**k) * l**k*X**(t-2*k)
            for k in range(t//2+1))
    return Z

def gperm(lens):
    P=[0];a=1
    for l in lens:
        P.extend(list(range(a+1,a+l))+[a]);a+=l
    return tuple(P)

def marked():
    for n in range(2,8):
        fs=list(valid_functions(n))
        B=0
        for P in partitions(n-1):
            g=gperm(P);c=fix_f(P)
            raw=[f for f in fs if all(f[g[i]]==g[f[i]] for i in range(1,n))]
            req(len(raw)==c,"declared fixed-function census")
            g2=sum((s-1)//2 for s in P)+sum(gcd(a,b) for a,b in combinations(P,2))
            centralizer=1
            for l,t in Counter(P).items():centralizer*=l**t*factorial(t)
            B+=factorial(n-1)//centralizer*c*(1<<g2)
        req(B%factorial(n-1)==0,"noninteger orbit count")
        print(f"n={n}: {len(fs)} labeled functional carriers (root fixed); {B//factorial(n-1)} marked tournament orbits")
        if n<=5:
            for f in fs:
                L=chart(f)
                for q in permutations(range(1,n)):
                    p=(0,)+q
                    g=[-1]*n
                    for v in range(1,n):g[p[v]]=p[f[v]]
                    K=chart(g)
                    for a,b in combinations(range(n),2):
                        if (a,b) not in E(f):
                            req(K(p[a],p[b])==tuple(map(p.__getitem__,L(a,b))),"not equivariant")
    # g cycles 1,2,3 in K4 on {1,2,3,4}, three isolated vertices incl root 0.
    g=(0,2,3,1,4,5,6)
    K=set(combinations((1,2,3,4),2))
    def orbits(X):
        X=set(X);c=0
        while X:
            x=X.pop(); y=tuple(sorted(map(g.__getitem__,x)));c+=1
            while y!=x:X.remove(y);y=tuple(sorted(map(g.__getitem__,y)))
        return c
    req(orbits(set(combinations(range(7),2))-K)==9 and orbits(combinations(range(1,7),2))==7,
        "first hostile: wrong fixed orientation counts")
    print("rooted K4 + 3 singletons: 512 versus 128 fixed orientations")

def hp(A):
    n=len(A);E=[[0]*n for _ in range(1<<n)]
    for i in range(n):E[1<<i][i]=1
    for s in range(1,1<<n):
        for i,x in enumerate(E[s]):
            if x:
                for j in range(n):
                    if not(s>>j&1) and A[i]>>j&1:E[s|1<<j][j]+=x
    return [1]+[sum(x) for x in E[1:]]
def profile(A):
    H=hp(A)
    return H[-1],sum(x*y for x,y in zip(H,H[::-1]))
def good(A,P):
    return all(A[x]>>y&1 for x,y in zip(P,P[1:]))
def L(A,v):
    out=0
    for P in permutations(x for x in range(len(A)) if x!=v):
        for k in range(len(P)+1):
            if good(A,P[:k]) and good(A,P[k:]) and (k==len(P) or A[v]>>P[k]&1):
                out+=1
    return out
def glue(A,r,j):
    m=len(A);n=m+j+1
    B=[row|sum(1<<k for k in range(m,n) if i!=r or k!=n-1) for i,row in enumerate(A)]
    B.extend(sum(1<<k for k in range(i+1,n))|((1<<r) if i==n-1 else 0) for i in range(m,n))
    return tuple(B)
def labeled(n):
    for bits in product((0,1),repeat=comb(n,2)):
        A=[0]*n
        for (a,b), x in zip(combinations(range(n),2),bits):
            if x:a,b=b,a
            A[a]|=1<<b
        yield A
def strong(A):
    for r in range(len(A)):
        V=1<<r
        for _ in A:
            for j,a in enumerate(A):
                if V>>j&1:V|=a
        if V!=(1<<len(A))-1:return False
    return True
def rotations():
    for n in range(1,5):
        for A in labeled(n):
            H=hp(A)[-1]
            for r in range(n):
                ell=L(A,r)
                req(not any(H+(ell<<j) in (45,69) for j in (1,2,4)),"earlier genuine Fermat source")
                req(69 not in [H+(ell<<j) for j in range(1,7)],"earlier near-transitive source for 69")
    A=(2,8,3,21,7)
    req(profile(A)==(9,84) and [L(A,i) for i in range(5)]==[10,9,16,15,25],
       "wrong source menu")
    for r,w in ((1,(45,1416)),(3,(69,1968))):
        B=glue(A,r,2)
        req(strong(B) and profile(B)==w,"literal seam mismatch")
    print("same five-vertex strong source Q=84; exact rooted strong outputs (45,1416),(69,1968)")

def square_switch_history():
    def cycle_count(n):
        G=[sum(1<<(y-1) for y in range(1,n+1) if x!=y and isqrt(x+y)**2==x+y)
           for x in range(1,n+1)]
        def visit(x,U):
            if not U:return int(bool(G[x]&1))
            Y=G[x]&U;out=0
            while Y:
                s=Y&-Y;Y-=s;out+=visit(s.bit_length()-1,U^s)
            return out
        return visit(0,((1<<n)-1)^1)//2
    def squares(n):
        N={x:{y for y in range(1,n+1) if x!=y and isqrt(x+y)**2==x+y}
           for x in range(1,n+1)}
        return sum(comb(len(N[a]&N[b]),2) for a,b in combinations(N,2))//2
    req(cycle_count(34)==11,"Q34 cycle count")
    req([squares(n) for n in (34,45,46,68,69)]==[0,0,1,2,2],"historical exchange threshold")
    print("square-sum Q34 has 11 Hamilton cycles but zero four-cycle moves; Q45/Q46/Q68/Q69 four-cycle counts 0,1,2,2")

def ternary_gate():
    expected={41,51,55,59,65,67,69}
    witness={}
    b2=[]
    for m in (5,6):
        H={}
        for u,v in product(range(1<<m),repeat=2):
            full=(1<<m)-1
            def anchor(x,j=1):
                L=(1<<x)-1;R=full^((1<<(x+j))-1)
                p,q=u&L,v&L;s,t=((~u)&R),((~v)&R)
                return p and q and ((p|q).bit_count()>=2) and s and t and ((s|t).bit_count()>=2)
            if not any(anchor(x) for x in range(m)):continue
            for d in (0,1):
                A=tuple((((1<<m)-(2<<i)) | ((not (u>>i&1))<<m) | ((not (v>>i&1))<<(m+1)))
                        for i in range(m))+(u|(d<<(m+1)),v|((1-d)<<m))
                h=hp(A)[-1]
                H[h]=H.get(h,0)+1
                if m==6 and any(anchor(x,2) for x in range(m-1)):b2.append(h)
                if h in expected:witness.setdefault(h,A)
        req({h for h in H if h<=69}==expected,"ternary gate spectrum mismatch")
    req(min(b2)==105,"wrong first exact two-ternary-digit threshold")
    print("first H admitting two ternary-eligible vertices is 105")
    print("absolute H≤69 ternary-eligible counts near 45/69:",sorted(witness))
    print("H67/H69 literal witnesses:",witness[67],witness[69])
    for h,expected_hist in ((67,{1:19,2:9,3:2,4:4,8:1}),(69,{1:19,2:12,3:2,4:5})):
        T=Counter()
        A=witness[h]
        for u,v in ((5,6),(6,5)):
            M=Counter()
            for bits in product(range(3),repeat=5):
                X,Y,Z=[[j for j,k in enumerate(bits) if k==l] for l in range(3)]
                if good(A,tuple(X+[u]+Y+[v]+Z)):
                    M[(X[-1] if X else -1,Y[0] if Y else -1,
                       Y[-1] if Y else -1,Z[0] if Z else -1)] += 1
            T.update(M.values())
        req(T==expected_hist,"wrong endpoint-fiber sum decomposition")
    print("H67/H69 endpoint-fiber sizes directly recovered")
    A=(62,28,56,112,32,2,55)
    F=Counter()
    for path in permutations(i for i in range(7) if i!=2):
        for use in product((0,1),repeat=7):
            q=tuple(x for i,b in enumerate(use) for x in ((2,) if b else ())
                    + ((path[i],) if i<len(path) else ()))
            if good(A,q):F[sum(use)]+=1
    F.pop(0,None)
    req(F=={1:41,2:32,3:8},"unexpected ternary-ray visit profile")
    print("strong all-length family H(b)=1+8·2^b+8·3^b certified by finite visit profile")

if __name__=="__main__":
    marked()
    rotations()
    square_switch_history()
    ternary_gate()
