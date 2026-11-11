"""unitops.py - which one-vertex extensions of the base preserve a unit set?
For every labelled tournament Q of order m <= MMAX and every unit set O (W(O) = H(Q) - 1, |O| in {2,3}),
try every new vertex y with out-set Y (subset of V(Q)), keeping O (y not in O), and test whether (Q + y, O) is unit.
Tally the successes by a description of Y relative to O:  (|Y meet O|, Y contains V-O ?, Y meet (V-O) empty ?, ...)."""
import sys, itertools, collections
def tables(A):
    n=len(A); full=1<<n
    out=[sum(1<<u for u in range(n) if A[v][u]) for v in range(n)]
    inn=[sum(1<<u for u in range(n) if A[u][v]) for v in range(n)]
    en=[[0]*n for _ in range(full)]; st=[[0]*n for _ in range(full)]
    for v in range(n): en[1<<v][v]=1; st[1<<v][v]=1
    for m in range(1,full):
        for v in range(n):
            c=en[m][v]
            if c:
                nx=out[v]&~m
                while nx:
                    u=(nx&-nx).bit_length()-1; nx&=nx-1; en[m|1<<u][u]+=c
            c=st[m][v]
            if c:
                nx=inn[v]&~m
                while nx:
                    u=(nx&-nx).bit_length()-1; nx&=nx-1; st[m|1<<u][u]+=c
    return en,st
def WUH(A):
    n=len(A); full=(1<<n)-1; en,st=tables(A); hs=[sum(e) for e in en]; hs[0]=1
    W=[[0]*n for _ in range(n)]
    for X in range(1,full):
        Y=full^X
        for a in range(n):
            if X>>a&1 and en[X][a]:
                for c in range(n):
                    if Y>>c&1 and st[Y][c]:
                        q=en[X][a]*st[Y][c]; W[a][c]+=q; W[c][a]+=q
    U=[sum(st[X][y]*hs[full^X] for X in range(1,full+1) if X>>y&1) for y in range(n)]
    return W,U,hs[full]
def strong(A):
    n=len(A); seen={0}; stack=[0]
    while stack:
        v=stack.pop()
        for u in range(n):
            if A[v][u] and u not in seen: seen.add(u); stack.append(u)
    if len(seen)<n: return False
    seen={0}; stack=[0]
    while stack:
        v=stack.pop()
        for u in range(n):
            if A[u][v] and u not in seen: seen.add(u); stack.append(u)
    return len(seen)==n
MMAX=int(sys.argv[1]) if len(sys.argv)>1 else 5
tally=collections.Counter(); tot=collections.Counter(); examples={}
for m in range(3,MMAX+1):
    pairs=[(i,j) for i in range(m) for j in range(i+1,m)]
    for bits in range(1<<len(pairs)):
        A=[[0]*m for _ in range(m)]
        for k,(i,j) in enumerate(pairs):
            if bits>>k&1: A[i][j]=1
            else: A[j][i]=1
        if not strong(A): continue
        W,U,H=WUH(A)
        for k in (2,3):
            for O in itertools.combinations(range(m),k):
                if sum(W[a][b] for i,a in enumerate(O) for b in O[i+1:])!=H-1: continue
                M=sum(U[y] for y in O)
                for Y in range(1<<m):
                    B=[row[:]+[0] for row in A]+[[0]*(m+1)]
                    for u in range(m):
                        if Y>>u&1: B[m][u]=1
                        else: B[u][m]=1
                    Yo=sum(1 for o in O if Y>>o&1); rest=[u for u in range(m) if u not in O]
                    Yr=sum(1 for u in rest if Y>>u&1)
                    desc=(k,Yo,"none" if Yr==0 else ("all" if Yr==len(rest) else "some"))
                    tot[desc]+=1
                    W2,U2,H2=WUH(B)
                    if sum(W2[a][b] for i,a in enumerate(O) for b in O[i+1:])==H2-1:
                        M2=sum(U2[y] for y in O); tally[desc]+=1
                        examples.setdefault(desc,[]).append((m,H,M,H2,M2))
print("pattern (|O|, |Y meet O|, Y on the rest): preserved / tried")
for d in sorted(tot): print(d,tally[d],"/",tot[d], examples.get(d,[])[:6])
