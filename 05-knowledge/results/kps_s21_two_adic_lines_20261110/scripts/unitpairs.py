"""unitpairs.py - structure of unit-offset roots: out-degree of the root and the identity
   delta' = sum over unordered pairs {e,f} of out-neighbours of v of W_ef,  W_ef = Q(e,f)+Q(f,e) in B = L - v,
   so unit offset  <=>  sum_{pairs in N+(v)} W_ef = H(B) - 1."""
import sys, collections
def decode(code,m):
    A=[[0]*m for _ in range(m)]; b=0
    for i in range(m):
        for j in range(i+1,m):
            if code>>b&1: A[i][j]=1
            else: A[j][i]=1
            b+=1
    return A
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
fn=sys.argv[1]; maxm=int(sys.argv[2])
deg=collections.Counter(); ok=True; n=0
for line in open(fn):
    m,code,v,H,d,E,U,M=map(int,line.split())
    A=decode(code,m); O=[j for j in range(m) if A[v][j]]
    deg[(m,len(O))]+=1
    if m>maxm: continue
    keep=[i for i in range(m) if i!=v]; B=[[A[i][j] for j in keep] for i in keep]; pos={x:i for i,x in enumerate(keep)}
    en,st=tables(B); nb=m-1; full=(1<<nb)-1
    def Q(a,b):
        t=0
        for X in range(1,full):
            if X>>a&1 and not X>>b&1 and en[X][a] and st[full^X][b]: t+=en[X][a]*st[full^X][b]
        return t
    dp=sum(Q(pos[e],pos[f])+Q(pos[f],pos[e]) for i,e in enumerate(O) for f in O[i+1:])
    n+=1
    if dp!=d-1: ok=False
print("out-degree of the root among unit-offset roots, by order:",{m:{k:deg[(m,k)] for k in range(0,9) if deg[(m,k)]} for m in range(2,10)})
print("identity sum_{pairs in N+(v)} W_ef = H(L-v) - 1 on %d roots of order <= %d:"%(n,maxm),"PASS" if ok else "FAIL")
