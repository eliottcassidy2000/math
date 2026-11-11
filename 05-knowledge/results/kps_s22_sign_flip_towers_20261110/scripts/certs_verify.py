"""certs_verify.py - literal recount of certificates "sign c v2 order m code Omask levels":
builds the pair (B+, O+) explicitly (start pair, then the tower levels) and recomputes its state (r, M, h, E) by the
subset dynamic programme pairstate.c; checks r = sign, M = c 2^v2, order."""
import subprocess, sys
def adj(m, code):
    A=[[0]*m for _ in range(m)]; b=0
    for i in range(m):
        for j in range(i+1,m):
            if code>>b&1: A[i][j]=1
            else: A[j][i]=1
            b+=1
    return A
def level_up(A,O,a):
    n=len(A); N=n+a+1; T=[[0]*N for _ in range(N)]
    for i in range(n):
        for j in range(n): T[i][j]=A[i][j]
        for j in range(a): T[i][n+j]=1
    for i in range(a):
        for j in range(i+1,a): T[n+i][n+j]=1
        T[n+i][N-1]=1
    for i in range(n):
        if i in O: T[N-1][i]=1
        else: T[i][N-1]=1
    return T
lines=[]; meta=[]
for l in open(sys.argv[1]):
    if l[0]=="#": continue
    p=l.split(); sign,c,v,order,m,code,Om=map(int,p[:7]); levels=[] if p[7]=="-" else [tuple(map(int,z.split(","))) for z in p[7].split(";")]
    A=adj(m,code); O=[i for i in range(m) if Om>>i&1]
    for (d,i,mm) in levels:
        a=i+d+mm; n=len(A); A=level_up(A,O,a); O=[n+i-1,n+i-1+d]
    n=len(A)
    if n>24: print("SKIP (order %d too large for the subset programme): c=%d"%(n,c)); continue
    lines.append(" ".join([str(n)]+[str(sum(1<<j for j in range(n) if A[i][j])) for i in range(n)]+[str(sum(1<<o for o in O))]))
    meta.append((sign,c,v,order))
out=subprocess.run(["./pairstate"],input="\n".join(lines)+"\n",capture_output=True,text=True).stdout.split("\n")
bad=0
for (sign,c,v,order),o in zip(meta,out):
    n,r,M,h,E=map(int,o.split())
    if not (n==order and r==sign and M==c*2**v): bad+=1; print("MISMATCH",sign,c,v,order,"->",n,r,M,h,E)
print("%s: %d certificates recounted, %d mismatches"%(sys.argv[1],len(meta),bad))
