# Check the R_1-gluing transfer at the three new roots against the brute-force root statistics.
import random, subprocess
def tobits(adj,n): return ''.join('1' if adj[i][j] else '0' for i in range(n) for j in range(i+1,n))
def stats(adj,n):
    out=subprocess.run(['./hl'],input=tobits(adj,n)+'\n',capture_output=True,text=True).stdout.split()
    H=int(out[1]); d={}
    for it in out[2:]:
        v,r=it.split(':'); h,L,E,S=map(int,r.split(',')); d[int(v)]=(h,L,E,S)
    return H,d
random.seed(7); bad=0; tests=0
for trial in range(300):
    m=random.randint(2,7)
    adj=[[False]*m for _ in range(m)]
    for i in range(m):
        for j in range(i+1,m):
            if random.random()<0.5: adj[i][j]=True
            else: adj[j][i]=True
    v=random.randrange(m)
    H,d=stats(adj,m); h,L,E,S=d[v]; D=L-H
    n=m+2; r1,r2=m,m+1
    B=[[False]*n for _ in range(n)]
    for i in range(m):
        for j in range(m): B[i][j]=adj[i][j]
    for a in range(m):
        if a!=v: B[a][r1]=True; B[a][r2]=True
    B[v][r1]=True; B[r1][r2]=True; B[r2][v]=True
    H2,d2=stats(B,n)
    def DS(w): hh,LL,EE,SS=d2[w]; return (LL-H2,hh,EE,SS)
    exp={v:(H-h+2*D,h,h,S), r1:(H-E,2*H+D-E,H+D,S), r2:(0,H,H,S)}
    tests+=1
    if H2!=3*H+2*D: bad+=1; print("H mismatch",H,D,H2)
    for w in (v,r1,r2):
        if DS(w)!=exp[w]: bad+=1; print("root",w,"got",DS(w),"exp",exp[w])
print("tests",tests,"mismatches",bad)
