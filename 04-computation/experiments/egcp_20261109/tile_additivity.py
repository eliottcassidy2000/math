# brute-force check of THM-5041 (2b)/(2c): tile-wise additivity of H on random digraphs, and H(U,L)=H(L,U), loops irrelevant
import random, itertools
def H(n, A):
    # endpoint DP
    E={}
    for v in range(n): E[(1<<v,v)]=1
    for size in range(2,n+1):
        for S in range(1<<n):
            if bin(S).count('1')!=size: continue
            for v in range(n):
                if not S>>v&1: continue
                Sp=S&~(1<<v); tot=0
                for u in range(n):
                    if Sp>>u&1 and A[u][v]: tot+=E.get((Sp,u),0)
                E[(S,v)]=tot
    full=(1<<n)-1
    return sum(E.get((full,v),0) for v in range(n))
rng=random.Random(7); bad=0; tests=0
for _ in range(300):
    n=rng.randint(2,6)
    A=[[rng.randint(0,1) if i!=j else rng.randint(0,1) for j in range(n)] for i in range(n)]
    i,j=sorted(rng.sample(range(n),2))
    vals={}
    for a,b in itertools.product((0,1),repeat=2):
        B=[row[:] for row in A]; B[i][j]=a; B[j][i]=b; vals[(a,b)]=H(n,B)
    tests+=1
    if vals[(1,0)]+vals[(0,1)]!=vals[(1,1)]+vals[(0,0)]: bad+=1
    At=[[A[j][i] for j in range(n)] for i in range(n)]
    Al=[[A[i][j] if i!=j else 1-A[i][j] for j in range(n)] for i in range(n)]
    if H(n,A)!=H(n,At) or H(n,A)!=H(n,Al): bad+=1
print("tile additivity / transpose / loop checks:",tests,"bad",bad)
