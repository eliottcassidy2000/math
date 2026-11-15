# Check: the double cone D(B) = (a -> B -> c -> a, a and c single vertices, B => c, c => a, a => B)
# has P_{D(B)}(q) = q [ P_B(q+1) + q P_B(q) - P_B(q-1) ] for q = 1, 2, 3, with P_X(q) = ordered q-slot path covers.
import random, itertools
def Pcounts(out, n, qmax):
    full=(1<<n)-1; dp=[[0]*n for _ in range(1<<n)]; H=[0]*(1<<n)
    for v in range(n): dp[1<<v][v]=1
    for S in range(1,1<<n):
        for v in range(n):
            x=dp[S][v]
            if not x: continue
            H[S]+=x; nb=out[v]&~S
            while nb:
                u=(nb&-nb).bit_length()-1; nb&=nb-1; dp[S|1<<u][u]+=x
    H[0]=1
    P={1:H[full]}; F=H[:]                      # F_q(R) = # ordered q-slot covers of R (empty slots allowed)
    for q in range(2,qmax+1):
        G=[0]*(1<<n)
        for R in range(1<<n):
            s=0; S=R
            while True:
                s+=H[S]*F[R^S]
                if S==0: break
                S=(S-1)&R
            G[R]=s
        F=G; P[q]=F[full]
    return P
rng=random.Random(4); bad=0
for t in range(60):
    m=rng.randint(1,6)
    outB=[0]*m
    for i in range(m):
        for j in range(i+1,m):
            if rng.random()<0.5: outB[i]|=1<<j
            else: outB[j]|=1<<i
    PB=Pcounts(outB,m,4); PB[0]=0
    n=m+2; A,C=m,m+1                            # vertices 0..m-1 = B, m = a, m+1 = c
    out=[o for o in outB]+[0,0]
    for v in range(m): out[A]|=1<<v; out[v]|=1<<C
    out[C]|=1<<A
    PD=Pcounts(out,n,3)
    for q in (1,2,3):
        pred=q*(PB[q+1]+q*PB[q]-PB[q-1])
        if PD[q]!=pred: bad+=1
print('double-cone operator identity on 60 random cores of order 1-6, q = 1,2,3: %d mismatches'%bad)
