# BFS over bottleneck-free unit words; minimal forbidden (rank-one) factors; exact via two primes
import pickle, numpy as np, math, sys, itertools
from sympy import ilcm
D=pickle.load(open('w30_hankel_13units.pkl','rb'))
U=D['units']; R=D['red']; n=len(U)
Lw=[len(u.split('.')) for u in U]
PR=[67108859, 67108837]
def modmat(M,p):
    A=np.zeros((8,8),dtype=np.int64)
    for i in range(8):
        for j in range(8):
            q=M[i,j]; A[i,j]=(int(q.p)*pow(int(q.q),-1,p))%p
    return A
MM=[[modmat(M,p) for M in R] for p in PR]
def rank_mod(A,p):
    A=A.copy()%p; r=0; rows,cols=A.shape
    for c in range(cols):
        piv=None
        for i in range(r,rows):
            if A[i,c]%p: piv=i;break
        if piv is None: continue
        A[[r,piv]]=A[[piv,r]]
        inv=pow(int(A[r,c]),-1,p)
        A[r]=(A[r]*inv)%p
        for i in range(rows):
            if i!=r and A[i,c]:
                A[i]=(A[i]-A[i,c]*A[r])%p
        r+=1
        if r==rows: break
    return r
def rank(Ps):
    return max(rank_mod(Ps[k],PR[k]) for k in range(2))
maxk=int(sys.argv[1])
level={(i,):[MM[0][i],MM[1][i]] for i in range(n)}
rk={w:rank(P) for w,P in level.items()}
alive={w:P for w,P in level.items() if rk[w]>=2}
print('k=1 alive',len(alive),'ranks',sorted(set(rk.values())))
minimal=[w for w in level if rk[w]<2]
for k in range(2,maxk+1):
    new={}; mins=[]
    for w,P in alive.items():
        for i in range(n):
            w2=w+(i,)
            if w2[1:] not in alive: continue   # suffix factor must be alive
            P2=[(MM[t][i]@P[t])%PR[t] for t in range(2)]
            r=rank(P2)
            if r>=2: new[w2]=P2
            else: mins.append(w2)
    alive=new
    from collections import Counter
    print('k',k,'alive',len(alive),'new minimal forbidden',len(mins))
    if len(mins)<=40:
        for m in mins: print('   ',' | '.join(U[i] for i in m))
    else:
        c=Counter(); 
        for m in mins: c[len(set(m))]+=1
        print('   distinct-unit counts among minimal:',dict(c))
    sys.stdout.flush()
    pickle.dump({'alive':list(alive.keys())},open(f'alive_k{k}.pkl','wb'))
