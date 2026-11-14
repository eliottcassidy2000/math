# Frame decomposition of two-hub words: H = sum over hub orders and frames of prod r(i) (r in {1,2,3}).
import itertools, sys
from collections import Counter
def frames(x,y,e):
    m=len(x); out=[]
    for order in (0,1):
        X,Y = (x,y) if order==0 else (y,x)       # hub h1 uses X, h2 uses Y
        ok12 = (e==1) if order==0 else (e==0)     # h1 -> h2 needed when S1 empty
        A=[None]+[i for i in range(m) if X[i]==0] # a = max S0 (beats h1)
        B=[i for i in range(m) if X[i]==1]        # b = min S1 (beaten by h1)
        C=[i for i in range(m) if Y[i]==0]        # c = max S1 (beats h2)
        D=[None]+[i for i in range(m) if Y[i]==1] # d = min S2 (beaten by h2)
        S1opts=[(b,c) for b in B for c in C if b<=c]+([None] if ok12 else [])
        for a in A:
            for d in D:
                for bc in S1opts:
                    ends=set(z for z in (a,d) if z is not None)
                    if bc is not None:
                        b,c=bc
                        if b in ends or c in ends: continue
                        ends|={b,c}
                    if a is not None and d is not None and a==d: continue
                    w=1; al=be=0
                    for i in range(m):
                        if i in ends: continue
                        r=(a is not None and i<a) + (bc is not None and bc[0]<i<bc[1]) + (d is not None and i>d)
                        if r==0: w=0; break
                        if r==2: al+=1
                        if r==3: be+=1
                        w*=r
                    # also S0 empty: a None means S0 empty -> vertices can't be in S0; S2 empty similarly
                    if w: out.append((order,a,bc,d,al,be))
    return out
def H_dp(x,y,e):
    return sum(2**f[4]*3**f[5] for f in frames(x,y,e))
if __name__=="__main__":
    import random
    # brute-force check of the frame identity against direct Hamiltonian path counting on small words
    def ham_count(x,y,e):
        m=len(x); n=m+2; u1=m; u2=m+1
        beats=[[False]*n for _ in range(n)]
        for i in range(m):
            for j in range(i+1,m): beats[i][j]=True
            if x[i]: beats[u1][i]=True
            else: beats[i][u1]=True
            if y[i]: beats[u2][i]=True
            else: beats[i][u2]=True
        if e: beats[u1][u2]=True
        else: beats[u2][u1]=True
        cnt=0
        for p in itertools.permutations(range(n)):
            if all(beats[p[k]][p[k+1]] for k in range(n-1)): cnt+=1
        return cnt
    random.seed(5); bad=0
    for _ in range(150):
        m=random.randint(1,6); x=[random.getrandbits(1) for _ in range(m)]; y=[random.getrandbits(1) for _ in range(m)]; e=random.getrandbits(1)
        if H_dp(x,y,e)!=ham_count(x,y,e): bad+=1; print("mismatch",x,y,e,H_dp(x,y,e),ham_count(x,y,e))
    print("frame identity vs brute force (150 random words, m<=6): mismatches",bad)
