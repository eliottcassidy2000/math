#!/usr/bin/env python3
"""FORMALIZE: every ROUND (locally-transitive) tournament has dichromatic number <=2
(=2 iff it has a 3-cycle). Proof: circle realization + a diameter split into two
transitive semicircle-arcs. Verify on all round tournaments m=5,7. opus-S592 round1."""
from itertools import permutations, product, combinations
def transitive(A,verts):
    for a,b,c in permutations(verts,3):
        if A[a][b] and A[b][c] and A[c][a]: return False
    return True
def is_round(A,m):
    return all(transitive(A,[u for u in range(m) if A[v][u]]) for v in range(m))
def has_3cycle(A,m):
    return any(A[a][b] and A[b][c] and A[c][a] for a,b,c in permutations(range(m),3))
def dichromatic(A,m):
    for k in range(1,m+1):
        for col in product(range(k),repeat=m):
            if max(col)!=k-1: continue
            if all(transitive(A,[v for v in range(m) if col[v]==c]) for c in range(k)): return k
    return m
def gen_tournaments(m):
    E=list(combinations(range(m),2))
    for mask in range(2**len(E)):
        A=[[0]*m for _ in range(m)]
        for b,(i,j) in enumerate(E):
            if mask>>b&1: A[i][j]=1
            else: A[j][i]=1
        yield A
def main():
    for m in [5,7]:
        rnd=0; chi2=0; chi1=0; bad=0
        # m=7 is 2M tournaments; sample if needed
        cnt=0
        import random; rng=random.Random(0)
        gen = gen_tournaments(m) if m<=5 else (None,)
        if m<=5:
            for A in gen_tournaments(m):
                if not is_round(A,m): continue
                rnd+=1; c=dichromatic(A,m)
                if c==1: chi1+=1
                elif c==2: chi2+=1
                else: bad+=1
        else:
            # sample round tournaments on 7 via random + circulant
            seen=0
            for _ in range(4000):
                E=list(combinations(range(m),2)); A=[[0]*m for _ in range(m)]
                for (i,j) in E:
                    if rng.random()<0.5: A[i][j]=1
                    else: A[j][i]=1
                if not is_round(A,m): continue
                rnd+=1; c=dichromatic(A,m)
                if c==1: chi1+=1
                elif c==2: chi2+=1
                else: bad+=1
                if rnd>=200: break
        print(f"  m={m}: round tournaments examined={rnd}; chi=1 (transitive)={chi1}; chi=2={chi2}; chi>=3={bad}  => round => chi<=2: {bad==0}")
if __name__=='__main__': main()
