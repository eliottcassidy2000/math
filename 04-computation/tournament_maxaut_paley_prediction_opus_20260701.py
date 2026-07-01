"""Paley-prime prediction (coding-theory). max|Aut| of circulant tournaments via the MULTIPLIER group:
Aut >= {x->ax+b : aS=S mod n}, |Aut|>=n*|mult(S)|. Peaks at Paley primes p=3 mod4 (|Aut|=p(p-1)/2).
Plus MFAS(Paley_p) exact(7)/heuristic(11) = covering radius at Paley primes."""
import numpy as np, itertools, random
def circulant_auts(n):
    # valid tournament connection sets S subset Z_n\{0}: exactly one of {d,-d} in S for each pair
    if n%2==0: return None
    half=(n-1)//2; pairs=[(d,n-d) for d in range(1,half+1)]
    best=(0,None)
    for choice in itertools.product(*[(d,nd) for d,nd in pairs]):
        S=set(choice)
        mult=[a for a in range(1,n) if set((a*s)%n for s in S)==S]  # multipliers fixing S
        aut=n*len(mult)
        if aut>best[0]: best=(aut,sorted(S),sorted(mult))
    return best
def mfas_exact_or_heur(S,n,exact_limit=8):
    A=np.zeros((n,n),np.int8)
    for i in range(n):
        for s in S: A[i,(i+s)%n]=1
    if n<=exact_limit:
        from itertools import permutations
        best=99
        for p in permutations(range(n)):
            B=A[np.ix_(p,p)]; best=min(best,int(np.tril(B,-1).sum()))
        return best,"exact"
    # heuristic: random restarts + adjacent-swap local search
    random.seed(0); best=99
    for _ in range(4000):
        p=list(range(n)); random.shuffle(p)
        improved=True
        while improved:
            improved=False
            for i in range(n-1):
                q=p[:]; q[i],q[i+1]=q[i+1],q[i]
                if np.tril(A[np.ix_(q,q)],-1).sum()<np.tril(A[np.ix_(p,p)],-1).sum():
                    p=q; improved=True
        best=min(best,int(np.tril(A[np.ix_(p,p)],-1).sum()))
    return best,"heur(UB)"
print("max|Aut| of circulant tournaments (lower bound on true max|Aut|), and MFAS = covering radius:")
print(f"{'n':>3} {'best|Aut|':>10} {'S':>18} {'Paley?':>7} {'MFAS(that S)':>14}")
exact_maxaut={3:3,4:3,5:5,6:9,7:21}   # exact max|Aut| n<=7 (enumeration)
for n in [3,5,7,9,11]:
    aut,S,mult=circulant_auts(n)
    isp = "yes" if (n in (7,11,19,23)) else ""
    mf,how=mfas_exact_or_heur(S,n)
    print(f"{n:>3} {aut:>10} {str(S):>18} {isp:>7} {mf:>10} {how}")
print(f"\nExact max|Aut| (all tournaments, n<=7): 3,3,5,9,21 (n=3..7). Circulant matches at odd n.")
print(f"Paley primes p=3 mod4: |Aut|=p(p-1)/2 = 7->21, 11->55, 19->171, 23->253 (KNOWN unique max).")
print(f"Between Paley primes: n=9 circulant |Aut|=27? (9*3), n=5 =5 -- much smaller peaks.")
print(f"=> max|Aut| (the covering-excess driver) PEAKS at Paley primes => predicted flip-rank k(n) jumps there.")
