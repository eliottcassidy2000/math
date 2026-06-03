"""The repo's OWN foundational theorem is the SOLVED GF(2) face of the (★) cancellation:
Rédei — every tournament has an ODD number of Hamiltonian paths. The path-count is a
permanent-shaped sum; its PARITY is forced (= a sign-reversing involution / determinant-
over-GF(2) cancellation), while the count itself varies wildly. This is what 'resolving the
cancellation' looks like: the unsolved LRC/Collatz are the real/arithmetic face of the same
shape. opus-2026-06-03-S599e."""
from itertools import combinations, product
import random
def ham_path_count(n, adj):
    # adj[u][v]=1 if u->v. dp over subsets: paths covering 'mask' ending at v.
    size=1<<n
    # dp[v] dict keyed by mask -> count; iterate masks ascending
    dp=[[0]*n for _ in range(size)]
    for v in range(n): dp[1<<v][v]=1
    for mask in range(size):
        for v in range(n):
            c=dp[mask][v]
            if not c: continue
            for w in range(n):
                if mask&(1<<w): continue
                if adj[v][w]:
                    dp[mask|(1<<w)][w]+=c
    full=size-1
    return sum(dp[full][v] for v in range(n))
def all_tournaments(n):
    pairs=list(combinations(range(n),2))
    for bits in product((0,1),repeat=len(pairs)):
        adj=[[0]*n for _ in range(n)]
        for (i,j),b in zip(pairs,bits):
            if b: adj[i][j]=1
            else: adj[j][i]=1
        yield adj
def rand_tournament(n,rng):
    adj=[[0]*n for _ in range(n)]
    for i,j in combinations(range(n),2):
        if rng.random()<0.5: adj[i][j]=1
        else: adj[j][i]=1
    return adj
def main():
    print("RÉDEI = the SOLVED GF(2)-cancellation face: #HamPaths(T) is ODD for every tournament T")
    for n in (3,4,5):
        odd=0; tot=0; counts=set()
        for adj in all_tournaments(n):
            c=ham_path_count(n,adj); tot+=1; counts.add(c)
            if c%2==1: odd+=1
        print(f"  n={n}: tournaments={tot}; all odd = {odd==tot}; distinct path-counts seen = {sorted(counts)}  (parity FIXED=1, count VARIES)")
    rng=random.Random(12345)
    for n in (6,7):
        odd=0; tot=0; counts=set()
        for _ in range(4000):
            adj=rand_tournament(n,rng); c=ham_path_count(n,adj); tot+=1; counts.add(c%2); 
            if c%2==1: odd+=1
        print(f"  n={n} (sample {tot}): all odd = {odd==tot}; parities seen = {sorted(counts)}")
    print("\nReading: the count = a permanent-shaped sum (varies a lot); its PARITY is an")
    print("invariant forced by cancellation (sign-reversing involution / det over GF(2)) = ⊕P.")
    print("LRC's p_0=Σ(−1)^|S|meas(∩) and Collatz's 2^E−3^k are the UNSOLVED real/arithmetic")
    print("face of the same (★): a signed sum whose sign must be controlled all-orders.")
if __name__=='__main__': main()
