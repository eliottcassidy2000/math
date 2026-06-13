"""PROPAGATION: do the H-gaps (7,21) propagate to the LRC worry-set? The worry-set = self-converse
ROUND tournaments (THM-402/407). If round tournaments AVOID H-values 7,21, then the LRC worry-set
INHERITS the impossibility — a concrete propagation channel. opus-2026-06-03-S599r."""
from itertools import combinations, product, permutations
def Hcount(n,adj):
    size=1<<n; dp=[[0]*n for _ in range(size)]
    for v in range(n): dp[1<<v][v]=1
    for mask in range(size):
        row=dp[mask]
        for v in range(n):
            c=row[v]
            if not c: continue
            av=adj[v]
            for w in range(n):
                if not(mask>>w&1) and av>>w&1: dp[mask|1<<w][w]+=c
    return sum(dp[size-1][v] for v in range(n))
def is_round(n,adj):  # locally transitive: out-nbhd has no 3-cycle
    for v in range(n):
        out=[w for w in range(n) if adj[v]>>w&1]
        for a,b,c in combinations(out,3):
            s=((adj[a]>>b&1)+(adj[b]>>c&1)+(adj[c]>>a&1))
            if s==3 or s==0: return False
    return True
def main():
    print("Do ROUND tournaments (= LRC worry-set carriers) avoid the H-gaps 7,21,...?")
    print(" n | round H-values | contains 7? 21? | all-tournament H-values (for contrast)")
    for n in range(3,8):
        rv=set(); allv=set()
        for bits in product((0,1),repeat=n*(n-1)//2):
            adj=[0]*n
            for (i,j),b in zip(combinations(range(n),2),bits):
                if b: adj[i]|=1<<j
                else: adj[j]|=1<<i
            h=Hcount(n,adj); allv.add(h)
            if is_round(n,adj): rv.add(h)
        print(f" {n} | {sorted(rv)} | 7:{7 in rv} 21:{21 in rv}")
    print("\n=> if round H-values exclude 7,21 at all n, the LRC worry-set provably cannot realize")
    print("   those tournament shapes: the H-impossibility PROPAGATES to LRC via THM-402/407.")
if __name__=='__main__': main()
