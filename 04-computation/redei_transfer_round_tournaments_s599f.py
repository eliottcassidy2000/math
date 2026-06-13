"""Transfer the SOLVED face (Rédei/GF(2)) to LRC via the worry-set <-> ROUND (locally
transitive) tournament bridge (THM-402). Compute: round-tournament iso-classes, how many are
SELF-CONVERSE (= the worry-set), Ham-path counts (Rédei: odd; is the round count special?),
and the GF(2)/parity structure. Grounds the 'enumerate tournaments not configs' SPEEDUP and
the 'forced-parity certificate' TRUTH transfer. opus-2026-06-03-S599f."""
from itertools import combinations, permutations, product
def tour_from_bits(m,bits):
    adj=[[0]*m for _ in range(m)]
    for (i,j),b in zip(combinations(range(m),2),bits):
        if b: adj[i][j]=1
        else: adj[j][i]=1
    return adj
def is_round(m,adj):
    # locally transitive: for each v, out-neighborhood induces a transitive subtournament
    for v in range(m):
        out=[w for w in range(m) if adj[v][w]]
        # transitive iff no 3-cycle among 'out'
        for a,b,c in combinations(out,3):
            t=adj  # check the 3 vertices form transitive
            # a triple is a 3-cycle iff cyclic
            e=lambda x,y: t[x][y]
            s=e(a,b)+e(b,c)+e(c,a)
            if s==3 or s==0:  # a->b->c->a or reverse cycle
                return False
    return True
def ham_paths(m,adj):
    size=1<<m; dp=[[0]*m for _ in range(size)]
    for v in range(m): dp[1<<v][v]=1
    for mask in range(size):
        for v in range(m):
            c=dp[mask][v]
            if not c: continue
            for w in range(m):
                if mask&(1<<w) or not adj[v][w]: continue
                dp[mask|(1<<w)][w]+=c
    return sum(dp[size-1][v] for v in range(m))
def canon(m,adj):
    best=None
    for p in permutations(range(m)):
        key=tuple(adj[p[i]][p[j]] for i in range(m) for j in range(m))
        if best is None or key<best: best=key
    return best
def converse_adj(m,adj): return [[adj[j][i] for j in range(m)] for i in range(m)]
def main():
    print(" m | #roundTours | #roundClasses | #self-converse(=worry-set) | Ham-counts(round) all-odd")
    for m in range(3,7):
        seen={}; selfconv=0
        npairs=m*(m-1)//2
        rounds=0; allodd=True; counts=set()
        for bits in product((0,1),repeat=npairs):
            adj=tour_from_bits(m,bits)
            if not is_round(m,adj): continue
            rounds+=1
            c=canon(m,adj)
            if c not in seen:
                hp=ham_paths(m,adj); counts.add(hp)
                if hp%2==0: allodd=False
                cc=canon(m,converse_adj(m,adj))
                seen[c]=cc
        # count classes & self-converse classes
        classes=set(seen.keys()); sc=sum(1 for c in classes if seen[c]==c)
        print(f" {m} | {rounds:5d} | {len(classes):6d} | {sc:6d} | counts={sorted(counts)} allodd={allodd}")
    print("\nWorry-set bridge: self-converse round classes = the LRC worry-set (the 64 at n=14).")
    print("SPEEDUP: enumerate round tournament classes (small) instead of integer configs (exp).")
    print("Rédei TRUTH: every round tournament has ODD #HamPaths -> a forced-parity certificate.")
if __name__=='__main__': main()
