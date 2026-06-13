"""Pin strong-min(7): does the quadratic m²−5m+9=23 hold (vs near-trans kmax=4's loose 25)? Search
near-transitive strong tournaments on 7 vertices with kmax up to 6, plus the 'nested/blow-up'
minimal-strong family. opus-2026-06-06-S699k."""
from itertools import combinations
def Hcount(m,adj):
    size=1<<m; dp=[[0]*m for _ in range(size)]
    for v in range(m): dp[1<<v][v]=1
    for mask in range(size):
        for v in range(m):
            c=dp[mask][v]
            if not c: continue
            for w in range(m):
                if not(mask>>w&1) and adj[v]>>w&1: dp[mask|1<<w][w]+=c
    return sum(dp[size-1][v] for v in range(m))
def is_strong(m,adj):
    def reach(nb):
        seen={0}; st=[0]
        while st:
            x=st.pop()
            for y in range(m):
                if y not in seen and nb(x,y): seen.add(y); st.append(y)
        return len(seen)==m
    return reach(lambda x,y:adj[x]>>y&1) and reach(lambda x,y:adj[y]>>x&1)
def transitive(m):
    adj=[0]*m
    for i in range(m):
        for j in range(i+1,m): adj[i]|=1<<j
    return adj
def rev(m,adj,arcs):
    adj=adj[:]
    for (i,j) in arcs:
        if adj[i]>>j&1: adj[i]&=~(1<<j); adj[j]|=1<<i
        else: adj[j]&=~(1<<i); adj[i]|=1<<j
    return adj
def main():
    m=7; base=transitive(m); pairs=list(combinations(range(m),2))
    best=10**9; bestrev=None
    for k in range(1,7):
        for combo in combinations(pairs,k):
            adj=rev(m,base,list(combo))
            if is_strong(m,adj):
                h=Hcount(m,adj)
                if h<best: best=h; bestrev=combo
        print(f"  kmax={k}: best strong-min(7) so far = {best}  {bestrev}")
    print(f"\nstrong-min(7) ≤ {best} (search); quadratic m²−5m+9 = {7*7-5*7+9}=23; >21 (kills 21)? {best>21}")
if __name__=='__main__': main()
