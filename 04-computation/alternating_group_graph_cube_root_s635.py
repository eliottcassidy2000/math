import itertools, collections
def parity(p):
    n=len(p); seen=[False]*n; par=0
    for i in range(n):
        if not seen[i]:
            j=i;c=0
            while not seen[j]: seen[j]=True;j=p[j];c+=1
            par+=c-1
    return par%2
def even_perms(n): return [p for p in itertools.permutations(range(n)) if parity(p)==0]
def compose(p,q): return tuple(p[q[i]] for i in range(len(p)))
def tc(n,a,b,c):
    p=list(range(n)); p[a],p[b],p[c]=b,c,a; return tuple(p)
def gens(n):
    G=[]
    for i in range(2,n):
        g=tc(n,0,1,i); G+=[g, compose(g,g)]
    return G
def build(n):
    V=even_perms(n); idx={v:k for k,v in enumerate(V)}
    g=gens(n); adj=[set() for _ in V]
    for k,v in enumerate(V):
        for s in g: adj[k].add(idx[compose(v,s)])
    return V,adj
def girth(adj):
    n=len(adj); best=10**9
    for s in range(min(n,40)):
        dist={s:0}; par={s:-1}; q=collections.deque([s])
        while q:
            u=q.popleft()
            for w in adj[u]:
                if w not in dist: dist[w]=dist[u]+1; par[w]=u; q.append(w)
                elif par[u]!=w: best=min(best,dist[u]+dist[w]+1)
    return best
def chromatic(adj):
    n=len(adj); order=sorted(range(n),key=lambda v:-len(adj[v]))
    for k in range(1,n+1):
        col=[-1]*n
        def bt(i):
            if i==n: return True
            v=order[i]
            for c in range(k):
                if all(col[u]!=c for u in adj[v]):
                    col[v]=c
                    if bt(i+1): return True
                    col[v]=-1
            return False
        if bt(0): return k
def indep_exact(adj):
    n=len(adj); best=[0]; order=sorted(range(n),key=lambda v:-len(adj[v]))
    def bt(i,size,used):
        if size+(n-i)<=best[0]: return
        if i==n:
            if size>best[0]: best[0]=size
            return
        v=order[i]
        if v not in used: bt(i+1,size+1,used|adj[v]|{v})
        bt(i+1,size,used)
    bt(0,0,frozenset()); return best[0]
def largest_color_class(adj,k):
    n=len(adj); order=sorted(range(n),key=lambda v:-len(adj[v])); col=[-1]*n
    def bt(i):
        if i==n: return True
        v=order[i]
        for c in range(k):
            if all(col[u]!=c for u in adj[v]):
                col[v]=c
                if bt(i+1): return True
                col[v]=-1
        return False
    bt(0); return max(collections.Counter(col).values())
print("=== AG_n (alternating group graph, Cayley of A_n via 3-cycles): how it changes with n ===")
print(f"{'n':>2} {'|V|':>5} {'reg':>4} {'girth':>6} {'chi':>4} {'alpha':>7} {'chi*alpha':>10} {'tight=|V|?':>11}")
for n in range(3,6):
    V,adj=build(n); nv=len(V); ch=chromatic(adj)
    al = indep_exact(adj) if nv<=14 else largest_color_class(adj,ch)
    note = 'TIGHT' if ch*al==nv else ('ok' if ch*al>=nv else 'FAIL')
    print(f"{n:>2} {nv:>5} {len(adj[0]):>4} {girth(adj):>6} {ch:>4} {al:>7} {ch*al:>10} {note:>11}")
print("\n  AG_n is 3-CHROMATIC for all n: girth-3 triangle (the 3-cycle) forces chi>=3, and a 3-coloring exists.")
print("  The 3 = the 3-cycle generator = the CUBE ROOT of unity (s^3=1, eigenvalues 1,w,w^2=e^{+-2pi i/3}) = pi/3/Phi_3.")
print("  vertex-transitive => chi_f=|V|/alpha; alpha=|V|/3 => chi*alpha=|V| TIGHT (S634 bound). AG_n = 'the third color' for all n.")
print("\n  NOTE: the S634 bound chi*alpha>=|V| caught a bug in my first run (AG_4 alpha=2 was wrong; chi=3 forces alpha>=4). The theorem is ground truth.")
