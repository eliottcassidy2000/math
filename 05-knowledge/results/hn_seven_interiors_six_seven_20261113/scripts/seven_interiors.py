# Seven private interiors.  New quotient cases: (i) 11-vertex 7-critical graphs, 4 port groups + 7 interiors;
# (ii) 10-vertex 7-critical graphs, 3 port groups + 7 interiors.  Same planar necessary conditions as six_interiors.py.
import subprocess, itertools, sys
def g6_to_adj(s):
    data=[ord(c)-63 for c in s]; n=data[0]; bits=[]
    for d in data[1:]:
        for k in range(5,-1,-1): bits.append((d>>k)&1)
    adj=[0]*n; idx=0
    for j in range(1,n):
        for i in range(j):
            if bits[idx]: adj[i]|=1<<j; adj[j]|=1<<i
            idx+=1
    return n,adj
def colourable(n,adj,k,mask):
    verts=[v for v in range(n) if mask>>v&1]
    verts.sort(key=lambda v:-bin(adj[v]&mask).count('1'))
    col={}
    def bt(i):
        if i==len(verts): return True
        v=verts[i]; used=0
        for u,c in col.items():
            if adj[v]>>u&1: used|=1<<c
        for c in range(k):
            if not used>>c&1:
                col[v]=c
                if bt(i+1): return True
                del col[v]
        return False
    return bt(0)
def critical_graphs(n,D):
    out=subprocess.run(['geng','-q',str(n),f'-D{D}'],capture_output=True,text=True).stdout.split()
    full=(1<<n)-1; res=[]
    for s in out:
        m,gam=g6_to_adj(s)
        adj=[(full^gam[v])&~(1<<v) for v in range(m)]   # complement
        if colourable(m,adj,6,full): continue
        if all(colourable(m,adj,6,full&~(1<<v)) for v in range(m)): res.append(adj)
    return res
def is_clique(adj,S): return all(adj[a]>>b&1 for a,b in itertools.combinations(S,2))
def popc(x): return bin(x).count('1')
def circle_ok(adj,T,Imask):
    # T a bitmask of interiors on one unit circle: induced unit graph is a linear forest or a 6-cycle
    vs=[v for v in range(len(adj)) if T>>v&1]
    if len(vs)<=1: return True
    degs=[popc(adj[v]&T) for v in vs]; e=sum(degs)//2
    if any(d>2 for d in degs): return False
    if len(vs)==6 and e==6 and all(d==2 for d in degs):   # must be one 6-cycle
        seen=1<<vs[0]; st=[vs[0]]
        while st:
            u=st.pop(); nb=adj[u]&T&~seen
            while nb:
                b=nb&-nb; w=b.bit_length()-1; seen|=b; st.append(w); nb^=b
        return seen==T
    # forest check: edges = vertices - components
    seen=0; comps=0
    for v in vs:
        if seen>>v&1: continue
        comps+=1; st=[v]; seen|=1<<v
        while st:
            u=st.pop(); nb=adj[u]&T&~seen
            while nb:
                b=nb&-nb; w=b.bit_length()-1; seen|=b; st.append(w); nb^=b
    return e==len(vs)-comps
def realizable(adj,P,I):
    n=len(adj); Imask=sum(1<<v for v in I)
    # interiors: no unit K4, no K_{2,3}
    for Q in itertools.combinations(I,4):
        if is_clique(adj,Q): return False
    for a,b in itertools.combinations(I,2):
        if popc(adj[a]&adj[b]&Imask)>=3: return False
    groups=[]
    for g in P:
        S=[x for x in I if adj[g]>>x&1]; opts=[]
        for assign in itertools.product((1,2,3),repeat=len(S)):
            T1=0;T2=0
            for x,a in zip(S,assign):
                if a&1: T1|=1<<x
                if a&2: T2|=1<<x
            if T2<T1: continue
            if popc(T1&T2)>2: continue
            if circle_ok(adj,T1,Imask) and circle_ok(adj,T2,Imask): opts.append((T1,T2))
        if not opts: return False
        groups.append(opts)
    # each port point vs interior x: common interior neighbours <= 2
    def point_ok(T):
        return all(popc(adj[x]&Imask&T)<=2 for x in I)
    groups=[[o for o in opts if point_ok(o[0]) and point_ok(o[1])] for opts in groups]
    if any(not g for g in groups): return False
    groups.sort(key=len)
    chosen=[]
    def bt(i):
        if i==len(groups):
            # two interiors: common nbrs among interiors + port points <= 2
            for a,b in itertools.combinations(I,2):
                c=popc(adj[a]&adj[b]&Imask)+sum(1 for T in chosen if T>>a&1 and T>>b&1)
                if c>2: return False
            return True
        for (T1,T2) in groups[i]:
            if any(popc(T1&T)>2 or popc(T2&T)>2 for T in chosen): continue
            chosen.extend([T1,T2])
            if bt(i+1): return True
            del chosen[-2:]
        return False
    return bt(0)
def _main():
  pass
results={}
if __name__=="__main__":
  for (n,D,k) in ((10,3,3),(11,4,4)):
     crit=critical_graphs(n,D)
     surv=[]
     for gi,adj in enumerate(crit):
         for P in itertools.combinations(range(n),k):
             if not is_clique(adj,P): continue
             I=[v for v in range(n) if v not in P]
             if realizable(adj,P,I): surv.append((gi,P))
     print(f"{n}-vertex 7-critical graphs: {len(crit)}; port groups k={k}, interiors={n-k}: surviving (graph,ports) pairs = {len(surv)}")
     for gi,P in surv[:20]:
         adj=crit[gi]; I=[v for v in range(n) if v not in P]
         print("   survivor graph",gi,"ports",P,"interior unit edges",sum(popc(adj[x]&sum(1<<y for y in I)) for x in I)//2)
     results[n]=(crit,surv)
