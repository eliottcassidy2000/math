"""Independent exhaustive fixed-path tilings; exact rational memory tests.
Run: python3 metagraph-memory/check.py
No repository imports, no numerical ranks, no external data.
"""
import itertools as it, json
from pathlib import Path
import sympy as s

def require(value):
    if not value: raise RuntimeError("Exact certificate check failed")

def census(n, linear=True):
    pairs=list(it.combinations(range(n),2)); free=[(i,j) for i,j in pairs if j>i+1]
    N=1<<len(free); full=N-1
    codes=[]; colors=[]
    for x in range(N):
        adj=[[0]*n for _ in range(n)]
        for i in range(n-1): adj[i][i+1]=1
        for k,(i,j) in enumerate(free):
            a,b=(i,j) if (x>>k)&1 else (j,i); adj[a][b]=1
        # Complete relabelling and converse minimization; degree sorting prunes only
        # permutations which cannot attain the lexicographically least degree row.
        # Canonical code includes sorted outdegrees, then adjacency bits.
        def canon(a):
            deg=[sum(row) for row in a]; blocks=[[i for i in range(n) if deg[i]==d] for d in sorted(set(deg))]
            best=None
            for chunks in it.product(*(it.permutations(b) for b in blocks)):
                p=sum((tuple(c) for c in chunks),()); code=tuple(a[p[i]][p[j]] for i,j in pairs)
                if best is None or code<best: best=code
            return tuple(sorted(deg)),best
        codes.append(min(canon(adj),canon(list(map(list,zip(*adj))))))
        refl=sum(((x>>free.index((n-1-j,n-1-i)))&1)<<k for k,(i,j) in enumerate(free))
        colors.append('blue' if refl==x else 'black')
    keys=sorted(set(codes)); node=[keys.index(c) for c in codes]; r=len(keys)
    cells={}; colored={}
    for x in range(N):
        key=(node[x],node[full^x]); cells.setdefault(key,[]).append(x)
        colored.setdefault(key,set()).add(colors[x])
    # Explicit optimal positive realization: point encodings and pair swap.
    pair_keys=sorted(cells); pair_index={p:i for i,p in enumerate(pair_keys)}
    swap=[pair_index[(j,i)] for i,j in pair_keys]
    require(all(swap[swap[a]]==a for a in range(len(swap))))
    for x in range(N):
        a=pair_index[(node[x],node[full^x])]
        for time in range(7):
            require(pair_keys[a][0]==node[x if time%2==0 else full^x])
            a=swap[a]
    mixed=[dict(pair=list(k),blue=next(x for x in cells[k] if colors[x]=='blue'),black=next(x for x in cells[k] if colors[x]=='black')) for k in cells if len(colored[k])==2]
    ans=dict(n=n,tilings=N,nodes=r,pair_states=len(cells),colored_pair_states=sum(map(len,colored.values())),mixed=mixed,positive_realization_checks=7*N)
    if not linear: return ans
    J=s.zeros(N,r); P=s.zeros(N,N)
    for x in range(N): J[x,node[x]]=1; P[x,full^x]=1
    F=J.T*J; A=J.T*P*J; Fi=F.inv(); T=Fi*A; E=J*Fi*J.T; Q=s.eye(N)-E
    defect=F-A*Fi*A; W=J.row_join(P*J)
    require(W.rank()==r+defect.rank())
    require(defect==(P*J-J*T).T*(P*J-J*T))
    # Rational-coordinate version of BD^kB* identity.
    C=Q*P*Q
    for k in range(5):
        require(J.T*P*Q*(C**k)*Q*P*J == (-1)**k*defect*(T**k))
    for t in range(7): require(Fi*J.T*(P**t)*J == (s.eye(r) if t%2==0 else T))
    adj={i:{j for j in range(r) if A[i,j]} for i in range(r)}; seen=set(); c=b=0
    for root in range(r):
        if root in seen: continue
        c+=1; shades={root:0}; todo=[root]; bip=True
        while todo:
            i=todo.pop(); seen.add(i)
            for j in adj[i]:
                if j in shades: bip &= shades[i]!=shades[j]
                else: shades[j]=1-shades[i]; todo.append(j)
        b+=bip
    require(defect.rank()==r-c-b)
    # Coarsest stable deterministic refinement from node observation.
    labels=node[:]
    while True:
        signatures=[(labels[x],labels[full^x]) for x in range(N)]
        dic={k:i for i,k in enumerate(sorted(set(signatures)))}; refined=[dic[k] for k in signatures]
        if len(set(refined))==len(set(labels)): break
        labels=refined
    require(len(set(refined))==len(cells))
    ans.update(components=c,bipartite_components=b,defect_rank=defect.rank(),linear_dimension=W.rank(),A=[list(map(int,A.row(i))) for i in range(r)],fibres=[int(F[i,i]) for i in range(r)],node=node,colors=colors,memory_checks=5,projected_time_checks=7,markov_failure=(T*T!=s.eye(r)))
    return ans

if __name__=='__main__':
    results=[census(n) for n in (3,4,5)]+[census(6,False)]
    out=Path(__file__).with_name('results.json'); out.write_text(json.dumps(results,indent=2)+'\n')
    for a in results: print({k:v for k,v in a.items() if k not in ('A','fibres','node','colors')})
