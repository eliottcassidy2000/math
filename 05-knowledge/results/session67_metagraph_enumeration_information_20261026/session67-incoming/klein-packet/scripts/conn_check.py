"""Connectivity lemma check: components of the 3-cycle hypergraph (singletons for isolated vertices) == strong components."""
import subprocess, itertools, sys
def comps_union(n, edges):
    parent=list(range(n))
    def f(x):
        while parent[x]!=x: parent[x]=parent[parent[x]]; x=parent[x]
        return x
    for e in edges:
        for a,b in zip(e,e[1:]): parent[f(a)]=f(b)
    return frozenset(frozenset(v for v in range(n) if f(v)==r) for r in set(f(v) for v in range(n)))
def scc(n, out):
    # Landau cuts
    sc=sorted(range(n), key=lambda v: bin(out[v]).count('1'))
    comps=[]; cum=0; cur=[]
    for i,v in enumerate(sc):
        cum+=bin(out[v]).count('1'); cur.append(v)
        if cum==(i+1)*i//2: comps.append(frozenset(cur)); cur=[]
    return frozenset(comps)
for n in range(3,9):
    lines=subprocess.run(['gentourng',str(n)],capture_output=True,text=True).stdout.split()
    bad=0
    for line in lines:
        out=[0]*n; p=0
        for i in range(n):
            for j in range(i+1,n):
                if line[p]=='1': out[i]|=1<<j
                else: out[j]|=1<<i
                p+=1
        cyc=[t for t in itertools.combinations(range(n),3) if sum((out[t[0]]>>t[1])&1 for _ in [0]) is not None and
             (((out[t[0]]>>t[1])&1 and (out[t[1]]>>t[2])&1 and (out[t[2]]>>t[0])&1) or (not (out[t[0]]>>t[1])&1 and not (out[t[1]]>>t[2])&1 and not (out[t[2]]>>t[0])&1))]
        if comps_union(n,cyc)!=scc(n,out): bad+=1
    print(f'n={n}: classes={len(lines)} violations of (hypergraph components == strong components): {bad}')
