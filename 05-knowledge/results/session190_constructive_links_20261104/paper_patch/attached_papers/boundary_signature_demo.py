"""Independent finite-state example transcribed from Fig. 2, Durias et al. (2026).
Uses Python standard library only; does not search global Hamiltonian cycles.
Terminals are outer corners 0,2,4,6 in cyclic order. Outer cycle is 0..7.
"""
from collections import defaultdict
import json

T=(0,2,4,6)
H=(10, [(i,(i+1)%8) for i in range(8)]+[(1,8),(8,7),(8,3),(7,9),(3,9),(9,5)])
H1=(12, [(i,(i+1)%8) for i in range(8)]+[(1,8),(8,7),(8,3),(8,9),(9,10),(9,11),(7,10),(3,11),(10,5),(11,5)])

def signatures(patch, deleted=None):
    n,ee=patch
    vertices=[v for v in range(n) if v!=deleted]
    terminals=tuple(v for v in T if v!=deleted)
    edges=[(a,b) for a,b in ee if deleted not in (a,b)]
    deg=[0]*n
    left=[0]*n
    for a,b in edges: left[a]+=1; left[b]+=1
    out={}
    chosen=[]
    def visit(i):
        if any(deg[v]>2 or (v not in terminals and deg[v]+left[v]<2) for v in vertices): return
        if i==len(edges):
            if any(deg[v]!=2 for v in vertices if v not in terminals): return
            adj=defaultdict(list)
            for a,b in chosen: adj[a].append(b); adj[b].append(a)
            unseen=set(vertices); pairs=[]
            while unseen:
                start=next(iter(unseen)); stack=[start]; comp=set()
                while stack:
                    v=stack.pop()
                    if v in comp: continue
                    comp.add(v); stack.extend(adj[v])
                unseen-=comp
                ends=sorted(v for v in comp if deg[v]==1)
                if len(comp)>1:
                    if len(ends)!=2: return # Reject cycle components.
                    pairs.append(tuple(ends))
                elif deg[start]!=0 or start not in terminals: return
            sig=(tuple(deg[t] for t in terminals),tuple(sorted(pairs)))
            out.setdefault(sig,tuple(chosen))
            return
        a,b=edges[i]; left[a]-=1; left[b]-=1
        visit(i+1)
        if deg[a]<2 and deg[b]<2:
            deg[a]+=1;deg[b]+=1;chosen.append((a,b));visit(i+1);chosen.pop();deg[a]-=1;deg[b]-=1
        left[a]+=1;left[b]+=1
    visit(0)
    return out

base=signatures(H);replacement=signatures(H1)
assert base.keys()==replacement.keys()
assert len(base)==14
terminal_checks={str(t):set(signatures(H,t))<=set(signatures(H1,t)) for t in T}
assert all(terminal_checks.values())
old={w:set(signatures(H,w)) for w in range(H[0]) if w not in T}
new={v:set(signatures(H1,v)) for v in range(H1[0]) if v not in T}
coverage={str(v):[w for w,s in old.items() if s<=ns] for v,ns in new.items()}
assert all(coverage.values())
result={
    'scope':'Standard-library check of two patches transcribed from Figure 2, not a verification of the full paper.',
    'terminals_in_cyclic_order':T,
    'H':{'order':H[0],'edges':H[1]},'H1':{'order':H1[0],'edges':H1[1]},
    'full_signature_count':len(base),'equal_full_signatures':True,
    'terminal_deletion_inclusions':terminal_checks,
    'internal_deletion_coverage_new_to_old':coverage,
    'old_deletion_signature_counts':{str(v):len(s) for v,s in old.items()},
    'new_deletion_signature_counts':{str(v):len(s) for v,s in new.items()},
    'full_signature_witnesses':[{'degree_vector':sig[0],'endpoint_pairs':sig[1],'H_edges':base[sig],'H1_edges':replacement[sig]} for sig in sorted(base)]
}
print(json.dumps(result,indent=2))
