#!/usr/bin/env python3
"""Independent literal controls; no repository imports or removable assertions."""
from collections import Counter, defaultdict
from itertools import combinations, product
from pathlib import Path
import json

HERE = Path(__file__).resolve().parent

def need(ok, message):
    if not ok:
        raise RuntimeError(message)

def tournament(k):
    return frozenset((j,i) if (i,j) in {(0,k),(k,2*k)} else (i,j)
                     for i in range(2*k+1) for j in range(i+1,2*k+1))

def vertices(k): return range(2*k+1)
def states(k): return range((1 << (k-1))+1)
def infinity(k): return 1 << (k-1)

def choice(k,x):
    return set(range(k)) if x == infinity(k) else {
        i+1 for i in range(k-1) if x & (1 << i)}

def path(k,state):
    x,y = state
    A, B = set(range(k)), set(range(k+1,2*k+1))
    I,K = choice(k,x), choice(k,y)
    J = {2*k-v for v in A-K}
    return tuple(sorted(I)+sorted(J)+[k]+sorted(A-I)+sorted(B-J))

def literal_paths(k):
    E = tournament(k)
    def visit(p, remaining):
        if not remaining:
            yield p
        for v in remaining:
            if not p or (p[-1],v) in E:
                yield from visit(p+(v,), remaining-{v})
    return set(visit((),set(vertices(k))))

def valid(k,p):
    E=tournament(k)
    return set(p)==set(vertices(k)) and len(p)==2*k+1 and all(
        (a,b) in E for a,b in zip(p,p[1:]))

def blue(k,p):
    E=tournament(k)
    r=dict(zip(p,reversed(p)))
    return all((r[b],r[a]) in E for a,b in E)

def project_state(k,s):
    """child k+1 -> parent k"""
    return tuple(infinity(k) if x==infinity(k+1) else x % infinity(k) for x in s)

def toggle(k,s,axis):
    t=list(s)
    if t[axis] != infinity(k+1):
        t[axis] ^= infinity(k)
    return tuple(t)

def act(k,s,g):
    for axis in range(2):
        if g >> axis & 1:
            s=toggle(k,s,axis)
    return s

def embed(k,v):return v if v<k else v+1 if v==k else v+2

def shrink(k,v):return v if v<k else v-1 if v==k+1 else v-2

def delete_path(k,p):
    return tuple(shrink(k,v) for v in p if v not in {k,k+2})

def pair(a,b):return tuple(sorted((a,b)))
def edges(p):return frozenset(pair(a,b) for a,b in zip(p,p[1:]))

def reverse_edges(E,F):
    return frozenset((b,a) if pair(a,b) in F else (a,b) for a,b in E)

def delete_tournament(k,E):
    return frozenset((shrink(k,a),shrink(k,b)) for a,b in E
                     if a not in {k,k+2} and b not in {k,k+2})

def cut_sidecar(k,p):
    D={k,k+2}
    retained=frozenset(pair(shrink(k,a),shrink(k,b)) for a,b in zip(p,p[1:])
                       if a not in D and b not in D)
    return edges(delete_path(k,p))-retained

def ports_direct(k,p):
    """Independent scan of maximal deleted blocks, not edge subtraction."""
    D={k,k+2};out=set();i=0
    while i<len(p):
        if p[i] not in D:
            i+=1;continue
        j=i
        while j<len(p) and p[j] in D:j+=1
        if i>0 and j<len(p):out.add(pair(shrink(k,p[i-1]),shrink(k,p[j])))
        i=j
    return frozenset(out)

def two_seam_sidecar(k,s):
    parent=project_state(k,s)
    A,B=set(range(k)),set(range(k+1,2*k+1))
    I,K=choice(k,parent[0]),choice(k,parent[1])
    J={2*k-v for v in A-K}
    Arest,Brest=A-I,B-J
    pre=(max(I) if I else None,min(J) if J else k)
    post=(max(Arest) if Arest else k,min(Brest) if Brest else None)
    e=1 if s[0]==infinity(k+1) else int(bool(s[0]&infinity(k)))
    h=1 if s[1]==infinity(k+1) else int(bool(s[1]&infinity(k)))
    out=set()
    if (e==1 or h==0) and None not in pre:out.add(pair(*pre))
    if (e==0 or h==1) and None not in post:out.add(pair(*post))
    return frozenset(out)

def scores(E,n):return sorted(sum((a,b) in E for b in range(n)) for a in range(n))

def run():
    reports=[];total_paths=0;total_actions=0
    for k in range(2,7):
        child_states=list(product(states(k+1),repeat=2))
        coded={s:path(k+1,s) for s in child_states}
        need(len(set(coded.values()))==len(coded),'decoder injectivity')
        if k<=4:
            need(set(coded.values())==literal_paths(k+1),'independent DFS versus decoder')
        fibres=defaultdict(set)
        cuts=Counter()
        for s,p in coded.items():
            need(valid(k+1,p),'source path legal')
            q=delete_path(k,p);parent=project_state(k,s)
            need(q==path(k,parent),'actual deletion product')
            need(valid(k,q),'parent path legal')
            fibres[parent].add(s)
            orbit={act(k,s,g) for g in range(4)}
            stab={g for g in range(4) if act(k,s,g)==s}
            expected={0}
            if s[0]==infinity(k+1):expected|={1}
            if s[1]==infinity(k+1):expected|={2}
            if s[0]==s[1]==infinity(k+1):expected|={3}
            need(stab==expected,'exact stabilizer')
            need(len(orbit)*len(stab)==4,'orbit stabilizer')
            for g,h in product(range(4),repeat=2):
                need(act(k,act(k,s,g),h)==act(k,s,g^h),'literal group law')
                total_actions+=1
            swapped=s[::-1]
            need(toggle(k,swapped,1)==toggle(k,s,0)[::-1],'reflection conjugates axes')
            need(path(k+1,swapped)==tuple(2*(k+1)-v for v in reversed(p)),
                 'actual reflected reversal')
            need(blue(k+1,p)==(s[0]==s[1]),'actual blue diagonal')
            C=cut_sidecar(k,p)
            need(C==ports_direct(k,p),'independent cut-block scan')
            need(C==two_seam_sidecar(k,s),'two-seam reconstruction from parent ports and bits')
            need(len(C)<=2,'two ports suffice')
            lhs=delete_tournament(k,reverse_edges(tournament(k+1),edges(p)))
            rhs=reverse_edges(tournament(k),edges(q)^C)
            need(lhs==rhs,'exact metagraph correction')
            cuts[len(C)]+=1
        types=Counter()
        for parent,fibre in fibres.items():
            s=next(iter(fibre))
            need(fibre=={act(k,s,g) for g in range(4)},'fibre exactly one orbit')
            types[len(fibre)]+=1
        q=2**(k-1)+1
        need(types==Counter({4:(q-1)**2,2:2*(q-1),1:1}),'fibre stratum counts')
        need(len(coded)==4*q*q-4*q+1,'H recurrence')
        need(cuts==Counter({0:2*q,1:2*q*q-4,2:2*q*q-6*q+5}),
             'all-order cut-defect census formula')
        need(cut_sidecar(k,tuple(vertices(k+1)))==frozenset({(k-1,k),(k,k+1)}),
             'singleton standard-path port obstruction')
        reps=list(product(states(k+1),repeat=2))
        # Eight distinct permutations: every translation followed by optional swap.
        funcs={tuple((act(k,s,g)[::-1] if r else act(k,s,g)) for s in reps)
               for g,r in product(range(4),range(2))}
        need(len(funcs)==8,'faithful D8 action')
        need(any(toggle(k,s,0)[0]!=toggle(k,s,0)[1]
                 for s in reps if s[0]==s[1]),'toggle does not preserve blue')
        reports.append(dict(parent_k=k,source_order=2*k+3,paths=len(coded),
                            direct_DFS_checked=k<=4,
                            fibre_orbit_counts=dict(sorted(types.items())),
                            cut_edge_count_distribution=dict(sorted(cuts.items()))))
        total_paths+=len(coded)
    k=2
    p0=(4,6,3,0,1,2,5);p1=(2,4,6,3,0,1,5)
    need(p0==path(3,(0,1)) and p1==path(3,(2,1)),'hostile coordinates')
    need(not blue(3,p0) and not blue(3,p1),'hostile black paths')
    need(delete_path(k,p0)==delete_path(k,p1)==(4,2,0,1,3),'hostile parent')
    a=delete_tournament(k,reverse_edges(tournament(3),edges(p0)))
    b=delete_tournament(k,reverse_edges(tournament(3),edges(p1)))
    need(scores(a,5)==[0,1,2,3,4],'transitive hostile target')
    need(scores(b,5)==[0,2,2,3,3],'nontransitive hostile target')
    need(scores(a,5)!=sorted(4-d for d in scores(b,5)),'not converse isomorphic')
    output=dict(status='FINITE-EXACT independent controls for separately proved all-order statements',
                source_path_checks=total_paths,group_composition_checks=total_actions,
                reports=reports,
                metagraph_hostile=dict(parent_path=list(delete_path(2,p0)),
                    source_paths=[list(p0),list(p1)],
                    source_states=[[0,1],[2,1]],
                    cut_sidecars=[sorted(cut_sidecar(2,p0)),sorted(cut_sidecar(2,p1))],
                    deleted_target_score_spectra=[scores(a,5),scores(b,5)]))
    (HERE/'balanced_group_ports.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps(output,indent=2))

if __name__=='__main__':run()
