"""Independent finite audit. No imports from the producer or its core.

Uses all vertex orders, all internal root-fixing permutations, alignment
enumeration, and subcube truth tables instead of the producer's DFS.
"""
from collections import Counter
from functools import lru_cache
from itertools import permutations, combinations
from pathlib import Path
import json

P=Path(__file__).resolve().parent
def check(b,msg):
    if not b:raise RuntimeError(msg)
def arc(t,i,j):return bool(t[i]&(1<<j))
def matrix(t):return [[arc(t,i,j) for j in range(len(t))] for i in range(len(t))]
@lru_cache(None)
def orders(n):return tuple(permutations(range(n)))
def is_path(t,p):return all(arc(t,a,b) for a,b in zip(p,p[1:]))
def target(t,p):
    kept=set(zip(p,p[1:]));n=len(t)
    return tuple(sum((arc(t,i,j) if (i,j) in kept or (j,i) in kept else not arc(t,i,j))<<j
                     for j in range(n) if j!=i) for i in range(n))
def symmetric(t,p):
    return all(arc(t,p[i],p[j])!=arc(t,p[-1-i],p[-1-j]) for i,j in combinations(range(len(t)),2))
@lru_cache(None)
def rooted(t,s,e):
    if len(t)==1:return 0
    pairs=list(combinations(range(len(t)),2));rest=[v for v in range(len(t)) if v not in [s,e]]
    best=None
    for mid in permutations(rest):
        p=(s,e)+mid;v=sum(arc(t,p[a],p[b])<<i for i,(a,b) in enumerate(pairs))
        if best is None or v<best:best=v
    return best

def aut_count(t,roots=()):
    return sum(all(p[r]==r for r in roots) and all(arc(t,i,j)==arc(t,p[i],p[j]) for i,j in combinations(range(len(t)),2)) for p in orders(len(t)))

def cover_alignment(a,b,roots,f):
    n=len(a);out=[[] for _ in a];ins=[0]*n
    for i,j in combinations(range(n),2):
        if arc(a,i,j)==arc(b,f[i],f[j]):
            u,v=(i,j) if arc(a,i,j) else (j,i);out[u].append(v);ins[v]+=1
            if len(out[u])>1 or ins[v]>1:return None
    inv=[0]*n
    for i,j in enumerate(f):inv[j]=i
    seen=set();whole=[]
    for r in range(0,len(roots),2):
        start,end=inv[roots[r]],inv[roots[r+1]]
        if ins[start]:return None
        walk=[];v=start
        while True:
            if v in seen:return None
            seen.add(v);walk.append(v)
            if not out[v]:break
            v=out[v][0]
        if v!=end:return None
        whole.extend(walk)
    return tuple(whole) if len(seen)==n else None

def guard_truth_table(n):
    free=list((i,j) for i in range(n) for j in range(i+2,n));lookup={e:k for k,e in enumerate(free)}
    gs=[0]*(1<<len(free));full=(1<<len(free))-1
    for pid,p in enumerate(orders(n)):
        need=0;fixed=0;possible=True
        for a,b in zip(p,p[1:]):
            i,j=sorted((a,b))
            if j==i+1:
                if a> b:possible=False;break
            else:
                bit=1<<lookup[i,j];fixed|=bit
                if a<b:need|=bit
        if not possible:continue
        optional=full^fixed;sub=optional;flag=1<<pid
        while True:
            gs[need|sub]|=flag
            if not sub:break
            sub=(sub-1)&optional
    return len(set(gs)),len(gs)

def main():
    atoms=json.loads((P/'strong-atlas.json').read_text(encoding='utf-8'))
    runs=json.loads((P/'run-atlas.json').read_text(encoding='utf-8'))
    summary=[];path_count=0;response_count=0
    for n in range(1,8):
        same=[a for a in atoms if a['n']==n]
        for a in same:
            t=tuple(a['source']);aut=aut_count(t);check(aut==a['Aut'],'source automorphisms')
            f=Counter();b=Counter();H=0
            for p in orders(n):
                if not is_path(t,p):continue
                H+=1;u=target(t,p);h=rooted(u,p[0],p[-1]);f[h]+=1;b[h]+=symmetric(t,p)
            check(H==a['H'],'independent H')
            expected_f=Counter();expected_b=Counter()
            for r in a['response']:
                key=rooted(tuple(r['rows']),r['roots'][0],r['roots'][1])
                expected_f[key]+=r['total']*aut;expected_b[key]+=r['blue']*aut
            check(f==expected_f and b==expected_b,'full strong joint-port response')
            path_count+=H;response_count+=len(f)
        cells,states=guard_truth_table(n)
        check(cells==states,'complete legal-guard truth table')
        summary.append({'n':n,'strong_classes':len(same),'guard_cells':cells,'states':states})
    # Every entry of the all-arity atlas checked by the independent agreement
    # forest criterion, with the automorphism denominator retained.
    cover_states=0;alignment_count=0
    for a in runs:
        t=tuple(a['source']);aut=aut_count(t)
        for level in a['responses']:
            total=0
            for h in level['states']:
                u=tuple(h['rows']);r=tuple(h['roots']);stabilizer=aut_count(u,tuple(set(r)))
                N=sum(cover_alignment(t,u,r,p) is not None for p in orders(len(t)))
                check(N==aut*stabilizer*h['multiplicity'],'all-arity double-coset coefficient')
                total+=h['multiplicity'];cover_states+=1;alignment_count+=N
            check(total==level['orbits'],'ordered-cover total')
    rec=json.loads((P/'atlas-normal.json').read_text(encoding='utf-8-sig'))
    w=rec['first_paired_port_guard_failure'];x,y=map(tuple,w['rows']);p=tuple(w['permutation']);n=len(x)
    check(rooted(x,0,n-1)==rooted(y,0,n-1),'hostile source roots')
    u=target(x,tuple(range(n)));v=target(y,tuple(range(n)))
    check(rooted(u,0,n-1)==rooted(v,0,n-1),'hostile target roots')
    check([is_path(x,p),is_path(y,p)]==w['legal'],'hostile guard mismatch')
    # All full permutations also certify the six-vertex target counts.
    r=rec['six_vertex_run_hostile'];t=tuple(r['source'])
    counts=[]
    for p,rows in zip(r['paths'],r['targets']):
        check(is_path(t,p) and target(t,p)==tuple(rows),'run target reconstruction')
        counts.append(sum(is_path(rows,q) for q in orders(len(rows))))
    check(counts==r['target_H']==[41,45],'six-vertex independent count')
    con=json.loads((P/'controls-normal.json').read_text(encoding='utf-8'))
    r=con['two_outer_ports_failure'];t=tuple(r['source']);outer=[];inner_counts=[]
    for z in r['records']:
        p=tuple(z['path']);u=target(t,p)
        check(is_path(t,p) and u==tuple(z['target']),'strong child run complement')
        inner_counts.append(sum(is_path(u,q) for q in orders(6)))
        outer.append(rooted(tuple(z['child_target']),*z['outer_roots']))
    check(outer[0]==outer[1] and inner_counts==[15,23],'same child two-root response, different assembly')
    for key,a,b in [('source_rooted_isomorphism',x,y),('target_rooted_isomorphism',target(x,tuple(range(7))),target(y,tuple(range(7))))]:
        f=con['alignment_amplification'][key]
        check(all(arc(a,i,j)==arc(b,f[i],f[j]) for i,j in combinations(range(7),2)),'amplification seed isomorphism')
    print(json.dumps({'full_strong_response_audit':summary,'raw_Hamilton_paths_checked':path_count,
                      'joint_response_states_checked':response_count,'all_arity_states_checked':cover_states,
                      'admissible_alignments_checked':alignment_count,'seven_vertex_guard_hostile':True,
                      'six_vertex_target_H':counts,'revisited_child_target_H':inner_counts,
                      'amplification_seed_alignment':True},sort_keys=True,indent=2))

if __name__=='__main__':main()
