#!/usr/bin/env python3
"""Exact, portable Session52 Fortunate / square-sum pattern audit.
Standard library only; all requirements remain active under python -O.
python3 check.py --write proposes the canonical results file; otherwise compare.
"""
import argparse
from collections import Counter
from functools import lru_cache
from itertools import combinations, groupby
import hashlib
import json
from math import gcd, isqrt, prod
from pathlib import Path

ROOT = Path(__file__).resolve().parent
PATTERN = [3, 2, 3, 5, 1, 1, 2]
EXPECTED_F = [3,5,7,13,23,17,19,23,37,61,67,61,71,47,107,59,61,109,89]
# Proper factors of the composite first wheel-surviving endpoint. Primality
# of these factors is unnecessary for their role as compositeness witnesses.
FACTORS = {4:13,5:23,9:127,10:107,11:433,12:1093,13:375569,15:941879,18:4129,19:92388407}

def require(ok, message):
    if not ok: raise ValueError(message)

def small_prime(n):
    return n >= 2 and all(n % d for d in range(2,isqrt(n)+1))

def mr_failure(n,a):
    if n<3 or n%2==0 or not 1<a<n:return False
    d,s=n-1,0
    while d%2==0:d//=2;s+=1
    x=pow(a,d,n)
    if x in (1,n-1):return False
    for _ in range(s-1):
        x=x*x%n
        if x==n-1:return False
    return True

def verify_arithmetic(doc):
    certs=doc['primality_certificates']; visiting=set(); verified=set()
    def prime(n):
        if n in verified:return
        require(n>=2,'invalid certified integer')
        require(n not in visiting,'cyclic prime certificate')
        require(str(n) in certs,'missing prime certificate')
        c=certs[str(n)]
        if n==2:
            require(c=={'base':1,'factors':[]},'bad certificate for 2')
            verified.add(2);return
        visiting.add(n);fs=c['factors'];a=c['base']
        require(bool(fs),'missing factorization')
        require(len({q for q,e in fs})==len(fs),'duplicate prime factor')
        require(all(isinstance(q,int) and isinstance(e,int) and 2<=q<n and e>=1 for q,e in fs),'bad factor domain')
        for q,e in fs:prime(q)
        require(prod(q**e for q,e in fs)==n-1,'factorization mismatch')
        require(1<a<n and pow(a,n-1,n)==1,'Lucas Fermat condition')
        require(all(gcd(pow(a,(n-1)//q,n)-1,n)==1 for q,e in fs),'Lucas order condition')
        visiting.remove(n);verified.add(n)
    def composite(n,w):
        require(not w.get('prime',False),'prime classified as a rejected offset')
        if 'divisor' in w:
            d=w['divisor'];require(1<d<n and n%d==0,'invalid proper divisor')
        elif 'mr_base' in w:require(mr_failure(n,w['mr_base']),'invalid Miller-Rabin failure')
        else:raise ValueError('unrecognized compositeness witness')
    ps=[p for p in range(2,74) if small_prime(p)]
    rows=[];P=1;endpoint_tests=0
    require(len(doc['rows'])==19,'wrong primorial universe')
    for n,r in enumerate(doc['rows'],1):
        p=ps[n-1];P*=p;q=ps[n]
        require((r['n'],r['p'],r['P'])==(n,p,P),'primorial indexing mismatch')
        F=r['plus']['offset'];ws=r['plus']['witnesses']
        require(len(ws)==F-1 and F==EXPECTED_F[n-1],'offset/minimality length mismatch')
        for d,w in enumerate(ws,2):
            endpoint_tests+=1
            if d<F:composite(P+d,w)
            else:
                require(w=={'prime':True},'terminal witness not prime');prime(P+d)
        # Wheel lower bound is also checked directly on every smaller d.
        require(F>=q,'wheel lower bound')
        for d in range(2,q):require(1<gcd(d,P)<P+d,'wheel rejection failure')
        Q=P+q;match=(F==q)
        if match:prime(Q);factor=None
        else:
            factor=FACTORS[n];require(1<factor<Q and Q%factor==0,'bad candidate factor')
        rows.append(dict(n=n,p_n=p,p_next=q,primorial=P,fortunate=F,
                         next_prime_gap=F-q,literal_prime_plus_one_gap=F-p-1,
                         candidate=Q,candidate_prime=match,
                         candidate_proper_factor=factor,
                         candidate_cofactor=Q//factor if factor else None))
    return rows,dict(positive_minimality_tests=endpoint_tests,used_lucas_primes=len(verified))

def graph(vertices):
    V=set(vertices)
    return {a:{b for b in V if a!=b and isqrt(a+b)**2==a+b} for a in V}

def components(G,removed=()):
    unseen=set(G)-set(removed);out=[]
    while unseen:
        todo=[min(unseen)];unseen.remove(todo[0]);C=[]
        while todo:
            v=todo.pop();C.append(v);new=G[v]&unseen;unseen-=new;todo+=sorted(new)
        out.append(sorted(C))
    return out

def paths(G):
    """Complete vertex-order DFS, quotienting reversal only."""
    if len(G)==1:return [(next(iter(G)),)]
    if len(components(G))!=1:return []
    leaves=[v for v in G if len(G[v])==1]
    starts=[min(leaves)] if leaves else sorted(G)
    out=[]
    def visit(path,unseen):
        if not unseen:
            if leaves or path[0]<path[-1]:out.append(tuple(path))
            return
        for w in sorted(G[path[-1]]&unseen):visit(path+[w],unseen-{w})
    for s in starts:visit([s],set(G)-{s})
    return out

def path_count_dp(G):
    """Independent bitmask dynamic program: all directed starts, divide by 2."""
    V=sorted(G);N=len(V);index={v:i for i,v in enumerate(V)}
    A=[sum(1<<index[w] for w in G[v]) for v in V];full=(1<<N)-1
    @lru_cache(None)
    def f(last,mask):
        if mask==full:return 1
        choices=A[last]&~mask;total=0
        while choices:
            bit=choices&-choices;choices-=bit
            total+=f(bit.bit_length()-1,mask|bit)
        return total
    if N==1:return 1
    total=sum(f(i,1<<i) for i in range(N))
    require(total%2==0,'odd reversal count')
    return total//2

def find_cover(G,k):
    """Find at most k disjoint spanning paths; only a positive-witness generator."""
    leaves=[v for v in G if len(G[v])<=1]
    starts=[min(leaves)] if leaves else sorted(G)
    def visit(pp,unseen):
        if not unseen:return pp
        for w in sorted(G[pp[-1][-1]]&unseen,key=lambda v:(len(G[v]&unseen),v)):
            out=visit(pp[:-1]+[pp[-1]+[w]],unseen-{w})
            if out:return out
        if len(pp)<k:
            for w in sorted(unseen,key=lambda v:(len(G[v]&unseen),v)):
                out=visit(pp+[[w]],unseen-{w})
                if out:return out
        return None
    for s in starts:
        out=visit([[s]],set(G)-{s})
        if out:return out
    return None

def verify_cover(G,cover):
    require(isinstance(cover,list) and all(p for p in cover),'empty path')
    entries=[v for p in cover for v in p]
    require(len(entries)==len(G) and set(entries)==set(G),'cover not a vertex partition')
    require(all(b in G[a] for p in cover for a,b in zip(p,p[1:])),'nonedge in cover')

def force_obstruction(G,s,t):
    """Sound endpoint-aware edge propagation, not assumed complete generally."""
    allowed={tuple(sorted((a,b))) for a in G for b in G[a]};forced=set()
    target={v:1 if v in (s,t) else 2 for v in G};steps=[]
    while True:
        old=(len(allowed),len(forced))
        for v in sorted(G):
            Av={e for e in allowed if v in e};Fv=Av&forced
            if len(Fv)>target[v] or len(Av)<target[v]:
                return dict(kind='degree',vertex=v,target=target[v],available=sorted(Av),forced=sorted(Fv),steps=steps)
            if len(Fv)==target[v]:
                deleted=Av-Fv
                if deleted:steps.append(dict(vertex=v,delete=sorted(deleted)));allowed-=deleted
            elif len(Av)==target[v]:
                added=Av-forced
                if added:steps.append(dict(vertex=v,force=sorted(added)));forced|=added
        H={v:set() for v in G}
        for a,b in forced:H[a].add(b);H[b].add(a)
        for C in components(H):
            if len(C)>1 and all(len(H[v])==2 for v in C):
                return dict(kind='forced_cycle',vertices=C,edges=sorted(e for e in forced if e[0] in C),steps=steps)
        if old==(len(allowed),len(forced)):return dict(kind='undecided',steps=steps)

def rle(seq):return [dict(value=value,length=sum(1 for _ in group)) for value,group in groupby(seq)]
def lengths(seq):return [r['length'] for r in rle(seq)]

def graph_audit():
    rows=[];all_paths={};negative={}
    for n in range(1,29):
        G=graph(range(1,n+1));C=components(G);P=paths(G)
        require(len(P)==path_count_dp(G),'DFS/DP count disagreement at '+str(n))
        cover=[];component_flags=[];component_counts=[]
        for cc in C:
            H={v:G[v]&set(cc) for v in cc};HP=paths(H)
            component_counts.append(len(HP));component_flags.append(bool(HP))
            if HP:part=[list(HP[0])]
            else:
                part=find_cover(H,2);require(part is not None,'two-path cover not found')
                # Independently reject every candidate endpoint pair, so the
                # zero-deficit impossibility does not rely on the path DFS.
                leaves=[v for v in sorted(H) if len(H[v])==1]
                ob=[]
                for s,t in combinations(sorted(H),2):
                    if not set(leaves)<={s,t}:continue
                    proof=force_obstruction(H,s,t)
                    require(proof['kind']!='undecided','forcing incomplete for claimed no-path component')
                    ob.append(dict(endpoints=[s,t],proof=proof))
                negative[f'{n}:{cc[0]}']=dict(leaves=leaves,
                    endpoint_rule='Every graph leaf must be a path endpoint.',
                    rejected_by_three_leaves=len(leaves)>2,
                    admissible_endpoint_obstructions=ob)
            verify_cover(H,part);cover+=part
        verify_cover(G,cover)
        old=graph(range(1,n));oldC=components(old)
        d=len(G[n]);touched=sum(bool(set(cc)&G[n]) for cc in oldC)
        E=sum(map(len,G.values()))//2;beta=E-n+len(C)
        require(d==isqrt(2*n-1)-isqrt(n),'birth formula')
        require(len(C)==len(oldC)+1-touched,'component recurrence')
        require(all(component_flags)==(len(cover)==len(C)),'path-cover deficit flag')
        row=dict(N=n,edges=E,components=len(C),component_vertex_sets=C,
                 component_hamilton_path_counts=component_counts,birth_neighbors=sorted(G[n]),
                 birth_degree=d,birth_touched_components=touched,cycle_rank=beta,
                 leaves=[v for v in sorted(G) if len(G[v])==1],
                 branch_vertices=[v for v in sorted(G) if len(G[v])>=3],
                 hamilton_path_count=len(P),global_traceable=bool(P),
                 componentwise_traceable=all(component_flags),minimum_path_cover=len(cover),
                 path_cover_deficit=len(cover)-len(C),path_cover=cover)
        rows.append(row);all_paths[str(n)]=[list(p) for p in P]
    return rows,all_paths,negative

def build():
    prime_rows,prime_checks=verify_arithmetic(json.loads((ROOT/'inputs/prime_offsets_certificates.json').read_text()))
    gr,all_paths,negative=graph_audit();F=[r['fortunate'] for r in prime_rows]
    ps=[r['p_n'] for r in prime_rows]+[71]
    flags=[r['candidate_prime'] for r in prime_rows]
    target=flags[:17];Gwindow=gr[9:26]
    require(lengths(target)==PATTERN,'prime run pattern')
    require(target==[r['componentwise_traceable'] for r in Gwindow],'componentwise finite-word identity')
    require(flags[17]!=gr[26]['componentwise_traceable'],'expected first extension failure absent')
    require([n+1 for n,b in enumerate(flags) if b]==[1,2,3,6,7,8,14,16,17],'equality-index audit')
    require([r['N'] for r in gr if r['global_traceable']]==[1,15,16,17,23,25,26,27,28],'small Hamilton support')
    require([r['N'] for r in gr if r['path_cover_deficit']]==[13,14,18,19,20,21,22,24],'path-cover-deficit support')
    require(all(r['path_cover_deficit'] in (0,1) for r in gr),'unexpected path-cover deficit')
    features={k:rle([r[k] for r in Gwindow]) for k in ['birth_degree','birth_touched_components','components','edges','cycle_rank','global_traceable','componentwise_traceable','minimum_path_cover','path_cover_deficit']}
    features.update({
        'leaf_count':rle([len(r['leaves']) for r in Gwindow]),
        'branch_vertex_count':rle([len(r['branch_vertices']) for r in Gwindow]),
        'forest':rle([r['cycle_rank']==0 for r in Gwindow]),
        'max_degree_at_most_two':rle([not r['branch_vertices'] for r in Gwindow]),
        'connected':rle([r['components']==1 for r in Gwindow]),
        'hamilton_cycle_exists':rle([False for r in Gwindow]),
    })
    # Every n<=28 has a leaf or is disconnected; hence no spanning cycle.
    require(all(r['components']>1 or r['leaves'] or r['N']<3 for r in gr),'cycle exclusion')
    cuts=[]
    for n,S in [(18,[7]),(19,[7]),(20,[5,7]),(21,[5,7]),(22,[5,7,14])]:
        C=components(graph(range(1,n+1)),S)
        require(len(C)>len(S)+1,'invalid cut certificate')
        cuts.append(dict(N=n,cut=S,remaining_components=C))
    birth_repairs=[]
    for n in [14,15,23,25]:
        require(gr[n-1]['minimum_path_cover']==gr[n-2]['minimum_path_cover']-1,'path-cover drop')
        old_cover=[]
        for path in gr[n-1]['path_cover']:
            if n not in path:old_cover.append(path);continue
            i=path.index(n)
            require(0<i<len(path)-1,'new vertex must be interior for path-cover drop')
            old_cover.extend([path[:i],path[i+1:]])
        verify_cover(graph(range(1,n)),old_cover)
        require(len(old_cover)==gr[n-2]['minimum_path_cover'],'nonminimum old cover')
        birth_repairs.append(dict(N=n,minimum_old_cover=old_cover,new_cover=gr[n-1]['path_cover']))
    for n in range(2,29):
        delta_pc=gr[n-1]['minimum_path_cover']-gr[n-2]['minimum_path_cover']
        require(-1<=delta_pc<=1,'path-cover Lipschitz bound')
        require(gr[n-1]['path_cover_deficit']-gr[n-2]['path_cover_deficit']==delta_pc-1+gr[n-1]['birth_touched_components'],'deficit update')
    alternatives={}
    for label,vertices in [('zero_to_N_minus_one',lambda n:range(n)),('zero_to_N',lambda n:range(n+1))]:
        vals=[]
        for n in range(10,27):
            G=graph(vertices(n));C=components(G)
            componentwise=all(find_cover({v:G[v]&set(cc) for v in cc},1) is not None for cc in C)
            vals.append(dict(N=n,global_traceable=(len(C)==1 and componentwise),componentwise_traceable=componentwise))
        alternatives[label]=dict(rows=vals,global_runs=rle([r['global_traceable'] for r in vals]),componentwise_runs=rle([r['componentwise_traceable'] for r in vals]))
    sorted_F=sorted(set(x for x in F if x<=71))
    equal_indices=[r['n'] for r in prime_rows if r['candidate_prime']]
    return dict(schema=1,status='FINITE-EXACT; not a universal arithmetic/graph bridge',
                prime_rows=prime_rows,prime_verification=prime_checks,graph_rows=gr,
                hamilton_paths=all_paths,no_path_propagation_certificates=negative,
                cut_certificates=cuts,path_cover_birth_repairs=birth_repairs,
                exact_word=target,run_lengths=PATTERN,
                prime_runs=rle(target),graph_componentwise_runs=rle([r['componentwise_traceable'] for r in Gwindow]),
                graph_global_runs=rle([r['global_traceable'] for r in Gwindow]),
                graph_full_componentwise_runs_through_28=rle([r['componentwise_traceable'] for r in gr]),
                graph_feature_runs_N10_to_26=features,
                shifted_graph_windows={str(start):{'global_runs':rle([r['global_traceable'] for r in gr[start-1:start+16]]),'componentwise_runs':rle([r['componentwise_traceable'] for r in gr[start-1:start+16]])} for start in [9,10,11]},
                alternative_vertex_conventions=alternatives,
                failed_prime_interpretations=dict(
                    literal_F_minus_prime_plus_one=[r['literal_prime_plus_one_gap'] for r in prime_rows],
                    literal_equality_runs=rle([r['literal_prime_plus_one_gap']==0 for r in prime_rows[:17]]),
                    raw_F_minus_next_prime=[r['next_prime_gap'] for r in prime_rows],
                    equality_index_gaps=[b-a for a,b in zip(equal_indices,equal_indices[1:])],
                    sorted_distinct_F_up_to_71=sorted_F,
                    sorted_distinct_minus_next_prime=[f-ps[i+1] for i,f in enumerate(sorted_F)]),
                first_extension_failure=dict(n=18,N=27,F=109,p_next=67,prime_gap=42,
                    candidate=117288381359406970983337,factor=4129,cofactor=28406001782370300553,
                    graph_hamilton_path_count=35,graph_path=all_paths['27'][0]))

def canonical(obj):return json.dumps(obj,indent=2,sort_keys=True)+'\n'

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');args=ap.parse_args()
    out=build();body=canonical(out);target=ROOT/'results.json'
    if args.write:target.write_text(body)
    else:require(target.read_text()==body,'results.json differs from recomputation')
    print(json.dumps(dict(status='PASS',run_lengths=out['run_lengths'],
        prime_equality_indices=[r['n'] for r in out['prime_rows'] if r['candidate_prime']],
        ordinary_graph_runs=out['graph_global_runs'],
        first_extension_failure=out['first_extension_failure'],
        results_sha256=hashlib.sha256(body.encode()).hexdigest()),indent=2,sort_keys=True))

if __name__=='__main__':main()
