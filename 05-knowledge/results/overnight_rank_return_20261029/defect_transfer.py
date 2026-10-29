"""Exact named path-cover queries using ordered backward-edge frontiers."""
from collections import defaultdict
from itertools import permutations
from random import Random
from math import factorial
import json

from complement_kernel import HERE, need, complement, paths
from run_recurrence import tribonacci_rows, joint_entry
from three_run import native


def defects(rows):
    return [(u,v) for v,row in enumerate(rows) for u in range(v) if row>>u&1]


def defect_width(n,edges):
    delta=[0]*(n+1)
    for u,v in edges: delta[u+1]+=1;delta[v+1]-=1
    sizes=[];current=0
    for t in range(n+1): current+=delta[t];sizes.append(current)
    return max(sizes,default=0),sizes


def query(n,edges,starts,ends,erase_pairing=False):
    """B_up xor edges; ports are named, paired and distinct within each list.

    Endpoint -1 is an anonymous open tail; -(i+2) is closed port i.
    Nonnegative endpoints are actual vertices still meeting future defects.
    """
    k=len(starts)
    need(len(ends)==k and len(set(starts))==len(set(ends))==k,'distinct ordered ports')
    need(all(0<=v<n for v in tuple(starts)+tuple(ends)), 'valid ports')
    edges=set(edges)
    need(all(0<=u<v<n for u,v in edges),'ordered defect edges')
    start_index={v:i for i,v in enumerate(starts)}
    end_index={v:i for i,v in enumerate(ends)}
    if any(start_index[v]!=end_index[v] for v in start_index.keys()&end_index.keys()):
        return 0,dict(peak=0,transitions=0,width=defect_width(n,edges)[0])
    last=[-1]*n
    backward=[set() for _ in range(n)]
    for u,v in edges: last[u]=max(last[u],v);backward[v].add(u)
    states={(0,()):1}
    peak=1;transitions=0
    for v in range(n):
        following=defaultdict(int)
        for (done,components),weight in states.items():
            incoming=[None] if v in start_index else [None]+[
                i for i,(head,tail) in enumerate(components)
                if tail==-1 or tail>=0 and tail not in backward[v]]
            outgoing=[None] if v in end_index else [None]+[
                i for i,(head,tail) in enumerate(components)
                if head>=0 and head in backward[v]]
            for left in incoming:
                for right in outgoing:
                    if left is not None and left==right: continue
                    head=components[left][0] if left is not None else -(start_index[v]+2) if v in start_index else v
                    tail=components[right][1] if right is not None else -(end_index[v]+2) if v in end_index else v
                    joined=[item for i,item in enumerate(components) if i!=left and i!=right]
                    new_done=done
                    if head<0 and tail<-1:
                        if head!=tail or done>>( -head-2)&1: continue
                        new_done|=1<<(-head-2)
                    else: joined.append((head,tail))
                    normalized=[]
                    valid=True
                    for a,b in joined:
                        if a>=0 and last[a]<=v:
                            valid=False;break
                        if b>=0 and last[b]<=v: b=-1
                        normalized.append((a,b))
                    if not valid: continue
                    if erase_pairing:
                        # Deliberately hostile quotient: preserve endpoint bags only.
                        normalized=list(zip(sorted(a for a,b in normalized),sorted(b for a,b in normalized)))
                    key=(new_done,tuple(sorted(normalized)))
                    following[key]+=weight;transitions+=1
        states=following;peak=max(peak,len(states))
    return states.get(((1<<k)-1,()),0),dict(peak=peak,transitions=transitions,width=defect_width(n,edges)[0])


def interval_target(lengths):
    need(all(a>=1 for a in lengths) and all(a!=1 or b!=1 for a,b in zip(lengths,lengths[1:])),
         'legal run composition')
    chunks=[];left=0
    for a in lengths:chunks.append(tuple(range(left,left+a)));left+=a
    path=tuple(v for block in reversed(chunks) for v in block)
    edges=set()
    for low,high in zip(chunks,chunks[1:]):
        edges.add((low[-1],high[0]));edges.add((low[0],high[-1]))
    return left,path,edges


def subset_query(table,starts,ends):
    """Independent subset convolution, with full vertex masks and no frontier."""
    full=len(table)-1
    states={0:1}
    for s,t in zip(starts,ends):
        following=defaultdict(int)
        for used,weight in states.items():
            remaining=full^used
            subset=remaining
            while subset:
                value=table[subset][s][t]
                if value:following[used|subset]+=weight*value
                subset=(subset-1)&remaining
        states=following
    return states.get(full,0)


def run():
    from two_run import subset_paths, two_runs
    from atom_precision import hostile
    from complement_verify import strong
    rng=Random(2026102912)
    checks=[];target_checks=[]
    records=json.loads((HERE/'strong_target.json').read_text())
    for number,record in enumerate(records,1):
        sources=[record['source']]
        if 'certificate' in record:sources.append(complement(sources[0],record['certificate']['path']))
        for rows in sources:
            n=len(rows);_,ports,M=native(rows)
            chosen=list(range(len(ports)**2)) if n<=4 else [rng.randrange(len(ports)**2) for _ in range(24)]
            for q in chosen:
                i,j=divmod(q,len(ports));value,stats=query(n,defects(rows),ports[i],ports[j])
                need(value==M[i][j],'all strong classes and their selected targets: native query')
            checks.append(dict(n=n,source=rows,queries=len(chosen)))
        if number%100==0:print('source queries',number,flush=True)
    for n in range(3,10):
        S=tribonacci_rows(n)
        for path in paths(S):
            U=complement(S,path);D=defects(U)
            need(defect_width(n,D)[0]<=2,'every interval-path target has width at most two')
            # Every three-port entry through five; independent mask queries thereafter.
            if n<=5:
                _,ports,M=native(U)
                for i,s in enumerate(ports):
                    for j,t in enumerate(ports):
                        value,_=query(n,D,s,t)
                        need(value==M[i][j],'complete small target matrices')
                count=len(ports)**2
            else:
                table=subset_paths(U);count=12
                for _ in range(count):
                    s=tuple(rng.sample(range(n),3));t=tuple(rng.sample(range(n),3))
                    value,_=query(n,D,s,t)
                    need(value==subset_query(table,s,t),'all interval complements: independent subset query')
            target_checks.append(dict(n=n,path=path,width=defect_width(n,D)[0],queries=count))
        print('interval complements',n,flush=True)
    controls=[]
    for flip in (False,True):
        rows=hostile(flip);ports,M,_=two_runs(rows)
        selected=[(i,j) for i in range(len(ports)) for j in range(len(ports)) if M[i][j]]
        rng.shuffle(selected)
        values=[]
        for i,j in selected[:64]:
            value,stats=query(12,defects(rows),ports[i],ports[j])
            need(value==M[i][j],'81/27 ordered two-port control')
            values.append(dict(starts=ports[i],ends=ports[j],value=value,stats=stats))
        from complement_kernel import rank_mod
        rank=rank_mod(M,3);mass=sum(map(sum,M))
        need(rank==(22 if flip else 24) and mass==9011898,'81/27 inherited separator')
        controls.append(dict(flip=flip,mass=mass,rank_mod3=rank,queries=values))
    long=[]
    for n in (20,100,1000):
        n,path,D=interval_target([2]*(n//2))
        s=(0,1,2);t=(n-3,n-2,n-1)
        value,stats=query(n,D,s,t)
        need(value>0 and stats['width']==2,'long actual-complement query')
        long.append(dict(n=n,composition=[2]*(n//2),starts=s,ends=t,value=value,stats=stats))
    previous=json.loads((HERE/'run_recurrence.json').read_text())
    interval_sources=[]
    for old in previous['long_queries']:
        n=old['n'];D={(i,i+1) for i in range(n-1)}
        s=tuple(n-1-v for v in old['starts']);t=tuple(n-1-v for v in old['ends'])
        value,stats=query(n,D,s,t)
        need(value==old['value'],'independent interval transfer in reversed vertex order')
        interval_sources.append(dict(n=n,value=value,stats=stats))
    pairing=dict(source=[12,29,24,16,1],path=[1,2,3,4,0],starts=[1,0,2],ends=[3,4,2])
    value,_=query(5,defects(pairing['source']),pairing['starts'],pairing['ends'])
    wrong,_=query(5,defects(pairing['source']),pairing['starts'],pairing['ends'],erase_pairing=True)
    need(value==0 and wrong==1,'component-pairing erasure invents a cover')
    need(complement(tribonacci_rows(5),pairing['path'])==pairing['source'],'hostile actual complement')
    pairing.update(exact=value,erased=wrong)
    cycle_source=[10,12,9,0]
    need(query(4,defects(cycle_source),(3,),(3,))[0]==0,
         'a closed C3 plus a singleton meets degree guards but is not a spanning path')
    paley=[sum(1<<((i+r)%7) for r in (1,2,4)) for i in range(7)]
    _,ports,M=native(paley);ix={p:i for i,p in enumerate(ports)}
    phase=lambda p:tuple(2*v%7 for v in p)
    need(all(M[i][j]==M[ix[phase(s)]][ix[phase(t)]]
             for i,s in enumerate(ports) for j,t in enumerate(ports)), 'simultaneous Frobenius transport')
    phase_hostile=None
    for i,s in enumerate(ports):
        for j,t in enumerate(ports):
            wrong=M[ix[phase(s)]][j]
            if M[i][j] and wrong and M[i][j]!=wrong:
                phase_hostile=dict(starts=s,ends=t,value=M[i][j],one_sided=wrong,transported_starts=phase(s))
                need(query(7,defects(paley),s,t)[0]==M[i][j], 'frontier keeps common phase')
                break
        if phase_hostile:break
    need(phase_hostile is not None,'nontrivial one-sided phase hostile')
    second=dict(n=4,first_path=[2,3,0,1],second_path=[0,2,1,3])
    U=complement(tribonacci_rows(4),second['first_path']);V=complement(U,second['second_path'])
    D=[(u,v) for u in range(4) for v in range(u+1,4) if V[u]>>v&1]
    need(strong(V) and defect_width(4,D)[0]==3,'second complement escapes width two in retained frame')
    second.update(first=U,second=V,defects=D,width=3)
    splice=[]
    for k in range(1,8):
        accepted=[]
        for pi in permutations(range(k)):
            seen={0};v=(pi[0]+1)%k
            while v not in seen:seen.add(v);v=(pi[v]+1)%k
            if len(seen)==k:accepted.append(pi)
        need(len(accepted)==factorial(k-1),'single-cycle splice retains factorially many pairings')
        splice.append(dict(k=k,accepted_count=len(accepted),examples=accepted[:3]))
    return dict(strong_source_queries=checks,interval_target_queries=target_checks,
                controls=controls,long_targets=long,interval_source_checks=interval_sources,
                pairing_hostile=pairing,cycle_hostile=dict(source=cycle_source,port=3,actual=0,degree_only=1),
                phase_hostile=phase_hostile,second_complement=second,splice_counts=splice)


def quick():
    rng=Random(2026102912)
    tests=0
    for n in range(3,8):
        sources=[tribonacci_rows(n)]
        for _ in range(6):
            rows=[0]*n
            for u in range(n):
                for v in range(u+1,n):
                    if rng.randrange(2):rows[u]|=1<<v
                    else:rows[v]|=1<<u
            sources.append(rows)
        for rows in sources:
            _,ports,M=native(rows)
            for _ in range(40):
                i,j=rng.randrange(len(ports)),rng.randrange(len(ports))
                value,_=query(n,defects(rows),ports[i],ports[j])
                need(value==M[i][j],f'native query n={n} rows={rows} ports={ports[i]},{ports[j]} {value}!={M[i][j]}')
                tests+=1
    print('PASS quick',tests,flush=True)


if __name__=='__main__':
    import sys
    if '--quick' in sys.argv:quick()
    else:
        receipt=json.loads(json.dumps(run()))
        path=HERE/'defect_transfer.json'
        if '--write' in sys.argv:path.write_text(json.dumps(receipt,separators=(',',':'))+'\n')
        else:need(receipt==json.loads(path.read_text()),'unchanged exact receipt')
        print('PASS bounded-defect joint-port transfer',flush=True)
