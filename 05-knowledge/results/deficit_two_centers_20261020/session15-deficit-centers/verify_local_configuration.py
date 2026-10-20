#!/usr/bin/env python3
"""Exhaustive, exact finite local lemma: 6-vertex strong => S4 + >=3 edge terminal."""
import itertools,json,hashlib
from collections import Counter
from pathlib import Path

def tournament(n,mask,fixed=None):
    a=[0]*n
    if fixed is not None:
        a[:len(fixed)]=fixed
    pairs=[(i,j) for i in range(n) for j in range(i+1,n)
           if fixed is None or j>=len(fixed)]
    for k,(i,j) in enumerate(pairs):
        u,v=(i,j) if mask>>k&1 else (j,i)
        a[u]|=1<<v
    return a

def strong_reach(a,S):
    S=set(S)
    for reverse in (False,True):
        seen={min(S)}
        while True:
            more={v for u in seen for v in S
                  if (a[v]>>u&1 if reverse else a[u]>>v&1)}
            if more<=seen:break
            seen|=more
        if seen!=S:return False
    return True

def strong_cuts(a,S):
    # Independent characterization: nonstrong iff a nontrivial dominating cut exists.
    S=list(S)
    for bits in range(1,(1<<len(S))-1):
        X=[S[i] for i in range(len(S)) if bits>>i&1]
        Y=[S[i] for i in range(len(S)) if not (bits>>i&1)]
        if all(a[x]>>y&1 for x in X for y in Y):return False
    return True

def witness(a,strong):
    for S in itertools.combinations(range(len(a)),4):
        if not strong(a,S):continue
        for v in range(len(a)):
            if v not in S and sum(a[u]>>v&1 for u in S)>=3:
                # Explicit directed four-cycle certificate, up to a chosen root.
                for tail in itertools.permutations(S[1:]):
                    C=(S[0],)+tail
                    if all(a[C[i]]>>C[(i+1)%4]&1 for i in range(4)):
                        ends=[u for u in C if a[u]>>v&1]
                        return dict(cycle=C,terminal=v,valid_ends=ends)
                raise AssertionError('strong S4 without cycle')
    return None

def run():
    results={}
    for mode in ('full_32768','fixed_s4_512'):
        bound=1<<15 if mode=='full_32768' else 1<<9
        fixed=None if mode=='full_32768' else [2,4,9,3]
        strong=strong_reach if fixed is None else strong_cuts
        rows=[];c=Counter()
        for mask in range(bound):
            a=tournament(6,mask,fixed)
            if not strong(a,range(6)):
                c['nonstrong']+=1;continue
            c['strong']+=1
            w=witness(a,strong)
            assert w is not None,(mode,mask,a)
            c[f"witness_indegree_{len(w['valid_ends'])}"]+=1
            rows.append(dict(mask=mask,**w))
        results[mode]=dict(counts=dict(c),witnesses=rows)
        print(mode,dict(c),flush=True)
    p=Path(__file__).with_name('local_configuration_verification.json')
    p.write_text(json.dumps(results,separators=(',',':'))+'\n')
    print('certificate_sha256',hashlib.sha256(p.read_bytes()).hexdigest())
if __name__=='__main__':run()
