#!/usr/bin/env python3
"""Independent brute-permutation audit; no producer canonicalizer in controls.
Run: python independent_audit.py --receipt receipt.json
The producer predictor is called with independent complete parent atlases.
A hard runtime guard forbids its canonicalizer from seeing a child-order object.
"""
import argparse
from collections import Counter
from fractions import Fraction
from functools import lru_cache
import hashlib
import importlib.util
from importlib.machinery import SourceFileLoader
import sys
sys.dont_write_bytecode = True
from itertools import combinations, permutations
import json
from pathlib import Path
import time

HERE = Path(__file__).resolve().parent
COUNTS = Counter()

def check(ok, name):
    COUNTS[name] += 1
    if not ok:
        raise RuntimeError(name)

@lru_cache(None)
def orders(n):
    return tuple(permutations(range(n)))

@lru_cache(None)
def pairs(n):
    return tuple(combinations(range(n), 2))

def code(rows, order):
    answer = 0
    for a,b in pairs(len(rows)):
        answer = (answer << 1) | ((rows[order[a]] >> order[b]) & 1)
    return answer

@lru_cache(None)
def canon(rows):
    return min(code(rows,p) for p in orders(len(rows)))

@lru_cache(None)
def rooted(rows, root):
    other=tuple(v for v in range(len(rows)) if v != root)
    return min(code(rows,p+(root,)) for p in permutations(other))

def matrix(n, word):
    rows=[0]*n
    for bit,(a,b) in enumerate(pairs(n)):
        if (word >> bit)&1:
            rows[a] |= 1<<b
        else:
            rows[b] |= 1<<a
    return tuple(rows)

@lru_cache(None)
def labelled_atlas(n):
    found={}
    for word in range(1 << len(pairs(n))):
        t=matrix(n,word)
        found.setdefault(canon(t),t)
    return tuple(found.values())

def path_atlas(n):
    """Every fixed directed Hamilton path presentation, independent constructor."""
    free=tuple((a,b) for a,b in pairs(n) if b != a+1)
    found={}; fibre=Counter()
    for word in range(1 << len(free)):
        rows=[0]*n
        for a in range(n-1):
            rows[a] |= 1<<(a+1)
        for bit,(a,b) in enumerate(free):
            if (word>>bit)&1:
                rows[a] |= 1<<b
            else:
                rows[b] |= 1<<a
        t=tuple(rows)
        found.setdefault(canon(t),t)
        fibre[rooted(t,n-1)] += 1
    return found,fibre

def reverse(rows):
    return tuple(sum(1<<w for w in range(len(rows)) if v!=w and ((rows[w]>>v)&1))
                 for v in range(len(rows)))

def brute_details(rows):
    n=len(rows); base=code(rows,tuple(range(n)))
    aut=0; aut_fix=0; ends=[0]*n
    for p in orders(n):
        if code(rows,p)==base:
            aut+=1
            aut_fix+=(p[-1]==n-1)
        if all((rows[p[i]]>>p[i+1])&1 for i in range(n-1)):
            ends[p[-1]]+=1
    return aut,aut_fix,ends

def product(row, kernel):
    answer={}
    for q,a in row.items():
        for r,b in kernel[q].items():
            answer[r]=answer.get(r,Fraction(0))+a*b
    return {q:a for q,a in answer.items() if a}

def run(producer):
    spec=importlib.util.spec_from_loader('pinned_producer',SourceFileLoader('pinned_producer',str(producer)))
    mod=importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    real_canonical=mod.canonical
    bound=0; observed=[]
    def guarded(rows):
        check(len(rows)<=bound,'hard parent-order canonicalization guard')
        observed.append(len(rows))
        return real_canonical(rows)
    guarded.cache_clear=real_canonical.cache_clear
    mod.canonical=guarded
    reports=[]
    for N in range(1,7):
        bound=N-1
        parents=list(labelled_atlas(bound))
        reps,report,data=mod.predict(parents)
        states=data['states']; children=data['children']; D=data['D']; C=data['C']
        keys=[rooted(t,N-1) for t in children]
        key_to_state={key:p for p,key in enumerate(keys)}
        check(len(key_to_state)==len(states),'pointed states distinct under full rooted permutations')
        child_keys=[canon(t) for t in children]
        native,native_fibres=path_atlas(N)
        check(set(native)==set(child_keys),'independent native cube covers exactly transfer classes')
        if N<=5:
            check({canon(t) for t in labelled_atlas(N)}==set(child_keys),
                  'all labelled child tournaments covered')
        independent_sc=sum(canon(reverse(t))==key for key,t in native.items())
        check(report['classes']==len(native),'independent unpointed count')
        check(report['self_converse']==independent_sc,'independent self-converse count')
        check(report['merged_classes']==(len(native)+independent_sc)//2,'independent merged count')
        end_support={}; actual_lifts={}; expected_D=[]; expected_C=[]
        for p,t in enumerate(children):
            row=Counter(key_to_state[rooted(t,v)] for v in range(N))
            c=key_to_state[rooted(reverse(t),N-1)]
            check(row==D[p],'full rooted-permutation transfer row')
            check(c==C[p],'full rooted-permutation converse map')
            aut,h,ends=brute_details(t)
            check(data['stab'][p]==h,'independent pointed stabilizer')
            check(h*row[p]==aut,'independent child isotropy reconstruction')
            check(ends[-1]==h*data['lifts'][p],'independent native fibre L over h')
            check(native_fibres.get(keys[p],0)==data['lifts'][p],
                  'all normalized word fibre multiplicities')
            end_support[p]={v for v,count in enumerate(ends) if count}
            actual_lifts[p]=ends[-1]//h
            expected_D.append(row); expected_C.append(c)
            block={q for q,k in enumerate(child_keys) if k==child_keys[p]}
            check(set(row)==block,'transfer support equals independently classified child block')
        E=[{q:Fraction(m,N) for q,m in row.items()} for row in expected_D]
        Pplus=[]; Pminus=[]
        for p in range(len(states)):
            c=expected_C[p]
            plus={q:(E[p].get(q,0)+E[c].get(q,0))/2 for q in set(E[p])|set(E[c])}
            minus={q:(E[p].get(q,0)-E[c].get(q,0))/2 for q in set(E[p])|set(E[c])}
            Pplus.append({q:a for q,a in plus.items() if a})
            Pminus.append({q:a for q,a in minus.items() if a})
            check(data['stab'][p]==data['stab'][c], 'converse preserves pointed stabilizer')
        for p in range(len(states)):
            check(product(Pplus[p],Pplus)==Pplus[p], 'weighted plus projector')
            check(product(Pminus[p],Pminus)==Pminus[p], 'weighted minus projector')
            check(not product(Pplus[p],Pminus), 'plus minus orthogonality')
            check(not product(Pminus[p],Pplus), 'minus plus orthogonality')
            for q,a in E[p].items():
                check(a/data['stab'][p]==E[q].get(p,0)/data['stab'][q],
                      'weighted self-adjoint all-root projector')
        check(sum(Pplus[p].get(p,0) for p in range(len(states)))==(len(native)+independent_sc)//2,
              'weighted plus projector trace')
        check(sum(Pminus[p].get(p,0) for p in range(len(states)))==(len(native)-independent_sc)//2,
              'weighted minus projector trace')
        good={p for p,l in actual_lifts.items() if l}
        check(good==data['good'],'independent endpoint admissibility')
        R={}; K={}
        for p in good:
            E=end_support[p]; rev=reverse(children[p])
            re=brute_details(rev)[2]
            Erev={v for v,c in enumerate(re) if c}
            R[p]={q:Fraction(m,len(E)) for q,m in expected_D[p].items() if q in good}
            kc=Counter(key_to_state[rooted(rev,v)] for v in Erev)
            K[p]={q:Fraction(m,len(Erev)) for q,m in kc.items()}
            check(R[p]==data['R'][p],'independent endpoint kernel')
            check(K[p]==data['K'][p],'independent anti-endpoint kernel')
        for p in good:
            check(product(R[p],R)==R[p],'independent R squared')
            check(product(K[p],K)==R[p],'independent K squared')
            check(product(R[p],K)==K[p],'independent R K')
            check(product(K[p],R)==K[p],'independent K R')
        trace=sum(expected_D[p][p] for p in range(len(states)))
        delta=sum(expected_D[p][p]-1 for p in range(len(states)))
        check(trace==N*len(native),'independent full trace')
        check(delta==N*len(native)-len(states),'independent symmetry defect')
        check(sum(R[p].get(p,0) for p in good)==len(native),'independent endpoint trace')
        check(sum(K[p].get(p,0) for p in good)==independent_sc,'independent anti-endpoint trace')
        if N==2:
            check(any(C[p] not in good for p in good),'naive endpoint converse restriction is invalid')
        if N==3:
            check(sorted(D[p][p] for p in range(len(states)))==[1,1,1,3],
                  'minimal pointed-division hostile')
        rec={'order':N,'parent_classes':len(parents),'pointed_states':len(states),
             'endpoint_states':len(good),'classes':len(native),'self_converse':independent_sc,
             'trace_D':trace,'symmetry_defect':delta,'max_predictor_canonical_order':max(observed)}
        reports.append(rec); print(json.dumps(rec,sort_keys=True),flush=True)
    return {'producer_sha256':hashlib.sha256(producer.read_bytes()).hexdigest(),
            'scope':'Independent full-permutation parent atlases through order 5; rooted and unrooted child comparisons through 6; all labelled children through 5; all native words through 6.',
            'reports':reports,'checks':dict(sorted(COUNTS.items())), 'total_checks':sum(COUNTS.values()),
            'known_count_prefix_used_as_input':False,'child_canonicalization_in_predictor':False}

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--producer',type=Path,default=HERE/'source'/'transfer.py.txt')
    parser.add_argument('--receipt',type=Path,default=HERE/'independent-receipt.json')
    args=parser.parse_args(); start=time.monotonic()
    result=run(args.producer)
    result['elapsed_seconds']=time.monotonic()-start
    args.receipt.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n')
    print('PASS',result['total_checks'],flush=True)
