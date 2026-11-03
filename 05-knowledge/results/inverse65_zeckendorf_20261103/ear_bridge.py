"""Lift the actual code-8 ear family into the fixed source-word bank."""
from itertools import combinations,permutations
from pathlib import Path
from collections import Counter
import json
from bank import *
HERE=Path(__file__).resolve().parent
LOOKUP={tuple(ATLAS[l]['source_rows']):l for l in LABELS}
A=(15,17,19,23,25,27,29,33,37,41)
B=(2,3,5,7,11,13,17,19,23,27)

def rows(n,code):
    out=[0]*n
    for k,(i,j) in enumerate(combinations(range(n),2)):
        u,v=(i,j) if code>>k&1 else (j,i)
        out[u]|=1<<v
    return tuple(out)

def add(old,mask):
    n=len(old)
    return tuple(r|((1<<n) if not(mask>>i&1) else 0) for i,r in enumerate(old))+(mask,)

def paths(r):
    return [p for p in permutations(range(len(r))) if all(r[u]>>v&1 for u,v in zip(p,p[1:]))]

def relabel(r,p):return tuple(sum(1<<j for j,v in enumerate(p) if r[u]>>v&1) for u in p)

def complement(r):
    out=list(r)
    for i,j in combinations(range(len(r)),2):
        if j>i+1:out[i]^=1<<j;out[j]^=1<<i
    return tuple(out)

def parse_source(r,lookup=LOOKUP):
    n=len(r);reach=[row|1<<i for i,row in enumerate(r)]
    for k in range(n):
        for i in range(n):
            if reach[i]>>k&1:reach[i]|=reach[k]
    groups=[];seen=set()
    for i in range(n):
        if i in seen:continue
        part=[j for j in range(n) if reach[i]>>j&1 and reach[j]>>i&1]
        seen.update(part);groups.append(part)
    letters=[]
    for g in groups:
        need(g==list(range(min(g),max(g)+1)),'SCC intervals along source path')
        atom=relabel(r,g);label=lookup.get(atom)
        if label is None:return None
        letters.append(label)
    return tuple(reversed(letters))

def main():
    parent=rows(5,8);need(len(paths(parent))==9,'actual H9 parent')
    receipts=[];hist=Counter();fail=[]
    for cut in range(1,31):
        child=add(parent,cut);hp=paths(child);hist[len(hp)]+=1;charts=[]
        for p in hp:
            source=complement(relabel(child,p));word=parse_source(source)
            if word is None:continue
            need(vector(word)[0]==len(hp),'matrix H agrees with literal ear')
            need(actual_target(word)==relabel(child,tuple(reversed(p))),'full source-preserving relabel')
            charts.append({'path_in_original_labels':p,'processing_word':word,'state':vector(word),
                           'port_in_original_labels':p[0]})
        if not charts:fail.append({'cut':cut,'H':len(hp)})
        receipts.append({'cut':cut,'H':len(hp),'rows':child,'accepted_chart_count':len(charts),
                         'canonical_lift':charts[0] if charts else None})
    need(tuple(sorted(hist))==A,'inherited ear support')
    intersection={h:hist[h] for h in B if h in hist}
    need(intersection=={17:4,19:2,23:6,27:2},'coefficient-preserving filter')
    full_lookup={tuple(e['source_rows']):l for l,e in ATLAS.items()}
    candidates={}
    for failed in fail:
        cut=failed['cut'];child=add(parent,cut);poss={}
        for p in paths(child):
            word=parse_source(complement(relabel(child,p)),full_lookup)
            need(word is not None and len(word)==1,'failed lift is a strong six-atom')
            poss[word[0]]=p
        candidates[cut]=poss
    common=set.intersection(*(set(p) for p in candidates.values())) if candidates else set()
    chosen=[min(common,key=lambda l:int(l[1:]))] if common else [min(c,key=lambda l:int(l[1:])) for c in candidates.values()]
    repaired=[]
    for cut,choices in candidates.items():
        label=next(l for l in chosen if l in choices);p=choices[label]
        need(actual_target((label,))==relabel(add(parent,cut),tuple(reversed(p))),'extra atom literal repair')
        repaired.append({'cut':cut,'new_atom':label,'path_in_original_labels':p,'port_in_original_labels':p[0]})
    seed_polynomial=Counter()
    phase_outputs={}
    for r in receipts:
        if r['H'] not in B:continue
        c=r['canonical_lift'];need(c is not None,'all filtered cuts lift')
        v=tuple(c['state']);seed_polynomial[(r['H'],v)]+=1
        phase_outputs.setdefault(str(r['H']),[]).append(step(v,MATS['a'])[0])
    result={'status':'PASS','parent_rows':parent,'histogram':dict(sorted(hist.items())),
            'filter_B':B,'filtered_coefficients':intersection,'cuts':receipts,'failed_lifts':fail,
            'repair_candidate_counts':{str(k):len(v) for k,v in candidates.items()},
            'shared_repair_atom_count':len(common),'minimum_added_atoms':chosen,'literal_repairs':repaired,
            'filtered_seed_state_polynomial':[{'initial_H':h,'state':v,'multiplicity':c} for (h,v),c in sorted(seed_polynomial.items())],
            'filtered_next_singleton_responses':phase_outputs}
    (HERE/'ear_bridge.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':'PASS','histogram':result['histogram'],'failed_lifts':fail,
                      'accepted_charts':sum(r['accepted_chart_count'] for r in receipts),
                      'minimum_added_atoms':chosen,'shared_repair_atom_count':len(common),
                      'filtered_next_singleton_responses':phase_outputs}))
if __name__=='__main__':main()
