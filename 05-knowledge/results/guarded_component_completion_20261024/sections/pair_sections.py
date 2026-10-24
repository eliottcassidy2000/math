"""Exact every-predicate census on each actual53helper pair at occupancy4."""
import json,sys
from pathlib import Path
from itertools import combinations
from collections import Counter
from hashlib import sha256
HERE=Path(__file__).resolve().parent;RESULTS=HERE.parents[1]
SUPPLY=RESULTS/'compatible_charged_helpers_20261023/supply/selected_supply.json'
def load():
    cert=json.loads(SUPPLY.read_text());bank=next(b for b in cert['banks']if len(b['source_ids'])==53)
    src=bank['source_ids'];by={r['source']:r for r in cert['helpers']};labels=[by[s]['label']for s in src]
    return cert,src,labels
def charge(labels,ids):
    x=0;y=0
    for i in ids:x+=labels[i][0];y^=sum(v<<j for j,v in enumerate(labels[i][1:]))
    return((x%128)<<4)|y
def table(labels):
    witnesses=[[]for _ in range(2048)];bits=[[0]*len(labels)for _ in range(2048)]
    for ids in combinations(range(len(labels)),4):
        q=charge(labels,ids);b=1<<len(witnesses[q]);m=sum(1<<i for i in ids)
        witnesses[q].append(m)
        for i in ids:bits[q][i]|=b
    return witnesses,bits
def main():
    cert,src,labels=load();ww,bits=table(labels);rows=[];available_hist=Counter();relation_safe=Counter();mixed=[]
    for i,j in combinations(range(53),2):
        hist=Counter();missing={r:[]for r in range(1,16)}
        for q in range(2048):
            full=(1<<len(ww[q]))-1;a,b=bits[q][i],bits[q][j]
            patterns=[full&~(a|b),b&~a,a&~b,a&b]
            available=sum(1<<k for k,v in enumerate(patterns)if v)
            hist[available]+=1;available_hist[available]+=1
            for r in range(1,16):
                if not r&available:missing[r].append(q)
        safe=[r for r in range(1,16)if not missing[r]]
        for r in safe:relation_safe[r]+=1
        for r in (2,4):
            if missing[r]:mixed.append(dict(pair=[src[i],src[j]],allowed_pattern='01'if r==2 else'10',missing=missing[r]))
        rows.append(dict(sources=[src[i],src[j]],available_pattern_histogram=dict(sorted(hist.items())),universally_safe_relation_masks=safe,
                         singleton11_covered=2048-len(missing[8])))
    out=dict(status='FINITE-EXACT full native53/k4 pair-predicate coverage',source_hash=sha256(SUPPLY.read_bytes().replace(b'\r\n',b'\n')).hexdigest(),
        sources=src,occupancy=4,subset_count=sum(map(len,ww)),pair_count=len(rows),
        pattern_bits={'00':1,'01':2,'10':4,'11':8},safe_pair_counts_by_relation=dict(sorted(relation_safe.items())),
        global_available_pattern_histogram=dict(sorted(available_hist.items())),mixed_singleton_failures=mixed,pair_rows=rows)
    out=json.loads(json.dumps(out));p=HERE/'pair_sections_receipt.json'
    if '--write'in sys.argv:p.write_text(json.dumps(out,indent=2))
    elif json.loads(p.read_text())!=out:raise RuntimeError('pair-predicate receipt mismatch')
    print('PASS',out['subset_count'],'quadruples;',len(rows),'pairs; safe masks',dict(relation_safe));print('available',dict(available_hist),'mixedfails',len(mixed))
if __name__=='__main__':main()
