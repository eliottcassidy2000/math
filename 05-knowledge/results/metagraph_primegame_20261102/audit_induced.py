#!/usr/bin/env python3
"""Independent simple-path enumeration, with no subset-count recurrence."""
from pathlib import Path
import json,sys
HERE=Path(__file__).resolve().parent

def need(ok,msg):
    if not ok:raise RuntimeError(msg)

def all_counts(rows):
    h=[0]*(1<<len(rows));h[0]=1
    def visit(mask,last):
        h[mask]+=1
        more=rows[last]&~mask
        while more:
            bit=more&-more;more-=bit;visit(mask|bit,bit.bit_length()-1)
    for v in range(len(rows)):visit(1<<v,v)
    return h

def main():
    receipt=json.loads((HERE/'INDUCED_BRANCHES.json').read_text())
    branches={(a['p'],a['q']):a for a in receipt['branches']}
    source_counts={};found=set();relation={};checks=0;visited=0
    for line in (HERE/'low_sources.txt').read_text().splitlines():
        values=list(map(int,line.split()));n,H=values[:2];rows=values[2:2+n]
        h=all_counts(rows);need(h[-1]==H,'full source count');visited+=sum(h)-1
        source_counts[H]=source_counts.get(H,0)+1
        for mask in range(1,len(h)-1):
            checks+=1;key=(h[mask],H)
            relation.setdefault(key,dict(source_count=h[mask],target_count=H,source_subset=mask,target_rows=rows))
            if key in branches:found.add(key)
    need(checks==receipt['subset_checks'],'declared subset universe')
    for key,b in branches.items():
        need(b['target_strong_sources']==source_counts.get(b['q'],0),'complete target count')
        need(b['embeds']==(key in found),'all branch truth values')
        if b['embeds']:
            h=all_counts(b['target_rows'])
            need(h[-1]==b['q'] and h[b['source_subset']]==b['p'],'positive witness')
    noncomposable=[]
    for p,r in sorted(set(branches)-found):
        if p==7 or not source_counts.get(r):continue
        mids=[q for q in range(p+2,r,2) if (p,q) in relation and (q,r) in relation]
        if mids:
            q=min(mids);noncomposable.append(dict(p=p,q=q,r=r,first=relation[p,q],second=relation[q,r]))
    wider=[]
    for p in range(3,125,2):
        if p==7 or any(p%d==0 for d in range(2,int(p**0.5)+1)):continue
        for r in range(p+2,min(3*p,125),2):
            if (p,r) in relation or (p,r) in branches:continue
            mids=[q for q in range(p+2,r,2) if (p,q) in relation and (q,r) in relation]
            if mids:
                q=min(mids);wider.append(dict(p=p,q=q,r=r,first=relation[p,q],second=relation[q,r]))
    out=dict(status='PASS',source_checks=sum(source_counts.values()),subset_checks=checks,enumerated_directed_simple_paths=visited,
             negative_branches=[[p,q] for p,q in sorted(set(branches)-found)],noncomposable_witnesses=noncomposable,wider_noncomposable_witnesses=wider,
             scope='Independent literal path enumeration; all-order interpretation uses strong-component reduction and inherited f(10)=125')
    target=Path(sys.argv[1]) if len(sys.argv)>1 else HERE/'AUDIT_INDUCED.json'
    target.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n')
    print('PASS independent induced branches; noncomposable triples',[(a['p'],a['q'],a['r']) for a in noncomposable+wider])
if __name__=='__main__':main()
