#!/usr/bin/env python3
"""Input preparation only: unmerge the frozen n=8 atlas and authenticate it.
No order-nine tournament is built here. Completeness uses distinct canonical
n=8 classes and orbit-stabilizer labelled mass, not a target count.
"""
import csv, hashlib, itertools, json, math
from pathlib import Path
from collections import Counter
HERE=Path(__file__).resolve().parent
COUNTS=Counter()
def require(c,m):
    COUNTS[m]+=1
    if not c: raise RuntimeError(m)
def colors(rows):
    n=len(rows); c=tuple(r.bit_count() for r in rows)
    while True:
        labels=sorted(set(c)); masks=[sum(1<<v for v in range(n) if c[v]==x) for x in labels]
        sig=[(c[v],tuple((rows[v]&m).bit_count() for m in masks)) for v in range(n)]
        names={s:i for i,s in enumerate(sorted(set(sig)))}; new=tuple(names[s] for s in sig)
        if len(set(new))==len(set(c)): return new
        c=new
def permute(rows,p):
    return tuple(sum(((rows[v]>>w)&1)<<j for j,w in enumerate(p)) for v in p)
def canonical_actions(rows):
    require(len(rows)<=8,'canonical order at most eight')
    c=colors(rows); cells=[tuple(v for v in range(len(rows)) if c[v]==x) for x in sorted(set(c))]
    best=None; orders=[]
    for pieces in itertools.product(*(itertools.permutations(cell) for cell in cells)):
        p=sum(pieces,()); value=permute(rows,p)
        if best is None or value<best: best=value;orders=[p]
        elif value==best: orders.append(p)
    q=orders[0]; inv={v:i for i,v in enumerate(q)}
    actions=sorted(tuple(inv[v] for v in p) for p in orders)
    require(all(permute(best,g)==best for g in actions),'full canonical actions preserve arcs')
    return best,actions
def decode(code,n=8):
    rows=[0]*n; pairs=list(itertools.combinations(range(n),2))
    for bit,(v,w) in zip(f'{code:0{len(pairs)}b}',pairs):
        if bit=='1':rows[v]|=1<<w
        else:rows[w]|=1<<v
    return tuple(rows)
def main():
    source=HERE/'archive/merged_nodes.tsv.txt'; raw=source.read_bytes(); merged=list(csv.DictReader(raw.decode().splitlines(),delimiter='\t'))
    codes={}
    for i,record in enumerate(merged):
        require(int(record['merged_id'])==i,'consecutive merged ids')
        for suffix in 'ab':
            index=int(record['class_'+suffix]);code=int(record['canon_'+suffix])
            require(index not in codes or codes[index]==code,'consistent repeated ordinary class ids')
            codes[index]=code
    require(sorted(codes)==list(range(len(codes))),'consecutive ordinary class ids')
    records=[]; seen=set(); mass=0; aut_hist=Counter(); almost=[]
    for i in range(len(codes)):
        rows=decode(codes[i]); canonical,actions=canonical_actions(rows)
        require(canonical not in seen,'all ordinary parent classes nonisomorphic');seen.add(canonical)
        size=len(actions); require(math.factorial(8)%size==0,'integral labelled class orbit')
        mass+=math.factorial(8)//size;aut_hist[size]+=1
        record={'rows':list(canonical),'actions':[list(g) for g in actions]};records.append(record)
        if set(r.bit_count() for r in canonical)=={3,4}: almost.append(i)
    require(mass==1<<28,'full labelled order-eight orbit mass')
    # Every tournament has eight vertices and one orientation for each of 28
    # pairs. Distinct class-orbits exhaust all 2^28 labelled tournaments.
    records.sort(key=lambda r:r['rows'])
    payload={'schema':'tournament-parent-atlas-actions-v1','maximum_parent_order':8,'parents':records}
    path=HERE/'parents-eight.json';path.write_text(json.dumps(payload,sort_keys=True,separators=(',',':'))+'\n')
    receipt={'status':'PASS','raw_source_sha256':hashlib.sha256(raw).hexdigest(),'merged_classes':len(merged),'ordinary_classes':len(records),'automorphism_order_histogram':sorted(aut_hist.items()),'almost_regular_parent_classes':len(almost),'labelled_mass':mass,'complete_labelled_universe':1<<28,'atlas_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'checks':dict(COUNTS)}
    (HERE/'input-receipt.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print(json.dumps(receipt,sort_keys=True))
    with (HERE/'parents-eight.txt').open('w') as out:
        out.write(str(len(records))+'\n')
        for record in records:
            out.write(' '.join(map(str,record['rows']))+' '+str(len(record['actions']))+'\n')
            for action in record['actions']:out.write(' '.join(map(str,action))+'\n')
if __name__=='__main__':main()
