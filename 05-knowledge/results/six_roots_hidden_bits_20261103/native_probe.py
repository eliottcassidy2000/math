#!/usr/bin/env python3
"""All 65 E8 path complements, then all one-birth cuts; exact bounded language."""
import json
from pathlib import Path
from verify import en,profile
from probe67 import strong
P=Path(__file__).resolve().parent
def count(rows):
    n=len(rows);d=[[0]*n for _ in range(1<<n)]
    for v in range(n):d[1<<v][v]=1
    for s in range(1,1<<n):
        for v in range(n):
            if not d[s][v]:continue
            todo=rows[v]&~s
            while todo:
                b=todo&-todo;todo-=b;d[s|b][b.bit_length()-1]+=d[s][v]
    return sum(d[-1])
def decode(t):
    return list(range(8)) if t is None else [i+1 for i in range(6) if t>>i&1]+[7,0]+[i+1 for i in range(6) if not t>>i&1]
def main():
    base=en(8);targets={195,197,199};census=[];witnesses={str(h):[] for h in sorted(targets)};best={};tested=0
    for t in [None]+list(range(64)):
        path=decode(t);protected={frozenset((a,b)) for a,b in zip(path,path[1:])};rows=base.copy()
        for i in range(8):
            for j in range(i+1,8):
                if frozenset((i,j)) not in protected:rows[i]^=1<<j;rows[j]^=1<<i
        h=count(rows);hit_counts={str(z):0 for z in sorted(targets)}
        if h<=max(targets):
            for mask in range(256):
                born=rows+[0]
                for i in range(8):
                    if mask>>i&1:born[i]|=256
                    else:born[8]|=1<<i
                out=count(born);tested+=1
                if out in targets and strong(born):
                    hit_counts[str(out)]+=1
                    got,k=profile(born)
                    record={'address':t,'path':path,'complement_H':h,'complement_rows':rows,'birth_mask':mask,'target_rows':born,'retained_root':0,'retained_root_J':k[0]-got,'J_all':[v-got for v in k]}
                    if len(witnesses[str(out)])<1:witnesses[str(out)].append(record)
                    if str(out) not in best or (min(record['J_all']),record['retained_root_J'])<(min(best[str(out)]['J_all']),best[str(out)]['retained_root_J']):best[str(out)]=record
        census.append({'address':t,'complement_H':h,'strong_target_hits':hit_counts})
    report={'source_rows':base,'source_H':65,'status':'FINITE-EXACT exhaustive one-native-complement then one-birth language; no repeatable return claim','all_paths':65,'births_tested':tested,'pruning':'Skip H(complement)>199 by induced-extension monotonicity','census':census,'first_witnesses':witnesses,'best_port_witnesses':best,'total_strong_hits':{h:sum(c['strong_target_hits'][h] for c in census) for h in witnesses}}
    (P/'native_probe.json').write_text(json.dumps(report,indent=2)+'\n');print('Native/birth census:',report['total_strong_hits'],'births tested:',tested)
if __name__=='__main__':main()
