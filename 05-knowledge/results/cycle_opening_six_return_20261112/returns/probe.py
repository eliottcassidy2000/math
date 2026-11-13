"""Exploratory exact same-parent indegree-one port collisions."""
from pathlib import Path
from collections import Counter
from itertools import combinations
from math import comb
import json,gzip,sys
HERE=Path(__file__).resolve().parent
RESULTS=HERE.parents[1]
sys.path.insert(0,str(RESULTS/'chanpak_five_core_bridges_20261110/alternatives'))
import check as a

def profiles(P):
    h,end=a.dp(P)
    full=(1<<len(P))-1
    rev=tuple(full^r^(1<<i) for i,r in enumerate(P))
    _,starts=a.dp(rev)
    e2=[sum(end[mask][v]*h[full^mask] for mask in range(1,full+1)) for v in range(len(P))]
    return h[-1],starts[-1],e2,[h[-1]+e2[v]-starts[-1][v] for v in range(len(P))]

def main():
    records=[json.loads(line) for line in gzip.open(RESULTS/'metagraph_aspect_sequences_20261111/records.jsonl.gz','rt')]
    need=lambda x,m: None if x else (_ for _ in ()).throw(RuntimeError(m))
    need(len(records)==7412,'complete inherited1..8 bank')
    hits=[];counts=Counter();parents=Counter()
    for rec in records:
        P=tuple(rec['rows']);n=len(P)
        h,s,e2,v=profiles(P)
        need(h==rec['H'],'freshHP')
        some=False
        for i,j in combinations(range(n),2):
            if v[i]==v[j] and s[i]!=s[j]:
                hit=dict(n=n,parent=P,parent_H=h,ports=[i,j],H=v[i],defects=[s[i],s[j]],e2=[e2[i],e2[j]],strong_parent=rec['strong'])
                hits.append(hit);counts[n]+=1;some=True
        if some:parents[n]+=1
    out=dict(status='FINITE-EXACT exploratory',scope='all 7412 ordinary parent classes of orders 1..8, every unordered pair of distinct entry ports',counts=dict(counts),parents=dict(parents),hits=hits)
    (HERE/'port_collisions.json').write_text(json.dumps(out,indent=2)+'\n')
    print('COUNTS',dict(counts),dict(parents))
    for hit in hits[:15]:print(hit)
if __name__=='__main__':main()
