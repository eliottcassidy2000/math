"""Independent unbounded-integer audit of every contained subset table."""
from collections import Counter
from pathlib import Path
import gzip
import hashlib
import json
import sys

HERE=Path(__file__).resolve().parent
RESULTS=HERE.parent.parent

def need(ok, msg):
    if not ok:
        raise RuntimeError(msg)

def all_counts(rows):
    n=len(rows)
    incoming=[sum(1<<u for u in range(n) if (rows[u]>>v)&1) for v in range(n)]
    dp=[[0]*n for _ in range(1<<n)]
    counts=[1]+[0]*((1<<n)-1)
    for mask in range(1,1<<n):
        vertices=mask
        while vertices:
            bit=vertices & -vertices
            vertices-=bit
            v=bit.bit_length()-1
            rest=mask ^ bit
            if not rest:
                dp[mask][v]=1
                continue
            prior=rest & incoming[v]
            while prior:
                b=prior & -prior
                prior-=b
                dp[mask][v]+=dp[rest][b.bit_length()-1]
        counts[mask]=sum(dp[mask])
    return counts

source=HERE.parent/'containment/receipt.json.gz'
head=RESULTS/'multicore_metagraph_20261113/next_target.json.gz'
d=json.loads(gzip.decompress(source.read_bytes()))
old=json.loads(gzip.decompress(head.read_bytes()))
targets=[1691,3311,6551,13031,25991]
expected=[3,3,4,2,2]
need([sum(r['target']==t for r in d['records']) for t in targets]==expected,
     'target multiplicities')
for t in targets:
    bound=next(x for x in old['bounds'] if x['H']==t)
    need(bound['complete'] and bound['strong_ready_classes']==expected[targets.index(t)],
         ('inherited completeness', t))
records=[]
for record in d['records']:
    rows=record['rows']
    counts=all_counts(rows)
    spectrum=sorted(Counter(counts).items())
    need([list(x) for x in spectrum]==record['spectrum'], ('full spectrum', record['target']))
    need(counts[-1]==record['target'], 'target')
    need(record['previous']==(record['target']+71)//2, 'predecessor')
    need(record['previous'] not in counts, 'forbidden predecessor')
    records.append({'target':record['target'],'parent_order':record['parent_order'],
                    'parent_rank':record['parent_rank'],'order':len(rows),
                    'subsets':len(counts),'distinct_counts':len(spectrum)})
    print(record['target'], record['parent_order'], record['parent_rank'], len(counts), flush=True)
need(sum(r['subsets'] for r in records)==321536, 'universe')
sources=[source,head,HERE.parent/'containment/count_spectrum.py',
         HERE.parent/'containment/PROOF.md',
         RESULTS/'multicore_metagraph_20261113/next_target.py',
         RESULTS/'two_entry_funnels_20261103/READY.md']
pins=[{'path':str(p.relative_to(RESULTS.parent.parent)).replace('\\','/'),
       'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in sources]
out={'status':'FINITE-EXACT independent all-subset recount; completeness inherited',
     'total_subsets':321536,'records':records,'source_pins':pins}
suffix='.optimized' if sys.flags.optimize else ''
(HERE/('containment_audit'+suffix+'.json')).write_text(json.dumps(out,indent=2)+'\n')
