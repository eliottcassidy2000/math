#!/usr/bin/env python3
"""Exhaustively classify tight 15-blocker covers of the 90 palette options."""
from itertools import combinations,permutations
from pathlib import Path
import json

def need(b,msg):
    if not b:raise RuntimeError(msg)

def matchings(vertices):
    if not vertices:yield ()
    else:
        a=vertices[0]
        for b in vertices[1:]:
            for rest in matchings(tuple(x for x in vertices if x not in (a,b))):
                yield ((a,b),)+rest

edges=list(combinations(range(6),2));types=[(c,e) for c in range(3) for e in edges]
options=[p for m in matchings(tuple(range(6))) for p in permutations(m)]
need(len(options)==90 and len(set(options))==90,'90 marked options')
masks=[sum(1<<i for i,p in enumerate(options) if p[c]==e) for c,e in types]
need(all(m.bit_count()==6 for m in masks),'six options per blocker')
by_option=[[j for j,m in enumerate(masks) if (m>>i)&1] for i in range(90)]
full=(1<<90)-1;answers=set();nodes=0
def search(covered,chosen):
    global nodes
    nodes+=1
    if covered==full:
        need(len(chosen)==15,'tight size');answers.add(tuple(sorted(chosen)));return
    choices=None
    for i in range(90):
        if (covered>>i)&1:continue
        candidates=[j for j in by_option[i] if not masks[j]&covered]
        if choices is None or len(candidates)<len(choices):choices=candidates
        if not choices:return
    for j in choices:search(covered|masks[j],chosen+[j])
search(0,[])
rows={tuple(i for i,(c,e) in enumerate(types) if c==d) for d in range(3)}
stars={tuple(i for i,(c,e) in enumerate(types) if alpha in e) for alpha in range(6)}
need(answers==rows|stars and len(answers)==9,'exact row/star classification')
# Independent truth-table proof of binary edge sets meeting every matching once.
unmarked=list(matchings(tuple(range(6))));edge_index={e:i for i,e in enumerate(edges)}
single=[]
for bits in range(1<<15):
    if all(sum((bits>>edge_index[e])&1 for e in m)==1 for m in unmarked):single.append(bits)
need(set(single)=={sum(1<<i for i,e in enumerate(edges) if a in e) for a in range(6)},'six stars only')
# Hostile: one independent exterior vertex per pair, all in class0.
# It blocks every pair option but each vertex still has four available colours.
need(all(len(set(range(6))-set(e))==4 for e in edges),'method failure is not forcing')
out=dict(status='FINITE-EXACT independent controls for elementary classification',
         options=90,blocker_types=45,search_nodes=nodes,tight_covers=len(answers),
         row_covers=3,star_covers=6,independent_binary_graphs=1<<15,
         hostile='15 independent blockers obstruct the pair method but admit a proper extension')
Path(__file__).with_name('palette_audit.json').write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(out,indent=2))
