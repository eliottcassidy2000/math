"""Freeze two discovered native four-source circuits and their shared-owner coupling."""
import json,copy
from collections import Counter
from itertools import product
import explore as e
blocks=json.loads((e.HERE/'candidate_blocks.json').read_text())
first=copy.deepcopy(blocks[0]);second=copy.deepcopy(blocks[4])
for s in second['slots']:
    s['before'],s['after']=s['after'],s['before'];s['sign']=-s['sign']
second['source_orientation']=[-1,-1,1,1]
slots=first['slots']+second['slots']
for s in slots:
    s['canon']=int(e.NODES[s['source']]['canon_a'])
    s['class_a']=int(e.NODES[s['source']]['class_a']);s['class_b']=int(e.NODES[s['source']]['class_b'])
    s['all_unchanged_clocks']=[[i,a] for i,(a,b) in enumerate(zip(s['before']['path'],s['after']['path'])) if a==b]
shared=[[[0,2],[4,5]],[[0,5],[4,2]]]
out={'targets':first['targets'],'first':first,'second_reversed':second,'slots':slots,'shared_owner_equalities':shared}
(e.HERE/'selected_blocks.json').write_text(json.dumps(out,indent=2))
survivors=[]
for bits in product(range(2),repeat=8):
    state=[s['after'] if b else s['before'] for s,b in zip(slots,bits)]
    if len({(s['source'],tuple(r['path'])) for s,r in zip(slots,state)})<8:continue
    if any(state[a]["path"][i]!=state[b]["path"][j] for (a,i),(b,j) in shared):continue
    if Counter(r['target'] for r in state)!=Counter(s['before']['target'] for s in slots):continue
    w=Counter()
    for b,s in zip(bits,slots):w[s['source']]+=b*s['sign']
    if sum(c*e.X[s] for s,c in w.items())%128:continue
    if any(sum(c*e.SIG[s][j] for s,c in w.items())%m for j,m in enumerate([4,2,2,2,2])):continue
    survivors.append(bits)
print('eight-slot admissible states',survivors)
