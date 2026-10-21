"""Bounded discovery only; production check.py never searches for a ROOT route."""
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent
data=json.loads((HERE.parent/'scale_credit_supply_20261021a/receipt.json').read_text())['case8192']
parents=sorted({r['parent'] for r in data['linear_rows']+data['frontier_supplement']})
rows=[]
for u in parents:
    x=u;w=[]
    for _ in range(10000):
        if x==1:break
        y=3*x+1;a=(y&-y).bit_length()-1;x=y>>a;w.append(a)
    if x!=1:raise ValueError('bounded discovery exhausted')
    rows.append({'parent':u,'mod15':u%15,'mod35':u%35,'ROOT_word':w})
(HERE/'inputs.json').write_text(json.dumps({'discovery_cap':10000,'parents':rows},indent=2)+'\n',encoding='utf-8')
print(json.dumps({'parents':len(rows),'max_ROOT_steps':max(len(r['ROOT_word']) for r in rows)}))
