"""One bounded discovery run; never used by the production proof checker."""
import json
from pathlib import Path

rows=[]
for e in range(14,81):
    n=(1<<e)-1;x=n;w=[]
    for _ in range(5000):
        z=3*x+1;a=(z & -z).bit_length()-1;x=z>>a;w.append(a)
        if x==1:break
    if x!=1:raise ValueError(('5000-step discovery cap',e))
    rows.append({'exponent':e,'source':n,'mod15':n%15,'mod35':n%35,'ROOT_word':w})
p=Path(__file__).resolve().parent/'inputs.json'
p.write_text(json.dumps({'scope':'exponents14..80 inclusive; strict odd ROOT discovery cap5000; no source omitted',
                         'rows':rows},indent=2)+'\n',encoding='utf-8')
print({'sources':len(rows),'max_ROOT_steps':max(len(r['ROOT_word']) for r in rows)})
