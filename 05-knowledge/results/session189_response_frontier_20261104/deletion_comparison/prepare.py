"""Extract only the three reported pairs from the exact census export."""
import json,pathlib
P=pathlib.Path(__file__).resolve().parent
pairs=[[20197,48562],[86814,89262],[23898,23984]]
ids=sum(pairs,[])
cases=[]
for s in (P/'data/tower_9.txt').read_text().splitlines()[1:]:
    if int(s.split()[0]) not in ids:continue
    code=int(s.split('|')[-1]);rows=[0]*9;k=35
    for q in range(1,9):
        for a in range(q):
            if code>>k&1:rows[a]|=1<<q
            else:rows[q]|=1<<a
            k-=1
    cases.append(dict(node=int(s.split()[0]),source_line=s,code=code,rows=rows))
cases.sort(key=lambda d:ids.index(d['node']))
response=[]
for g in cases:
    for op in [False,True]:
        h=dict(g);h['converse']=op
        h['rows']=[511^x^(1<<i)for i,x in enumerate(g['rows'])]if op else g['rows']
        response.append(h)
(P/'representatives.json').write_text(json.dumps(dict(pairs=pairs,cases=cases,response_cases=response),indent=2)+'\n')
(P/'graphs.txt').write_text(str(len(response))+'\n'+'\n'.join('9 '+' '.join(map(str,g['rows']))for g in response)+'\n')
