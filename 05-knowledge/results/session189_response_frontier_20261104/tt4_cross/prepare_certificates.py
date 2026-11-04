import json,pathlib
from exterior import capacities
from prepare import ORDERS
P=pathlib.Path(__file__).resolve().parent
bank=json.loads((P/'bank.json').read_text());control=json.loads((P/'controls.json').read_text());inner=json.loads((P/'inner_kernels.json').read_text());keys=[(x['p'],x['q'],x['r'])for x in inner[0]['coefficients']]
def coef(g,p,q,r,i):
    b=bank[g];index={v:j for j,v in enumerate(b['masks'])};return b['coefficients'][index[(1<<p)|(1<<q)]*len(index)+index[(1<<r)|(1<<i)]]
individual=[];cap=[]
for caseid,c in enumerate(control['cases']):
    for h in range(4):
        for p,q,r in keys:
            for i in range(4,len(bank[c['source']]['rows'])):individual.append([caseid,h,p,q,r,i,coef(c['source'],p,q,r,i)-coef(c['local_replacements'][h],p,q,r,i)])
    V=c['virtual'][:]
    if c['ambient']:
        n=len(V);V=[x|(1<<n)for x in V]+[1<<(n+1),(1<<n)-1]
    C=capacities(V);G=[i-3 for i in c['group']];cap.append(dict(case=caseid,rows=V,matrix=C,group=G,symmetric_capacity=sum(C[0][i]+C[i][0]for i in G)))
(P/'individual_quartics.json').write_text(json.dumps(individual,separators=(',',':'))+'\n')
(P/'capacity_controls.json').write_text(json.dumps(cap,indent=2)+'\n')
print(len(individual),[x['symmetric_capacity']for x in cap])
