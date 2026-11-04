import json,pathlib
P=pathlib.Path(__file__).resolve().parent
bank=json.loads((P/'bank.json').read_text());control=json.loads((P/'controls.json').read_text())
def augmented(rows,a,c):
    n=len(rows);r=rows[:]+[1<<(n+1),0]
    for i in range(n):
        if a>>i&1:r[i]|=1<<n
        else:r[n]|=1<<i
        if c>>i&1:r[i]|=1<<(n+1)
        else:r[n+1]|=1<<i
    return r
with(P/'independent_input.txt').open('w')as f:
    for b in bank:
        for a in b['masks']:
            for c in b['masks']:
                r=augmented(b['rows'],a,c);print(len(r),*r,file=f)
queries=[];index={};cases=[]
def put(g,a,c):
    key=(g,a,c)
    if key not in index:
        index[key]=len(queries);queries.append(dict(graph=g,a=a,c=c,rows=augmented(bank[g]['rows'],a,c)))
    return index[key]
for caseid,case in enumerate(control['cases']):
    groupmask=sum(1<<i for i in case['group']);base=(1<<(case['core_order']+1))if case['ambient']else 0
    for h,combo in enumerate(control['combinations']):
        item=dict(case=caseid,target=h,constant=control['constants'][h],terms=[])
        for p,q,r,w in combo:
            term=dict(p=p,q=q,r=r,weight=w,queries=[])
            for aa in range(4):
                for cc in range(4):
                    a=base|((1<<p)if aa&1 else 0)|((1<<q)if aa&2 else 0)
                    c=base|((1<<r)if cc&1 else 0)|(groupmask if cc&2 else 0)
                    term['queries'].append(dict(sign=(-1)**(4-aa.bit_count()-cc.bit_count()),source=put(case['source'],a,c),target=put(case['targets'][h],a,c)))
            item['terms'].append(term)
        cases.append(item)
(P/'physical_queries.json').write_text(json.dumps(dict(cases=cases,queries=queries),separators=(',',':'))+'\n')
with(P/'physical_input.txt').open('w')as f:
    for q in queries:print(len(q['rows']),*q['rows'],file=f)
print(len(bank),len(queries))
