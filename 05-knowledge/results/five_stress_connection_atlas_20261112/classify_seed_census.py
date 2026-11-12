from pathlib import Path
import json
from collections import Counter

HERE=Path(__file__).resolve().parent

def run():
    data=json.loads((HERE/'tournament_receipt.json').read_text())
    records=[];levels={}
    for q in data['seed_records']:
        w=q['word'];rows=q['rows'];O=q['selected'];n=len(w)
        runs=[]
        for i,v in enumerate(w):
            if not runs or w[runs[-1][-1]]!=v:runs.append([i])
            else:runs[-1].append(i)
        if len(runs)%2==0:runs=[runs[-1]+runs[0]]+runs[1:-1]
        label='other';pars={}
        if q['H']==1:label='transitive singleton mark'
        elif len(runs)==3:
            # Follow the actual directed C3, not merely projective run order.
            start=0;nexti=next(j for j in (1,2) if rows[runs[0][0]]>>runs[j][0]&1)
            order=[0,nexti,3-nexti];runs=[runs[j] for j in order]
            selected=[[i for i,v in enumerate(r) if O>>v&1] for r in runs]
            marked=[j for j,s in enumerate(selected) if s]
            if len(marked)==2 and all(len(runs[j])==1 for j in marked):
                label='book';pars={'L':n-2}
            elif len(marked)==1:
                j=marked[0];pos=selected[j]
                if len(pos)==2 and sorted(len(runs[k]) for k in range(3) if k!=j)==[1,2] and pos[1]-pos[0] in (1,2):
                    a=pos[0];b=pos[1]-pos[0]-1;c=len(runs[j])-pos[1]-1
                    orient=[len(runs[(j+1)%3]),len(runs[(j+2)%3])]
                    label='marked pair family';pars={'a':a,'b':b,'c':c,'other_blocks':orient}
            if label=='other':label='C3 exception in this census'
        elif len(runs)==5:label='C5 exception in this census'
        item={**q,'class':label,'parameters':pars,'marked_blocks':[''.join('1' if O>>v&1 else '0' for v in r) for r in runs]}
        records.append(item)
        levels.setdefault(n,Counter())[label]+=1
    out={'scope':data['scope'],'by_size':{n:dict(v) for n,v in levels.items()},'totals':dict(Counter(r['class'] for r in records)),'exceptions':[r for r in records if 'exception' in r['class'] or r['class']=='other']}
    for n,v in levels.items():
        if n>=5:
            if v['transitive singleton mark']+v['book']+v['marked pair family']!=5*n-17:raise RuntimeError(('family count',n))
    (HERE/'seed_classification.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf8')
    print(json.dumps(out['totals']))
    print(json.dumps([(r['n'],r['H'],r['U'],r['marked_blocks']) for r in out['exceptions']]))

if __name__=='__main__':run()
