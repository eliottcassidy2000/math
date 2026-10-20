"""Hostile probes for inflating the context of the six order-five certificates."""
import json
from pathlib import Path
from inventory import presentations,join,blue,move,need,arc
HERE=Path(__file__).resolve().parent

def inflate(p,m):
    q=[];i=0
    while i<len(p):
        if p[i]<2:
            need(i+1<len(p) and p[i+1]==1-p[i],'split context')
            q.extend(range(m) if p[i]==0 else reversed(range(m)));i+=2
        else:q.append(p[i]+m-2);i+=1
    return tuple(q)

def check(a,row):
    m=len(a);state=join(a,row['atom']);flags=[];scores=[]
    for p in row['paths']:
        p=inflate(p,m)
        missing=[(u,v) for u,v in zip(p,p[1:]) if not arc(state,u,v)]
        if missing:return dict(failed='directedness',step=len(flags)+1,A=a,missing=missing)
        flags.append(blue(state,p))
        scores.append([state[v].bit_count() for v in p])
        state=move(state,p)
    local=tuple(row['target'][i]>>2 for i in range(2,7))
    need(state==join(a,local),'wrong inflated target')
    if any(flags):return dict(failed='blackness',A=a,flags=flags,scores=scores)
    return None

def main():
    rows=json.loads(HERE.joinpath('five_atom_probe.json').read_text());report=[]
    for idx,row in enumerate(rows):
        counts=[];hostiles=[]
        for m in range(1,8):
            count=fail=0
            for a in presentations(m):
                count+=1; result=check(a,row)
                if result:
                    fail+=1
                    if len(hostiles)<4:hostiles.append(result)
            counts.append(dict(m=m,count=count,fail=fail))
        result=dict(atom=idx,counts=counts,hostiles=hostiles)
        report.append(result);print(json.dumps(result),flush=True)
    HERE.joinpath('five_atom_lift.json').write_text(json.dumps(report,indent=2)+'\n')

if __name__=='__main__':main()
