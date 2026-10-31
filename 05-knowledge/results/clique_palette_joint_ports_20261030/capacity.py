"""Exact finite capacity recurrence for the hereditary join/C3/Paley7 substitution class."""
from itertools import combinations
from pathlib import Path
import argparse
import json

def templates():
    p=[sum(1 << j for j in range(7) if (j-i)%7 in (1,2,4)) for i in range(7)]
    tri=[]
    for t in combinations(range(7),3):
        if not all(sum(p[i]>>j&1 for j in t)==1 for i in t):
            tri.append(t)
    return [('join2',2,[(0,1)]),('C3',3,[(0,1),(0,2),(1,2)]),('P7',7,tri)]

def capacity(bound):
    values=[0,1]; witnesses=[None,{'leaf':True}]; counts=[]
    for k in range(2,bound+1):
        best=values[-1]; witness={'reuse_budget':k-1}
        seen=0
        for name,n,sets in templates():
            # The set system is vertex-transitive, so one maximum demand may be at port0.
            assigned=[0]*n
            touching=[[s for s in sets if i in s] for i in range(n)]
            def visit(i,limit,score):
                nonlocal best,witness,seen
                if i==n:
                    if sum(x>0 for x in assigned)<2:
                        return
                    seen+=1
                    if score>best:
                        best=score; witness={'template':name,'demands':assigned[:]}
                    return
                upper=limit
                for s in touching[i]:
                    upper=min(upper,k-sum(assigned[j] for j in s if j<i))
                for a in range(upper+1):
                    assigned[i]=a
                    visit(i+1,limit,score+values[a])
                assigned[i]=0
            for first in range(1,k):
                assigned[0]=first
                visit(1,first,values[first])
        values.append(best); witnesses.append(witness); counts.append({'budget':k,'candidates':seen})
    return {'values':values,'witnesses':witnesses,'enumeration':counts,
            'scope':'Exact recurrence for the declared hereditary substitution palette; not unrestricted Ramsey numbers.'}

if __name__=='__main__':
    ap=argparse.ArgumentParser(); ap.add_argument('--bound',type=int,default=12); ap.add_argument('--output',type=Path)
    args=ap.parse_args(); result=capacity(args.bound)
    text=json.dumps(result,indent=2,sort_keys=True)+'\n'
    if args.output: args.output.write_text(text,encoding='utf-8',newline='\n')
    else: print(text,end='')
