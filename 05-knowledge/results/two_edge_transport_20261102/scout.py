"""Exhaustive two-edge escape and restricted marked-carrier probes."""
from collections import Counter, defaultdict
from itertools import combinations
from pathlib import Path
import importlib.util
import json
import sys
HERE=Path(__file__).resolve().parent
PREV=HERE.parent/'boolean_realization_71_17_20261102'
sys.path.insert(0,str(PREV))
spec=importlib.util.spec_from_file_location('previous_counts',PREV/'scout.py')
old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
rows,h_count=old.rows,old.h_count

def main():
    n=6;edges=list(combinations(range(n),2));d=len(edges);target=17
    counts=[h_count(rows(n,x)) for x in range(1<<d)]
    traps=[x for x,h in enumerate(counts) if (h-target)**2==4 and
           all((counts[x^(1<<e)]-target)**2>4 for e in range(d))]
    records=[];types=defaultdict(list)
    for x in traps:
        escapes=[]
        for a,b in combinations(range(d),2):
            y=x^(1<<a)^(1<<b)
            if counts[y]!=target:continue
            mids=[counts[x^(1<<a)],counts[x^(1<<b)]]
            barrier=min((v-target)**2 for v in mids)
            escapes.append({'edges':[a,b],'target_mask':y,'intermediate_counts':mids,'barrier':barrier,
                            'adjacent_edges':bool(set(edges[a])&set(edges[b]))})
        out=rows(n,x)
        universal=[(i,'source' if r.bit_count()==n-1 else 'sink') for i,r in enumerate(out) if r.bit_count() in (0,n-1)]
        key=(counts[x],len(escapes),min((z['barrier'] for z in escapes),default=None),
             sum(z['adjacent_edges'] for z in escapes),tuple(t for _,t in universal))
        types[key].append(x)
        records.append({'source':x,'H':counts[x],'universal':universal,'escapes':escapes})
    result={'trap_count':len(traps),'escape_count_histogram':dict(sorted(Counter(len(r['escapes']) for r in records).items())),
            'types':[{'key':key,'count':len(xs),'representative':xs[0],'out_rows':rows(n,xs[0])} for key,xs in types.items()],
            'records':records}
    (HERE/'escape_scout.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k!='records'},indent=2))

if __name__=='__main__':main()
