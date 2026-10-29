from observable_quotient import *
from verify_surplus_jets import jet
from math import lcm,gcd
import re
from pathlib import Path
banks={}
for line in Path(__file__).with_name('prime_types12.inherited.txt').read_text().splitlines():
    count,aut=map(int,re.search(r'count=(\d+).*aut=(\d+)',line).groups())
    terms={}
    for token in line.split('AUT: ')[1].split(' SC=')[0].split():
        cycles,mult=token.split(':');part=[]
        for pair in cycles.rstrip(',').split(','):
            d,c=map(int,pair.split('^'));part += [d]*c
        terms[tuple(part)]=F(int(mult),aut)
    key=tuple(sorted(terms.items()))
    banks[key]=banks.get(key,0)+count
indices=[dict(k) for k in banks];counts=list(banks.values())
for h in range(3,7):
    jets=[jet(z,h) for z in indices]
    observed=[]
    for j in jets:
        out=defaultdict(F)
        for p,v in j.items():
            for r,w in grid_column(column(p,h)).items():out[r]+=v*w
        observed.append(dict(out))
    rank,kernel=echelon(observed)
    for rel in kernel:
        delta=defaultdict(F)
        for i,v in rel.items():
            for p,w in jets[i].items():delta[p]+=v*w
        delta={p:v for p,v in delta.items() if v}
        if not delta:continue
        denominator=lcm(*(v.denominator for v in rel.values()))
        rel={i:int(v*denominator) for i,v in rel.items()};g=gcd(*rel.values());rel={i:v//g for i,v in rel.items()}
        if not all(abs(v)<=counts[i] for i,v in rel.items()):continue
        delta=defaultdict(F)
        for i,v in rel.items():
            for p,w in jets[i].items():delta[p]+=v*w
        delta={p:v for p,v in delta.items() if v}
        observed_delta=defaultdict(F)
        for p,v in delta.items():
            for row,w in grid_column(column(p,h)).items():observed_delta[row]+=v*w
        require(all(not v for v in observed_delta.values()),'normalized observed difference not zero')
        require(h==3 and delta=={(1,1,1):F(-24),(3,):F(24)},'unexpected witness')
        next_delta=defaultdict(F)
        for i,v in rel.items():
            for p,w in jet(indices[i],h+1).items():next_delta[p]+=v*w
        # E=x for the actual palette TT1,TT2.
        next_count=sum(v for p,v in next_delta.items() if sum(p)==h+1)
        require(next_count==-18,'next horizon separation failed')
        print(json.dumps({'h':h,'next_horizon_TT1_TT2_difference':str(next_count),'inventory_relation':rel,'available':counts,
          'fits_distinct':all(abs(v)<=counts[i] for i,v in rel.items()),
          'delta_jet':{str(p):str(v) for p,v in delta.items()},
          'indices':[{str(k):str(v) for k,v in z.items()} for z in indices]},indent=2))
        raise SystemExit
print('NO WITNESS')
