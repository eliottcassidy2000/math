#!/usr/bin/env python3
from pathlib import Path
from itertools import combinations,permutations
import json
from independent_carrier import packet,check
P=Path(__file__).resolve().parents[1]
source=json.loads((P/'renewal/RESULTS.json').read_text())
checked=[]
for w in source['literal_nine'][:8]:
    p=packet(tuple(w['rows']))
    check(p==w['profile'],'independent full bank packet')
    checked.append(p)
check([p['H']for p in checked]==list(range(54365,54380,2)),'eight raw scalar interval')
check({(p['h'],p['V'])for p in checked}=={(15,139878)},'eight common raw receiver')
# Independent direct Hamilton permutations + reachability, no inherited DP.
S4=(2,4,9,3);census=[]
for n in range(1,6):
    pairs=list(combinations(range(n),2));five=0
    for code in range(1<<len(pairs)):
        T=[0]*n
        for j,(u,v)in enumerate(pairs):
            if code>>j&1:T[u]|=1<<v
            else:T[v]|=1<<u
        strong=True
        for start in range(n):
            todo=[start];seen={start}
            while todo:
                u=todo.pop()
                for v in range(n):
                    if T[u]>>v&1 and v not in seen:seen.add(v);todo.append(v)
            strong &= len(seen)==n
        if not strong:continue
        H=sum(all(T[p[i]]>>p[i+1]&1 for i in range(n-1))for p in permutations(range(n)))
        if H!=5:continue
        five+=1;check(n==4,'strong H5 order bound finite base')
        check(any(all((T[i]>>j&1)==(S4[p[i]]>>p[j]&1)for i in range(n)for j in range(n))for p in permutations(range(n))),'strong H5 unique core')
    census.append(five)
check(census==[0,0,0,24,0],'exact H5 core census')
print(json.dumps(dict(status='PASS',raw_bank_profiles=checked,strong_H5_labelled_by_order=census,method='independent backward start/end tables and root splitting; direct permutation core census'),indent=2,sort_keys=True))
