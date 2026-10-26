#!/usr/bin/env python3
"""Optional Z3 witness search: repair first word by at most two position swaps."""
from itertools import combinations
import json,sys
from pathlib import Path
import z3
from sat_repairs import solve
from search_repairs import directed,blue

def main():
    here=Path(__file__).resolve().parent;index=int(sys.argv[1])
    rows=json.loads((here.parent/'port_residual_integral_20261025/low_orders.normal.json').read_text())['failed_ledger_hostiles']
    row=[r for r in rows if r['n']==14][index];T=tuple(row['source_rows'])
    q=tuple(json.loads((here/'sat_repair_0.json').read_text())['words'][0])
    seen={q};layer=[q];bank=[]
    for depth in range(3):
        bank.extend(r for r in sorted(layer) if directed(T,r) and not blue(T,r))
        following=[]
        for r in layer:
            for a,b in combinations(range(14),2):
                s=list(r);s[a],s[b]=s[b],s[a];s=tuple(s)
                if s not in seen:seen.add(s);following.append(s)
        layer=following
    trials=[];result={'solver_result':'unknown','reason':'bounded bank exhausted'}
    for first in bank[:40]:
        r=solve(T,10000,first);trials.append({'first':first,'result':r['solver_result']})
        if r['solver_result']=='sat':result=r;break
    result.update({'n':14,'ports':row['ports'],'source_rows':T,'z3_version':z3.get_version_string(),
        'scope':'literal three-black detour; first word within two swaps of stored seed; at most 40 trials',
        'candidate_first_words':len(bank),'trials':trials})
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
