"""Complete source-specific admission from each off-hub native configuration."""
from collections import defaultdict,Counter
from itertools import combinations
from pathlib import Path
import json,sys
HERE=Path(__file__).resolve().parent;RESULTS=HERE.parents[1];N=0
def need(v,m):
    global N
    N+=1
    if not v:raise RuntimeError(m)
def main():
    d=json.loads((RESULTS/'compatible_charged_helpers_20261023/supply/selected_supply.json').read_text())
    new=json.loads((HERE/'native_2125_receipt.json').read_text());by={z['source']:z for z in d['helpers']}
    ids=sorted(by)+[2125];labels={s:z['label']for s,z in by.items()};labels[2125]=new['label']
    costs={s:z['selected']['cost']for s,z in by.items()};costs[2125]=8;free=[s for s in ids if s not in(368,776,3139)]
    def charge(ss):
        y=0
        for s in ss:y^=sum(b<<j for j,b in enumerate(labels[s][1:]))
        return((sum(labels[s][0]for s in ss)%128)<<4)|y
    hub=defaultdict(list)
    for t in combinations(free,3):hub[charge((368,)+t)].append(frozenset((368,)+t))
    need(set(hub)==set(range(2048))and sum(map(len,hub.values()))==20825,'complete20825-state native hub')
    rows=[];hist=Counter();off_by_charge=defaultdict(list)
    for p in combinations(free,2):
        source=frozenset((776,3139)+p);q=charge(source);options=[]
        for target in hub[q]:
            changed=source^target;options.append((sum(costs[s]for s in changed),len(changed),tuple(sorted(target))))
        cost,support,target=min(options);need(charge(target)==q and support<=8,'minimum admission retains same actual charge')
        need(cost==min(v[0]for v in options),'full source-specific minimum over complete target fibre')
        rows.append(dict(source=sorted(source),charge=q,target=target,minimum_clock_cost=cost,support=support,
                         target_fibre_count=len(options),cost_minimizer_count=sum(v[0]==cost for v in options)))
        hist[cost]+=1
        off_by_charge[q].append(source)
    worst=max(hist);hard=[r for r in rows if r['minimum_clock_cost']==worst]
    need(len(rows)==1275,'all actual off-hub configurations, not one per charge')
    common=[]
    for q,sources in off_by_charge.items():
        for a,b in combinations(sources,2):
            choices=[(max(sum(costs[s]for s in a^t),sum(costs[s]for s in b^t)),tuple(sorted(t)))for t in hub[q]]
            bound,target=min(choices)
            need(charge(target)==q,'actual common-hub two-endpoint witness')
            common.append(dict(first=sorted(a),second=sorted(b),charge=q,common_hub=target,minimum_peak=bound))
    out=dict(status='FINITE-EXACT complete native off-hub admission minima',source_count=len(rows),
             hub_state_count=20825,minimum_clock_histogram=dict(sorted(hist.items())),worst_minimum_clock_cost=worst,
             max_support=max(r['support']for r in rows),worst_rows=hard,all_minimum_witnesses=rows,
             same_charge_off_pair_count=len(common),worst_common_hub_peak=max(r['minimum_peak']for r in common),
             all_common_hub_witnesses=common,checks=N)
    out=json.loads(json.dumps(out));path=HERE/'admission_receipt.json'
    if '--write'in sys.argv:path.write_text(json.dumps(out,indent=2,sort_keys=True)+'\n',newline='\n')
    else:need(out==json.loads(path.read_text()),'frozen complete admission minima')
    print(f'PASS {N}; worst optimal admission={worst}; same-charge pairs={len(common)}, common peak={out["worst_common_hub_peak"]}; histogram={dict(sorted(hist.items()))}')
if __name__=='__main__':main()
