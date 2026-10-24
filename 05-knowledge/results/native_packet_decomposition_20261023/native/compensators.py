"""Exact compensator occupancy and source-supply cuts inside integer-current fibre."""
import json
import sys
from collections import deque,Counter
from itertools import combinations
from pathlib import Path
from fibre import DSU
HERE=Path(__file__).resolve().parent
def route(states,a,b,budget=6):
    prev={a:None};todo=deque([a])
    while todo and b not in prev:
        x=todo.popleft()
        for y in states:
            if y not in prev and (x^y).bit_count()<=budget:prev[y]=x;todo.append(y)
    if b not in prev:return None
    out=[]
    while b is not None:out.append(b);b=prev[b]
    return list(reversed(out))
def comps(states,budget):
    ds=DSU(len(states))
    for i,x in enumerate(states):
        for j in range(i):
            if (x^states[j]).bit_count()<=budget:ds.join(i,j)
    return sorted(Counter(ds.root(i) for i in range(len(states))).values())
def main():
    f=json.loads((HERE/'fibre_exact_receipt.json').read_text());dom=json.loads((HERE/'domains.json').read_text());by={z['source']:z for z in dom['sources']}
    sources=f['sources'];outside=sum(1<<sources.index(s) for s in f['outside_sources']);a,b=f['initial'],f['terminal'];states=f['states']
    one=[m for m in states if (m&outside).bit_count()<=1];path=route(one,a,b)
    if path is None:raise RuntimeError('one-active route absent')
    singleton=[]
    for s in f['outside_sources']:
        keep=1<<sources.index(s);ss=[m for m in states if not m&outside&~keep]
        singleton.append(dict(source=s,states=ss,components_at6=comps(ss,6),endpoints_connected=route(ss,a,b) is not None))
    pairrows=[]
    for pair in combinations(f['outside_sources'],2):
        keep=sum(1<<sources.index(s) for s in pair);ss=[m for m in states if not m&outside&~keep]
        pp=route(ss,a,b)
        if pp is None:continue
        restrictedone=[m for m in ss if (m&outside).bit_count()<=1]
        pairrows.append(dict(sources=pair,states=ss,components_at6=comps(ss,6),one_active_route=route(restrictedone,a,b)))
    candidates=[p for p in pairrows if len(p['components_at6'])==1 and p['one_active_route'] is not None]
    selected=min(candidates,key=lambda p:(len(p['one_active_route']),len(p['states']),p['sources'])) if candidates else next(p for p in pairrows if list(p['sources'])==[2176,3015])
    path=selected['one_active_route'];ledger=[]
    sys.path.insert(0,str(HERE.parents[1]/'injective_native_blocks_20261023/search'))
    import explore as e
    target_binary=sum((e.Y[j][1965]%2)<<(j-1) for j in range(1,5))
    def binary(mask):
        y=0
        for i,s in enumerate(sources):
            if mask>>i&1:y^=sum(b<<j for j,b in enumerate(by[s]['code'][1:]))
        return y
    for m in path:
        old=sum(by[s]['X'] for i,s in enumerate(sources) if (m&~outside)>>i&1);extra=sum(by[s]['X'] for i,s in enumerate(sources) if (m&outside)>>i&1)
        oldcount=(m&~outside).bit_count();extracount=(m&outside).bit_count();dc=oldcount-8
        yp=binary(m&~outside)^binary(a);yq=binary(m&outside)
        ledger.append(dict(mask=m,V_sources=[s for i,s in enumerate(sources) if m>>i&1],
                           outside_V_sources=[s for i,s in enumerate(sources) if (m&outside)>>i&1],
                           old16_V_sum=old,outside_V_sum=extra,old16_V_count=oldcount,outside_V_count=extracount,
                           source_charge_X_change_old16=2*(old-1016),source_charge_X_change_outside=2*extra,
                           target_U_degree_delta_old16=-2*dc,target_U_degree_delta_outside=-2*extracount,
                           U_debt_old16=2*(old-1016-146*dc),U_debt_outside=2*(extra-146*extracount),U_debt_total=2*(old+extra-1016),
                           old16_binary_source_charge=yp,outside_binary_source_charge=yq,
                           old16_binary_residual=yp^(target_binary if dc%2 else 0),outside_binary_residual=yq^(target_binary if extracount%2 else 0),
                           V_debts='negatives of U debts'))
    transitions=[]
    for x,y in zip(path,path[1:]):
        changed=[s for i,s in enumerate(sources) if (x^y)>>i&1]
        transitions.append(dict(before=x,after=y,changed_sources=changed,slot_cost=len(changed),
                                 ordinary_owner_clock_occurrences=sum(sum(u!=v for u,v in zip(by[s]['old_before'],by[s]['old_after'])) for s in changed)))
    out=dict(one_active_full28_states=one,single_extra_source_controls=singleton,eligible_two_source_helper_pairs=pairrows,
             selected_two_source_helper=selected,selected_route=path,selected_route_debt_ledger=ledger,selected_route_events=transitions,
             selected_full_restricted_thresholds=[dict(budget=k,component_sizes=comps(selected['states'],k)) for k in (2,4,6,8)],
             selected_one_active_states=[m for m in selected['states'] if (m&outside).bit_count()<=1],
             selected_one_active_components_at6=comps([m for m in selected['states'] if (m&outside).bit_count()<=1],6),
             remaining_fixed_U_sources=[s for s in f['outside_sources'] if s not in selected['sources']])
    (HERE/'compensator_receipt.json').write_text(json.dumps(out,indent=2))
    print('one active states',len(one),'singleton successes',sum(z['endpoints_connected'] for z in singleton),'pair successes',len(pairrows),'fullyconnected oneactive pairs',len(candidates))
    print('selected',selected['sources'],'states',len(selected['states']),'components',selected['components_at6'],'route',path,'events',transitions)
    print(ledger)
if __name__=='__main__':main()
