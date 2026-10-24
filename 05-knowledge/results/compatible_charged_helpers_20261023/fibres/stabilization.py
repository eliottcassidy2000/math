"""Occupied-context stabilization: all hostile margin4 pairs and a native repair."""
import json,sys
from collections import Counter,deque
from itertools import combinations
from pathlib import Path
from census import load,groups,components
HERE=Path(__file__).resolve().parent
def route(states,a,b):
    prev={a:None};todo=deque([a])
    while todo and b not in prev:
        x=todo.popleft()
        for y in states:
            if y not in prev and (x^y).bit_count()<=6:prev[y]=x;todo.append(y)
    if b not in prev:return None
    out=[]
    while b is not None:out.append(b);b=prev[b]
    return out[::-1]
def main():
    rows=load();sources=[r['source']for r in rows];g5=groups(5,rows);labels={}
    for key,ss in g5.items():
        for i,cc in enumerate(components(ss)):
            for m in cc:labels[m]=(key,i)
    old=json.loads((HERE/'census_k4.json').read_text());pairs=[];hist=Counter()
    for z in old['disconnected']:
        for i,cc in enumerate(z['components']):
            for dd in z['components'][i+1:]:
                for a in cc:
                    for b in dd:
                        eligible=[j for j in range(28)if not(a|b)>>j&1 and labels[a|1<<j]==labels[b|1<<j]]
                        if not eligible:raise RuntimeError('unrepaired margin4 pair')
                        hist[len(eligible)]+=1;pairs.append(dict(before=a,after=b,eligible_sources=[sources[j]for j in eligible]))
    bank=[176,1394,1477,2071]
    if not all(set(bank)&set(p['eligible_sources'])for p in pairs):raise RuntimeError('four-choice bank')
    # Exhaustive finite minimality: every three-source choice bank has a hostile pair.
    bank_hostiles=[]
    for triple in combinations(sources,3):
        i=next((i for i,p in enumerate(pairs)if not set(triple)&set(p['eligible_sources'])),None)
        if i is None:raise RuntimeError('three-choice bank succeeds')
        bank_hostiles.append(dict(bank=triple,pair_index=i))
    z=next(z for z in old['disconnected']if z['charge']==[15,3]);a,b=z['states'];core=a|b;outside=[j for j in range(28)if not core>>j&1]
    singleton=[];two=[];integer_two=[]
    for j in outside:
        key=((15+rows[j]['code'][0])%128,3^sum(v<<i for i,v in enumerate(rows[j]['code'][1:])))
        one=[m for m in g5[key]if not m&~(core|1<<j)]
        if route(one,a|1<<j,b|1<<j) is not None:raise RuntimeError('one new source unexpectedly repairs core')
        singleton.append(dict(occupied_source=sources[j],states=sorted(one)))
        for t in outside:
            if t==j:continue
            states=sorted(m for m in g5[key]if not m&~(core|1<<j|1<<t))
            rr=route(states,a|1<<j,b|1<<j)
            if rr is not None:two.append(dict(occupied_source=sources[j],unoccupied_source=sources[t],states=states,route=rr))
            exact=[m for m in states if sum(r['X']for i,r in enumerate(rows)if m>>i&1)==527+rows[j]['X']]
            ee=route(exact,a|1<<j,b|1<<j)
            if ee is not None:integer_two.append(dict(occupied_source=sources[j],unoccupied_source=sources[t],states=exact,route=ee))
    selected=next(z for z in integer_two if z['occupied_source']==1480 and z['unoccupied_source']==3000)
    rr=selected['route'];selected['route_source_sets']=[[sources[i]for i in range(28)if m>>i&1]for m in rr]
    selected['route_X_sums']=[sum(r['X']for i,r in enumerate(rows)if m>>i&1)for m in rr]
    selected['events']=[]
    for m,n in zip(rr,rr[1:]):
        ids=[i for i in range(28)if(m^n)>>i&1]
        selected['events'].append(dict(changed_sources=[sources[i]for i in ids],changed_slots=len(ids),
            ordinary_owner_clock_work=sum(sum(x!=y for x,y in zip(rows[i]['old_before'],rows[i]['old_after']))for i in ids)))
    result=dict(scope='Complete28binary menu; occupied-helper stabilization changes endpoint context, not original margin4 endpoints',
        all_disconnected_margin4_state_pairs=pairs,eligible_helper_histogram=dict(sorted(hist.items())),minimum_choice_bank=bank,
        all_three_choice_bank_hostiles=bank_hostiles,minimal_hostile=next(z for z in old['disconnected']if z['charge']==[15,3]),
        all_twenty_single_source_hostiles=singleton,successful_two_source_contexts=two,successful_exact_integer_contexts=integer_two,
        selected_ten_source_exact_context=selected)
    result=json.loads(json.dumps(result));out=HERE/'stabilization_receipt.json'
    if '--write' in sys.argv:out.write_text(json.dumps(result,indent=2))
    elif json.loads(out.read_text())!=result:raise RuntimeError('stabilization receipt')
    print('PASS',len(pairs),'hostile statepairs; helper histogram',dict(hist),'; choicebank4 minimal;',len(two),'two-helper contexts;',len(integer_two),'integer contexts')
    print('selected',selected)
if __name__=='__main__':main()
