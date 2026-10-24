#!/usr/bin/env python3
"""Complete one-active-helper layers for the NEW primitive18 circuit, not the old16."""
from pathlib import Path
from collections import defaultdict,deque
from itertools import combinations
import json,subprocess
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[2];CHECKS=0
def need(x,m):
    global CHECKS
    CHECKS+=1
    if not x:raise RuntimeError(m)
def load(path,rev='0f9dec9c3e1d'):
    return json.loads(subprocess.check_output(['git','show',rev+':05-knowledge/results/'+path],cwd=ROOT))
def main():
    inv=load('session45_injective_multisource_moves_20261023/multisource-trades/target_pair_inventory.json')
    rows={r['source']:r for r in inv['records']};sources=sorted(rows);index={s:i for i,s in enumerate(sources)}
    labels={s:(sum(x<<j for j,x in enumerate(rows[s]['source_label'][:4])),rows[s]['source_label'][4]) for s in sources}
    plus=[176,1394,1457,1477,1532,1923,2087,2176,3369];minus=[1425,1549,1841,1921,2185,2553,2835,2905,3000]
    old=sorted(plus+minus);helpers=sorted(set(sources)-set(old))
    def mask(ss):return sum(1<<index[s] for s in ss)
    def charge(ss):
        p=x=0
        for s in ss:p^=labels[s][0];x+=labels[s][1]
        return p,x%128
    wanted=charge(minus);need(charge(plus)==wanted,'actual endpoints share full paired charge')
    # All old-source subsets of size8 or9: independent of the primitive18 bit-DP.
    layers={0:[]};bycharge=defaultdict(list)
    for I in combinations(old,9):
        if charge(I)==wanted:layers[0].append(mask(I))
    for I in combinations(old,8):bycharge[charge(I)].append(mask(I))
    for h in helpers:
        q=(wanted[0]^labels[h][0],(wanted[1]-labels[h][1])%128)
        layers[h]=[m|mask([h]) for m in bycharge[q]]
    initial=mask(minus);terminal=mask(plus)
    need(set(layers[0])=={initial,terminal},'complete no-helper layer agrees with independent primitive audit')
    states=sorted(m for layer in layers.values() for m in layer);need(len(states)==len(set(states)),'disjoint actual helper layers')
    layerof={m:h for h,layer in layers.items() for m in layer}
    adj={b:{m:[] for m in states} for b in (4,6,8)}
    for a,c in combinations(states,2):
        d=(a^c).bit_count()
        for b in adj:
            if d<=b:adj[b][a].append(c);adj[b][c].append(a)
    def route(H,b):
        allowed=set(H)|{0};seen={initial:None};todo=deque([initial])
        while todo:
            m=todo.popleft()
            if m==terminal:
                path=[]
                while m is not None:path.append(m);m=seen[m]
                return path[::-1],len(seen)
            for q in adj[b][m]:
                if q not in seen and layerof[q] in allowed:seen[q]=m;todo.append(q)
        return None,len(seen)
    reports=[];best=None
    for b in (4,6,8):
        found=[];counts=[]
        for k in range(1,4):
            ok=0
            for H in combinations(helpers,k):
                path,n=route(H,b)
                if path is not None:
                    ok+=1
                    if not found:found=[H,path]
            counts.append({'helper_count':k,'successful_subsets':ok,'subsets_tested':sum(1 for _ in combinations(helpers,k))})
            if ok:break
        reports.append({'budget':b,'counts':counts,'first_helpers':found[0] if found else None})
        if best is None and found:best=(b,*found)
    need(best is not None,'a supplied charged-helper route exists')
    b,H,path=best
    component_reports=[]
    for choice in [tuple(H),tuple(helpers)]:
        allowed=[m for m in states if layerof[m] in set(choice)|{0}]
        for budget in (4,6):
            unseen=set(allowed);components=[]
            while unseen:
                q=deque([min(unseen)]);component={q[0]};unseen.remove(q[0])
                while q:
                    for x in adj[budget][q.popleft()]:
                        if x in unseen:unseen.remove(x);component.add(x);q.append(x)
                components.append(component)
            component_reports.append({'helpers':choice,'budget':budget,'states':len(allowed),
              'component_sizes':sorted(map(len,components)),
              'endpoints_connected':any(initial in c and terminal in c for c in components)})
    pot=load('metagraph_order8_electrical_20261023/receipt.json')['order8'];X=pot['higher_layers'][-1]['lifts'][0];Y=pot['radical_lifts_mod4']
    rawminus=sum(X[s] for s in minus);target_X=X[1965];ledgers=[];work=[];ports=[]
    domains=load('native_packet_decomposition_20261023/native/domains.json','a793424e20')
    # Native domain evidence is inherited, with complete raw-path exception retained.
    needed=set(old)|set(H);domainmap={r['source']:r for r in domains['sources']}
    need(needed<=set(domainmap),'every helper has a complete authenticated native domain')
    full=(1<<21)-1;chords=[(i,j) for i,j in combinations(range(8),2) if j>i+1];chordindex={z:i for i,z in enumerate(chords)}
    def orbit(w):
        R=sum(((w>>k)&1)<<chordindex[7-j,7-i] for k,(i,j) in enumerate(chords))
        return {w,full^w,R,full^R}
    for m in path:
        V={s for s in needed if m>>index[s]&1};need(len(V)==9 and charge(V)==wanted,'margin and full current at each integer arrival')
        occupied=set()
        for s in needed:
            rec=rows[s]['after' if s in V else 'before'];o=orbit(rec['mask'])
            need(len(o)==4 and not occupied&o,'full normalized resource injection at every arrival');occupied|=o
        active=V-set(old);need(len(active)<=1,'at most one charged helper active')
        oldV=V&set(old);dcount=len(oldV)-9;oldsum=sum(X[s] for s in oldV);hsum=sum(X[s] for s in active)
        oldcurrent=2*(oldsum-rawminus-target_X*dcount);helpercurrent=2*(hsum-target_X*len(active))
        need((oldcurrent+helpercurrent)%256==0,'physical owned X current includes the diagonal degree term')
        binary=[]
        for j in range(4):
            z=sum(Y[j][s] for s in V)-sum(Y[j][s] for s in minus)
            need(z%2==0,'every target-owned Y charge cancels');binary.append(z%2)
        ledgers.append({'mask':m,'V_sources':sorted(V),'active_helper':sorted(active),'old_U_current_X':oldcurrent,
          'helper_U_current_X':helpercurrent,'total_U_current_X':oldcurrent+helpercurrent,'V_currents':'negatives of U currents','Y_residuals':binary})
        ports.append(len(occupied))
    for a,c in zip(path,path[1:]):
        changed=[s for s in needed if (a^c)>>index[s]&1];need(len(changed)<=b,'actual atomic slot bound')
        work.append(sum(rows[s]['cost'] for s in changed))
    need(all(not initial>>index[h]&1 and not terminal>>index[h]&1 for h in H),'all helper endpoints restored')
    return {'status':'FINITE-EXACT native one-active-helper route for new primitive18',
      'sources':sources,'old_sources':old,'helper_candidates':helpers,'layer_sizes':{str(k):len(v) for k,v in layers.items()},
      'complete_one_active_states':len(states),'search':reports,'components':component_reports,
      'selected_budget':b,'selected_helpers':H,
      'route':path,'route_owner_clock_costs':work,'total_owner_clock_work':sum(work),'ledger':ledgers,
      'closed_ports_per_arrival':ports,'current_scope':'full original dyadic bank, not exact unscaled integer-X conservation',
      'checks':CHECKS}
if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
