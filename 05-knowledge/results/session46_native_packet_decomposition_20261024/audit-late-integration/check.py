#!/usr/bin/env python3
"""Bounded independent controls for the conditional multi-target decomposition.
This does not reconstruct incoming native paths or replay incoming verifiers.
"""
from collections import Counter, defaultdict
from itertools import combinations, product
import argparse,json
from pathlib import Path

CHECKS=0

def need(ok,why):
    global CHECKS
    CHECKS+=1
    if not ok: raise RuntimeError(why)


def cycles(old,new):
    remaining={i for i,(a,b) in enumerate(zip(old,new)) if a!=b}
    out=[]
    while remaining:
        i=min(remaining)
        vertices=[old[i]]; edges=[];at=old[i]
        while True:
            choices=[j for j in sorted(remaining) if old[j]==at and j not in edges]
            if not choices: raise RuntimeError('balanced directed remainder dead end')
            e=choices[0];edges.append(e);at=new[e]
            if at in vertices:
                start=vertices.index(at);cycle=edges[start:]
                out.append(cycle);remaining.difference_update(cycle);break
            vertices.append(at)
    return out


def target_charge(state,labels,t,m):
    return tuple(sum(labels[i] for i,v in enumerate(state) if v==j)%m for j in range(t))


def route(old,new,labels,t,m):
    inventory=cycles(old,new)
    need(sorted(i for c in inventory for i in c)==[i for i in range(len(old)) if old[i]!=new[i]],
         'cycle decomposition lost or repeated a slot')
    debts=[]
    for c in inventory:
        need(len(c)<=t,'simple cycle too long')
        need(Counter(old[i] for i in c)==Counter(new[i] for i in c),'cycle lost exact margins')
        d=tuple((sum(labels[i] for i in c if new[i]==v)-
                 sum(labels[i] for i in c if old[i]==v))%m for v in range(t))
        need(sum(d)%m==0,'cycle debt not in sum-zero group')
        debts.append(d)
    state=list(old);rest=list(range(len(inventory)));maximum=0
    bound=t*(1+(t-1)*(m-1)) # cyclic p-group C2 or C3 used below
    while rest:
        selected=None
        for k in range(1,len(rest)+1):
            for chosen in combinations(rest,k):
                if all(sum(debts[j][v] for j in chosen)%m==0 for v in range(t)):
                    selected=chosen;break
            if selected is not None:break
        need(selected is not None,'zero-sum total has no packet')
        need(len(selected)<=1+(t-1)*(m-1),'cycle packet exceeds Davenport bound')
        event=[i for j in selected for i in inventory[j]]
        maximum=max(maximum,len(event))
        need(len(event)<=bound,'slot packet exceeds tD bound')
        for i in event:
            need(state[i]==old[i],'slot changed more than once')
            state[i]=new[i]
        need(Counter(state)==Counter(old),'arrival lost exact target margin')
        need(target_charge(state,labels,t,m)==target_charge(old,labels,t,m),
             'arrival lost a target-owned charge')
        rest=[j for j in rest if j not in selected]
    need(tuple(state)==tuple(new),'route misses endpoint')
    return len(inventory),maximum


def exhaustive_controls():
    reports=[]
    for m in (2,3):
        for t in (2,3):
            pairs=0;nontrivial=0;largest=0
            for n in range(1,5):
                for labels in product(range(m),repeat=n):
                    fibers=defaultdict(list)
                    for state in product(range(t),repeat=n):
                        key=(tuple(state.count(v) for v in range(t)),target_charge(state,labels,t,m))
                        fibers[key].append(state)
                    for states in fibers.values():
                        for old,new in product(states,repeat=2):
                            count,maximum=route(old,new,labels,t,m)
                            pairs+=1;nontrivial+=old!=new;largest=max(largest,maximum)
            reports.append({'cyclic_group_order':m,'targets':t,'maximum_slots':4,
                            'accepted_ordered_endpoint_pairs':pairs,
                            'nontrivial_pairs':nontrivial,'largest_event_observed':largest,
                            'proved_bound':t*(1+(t-1)*(m-1))})
    return reports


def longer_cycle_packet():
    # Three cycles in a prescribed decomposition with a minimal sequence in C2^2.
    # This does not assert the physical endpoint trade itself is primitive.
    old=(0,1,2)*3;new=(1,2,0)*3
    labels=(0,0,1, 0,1,0, 0,1,1)
    explicit=[(0,1,2),(3,4,5),(6,7,8)]
    debts=[]
    for c in explicit:
        debts.append(tuple((sum(labels[i] for i in c if new[i]==v)-
                            sum(labels[i] for i in c if old[i]==v))%2 for v in range(3)))
    zero_sets=[]
    for bits in product((0,1),repeat=3):
        if all(sum(bits[j]*debts[j][v] for j in range(3))%2==0 for v in range(3)):
            zero_sets.append(bits)
    need(zero_sets==[(0,0,0),(1,1,1)],'prescribed cycle sequence not minimal zero-sum')
    count,maximum=route(old,new,labels,3,2)
    need(maximum<=9,'longer example bound')
    need(3*(1+2*(132-1))==789 and 4*(1+3*(132-1))==1576,'original-code bound arithmetic')
    return {'slots':9,'prescribed_cycle_debts':debts,
            'proper_nonempty_zero_sum_cycle_subsets':0,'algorithm_cycles':count,
            'algorithm_maximum_event':maximum,
            'native_sharpness_claimed':False}


def main():
    p=argparse.ArgumentParser();p.add_argument('--out',required=True);args=p.parse_args()
    report={'exhaustive_models':exhaustive_controls(),'longer_cycle_control':longer_cycle_packet()}
    report.update(status='PASS',checks=CHECKS,
                  scope='Independent abstract multi-target algorithm controls; no incoming native replay')
    Path(args.out).write_text(json.dumps(report,sort_keys=True,indent=2)+'\n')
    print(json.dumps({'status':'PASS','checks':CHECKS,'out':args.out}))

if __name__=='__main__':main()
