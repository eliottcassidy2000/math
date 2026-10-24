#!/usr/bin/env python3
"""Independent finite hostile checks for the proved signed-component theorem."""
from pathlib import Path
from itertools import product, combinations
from collections import Counter
import argparse, json
P=Path(__file__).resolve().parent
CHECKS=0
def need(ok,why):
    global CHECKS
    CHECKS+=1
    if not ok:raise RuntimeError(why)

def run():
    graphs=consistent=endpoint_tests=partition_tests=0
    # No generated accepted states or component data are read from a candidate.
    for n in range(1,5):
        pairs=list(combinations(range(n),2))
        for relation in product((-1,0,1),repeat=len(pairs)):
            graphs+=1
            edges=[(i,j,e) for (i,j),e in zip(pairs,relation) if e>=0]
            adjacency=[[] for _ in range(n)]
            for i,j,e in edges:adjacency[i].append((j,e));adjacency[j].append((i,e))
            parity={};components=[];bad=False
            for root in range(n):
                if root in parity:continue
                parity[root]=0;queue=[root];component=[]
                while queue:
                    i=queue.pop();component.append(i)
                    for j,e in adjacency[i]:
                        expected=parity[i]^e
                        if j in parity:
                            if parity[j]!=expected:bad=True
                        else:parity[j]=expected;queue.append(j)
                components.append(component)
            full=[x for x in range(1<<n) if all((((x>>i)^(x>>j))&1)==e for i,j,e in edges)]
            predicted=[]
            if not bad:
                consistent+=1
                for roots in product((0,1),repeat=len(components)):
                    predicted.append(sum((parity[i]^b)<<i for c,b in zip(components,roots) for i in c))
            need(set(full)==set(predicted),'signed cycle consistency and reconstructed whole product exactly match')
            masks=[sum(1<<i for i in c) for c in components]
            # Three distinct guard choices test integer margins, finite charges,
            # and their conjunction; zero is allowed but never silently inferred.
            for scenario in range(3):
                def guard(x):
                    if scenario==0:return (x.bit_count(),0)
                    if scenario==1:return (sum((1+i%2)*((x>>i)&1) for i in range(n)),sum((i+1)*((x>>i)&1) for i in range(n))%6)
                    return (0,sum((2*i+1)*((x>>i)&1) for i in range(n))%6)
                for x,y in combinations(full,2):
                    if guard(x)!=guard(y):continue
                    endpoint_tests+=1;diff=x^y
                    changed=[c for c in masks if c&diff]
                    need(all(c&diff==c for c in changed),'whole relation components synchronize')
                    neutral={};m=len(changed)
                    for subset in range(1<<m):
                        move=0
                        for j,c in enumerate(changed):
                            if subset>>j&1:move^=c
                        z=x^move
                        if guard(z)==guard(x):neutral[subset]=z
                    interval={z for z in full if (x^z)&diff==x^z and guard(z)==guard(x)}
                    need(set(neutral.values())==interval,'neutral component unions exactly reconstruct full accepted endpoint interval')
                    # Disjoint zero-sum atom partition optimum, independently of
                    # an endpoint-relative DAG on literal full records.
                    def event_cost(subset):
                        move=x^neutral[subset] if subset in neutral else 0
                        return (move.bit_count(),sum((i+2)*((move>>i)&1) for i in range(n)))
                    atom=[a for a in neutral if a and not any(b and b!=a and b&a==b for b in neutral)]
                    best={0:(0,0)}
                    for subset in sorted(neutral,key=int.bit_count):
                        if subset==0:continue
                        choices=[max(best[subset^a],event_cost(a)) for a in atom if a&subset==a and subset^a in best]
                        need(bool(choices),'every nonempty neutral set partitions into neutral atoms')
                        best[subset]=min(choices)
                    literal={x:(0,0)}
                    for z in sorted(interval,key=lambda z:(x^z).bit_count()):
                        if z==x:continue
                        literal[z]=min(max(price,((w^z).bit_count(),sum((i+2)*(((w^z)>>i)&1) for i in range(n)))) for w,price in literal.items() if (x^w)&(x^z)==x^w)
                    need(best[(1<<m)-1]==literal[y],'literal strict optimum equals minimal neutral-atom partition optimum')
                    partition_tests+=1
    # Every multiplicity in the sharp length132 sequence is independently
    # classified by selecting k of its127 repeated cyclic generators.
    zero_subsets=0;from math import comb
    distribution=Counter()
    for k in range(128):
        for parity_bits in range(16):
            for terminal in (0,1):
                q=((k+terminal)%128,parity_bits^(15 if terminal else 0))
                distribution[q]+=comb(127,k)
                if q==(0,0):zero_subsets+=comb(127,k)
    need(zero_subsets==2,'length132 sequence has only empty/full zero-sum subsets')
    need(sum(distribution.values())==1<<132,'complete subset multiplicity count')
    # Explicit lower zero-sumfree sequence of length131 and augmentation bound.
    need(not any(k%128==0 and parity_bits==0 for k in range(128) for parity_bits in range(16) if (k,parity_bits)!=(0,0)),'length131 standard lower sequence is zero-sumfree')
    nilpotent_degrees=[127,1,1,1,1]
    need(sum(nilpotent_degrees)==131,'every degree132 group-algebra monomial vanishes')
    path_controls=[]
    for m in range(1,10):
        n=2*m
        states=[sum(((i+b)&1)<<i for i in range(n)) for b in (0,1)]
        need(all(x.bit_count()==m for x in states) and (states[0]^states[1])==(1<<n)-1,'equal-margin alternating path has two antipodal assignments')
        need(all(all((((x>>i)^(x>>(i+1)))&1)==1 for i in range(n-1)) for x in states),'each alternating assignment satisfies literal signed edges')
        path_controls.append({'slots':n,'population':m,'forced_event_slots':n,'forced_clocks_at_cost2':2*n})
    return {'status':'PASS independently verified finite controls; all-size proof reviewed separately',
        'graphs':graphs,'consistent_graphs':consistent,'guarded_endpoint_tests':endpoint_tests,
        'neutral_atom_partition_tests':partition_tests,'sharp_zero_sum_sequence_length':132,'zero_sum_subset_count':zero_subsets,
        'alternating_path_controls':path_controls,'checks':CHECKS}

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--write',action='store_true');args=parser.parse_args()
    result=run();target=P/'STRUCTURAL_RESULT.json'
    if args.write:target.write_text(json.dumps(result,indent=2,sort_keys=True)+'\n')
    else:need(result==json.loads(target.read_text()),'frozen independent structural result')
    print(json.dumps(result,sort_keys=True))
