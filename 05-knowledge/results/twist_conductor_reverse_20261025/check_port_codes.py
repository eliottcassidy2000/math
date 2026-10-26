#!/usr/bin/env python3
from collections import Counter
from functools import reduce
from itertools import product,permutations
from operator import xor
import json

def need(ok,msg):
    if not ok:raise RuntimeError(msg)
def canonical(word):
    return min(tuple(p[a-1] for a in word) for p in permutations((1,2,3)))
def main():
    results=[]
    for b in (5,7,9):
        count=0;orbits=set();patterns=Counter()
        for word in product((1,2,3),repeat=b):
            if reduce(xor,word,0):continue
            mult=tuple(word.count(i) for i in (1,2,3))
            need(all(m%2 for m in mult),'three odd colour classes')
            count+=1;orbits.add(canonical(word));patterns[tuple(sorted(mult,reverse=True))]+=1
        need(count==6*len(orbits) and len(orbits)==(3**(b-1)-1)//8,'free kernel count')
        results.append({'ports':b,'ordered_assignments':count,'free_kernel_count':len(orbits),
                        'patterns':{str(k):v//6 for k,v in sorted(patterns.items())}})
    hostiles=[]
    for b in (5,6,7):
        ports=tuple(1<<i for i in range(b-1))+((1<<(b-1))-1,)
        h=ports[0]^ports[1]^ports[2]^ports[3]
        sums={a^c for i,a in enumerate(ports) for c in ports[i+1:]}
        need((h in ports)==(b==5) and (h in sums)==(b==6),'Petersen five-port boundary')
        hostiles.append({'ports':b,'disjoint_pair_sum_is_inertia':h in ports,'is_pair_sum':h in sums})
    print(json.dumps({'status':'FINITE-EXACT PASS','counts':results,'disjoint_pair_hostiles':hostiles},indent=2,sort_keys=True))
if __name__=='__main__':main()
