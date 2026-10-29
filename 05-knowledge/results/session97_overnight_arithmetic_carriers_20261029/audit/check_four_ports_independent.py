"""Independent explicit-generator verification of the eight-code hostile."""
from itertools import combinations
from collections import Counter
import json

def need(b,m):
    if not b:raise ArithmeticError(m)

def main():
    gs=[75035,110921,73235,83903,92135,102641,118733,129989]
    banks={};tables={};profiles={}
    for g in gs:
        W={0}
        for k in range(15):W|={x^(g<<k) for x in W}
        need(len(W)==32768,'code dimension');banks[g]=W
        extended=W|{x^((1<<32)-1) for x in W}
        octads=[x for x in extended if x.bit_count()==8]
        need(len(octads)==620,'octad count')
        counters={k:Counter() for k in range(1,5)}
        for x in octads:
            S=[i for i in range(32) if x>>i&1]
            for k in range(1,5):counters[k].update(combinations(S,k))
        for k,num,value in [(1,32,155),(2,496,35),(3,4960,7)]:
            need(len(counters[k])==num and set(counters[k].values())=={value},'all named low ports')
        hist=Counter(counters[4].values());hist[0]=35960-len(counters[4])
        profiles[g]=dict(hist);tables[g]=frozenset(counters[4].items())
        need((counters[4][(0,1,2,3)]==0)==(g in [75035,110921]),'zero witness')
    need(len(set(tables.values()))==8,'eight distinct full profiles')
    orbits=[]
    for g in [75035,73235]:
        orbit=set()
        for a in range(1,31):
            moved=[sum(1<<((a*i)%31) for i in range(31) if (g<<k)>>i&1) for k in range(15)]
            hits=[h for h in gs if all(v in banks[h] for v in moved)]
            need(len(hits)==1,'unique multiplier image');orbit.add(hits[0])
        orbits.append(sorted(orbit))
    need(sorted(map(len,orbits))==[2,6] and set(orbits[0]).isdisjoint(orbits[1]),'two multiplier orbits')
    need(profiles[75035]=={0:11160,1:12400,2:7440,3:3720,4:1240},'profileA')
    need(profiles[73235]=={0:0,1:34720,7:1240},'profileB')
    print(json.dumps({'multiplier_orbits':orbits,'full_four_port_tables':8,
        'low_port_values':[155,35,7],'four_port_profiles':[profiles[75035],profiles[73235]],
        'zero_witness':[0,1,2,3]},sort_keys=True))
if __name__=='__main__':main()
