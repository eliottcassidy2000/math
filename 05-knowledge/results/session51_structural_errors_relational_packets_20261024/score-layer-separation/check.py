#!/usr/bin/env python3
"""Finite arithmetic/actual-port controls for the all-size score-layer proof."""
from collections import Counter
from itertools import combinations
from math import comb
from pathlib import Path
import argparse,hashlib,json
HERE=Path(__file__).resolve().parent
def need(ok,msg):
    if not ok:raise RuntimeError(msg)
def c3_from_scores(p,a):
    b=2*p-1-a
    return comb(2*p,3)-p*(comb(a,2)+comb(b,2))
def triples(T):
    return sum(((T[a]>>b)&1)==((T[b]>>c)&1)==((T[c]>>a)&1) for a,b,c in combinations(range(len(T)),3))
def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',default=str(HERE/'normal.json'));args=ap.parse_args()
    pin=json.loads((HERE/'SOURCE_PINS.json').read_text())
    need(hashlib.sha256((HERE/'sources/PATH_TRANSPORT_AND_FOREST_RULE.md').read_bytes()).hexdigest()==pin['sha256'],'inherited exact transport source pin')
    for source in pin.get('additional_sources',[]):need(hashlib.sha256((HERE/source['path']).read_bytes()).hexdigest()==source['sha256'],'inherited leg-law source pin')
    arithmetic=0
    for p in range(1,501):
        allowed=range(p//2,p) # ceil((p-1)/2)=floor(p/2)
        if p<=2:need(list(allowed)==[p-1],'small p has only balanced feasible score layer')
        for a in allowed:
            b=2*p-1-a;k=b-a;c=c3_from_scores(p,a);arithmetic+=1
            need(k>0 and k%2==1,'score gap positive odd')
            shifts=(-1,-k-1,k-1)
            need((0 in shifts)==(k==1),'same-layer cycle preservation only balanced')
            if a+1<p:
                cnext=c3_from_scores(p,a+1);knext=k-2
                need(cnext-c==p*(k-1),'adjacent source cycle gap')
                gap=(cnext-knext-1)-(c+k-1)
                need(gap==(p-2)*(k-1)>0,'strict separation of entire neighboring intervals')
    counts=Counter();shifts=Counter();negative_hostile=False
    for n in (2,4,6):
        p=n//2;pairs=list(combinations(range(n),2));full=(1<<n)-1
        for bits in range(1<<len(pairs)):
            T=[0]*n
            for k,(i,j) in enumerate(pairs):
                a,b=(i,j) if bits>>k&1 else (j,i);T[a]|=1<<b
            degrees=[r.bit_count() for r in T];hist=Counter(degrees)
            if len(hist)!=2 or sorted(hist.values())!=[p,p]:continue
            a,b=sorted(hist);gap=b-a;c=triples(T);counts[f'n{n}_labelled_sources']+=1
            need(a+b==n-1 and c==c3_from_scores(p,a),'actual two-level source cycle value')
            target_cycles=set()
            def walk(path,remaining):
                nonlocal negative_hostile
                if remaining:
                    for v in range(n):
                        if remaining>>v&1 and (not path or T[path[-1]]>>v&1):walk(path+(v,),remaining^(1<<v))
                    return
                counts[f'n{n}_actual_ports']+=1
                U=[full^r^(1<<v) for v,r in enumerate(T)]
                for x,y in zip(path,path[1:]):U[x]|=1<<y;U[y]&=~(1<<x)
                ct=triples(U);target_cycles.add(ct)
                ds,dt=degrees[path[0]],degrees[path[-1]]
                need(ct-c==ds-dt-1 and ct-c in (-1,-gap-1,gap-1),'actual retained-path cycle jump')
                q={path[i]:path[n-1-i] for i in range(n)}
                blue=all(((T[i]>>j)&1)==((T[q[j]]>>q[i])&1) for i,j in pairs)
                if blue:
                    counts[f'n{n}_blue_ports']+=1;need(ds+dt==n-1 and ds!=dt,'blue endpoint scores complementary')
                    need(ct-c in (-gap-1,gap-1),'blue two-value necessary list')
                targethist=Counter(r.bit_count() for r in U)
                if len(targethist)==2 and sorted(targethist.values())==[p,p]:
                    need(min(targethist)==a and gap==1 and ds==b and dt==a,'actual union-of-layers adjacency only balanced within layer')
                shifts[n,a,ct-c,blue]+=1
            walk((),full)
            if n==6 and a==1:
                need(c==2 and c-gap-1==-2 and -2 not in target_cycles,'formal allowed negative target is not claimed attained')
                negative_hostile=True
    need(negative_hostile,'unattained-value hostile encountered')
    result={'status':'PASS','scope':'arithmetic p<=500 and every actual retained-path port on labelled two-score equal-multiplicity tournaments at n2,4,6; all-size statement is proved separately','arithmetic_layers_checked':arithmetic,'actual_universes':dict(counts),'actual_shift_histogram':[[list(k),v] for k,v in sorted(shifts.items())],'order_ten_layers':[{'a':a,'b':9-a,'source_c3':c3_from_scores(5,a),'possible_blue_target_c3':[c3_from_scores(5,a)-(9-2*a)-1,c3_from_scores(5,a)+(9-2*a)-1]} for a in (3,4)],'negative_formal_value_hostile':'n6 scores1^3,4^3 have c3=2; formal allowed value -2 is impossible and never claimed attained'}
    Path(args.output).write_text(json.dumps(result,indent=2,sort_keys=True)+'\n');print(json.dumps(result,indent=2,sort_keys=True))
if __name__=='__main__':main()
