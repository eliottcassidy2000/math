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
def ordered(blocks):return tuple(sorted((tuple(sorted(b)) for b in blocks),key=lambda b:b[0]))
def partition(word):return ordered([{i for i,c in enumerate(word) if c==a} for a in (1,2,3)])
def extend(P,i,j):
    b=sum(map(len,P));B=[set(s) for s in P]
    if i!=j:
        a=min(B[i]|B[j]);source=i if a in B[i] else j;target=j if source==i else i
        B[source].remove(a);B[target].add(a)
    B[i].add(b);B[j].add(b+1)
    return ordered(B)
def decode(P):
    b=sum(map(len,P));r,s=b-2,b-1;B=[set(a) for a in P]
    i=next(j for j,a in enumerate(B) if r in a);j=next(j for j,a in enumerate(B) if s in a)
    B[i].remove(r);B[j].remove(s)
    if i!=j and not(B[i]|B[j]):return None
    if i!=j:
        a=min(B[i]|B[j]);source=i if a in B[i] else j;target=j if source==i else i
        B[source].remove(a);B[target].add(a)
    Q=ordered(B);return Q,Q.index(tuple(sorted(B[i]))),Q.index(tuple(sorted(B[j])))
def encode(P):
    m=(sum(map(len,P))-1)//2
    if m==1:return (1,)
    data=decode(P)
    if data is None:return (1,)+(0,)*(m-1)
    Q,i,j=data;return (i+3*j,)+encode(Q)
def from_code(code):
    if len(code)==1:
        need(code==(1,),'normalized projective base');return ((0,),(1,),(2,))
    if not any(code[1:]):
        need(code[0]==1,'exception normalization');b=2*len(code)-1
        return (tuple(range(b)),(b,),(b+1,))
    Q=from_code(code[1:]);a=code[0];return extend(Q,a%3,a//3)
def main():
    results=[];previous={((0,),(1,),(2,))};codec_tests=0
    for b in (5,7,9):
        count=0;orbits=set();patterns=Counter()
        for word in product((1,2,3),repeat=b):
            if reduce(xor,word,0):continue
            mult=tuple(word.count(i) for i in (1,2,3))
            need(all(m%2 for m in mult),'three odd colour classes')
            count+=1;orbits.add(canonical(word));patterns[tuple(sorted(mult,reverse=True))]+=1
        need(count==6*len(orbits) and len(orbits)==(3**(b-1)-1)//8,'free kernel count')
        exact={partition(w) for w in orbits};generated=set()
        for P in previous:
            for i,j in product(range(3),repeat=2):
                Q=extend(P,i,j)
                need(decode(Q)==(P,i,j),'inverse extension with actual old port')
                need(Q not in generated,'disjoint nine extensions');generated.add(Q)
        exceptional=(tuple(range(b-2)),(b-2,),(b-1,))
        need(exceptional not in generated and decode(exceptional) is None,'unique exceptional partition')
        generated.add(exceptional);need(generated==exact,'complete extension coverage')
        codes=set()
        for P in exact:
            c=encode(P);need(from_code(c)==P,'projective address decoder');codes.add(c);codec_tests+=1
        m=(b-1)//2
        expected={v for v in product(range(9),repeat=m) if any(v) and next(x for x in reversed(v) if x)==1}
        need(codes==expected,'all normalized F9 projective points')
        previous=exact
        results.append({'ports':b,'ordered_assignments':count,'free_kernel_count':len(orbits),
                        'patterns':{str(k):v//6 for k,v in sorted(patterns.items())}})
    hostiles=[]
    for b in (5,6,7):
        ports=tuple(1<<i for i in range(b-1))+((1<<(b-1))-1,)
        h=ports[0]^ports[1]^ports[2]^ports[3]
        sums={a^c for i,a in enumerate(ports) for c in ports[i+1:]}
        need((h in ports)==(b==5) and (h in sums)==(b==6),'Petersen five-port boundary')
        hostiles.append({'ports':b,'disjoint_pair_sum_is_inertia':h in ports,'is_pair_sum':h in sums})
    print(json.dumps({'status':'FINITE-EXACT PASS','counts':results,'disjoint_pair_hostiles':hostiles,
        'exact_codec_roundtrips':codec_tests,'recursive_kernel_rule':'A_(m+1)=9 A_m+1; named least-port transfer',
        'projective_incidence_claimed':False},indent=2,sort_keys=True))
if __name__=='__main__':main()
