#!/usr/bin/env python3
"""Small hostile controls for the separately proved all-mechanical corollary."""
from collections import Counter
from hashlib import sha256
from math import gcd
from pathlib import Path
import json

checks=0
def need(value,reason):
    global checks
    checks+=1
    if not value:
        raise RuntimeError(reason)

def carry(bits):
    B=0
    for i,bit in enumerate(bits):
        B=(3*B+(1<<i)) if bit else B
    return B

def main():
    shapes=0
    counts=Counter()
    hits=[]
    for q in range(2,65):
        for p in range(1,q):
            shapes+=1
            source=[((i+1)*p)//q-(i*p)//q for i in range(q)]
            labels=[i if source[i] else None for i in range(q)]
            D=(1<<q)-3**p
            d=gcd(p,q)
            for i in range(q):
                a,b=i,(i+1)%q
                if source[a]==source[b]:
                    continue
                first_set={a,b}
                first_label=a if source[a] else b
                middle=labels[:]
                middle[a],middle[b]=middle[b],middle[a]
                for j in sorted({i,(i-1)%q,(i+1)%q}):
                    c,e=j,(j+1)%q
                    if (middle[c] is None)==(middle[e] is None):
                        continue
                    second_label=middle[c] if middle[c] is not None else middle[e]
                    final=middle[:]
                    final[c],final[e]=final[e],final[c]
                    bits=[int(x is not None) for x in final]
                    if first_set=={c,e}:
                        need(bits==source,'same physical edge not undo')
                        counts['undo']+=1
                        continue
                    need(bits!=source,'unclassified nonedge undo')
                    if first_label==second_label:
                        counts['single_one']+=1
                        need(q>2*p,'single-one gap does not force q>2p')
                    else:
                        counts['two_ones']+=1
                        need(q<2*p,'two-one gap does not force q<2p')
                    if d>1:
                        counts['nonprimitive_nonundo']+=1
                    B=carry(bits)
                    if B%abs(D)==0:
                        need((p,q)==(2,3),'unexpected integral overlapping output')
                        hits.append({'p':p,'q':q,'edges':[i,j],'word':''.join(map(str,bits)),
                                     'fixed_point':B//D})
    need(len(hits)==2,'exception count')
    need(shapes==2016,'shape coverage')
    print(json.dumps({'status':'FINITE-EXACT structural hostile controls only',
        'script_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
        'q_range_inclusive':[2,64],'shape_count':shapes,'checks':checks,
        'counts':dict(counts),'integral_nonundo_outputs':hits},indent=2,sort_keys=True))
if __name__=='__main__':
    main()
