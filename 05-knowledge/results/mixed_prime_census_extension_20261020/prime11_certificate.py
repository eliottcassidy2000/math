#!/usr/bin/env python3
"""Independent bounded certificate for the elementary degree-eleven proof.
Exhausts 32 skew connection sets and only 5! permutations of each out-neighborhood.
The full vertex stabilizer is reconstructed from distinct in-neighbor signatures.
"""
from itertools import product,permutations
from collections import Counter
import json
CHECKS=0

def need(x,msg):
    global CHECKS
    CHECKS+=1
    if not x:raise RuntimeError(msg)

def main():
    p=11;sets=[frozenset(i if b else p-i for i,b in zip(range(1,6),bits)) for bits in product((0,1),repeat=5)]
    rows=[]
    for S in sets:
        out=sorted(S);inside=sorted(set(range(1,p))-S)
        A=[[int((v-u)%p in S) for v in range(p)] for u in range(p)]
        signatures={tuple(A[u][v] for v in out):u for u in inside}
        need(len(signatures)==5,'in-neighbor signatures not distinct')
        stabilizer=[];neighbor_automorphisms=[]
        for images in permutations(out):
            q=dict(zip(out,images))
            if not all(A[u][v]==A[q[u]][q[v]] for u in out for v in out):continue
            neighbor_automorphisms.append(images)
            f={0:0,**q}
            for u in inside:
                # Required target signature is A[u,v] at target q(v).
                sig=tuple(A[u][out[images.index(v)]] for v in out)
                if sig not in signatures:break
                f[u]=signatures[sig]
            if len(f)==p and len(set(f.values()))==p and all(A[u][v]==A[f[u]][f[v]] for u in range(p) for v in range(p)):
                stabilizer.append(tuple(f[i] for i in range(p)))
        expected=5 if S in (frozenset((1,3,4,5,9)),frozenset((2,6,7,8,10))) else 1
        need(len(stabilizer)==expected,'wrong degree-eleven stabilizer')
        need(len(neighbor_automorphisms)==expected,'neighborhood bound not sharp')
        for q in stabilizer:need(all(q[i]==q[1]*i%p for i in range(p)),'non-multiplier stabilizer')
        orbit={frozenset(a*x%p for x in S) for a in range(1,p)}
        rows.append({'connection_set':sorted(S),'multiplier_orbit_representative':min(map(sorted,orbit)),'multiplier_orbit_size':len(orbit),'out_neighbor_score_sequence':sorted(sum(A[u][v] for v in out) for u in out),'out_neighbor_matrix':[[A[u][v] for v in out] for u in out],'in_neighbor_signatures':[[A[u][v] for v in out] for u in inside],'vertex_stabilizer':stabilizer,'full_automorphism_order':11*len(stabilizer)})
    reps={tuple(r['multiplier_orbit_representative']) for r in rows}
    need(len(reps)==4,'wrong representative count')
    need(sorted(Counter(r['full_automorphism_order'] for r in rows).items())==[(11,30),(55,2)],'wrong order histogram')
    print(json.dumps({'status':'PASS','checks':CHECKS,'all_skew_connection_sets':32,'multiplier_representatives':[list(x) for x in sorted(reps)],'automorphism_order_histogram':dict(sorted(Counter(r['full_automorphism_order'] for r in rows).items())),'rows':rows},sort_keys=True,indent=2))
if __name__=='__main__':main()
