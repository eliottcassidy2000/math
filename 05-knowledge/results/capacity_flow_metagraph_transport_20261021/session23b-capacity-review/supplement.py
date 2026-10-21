#!/usr/bin/env python3
"""Independent checks added for the producer's final frozen refinements."""
import json
from collections import Counter
from itertools import product
from audit import check, pairs, normalized_regular, tree_flows, score, hamiltonian_paths, automorphisms


def main():
    hists=sorted([sorted(Counter(map(sum, tree_flows(5,A,1,zero_path=True))).items()) for A in normalized_regular(5)])
    check(hists==sorted([[(0,1),(3,1),(5,1)],[(0,1),(4,1),(5,1)],[(0,1),(3,1),(4,1)]]),'three n5 use histograms')
    scores={}
    for n in (5,7):
        A=pairs(n); free=[i for i,(u,v) in enumerate(A) if v!=u+1]
        r=(n-1)//2; target=(r-1,)+(r,)*(n-2)+(r+1,); count=0
        for bits in product((0,1),repeat=len(free)):
            T=A[:]
            for i,bit in zip(free,bits):
                if bit:T[i]=T[i][::-1]
            count+=score(n,T)==target
        scores[str(n)]=count
    check(scores=={'5':1,'7':34},'fixed-path endpoint-shift score fibre')
    Hhist=Counter(); reps={}
    for A in normalized_regular(7):
        H=hamiltonian_paths(7,A); Hhist[H]+=1; reps.setdefault(H,A)
    check(Hhist==Counter({171:57,175:25,189:9}),'H distribution')
    auts={H:automorphisms(7,A) for H,A in reps.items()}
    check(auts=={171:3,175:7,189:21},'aut orders')
    check(sum(H//a for H,a in auts.items())==91,'H/Aut sum')
    return {'status':'PASS','n5_M1_use_histograms':hists,'endpoint_shift_score_fibres':scores,'n7_H_histogram':dict(Hhist),'representative_automorphism_orders':auts,'H_over_Aut_sum':91}

if __name__=='__main__':
    print(json.dumps(main(),indent=2,sort_keys=True))
