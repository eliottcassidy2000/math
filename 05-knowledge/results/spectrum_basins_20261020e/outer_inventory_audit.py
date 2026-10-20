"""Independent audit of the parent's cross-feedback four-R identity.

All fixed-path A presentations through order six. The general adjacency
identity and middle score sums are checked for every A; four-blackness is
claimed/tested only for A whose first path vertex is a source.
"""
from collections import Counter
from pathlib import Path
import json
from outer_probe import G, fixed_paths


def reverse_path(a,p):
    G.need(G.is_path(a,p),'not a Hamiltonian path')
    b=list(a)
    for u,v in zip(p,p[1:]):
        b[u]^=1<<v;b[v]^=1<<u
    return tuple(b)


def main():
    counts=Counter();controls=[]
    for m in range(1,7):
        for A in fixed_paths(m):
            p=tuple(range(m));x,y,z,a,b,c=range(m,m+6)
            source=G.join([A,(0,),(0,),(0,),(2,4,1)])
            words=[p+(x,y,z,b,c,a), (a,b,z,y,x)+p[::-1]+(c,),
                   (x,b,a,c)+p+(y,z), (c,a,b,x,z,y)+p[::-1]]
            current=source;sums=[];colours=[];pair_sums=[]
            for q in words:
                G.need(G.is_path(current,q),'literal directedness')
                sums.append(current[q[0]].bit_count()+current[q[-1]].bit_count())
                pair_sums.append(current[q[1]].bit_count()+current[q[-2]].bit_count())
                colours.append(G.blue(current,q))
                current=reverse_path(current,q)
            G.need(current==G.join([A,(2,4,1),(4,1,2)]),'literal target')
            G.need(sums[1:3]==[3,8],'middle endpoint score sums')
            G.need(pair_sums[2]==5,'third move interior score witness')
            counts['all_presentations']+=1
            if A[0].bit_count()==m-1:
                G.need(not any(colours),'source-guard blackness')
                counts['first_vertex_source']+=1
            if m==3 and A==(6,4,0):
                controls.append(dict(A=A,paths=words,colours=colours,
                                     endpoint_sums=sums,second_penultimate_sums=pair_sums))
    result=dict(universe='all fixed-path A presentations of orders 1..6',
                counts=dict(counts),special_order_three_control=controls)
    Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2)+'\n')
    print(result)


if __name__=='__main__':main()
