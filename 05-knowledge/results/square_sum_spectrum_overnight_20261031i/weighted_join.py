#!/usr/bin/env python3
"""Actual signed Hamilton-polynomial product on a sixteen-cycle fibre."""
from itertools import product
import json
from elliptic_ports import Atlas,from_cube,join
from one_copy_ports import edge,need

def weight(es,weights):
    ans=1
    for e in es:ans*=weights.get(e,1)
    return ans

def main():
    atlas=Atlas();A=from_cube(atlas.cycle(2812,'A'));B=from_cube(atlas.cycle(2861,'A'),(196,2304))
    C,r=join(A,B);a,b=r['source_scales'];wa={};wb={}
    for graph,weights in [(A,wa),(B,wb)]:
        for R,S in graph.moves[:2]:
            for e in R|S:weights[e]=1+(e[0]+3*e[1])%4
    source_sums=[];codes=[''.join(z)+'0'*23 for z in product('01',repeat=2)]
    for graph,weights in [(A,wa),(B,wb)]:
        source_sums.append(sum(weight(graph.edges(code)-{graph.port},weights) for code in codes))
    weights={edge(a*u,a*v):w for (u,v),w in wa.items()}
    weights.update({edge(b*u,b*v):w for (u,v),w in wb.items()})
    seams=[tuple(e) for e in r['seams']];weights[seams[0]]=2;weights[seams[1]]=-3
    total=sum(weight(C.edges(left+right),weights) for left in codes for right in codes)
    need(total==-6*source_sums[0]*source_sums[1],'exact signed source-weight product')
    print(json.dumps({'status':'passed','complete_cycle_fibre':16,'source_path_sums':source_sums,
          'seam_weights':[2,-3],'joined_signed_sum':total,
          'scope':'The sixteen actual edge-set products agree with the source-preserving polynomial identity; no positivity or averaging assumption is used.'},indent=2,sort_keys=True))

if __name__=='__main__':main()
