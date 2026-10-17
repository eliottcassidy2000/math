#!/usr/bin/env python3
"""Hostile probe: does a high-degree spindle product actually mix addresses?"""
import json
from fractions import Fraction as Q
from itertools import product
import numpy as np
from escaped_ports_20261017d import make,need,P


def main():
    f,z,j,alpha,ports,units,eta,conj,escape,reduce=make('sqrt3')
    rho=f.scale(f.add(f.scalar(5),j),Q(1,6));a=f.add(f.one,eta)
    M=[f.zero,f.one,eta,a,rho,f.mul(rho,eta),f.mul(rho,a)]
    out=[]
    for name,U,V in [('mixed_golden',units[0],units[1]),('cubic_square',units[0],f.power(units[0],2))]:
        address=list(product(range(7),repeat=3))
        points=[f.add(M[a],f.add(f.mul(U,M[b]),f.mul(V,M[c]))) for a,b,c in address]
        need(len(set(points))==343,'all three source coordinates retained injectively')
        r=np.array([reduce(p) for p in points],dtype=np.int64)
        counts=[0,0,0,0];mixed=[];edges=[]
        for a in range(343):
            delta=r[a+1:]-r[a]
            for off in np.flatnonzero((delta[:,0]**2+delta[:,1]**2)%P==1):
                b=a+1+int(off);v=f.sub(points[a],points[b])
                if f.mul(v,f.apply(v,conj))==f.one:
                    h=sum(x!=y for x,y in zip(address[a],address[b]));counts[h]+=1;edges.append([a,b])
                    if h>1:mixed.append([address[a],address[b]])
        out.append(dict(name=name,vertices=343,edge_count=len(edges),edges_by_changed_coordinates=counts,
                        mixed_edges=mixed,edges=edges))
    return dict(status='FINITE-EXACT induced origin-block probes; not a boundary forcing test',runs=out)


if __name__=='__main__':print(json.dumps(main(),sort_keys=True,separators=(',',':')))
