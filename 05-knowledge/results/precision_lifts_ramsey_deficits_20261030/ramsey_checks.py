#!/usr/bin/env python3
"""Exact local-port census and an explicit obstruction to the literal R366 target."""
from pathlib import Path
from collections import Counter
from itertools import combinations
import gzip,hashlib,json
P=Path(__file__).parent
def need(x,m):
    if not x:raise RuntimeError(m)
def decode(line):
    n=line[0]-63;need(0<=n<63,'small graph6')
    a=[0]*n;k=0
    for j in range(n):
        for i in range(j):
            if ((line[1+k//6]-63)>>(5-k%6))&1:a[i]|=1<<j;a[j]|=1<<i
            k+=1
    return a
def complement(a):
    allbits=(1<<len(a))-1
    return [allbits^(1<<i)^s for i,s in enumerate(a)]
def clique_count(a,k,candidates=None):
    if candidates is None:candidates=(1<<len(a))-1
    if k==0:return 1
    if candidates.bit_count()<k:return 0
    if k==1:return candidates.bit_count()
    count=0
    while candidates:
        bit=candidates&-candidates;candidates-=bit;i=bit.bit_length()-1
        count+=clique_count(a,k-1,candidates&a[i])
    return count
def no_four(a,subset=None):
    s=(1<<len(a))-1 if subset is None else subset
    while s:
        b=s&-s;s-=b;i=b.bit_length()-1;js=s&a[i]
        while js:
            b=js&-js;js-=b;j=b.bit_length()-1;rest=js&a[j]
            while rest:
                b=rest&-rest;rest-=b
                if rest&a[b.bit_length()-1]:return False
    return True
def independent_no_six(a):
    for i in range(len(a)):
        for j in range(i):
            if a[i]>>j&1:
                if not no_four(a,a[i]&a[j]):return False
    return True
def run():
    # Independent clique procedures and self-complementing symmetry.
    p=101;qr={x*x%p for x in range(1,p)}
    a=[sum(1<<j for j in range(p) if (j-i)%p in qr) for i in range(p)]
    b=complement(a);nonres=next(i for i in range(1,p) if i not in qr)
    need(all(bool(a[i]>>j&1)!=bool(a[(nonres*i)%p]>>((nonres*j)%p)&1) for i in range(p) for j in range(i)),'complement isomorphism')
    control=[]
    for graph in (a,b):
        c5=clique_count(graph,5);c6=clique_count(graph,6)
        need(c5>0 and c6==0 and independent_no_six(graph),'Paley101 clique audit')
        control.append(dict(K5=c5,K6=c6))
    # R344: with n=73, red degrees are14..17. Check every supplied R44 port.
    banks=[]
    for n,expected in ((14,130816),(15,640),(16,2),(17,1)):
        f=P/'data'/('r44_%d.g6'%n+('.gz' if n==14 else ''));raw=f.read_bytes()
        data=gzip.decompress(raw) if n==14 else raw
        lines=data.splitlines();need(len(lines)==expected,'catalog count')
        profiles=Counter();edges=Counter();triangles=Counter()
        for line in lines:
            graph=decode(line);comp=complement(graph)
            need(len(graph)==n and no_four(graph) and no_four(comp),'catalog R44 predicate')
            ds=tuple(sorted(x.bit_count() for x in graph));ec=sum(ds)//2
            profiles[ds]+=1;edges[ec]+=1
            # The number of monochromatic triangles is Goodman's exact identity.
            mono=(n*(n-1)*(n-2)//6)-sum(d*(n-1-d) for d in ds)//2
            need(mono>=0,'triangle inventory')
            triangles[mono]+=1
        banks.append(dict(order=n,graphs=expected,sha256=hashlib.sha256(raw).hexdigest(),
                          source='https://users.cecs.anu.edu.au/~bdm/data/'+f.name,
                          edge_counts=sorted(edges.items()),degree_profiles=[dict(degrees=k,count=v) for k,v in sorted(profiles.items())],
                          monochromatic_triangle_counts=sorted(triangles.items())))
        print('R44 port',n,expected,'profiles',len(profiles),'edge range',min(edges),max(edges),flush=True)
    # Exact inventory relaxation survivor at n=73. This is NOT a coloring.
    n=73;degrees=(16,28,28);port_edges=((0,60,60),(96,102,180),(96,180,102))
    need(sum(degrees)==n-1 and all(n*d%2==0 for d in degrees),'degree parity')
    need(all(sum(port_edges[i])==degrees[i]*(degrees[i]-1)//2 for i in range(3)),'port total')
    need(all(n*port_edges[i][i]%3==0 for i in range(3)),'triangle divisibility')
    for i in (1,2):
        # Internal R334 cap vector on the blue/green port, in global color order.
        cap=(8,8,16) if i==1 else (8,16,8)
        need(all(2*port_edges[i][j]<=degrees[i]*cap[j] for j in range(3)),'internal degree cap')
    mono=sum(n*port_edges[i][i]//3 for i in range(3))
    repeated=sum(n*port_edges[i][j] for i in range(3) for j in range(3) if i!=j)
    total=n*(n-1)*(n-2)//6;rainbow=total-mono-repeated
    mixed=n*sum(degrees[i]*degrees[j] for i in range(3) for j in range(i))
    need(rainbow>=0 and 2*repeated+3*rainbow==mixed,'global three-color triangle identity')
    out=dict(status='FINITE-EXACT; catalog completeness CITED; requested upper bounds OPEN or REFUTED as stated',
             paley101=dict(quadratic_residues=sorted(qr),nonresidue=nonres,color_cliques=control,
                           doubled_coloring_vertices=202,proved_lower_bound_R366=203),
             R344_red_port_banks=banks,
             inventory_relaxation_survivor=dict(n=n,degrees=degrees,port_edge_rows=port_edges,
                                               mono_triangles=mono,two_color_triangles=repeated,rainbow_triangles=rainbow))
    (P/'ramsey_checks.json').write_text(json.dumps(out,indent=2)+'\n')
    print('Paley101',control,'inventory relaxation',out['inventory_relaxation_survivor'])
if __name__=='__main__':run()
