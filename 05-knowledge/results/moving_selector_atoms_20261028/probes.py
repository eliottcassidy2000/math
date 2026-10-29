"""Moving vertex probes and ordered-run attachment controls."""
from pathlib import Path
from collections import defaultdict,Counter
from itertools import permutations,product
import importlib.util,sys,json
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[2]
sys.path.insert(0,str(HERE.parent/'minimal_blue_selectors_20261028'))
import blue_minimum as m
r=m.r;core=m.core;need=m.need
controls=r.old.source('cursor_controls','05-knowledge/results/root_ports_closed_carriers_20261028/check_carriers.py')


def cycle3():return (2,4,1)


def transitive(n):return tuple(sum(1<<j for j in range(i+1,n)) for i in range(n))


def rooted_maps(t,u,roots):
    return [p for p in r.oracle.isomorphisms(t,u,True) if all(p[a]==b for a,b in roots)]


def padded_cursor():
    out=[]
    for k in range(5):
        a,b=(core.presentation(7,x) for x in (1072,4176))
        for side in (True,False):
            for _ in range(k):
                a=core.complement(r.adjoin(core.complement(a),side,False))
                b=core.complement(r.adjoin(core.complement(b),side,False))
        n=len(a);maps=[r.oracle.isomorphisms(t,u,True) for t,u in ((a,b),(core.complement(a),core.complement(b)))]
        scans=[]
        for i in range(n):
            scans.append([not any(p[0]==0 and p[-1]==n-1 and p[i]==i for p in G) for G in maps])
        first=next(i for i,v in enumerate(scans) if any(v))
        need(first==k+1,'moving cursor reaches hidden old second vertex')
        out.append(dict(order=n,padding=k,first_separating_position=first,scans=scans,rows=[a,b]))
    return out


def singleton_run_hostile():
    A,_=core.substitute(cycle3(),(cycle3(),(0,),(0,)))
    words=((3,0,1,2,4),(3,0,2,1,4));G=r.oracle.isomorphisms(A,A,True)
    need(core.is_strong(A) and len(G)==3,'five-vertex strong source')
    for i in range(5):
        need(any(g[3]==3 and g[4]==4 and g[words[0][i]]==words[1][i] for g in G),'all one-cursor rooted probes agree')
    need(not any(all(g[x]==y for x,y in zip(*words)) for g in G),'joint ordered ports differ')
    source,_=core.substitute(cycle3(),(A,transitive(4),transitive(4)))
    paths=[];partners=[]
    for word in words:
        p=[]
        for i,a in enumerate(word):
            p.append(a)
            if i<4:p.extend((5+i,9+i))
        need(all(core.edge(source,x,y) for x,y in zip(p,p[1:])),'actual quotient walk Hamilton path')
        paths.append(tuple(p));partners.append(core.complement(source,p))
    H=[r.oracle.hamilton(t) for t in partners]
    need(H==[core.hamilton(t) for t in partners],'independent forward/backward Hamilton DP')
    need(H==[2483993,2491877],'actual run-order separation')
    canon=controls.Canon();codes=[canon.maps(t)[0] for t in partners]
    strong_source,_=core.substitute(cycle3(),(A,(2,4,9,3),(2,4,9,3)))
    strong_partners=[core.complement(strong_source,p) for p in paths]
    strong_H=[r.oracle.hamilton(t) for t in strong_partners]
    need(strong_H==[core.hamilton(t) for t in strong_partners]==[2899675,2906091],'all-strong children context')
    need(all(core.is_strong(t) for t in (A,(2,4,9,3),strong_source,*strong_partners)),'strong sources and partners')
    return dict(child=A,child_Aut=len(G),ordered_singleton_runs=words,source=source,paths=paths,partners=partners,
                H=H,ordinary_codes=codes,ordinary_classes_differ=codes[0]!=codes[1],
                blue=[r.blue(core.relabel(source,p)) for p in paths],
                single_cursor_profiles_equal=True,full_joint_port_orbits_differ=True,
                all_strong_children=dict(source=strong_source,partners=strong_partners,H=strong_H))


def profile(t):
    canon=controls.Canon();code,orders=canon.maps(t)
    inverses=[]
    for order in orders:
        inv=[0]*len(t)
        for i,v in enumerate(order):inv[v]=i
        inverses.append(inv)
    return code,tuple(min((p[0],p[-1],p[i]) for p in inverses) for i in range(len(t)))


def wreath_presentations():
    A,_=core.substitute(cycle3(),(cycle3(),)*3)
    masks=sorted({core.free_mask(core.relabel(A,p)) for p in core.paths(A)})
    cells=defaultdict(list)
    for x in masks:
        t=core.presentation(9,x)
        cells[profile(t),profile(core.complement(t)),m.state(t)].append(x)
    return dict(order=9,source='C3[C3,C3,C3]',Hamilton_paths=r.oracle.hamilton(A),
                presentations=len(masks),cursor_profile_cells=len(cells),
                collisions=[v for v in cells.values() if len(v)>1])


def main():
    result=dict(status='FINITE-EXACT',padded=padded_cursor(),run_hostile=singleton_run_hostile(),wreath_scan=wreath_presentations())
    (HERE/'probe-receipt.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(cursor=[(v['order'],v['first_separating_position']) for v in result['padded']],
                         run_partner_H=result['run_hostile']['H'],run_class_difference=result['run_hostile']['ordinary_classes_differ'],
                         wreath=result['wreath_scan'])));return result


if __name__=='__main__':main()
