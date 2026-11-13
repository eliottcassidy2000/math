"""Count-decreasing triangle grammar, with complete small-parent source audit."""
from pathlib import Path
from functools import lru_cache
import json,gzip,sys
from triangle_probe import bank,need,save,HERE,scan,attach
import numpy as np
sys.path.insert(0,str(HERE.parent/'modular_centers_collatz_pell_20261111'))
import content_checks as cc

def induced(T,vertices):
    return tuple(sum(((T[v]>>w)&1)<<j for j,w in enumerate(vertices)) for v in vertices)

def strong(T):
    n=len(T);reach=[T[v]|1<<v for v in range(n)]
    for k in range(n):
        for v in range(n):
            if reach[v]>>k&1:reach[v]|=reach[k]
    return all(x==(1<<n)-1 for x in reach)

def near(T):
    n=len(T)
    for a in range(n):
        for b in range(n):
            if not(T[a]>>b&1):continue
            U=list(T);U[a]^=1<<b;U[b]^=1<<a
            if sorted(v.bit_count() for v in U)==list(range(n)):return (a,b)
    return None

@lru_cache(None)
def parse(n,code):
    T=tuple(cc.decode(n,code));edge=near(T)
    if edge is not None:return dict(kind='near',n=n,code=code,back_edge=edge)
    for a in range(n):
        for b in range(a+1,n):
            keep=[v for v in range(n) if v not in (a,b)];U=induced(T,keep)
            if not strong(U):continue
            parent=parse(n-2,cc.canonical(U))
            if parent is None:continue
            for v in keep:
                if (T[a]>>b&1 and T[b]>>v&1 and T[v]>>a&1) or (T[b]>>a&1 and T[a]>>v&1 and T[v]>>b&1):
                    return dict(kind='triangle',n=n,code=code,deleted=[a,b],shared=v,parent=parent)
    return None

def run():
    data=json.loads(gzip.decompress((HERE/'triangle_probe.json.gz').read_bytes()))
    ordinary={};fibres={};witnesses={};rejected={}
    for rec in data['records']:
        n=rec['n'];code=rec['ordinary_code'];h=rec['core_H'];tree=parse(n,code)
        ordinary[n,code]=dict(n=n,code=code,H=h,accepted=tree is not None,tree=tree)
        if tree is None:rejected.setdefault(h,[]).append((n,code));continue
        for value in map(int,rec['spectrum']):
            item=dict(parent_n=n,parent_code=code,parent_rows=rec['rows'],parent_H=h,target=value,excess=value-3*h)
            if value not in witnesses or (h,n,code)<(witnesses[value]['parent_H'],witnesses[value]['parent_n'],witnesses[value]['parent_code']):witnesses[value]=item
    # The complete strong floor at order8 is45; hence all parents of a
    # triangle output below135 occur in this order<=7 parent bank.
    limit=134
    missing=[v for v in range(9,limit+1,2) if v not in witnesses]
    levels=[]
    for n in range(3,8):
        rr=[x for x in ordinary.values() if x['n']==n]
        levels.append(dict(n=n,strong_classes=len(rr),accepted=sum(x['accepted'] for x in rr),rejected=sum(not x['accepted'] for x in rr)))
    # A wider finite positive list; no all-order exclusion is claimed there.
    first_gap=next(v for v in range(9,max(witnesses)+3,2) if v!=21 and v not in witnesses)
    selected={str(v):witnesses[v] for v in sorted(witnesses) if v<=2881}
    masks={}
    for item in selected.values():
        T=tuple(item['parent_rows'])
        if T not in masks:
            values=scan(np.array(T,dtype=np.int64),len(T)-1);first={}
            for mask,h in enumerate(values):first.setdefault(int(h),mask)
            masks[T]=first
        item['mask']=masks[T][item['target']]
        item['child_rows']=list(map(int,attach(np.array(T,dtype=np.int64),len(T)-1,item['mask'])))
    save('peeling',dict(scope='Exact recursive parser on all strong classes through7; all triangle exits of accepted parents; complete negative scope9..134 by the Busch bound',levels=levels,classes=list(ordinary.values()),all_order_missing_9_134=missing,first_finite_gap_except21=first_gap,witnesses=selected))
    print('LEVELS',levels)
    print('ALL-ORDER MISSING9..134',missing,'first finite gap except21',first_gap)
    print('SAMPLE PARENTS',[(v,witnesses[v]['parent_H'],witnesses[v]['parent_n'],witnesses[v]['excess']) for v in (23,25,33,35,39,47,49,55,63,65,67,881) if v in witnesses])

if __name__=='__main__':run()
