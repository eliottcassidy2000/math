"""All labelled rooted tournaments through order5; explicit checks survive -O."""
from itertools import combinations
from pathlib import Path
import hashlib
import json
import sys

HERE=Path(__file__).resolve().parent
RESULTS=HERE.parent.parent

def need(ok,msg):
    if not ok: raise RuntimeError(msg)

def decode(n,code):
    rows=[0]*n
    for k,(i,j) in enumerate(combinations(range(n),2)):
        u,v=(i,j) if code>>k&1 else(j,i)
        rows[u]|=1<<v
    return rows

def path_table(rows):
    n=len(rows)
    end=[[0]*n for _ in range(1<<n)]
    h=[1]+[0]*((1<<n)-1)
    for mask in range(1,1<<n):
        for v in range(n):
            if not mask>>v&1: continue
            rest=mask^(1<<v)
            end[mask][v]=sum(end[rest][u] for u in range(n) if rows[u]>>v&1) if rest else 1
        h[mask]=sum(end[mask])
    return end,h

def ending_cover(end,h,mask,r):
    value=0
    part=mask
    while part:
        if part>>r&1: value+=end[part][r]*h[mask^part]
        part=(part-1)&mask
    return value

totals=dict(tournaments=0,rooted_tournaments=0,positive_start_restrictions=0,
            unit_ready_roots=0,unit_ready_restrictions=0,subdoubling_subsets=0)
for n in range(1,6):
    full=(1<<n)-1
    for code in range(1<<(n*(n-1)//2)):
        rows=decode(n,code)
        end,h=path_table(rows)
        start,_=path_table([full^(1<<v)^rows[v] for v in range(n)])
        totals['tournaments']+=1
        for r in range(n):
            totals['rooted_tournaments']+=1
            if start[full][r]>0:
                for mask in range(1,full+1):
                    if mask>>r&1:
                        need(start[mask][r]<=start[full][r],('start monotone',n,code,r,mask))
                        totals['positive_start_restrictions']+=1
            F=ending_cover(end,h,full,r)
            if F-h[full]!=1: continue
            totals['unit_ready_roots']+=1
            pred=[v for v in range(n) if rows[v]>>r&1]
            need(len(pred)==1,'unit indegree')
            x=pred[0]
            need(start[full^(1<<r)][x]==1,'unique-start inverse')
            for mask in range(1,full+1):
                if mask>>r&1:
                    f=ending_cover(end,h,mask,r)
                    b=f-h[mask]
                    need(b in (0,1),('two-type menu',n,code,r,mask,b))
                    need(h[mask^(1<<r)]-start[mask][r]==b,'deletion-start response')
                    totals['unit_ready_restrictions']+=1
                if h[full]<2*h[mask]-1:
                    need(mask>>r&1 and mask>>x&1,'subdoubling inherited ports')
                    need(h[mask^(1<<r)]<=h[full^(1<<r)]<h[mask],'subdoubling deletion budget')
                    totals['subdoubling_subsets']+=1
end,h=path_table([2,0])
start,_=path_table([0,1])
need(start[3][1]==0 and start[2][1]==1,'positivity guard hostile')
source_menus=[]
for name,rows in [('A2',[518,388,456,1011,487,455,387,257,1,502]),
                  ('B9',[354,228,425,179,231,192,140,1,250])]:
    n=len(rows); full=(1<<n)-1
    end,h=path_table(rows)
    start,_=path_table([full^(1<<v)^rows[v] for v in range(n)])
    need(h[full]==881,'literal source count')
    for origin in range(n):
        reach=1<<origin
        while True:
            enlarged=reach
            for v in range(n):
                if reach>>v&1: enlarged|=rows[v]
            if enlarged==reach: break
            reach=enlarged
        need(reach==full,('strong literal source',name,origin))
    states=[]
    for r in range(n):
        f=ending_cover(end,h,full,r)
        need(f-h[full] not in (0,1),('forbidden source menu',name,r))
        states.append({'root':r,'state':[h[full],start[full][r],h[full^(1<<r)],f],
                       'readiness_defect':f-h[full]})
    source_menus.append({'name':name,'rows':rows,'states':states,
                         'strong':True,
                         'minimum_starts':min(s['state'][1] for s in states),
                         'minimum_readiness_defect':min(s['readiness_defect'] for s in states),
                         'forbidden_ready_interval':[881,1761],
                         'interval_upper_endpoint_is_open':True})
    need(min(s['state'][1] for s in states)>=2,'all-additions start barrier')
    need(min(s['readiness_defect'] for s in states)>=2,'all-additions old-root barrier')
sources=['chanpak_five_core_bridges_20261110/alternatives/NOTE.md',
         'lambda_emergent_20261103/READY_REPAIR.md','two_entry_funnels_20261103/READY.md',
         'changing_core_returns_20261113/preparation/NEUTRAL_881.md']
pins=[{'path':'05-knowledge/results/'+p,'sha256':hashlib.sha256((RESULTS/p).read_bytes()).hexdigest()} for p in sources]
out={'status':'FINITE-EXACT all labelled rooted tournaments of orders1..5',
     'counts':totals,'source_pins':pins,'source881_menus':source_menus,
     'positivity_guard_hostile':{'rows':[2,0],'root':1,'full_starts':0,'singleton_starts':1}}
suffix='.optimized' if sys.flags.optimize else ''
(HERE/('rooted_heredity'+suffix+'.json')).write_text(json.dumps(out,indent=2)+'\n')
print(json.dumps(totals))
