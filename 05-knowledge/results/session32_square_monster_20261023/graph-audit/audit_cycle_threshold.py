#!/usr/bin/env python3
"""Independent degree-constraint audit of G31/G32 cycle covers.
No subset-size loop, producer imports, external calls, or Git operations.
"""
from collections import Counter
from math import isqrt
from pathlib import Path
import hashlib,json

checks=0
def need(ok,why):
    global checks
    checks+=1
    if not ok:raise RuntimeError(why)

def edges(n):
    return sorted((x,q*q-x) for q in range(2,isqrt(2*n-1)+1)
                  for x in range(max(1,q*q-n),(q*q-1)//2+1)
                  if q*q-x<=n)

def components(n,E):
    p=list(range(n+1))
    def find(x):
        while x!=p[x]:p[x]=p[p[x]];x=p[x]
        return x
    for a,b in E:p[find(a)]=find(b)
    blocks={}
    for v in range(1,n+1):blocks.setdefault(find(v),[]).append(v)
    return sorted(blocks.values())

def initial_closure(n,E):
    state={e:-1 for e in E}
    while True:
        changes={}
        for v in range(1,n+1):
            incident=[e for e in E if v in e]
            pending=[e for e in incident if state[e]<0]
            deficit=2-sum(state[e]==1 for e in incident)
            need(0<=deficit<=len(pending),'root closure feasible')
            if deficit in (0,len(pending)):
                value=int(deficit>0)
                for e in pending:
                    need(e not in changes or changes[e]==value,'root closure consistency')
                    changes[e]=value
        if not changes:break
        state.update(changes)
    return {'forced':[e for e in E if state[e]==1],
            'optional':[e for e in E if state[e]<0],
            'excluded':[e for e in E if state[e]==0]}

def constrained_sets(n,E,targets):
    incident={v:[i for i,e in enumerate(E) if v in e] for v in range(1,n+1)}
    found=[];nodes=0
    def search(state):
        nonlocal nodes
        nodes+=1
        while True:
            changes={}
            for v in range(1,n+1):
                pending=[i for i in incident[v] if state[i]<0]
                deficit=targets[v]-sum(state[i]==1 for i in incident[v])
                if deficit<0 or deficit>len(pending):return
                if deficit in (0,len(pending)):
                    value=int(deficit>0)
                    for i in pending:
                        if i in changes and changes[i]!=value:return
                        changes[i]=value
            if not changes:break
            state=list(state)
            for i,value in changes.items():state[i]=value
        if all(x>=0 for x in state):
            found.append([E[i] for i,value in enumerate(state) if value==1]);return
        _,v=min((sum(state[i]<0 for i in incident[v]),v) for v in incident
                if any(state[i]<0 for i in incident[v]))
        i=next(i for i in incident[v] if state[i]<0)
        for value in (0,1):
            new=list(state);new[i]=value;search(new)
    search([-1]*len(E))
    return found,nodes

obstructions=[]
for n in range(1,31):
    E=edges(n);d=Counter(v for e in E for v in e)
    witness=min(v for v in range(1,n+1) if d[v]<2)
    obstructions.append({'n':n,'vertex':witness,'degree':d[witness]})

rows=[]
for n in (31,32):
    E=edges(n);covers,nodes=constrained_sets(n,E,{v:2 for v in range(1,n+1)})
    need(len({tuple(c) for c in covers})==len(covers),'distinct cover edge sets')
    records=[]
    for cover in covers:
        degree=Counter(v for e in cover for v in e)
        need(all(degree[v]==2 for v in range(1,n+1)),'cover degrees')
        parts=components(n,cover)
        records.append({'edges':cover,'components':parts,'cycle_sizes':sorted(map(len,parts))})
    records.sort(key=lambda x:(x['cycle_sizes'],x['edges']))
    closure=initial_closure(n,E)
    rows.append({'n':n,'graph_edges':len(E),'search_nodes':nodes,'covers':records,
                 'closure':closure,
                 'cover_count':len(records),'Hamiltonian_cycle_count':sum(len(x['components'])==1 for x in records)})
need([r['cycle_sizes'] for r in rows[0]['covers']]==[[11,20],[14,17]],'G31 cover types')
need([r['cycle_sizes'] for r in rows[1]['covers']]==[[3,6,23],[3,29],[5,27],[9,23],[32]],'G32 cover types')
need([r['Hamiltonian_cycle_count'] for r in rows]==[0,1],'first-cycle threshold')
need([(len(r['closure']['forced']),len(r['closure']['optional'])) for r in rows]==[(23,15),(23,19)],
     'propagation inventories')
# Certify the six-chain K4 proof directly against the independently built graph.
chains=[[1,3],[1,15],[1,24,12,4],[3,6,10,15],[3,22,14,11,5,4],[4,21,15]]
chain_edges=[tuple(sorted((a,b))) for p in chains for a,b in zip(p,p[1:])]
need(set(chain_edges)==set(rows[0]['closure']['optional']) and len(chain_edges)==15,'K4 optional chains')
expanded=[]
for u in (0,1):
    values=[0,1,1,u,1-u,u]
    selected=list(rows[0]['closure']['forced'])
    for chain,value in zip(chains,values):
        for i,(a,b) in enumerate(zip(chain,chain[1:])):
            if (value if i%2==0 else 1-value):selected.append(tuple(sorted((a,b))))
    need(all(Counter(v for e in selected for v in e)[v]==2 for v in range(1,32)),
         'Boolean K4 expansion degrees')
    expanded.append(tuple(sorted(selected)))
need(set(expanded)=={tuple(c['edges']) for c in rows[0]['covers']},'two-case proof exhausts covers')
# Removing the degree-two arrival32 gives an endpoint-conditioned G31 path.
paths,nodes=constrained_sets(31,edges(31),{v:1 if v in (4,17) else 2 for v in range(1,32)})
connected=[E for E in paths if len(components(31,E))==1]
need(len(connected)==1,'unique endpoint-conditioned G31 path')
cycle=next(c['edges'] for c in rows[1]['covers'] if len(c['components'])==1)
need(set(cycle)-{(4,32),(17,32)}==set(connected[0]),'cycle/path arrival bijection')

root=Path(__file__).resolve().parent.parent
sources=['cycle-threshold/README.md','cycle-threshold/check_cycle_threshold.py',
         'cycle-threshold/certificate.json','cycle-threshold/certificate.optimized.json']
source_pins={p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in sources}
producer=json.loads((root/'cycle-threshold/certificate.json').read_text())
need(producer['source_sha256']==source_pins['cycle-threshold/check_cycle_threshold.py'],
     'producer checker hash')
def path_edges(path,close=False):
    pairs=list(zip(path,path[1:]))
    if close:pairs.append((path[-1],path[0]))
    return {tuple(sorted(e)) for e in pairs}
for row in rows:
    n=row['n'];pr=producer['g'+str(n)];E=set(edges(n))
    need(E=={tuple(e) for e in pr['edges']},'complete producer graph')
    need(set(row['closure']['forced'])=={tuple(e) for e in pr['forced_edges']},'producer forced edges')
    need(set(row['closure']['optional'])|set(row['closure']['forced'])==
         {tuple(e) for e in pr['allowed_edges']},'producer allowed edges')
    actual={tuple(c['edges']) for c in row['covers']}
    listed=set()
    for cover in pr['two_factors']:
        combined=set().union(*(path_edges(c,True) for c in cover['cycles']))
        need(len(combined)==n and combined<=E,'displayed producer cycle cover')
        need(sorted(len(c) for c in cover['cycles'])==cover['component_lengths'],'producer component lengths')
        need(set(v for c in cover['cycles'] for v in c)==set(range(1,n+1)),'producer cover spans')
        listed.add(tuple(sorted(combined)))
    need(actual==listed,'complete producer cover sets')
    # Replay each claimed inference using only its currently allowed/forced sets.
    allowed=set(E);forced=set()
    for item in pr['forcing_trace']:
        v=item['vertex'];inc={e for e in allowed if v in e};add={tuple(e) for e in item['edges']}
        if item['action']=='force':
            need(len(inc)==2 and add==inc-forced,'valid forcing trace')
            forced|=add
        else:
            need(item['action']=='delete' and len({e for e in forced if v in e})==2,
                 'valid exclusion premise')
            need(add==inc-forced,'valid exclusion trace')
            allowed-=add
    need(forced==set(row['closure']['forced']) and allowed==E-set(row['closure']['excluded']),
         'forcing trace final closure')
need(path_edges(producer['ear_bijection']['g31_prescribed_endpoint_path'])==set(connected[0]),
     'producer endpoint-conditioned path')
need(path_edges(producer['g32']['hamiltonian_cycle'],True)==set(cycle),'producer unique cycle')
cover1120=next(c['edges'] for c in rows[0]['covers'] if c['cycle_sizes']==[11,20])
removed=set(cover1120)-set(cycle);added=set(cycle)-set(cover1120)
need(removed=={(1,24),(4,12),(8,17)} and added=={(1,8),(12,24),(4,32),(17,32)},
     'exact local surgery')
need(removed=={tuple(e) for e in producer['surgery_from_11_20_cover']['delete']} and
     added=={tuple(e) for e in producer['surgery_from_11_20_cover']['add']},'producer surgery')
chains=producer['g31']['forced_chains']
need(set().union(*(path_edges(c) for c in chains))==set(rows[0]['closure']['forced']),
     'forced-chain expansion')
owner={v:i for i,c in enumerate(chains) for v in c}
bare={tuple(sorted((owner[a],owner[b]))) for a,b in rows[0]['closure']['optional']}
reps=producer['g31']['bare_kernel_hostile_representatives'];q=[owner[v] for v in reps]
need(len(q)==len(set(q))==8 and path_edges(q,True)<=bare,'spurious bare quotient cycle')
# The specific two wrong port allocations in that quotient cycle.
need(owner[5]==owner[22] and owner[11]==owner[24],'hostile partner ports')
need({(4,5),(5,11),(11,14)}<=set(rows[0]['closure']['optional']),'hostile duplicate port edges')
out={'status':'PASS independent G31/G32 degree-constraint audit','explicit_checks':checks,
     'scope':'Complete spanning 2-factor edge sets for n31 and n32; all smaller orders excluded by degree<2 witnesses',
     'small_order_degree_obstructions':obstructions,'rows':rows,
     'G31_prescribed_4_17_path':{'edge_set':connected[0],'degree_feasible_count':len(paths),'search_nodes':nodes},
     'source_sha256':source_pins}
print(json.dumps(out,indent=2))
