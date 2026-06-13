#!/usr/bin/env python3
"""
gn_triple_view_s224.py — MANDATORY triple analysis: Blue, Black, Combined
opus-2026-03-23-S224

Every metric computed THREE ways:
  BLUE subgraph  = SC-preserving edges (SC↔SC and NS↔NS)
  BLACK subgraph = SC-changing edges (SC↔NS only)
  COMBINED       = all edges

This is the correct way to analyze G_n/Z_2: never look at just one view.
"""

import sys
import subprocess
from math import comb, factorial
from itertools import permutations
from collections import defaultdict, Counter, deque
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  TRIPLE-VIEW ANALYSIS: BLUE / BLACK / COMBINED")
print("  opus-2026-03-23-S224")
print("=" * 80)

# Build function (same as S222/S223)
def build(n):
    m=comb(n,2); P=[(i,j) for i in range(n) for j in range(i+1,n)]
    if n<=6:
        AP=list(permutations(range(n)))
        def b2a(b):
            a=[[0]*n for _ in range(n)]
            for k,(i,j) in enumerate(P):
                if b&(1<<k): a[i][j]=1
                else: a[j][i]=1
            return a
        def cn(a):
            best=None
            for p in AP:
                f=tuple(a[p[i]][p[j]] for i in range(n) for j in range(n))
                if best is None or f<best: best=f
            return best
        cm={}
        for b in range(1<<m):
            a=b2a(b); c=cn(a)
            if c not in cm: cm[c]={'bits':b,'adj':a,'size':0}
            cm[c]['size']+=1
        cls=list(cm.values())
        for i,c in enumerate(cls): c['cid']=i
        ny=False
    else:
        r=subprocess.run(['gentourng',str(n)],capture_output=True,text=True)
        ls=[l.strip() for l in (r.stdout or '').split('\n') if len(l.strip())==m and all(c in '01' for c in l.strip())]
        def b2a_n(b):
            a=[[0]*n for _ in range(n)]
            for k,(i,j) in enumerate(P):
                if b&(1<<(m-1-k)): a[i][j]=1
                else: a[j][i]=1
            return a
        cls=[{'bits':int(l,2),'adj':b2a_n(int(l,2)),'cid':i} for i,l in enumerate(ls)]
        ny=True

    def ss(a): return tuple(sorted(sum(a[i][j] for j in range(n)) for i in range(n)))
    def c3(a):
        c=0
        for i in range(n):
            for j in range(i+1,n):
                for k in range(j+1,n):
                    if a[i][j] and a[j][k] and a[k][i]: c+=1
                    if a[i][k] and a[k][j] and a[j][i]: c+=1
        return c
    def Hf(a):
        dp={}
        for v in range(n): dp[(1<<v,v)]=1
        for S in range(1,1<<n):
            for v in range(n):
                if not(S&(1<<v)): continue
                val=dp.get((S,v),0)
                if val==0: continue
                for u in range(n):
                    if S&(1<<u): continue
                    if a[v][u]: key=(S|(1<<u),u); dp[key]=dp.get(key,0)+val
        return sum(dp.get(((1<<n)-1,v),0) for v in range(n))
    def cp(a):
        scores=[sum(a[i][j] for j in range(n)) for i in range(n)]
        sg=defaultdict(list)
        for v in range(n): sg[scores[v]].append(v)
        gs=[sg[s] for s in sorted(set(scores))]
        if all(len(g)==1 for g in gs):
            p=[g[0] for g in gs]
            return tuple(a[p[i]][p[j]] for i in range(n) for j in range(n))
        best=None
        def gen(gs2):
            if not gs2: yield []; return
            for p in permutations(gs2[0]):
                for rest in gen(gs2[1:]): yield list(p)+rest
        for p in gen(gs):
            f=tuple(a[p[i]][p[j]] for i in range(n) for j in range(n))
            if best is None or f<best: best=f
        return best

    hm=defaultdict(list); cm2={}
    for cl in cls:
        a=cl['adj']; cl['score']=ss(a); cl['c3']=c3(a); cl['H']=Hf(a)
        cl['own_canon']=cp(a)
        comp=[[1-a[i][j] if i!=j else 0 for j in range(n)] for i in range(n)]
        cl['comp_canon']=cp(comp); cl['sc']=(cl['own_canon']==cl['comp_canon'])
        hm[(cl['score'],cl['c3'],cl['H'])].append(cl['cid'])
        cm2[cl['own_canon']]=cl['cid']
    for cl in cls: cl['comp_cid']=cm2.get(cl['comp_canon'],-1)

    def lk(a):
        s=ss(a);c=c3(a);h=Hf(a)
        cids=hm.get((s,c,h))
        if not cids: return -1
        if len(cids)==1: return cids[0]
        return cm2.get(cp(a),-1)

    edges=set(); al=defaultdict(set)
    for cl in cls:
        b=cl['bits']
        for arc in range(m):
            if ny: fb=b^(1<<(m-1-arc))
            else: fb=b^(1<<arc)
            if ny:
                fa=[[0]*n for _ in range(n)]
                for k,(i,j) in enumerate(P):
                    if fb&(1<<(m-1-k)): fa[i][j]=1
                    else: fa[j][i]=1
            else:
                fa=[[0]*n for _ in range(n)]
                for k,(i,j) in enumerate(P):
                    if fb&(1<<k): fa[i][j]=1
                    else: fa[j][i]=1
            nb=lk(fa)
            if nb>=0 and nb!=cl['cid']:
                edges.add((min(cl['cid'],nb),max(cl['cid'],nb)))
                al[cl['cid']].add(nb); al[nb].add(cl['cid'])

    # Merged
    mn={}; mid=0
    for cl in cls:
        ci=cl['cid']
        if ci in mn: continue
        if cl['sc']: mn[ci]=mid; mid+=1
        else:
            mn[ci]=mid
            comp=cl['comp_cid']
            if comp>=0 and comp!=ci: mn[comp]=mid
            mid+=1
    me=set(); ma=defaultdict(set)
    for(a,b) in edges:
        x,y=mn[a],mn[b]
        if x!=y: me.add((min(x,y),max(x,y))); ma[x].add(y); ma[y].add(x)

    mi={}
    for mm in range(mid):
        cids=[cl['cid'] for cl in cls if mn[cl['cid']]==mm]
        mi[mm]={'H':cls[cids[0]]['H'],'sc':cls[cids[0]]['sc'],
                'score':cls[cids[0]]['score'],'c3':cls[cids[0]]['c3']}

    # Classify edges as blue or black
    blue=set(); black=set()
    blue_adj=defaultdict(set); black_adj=defaultdict(set)
    for(a,b) in me:
        if mi[a]['sc']==mi[b]['sc']:
            blue.add((a,b)); blue_adj[a].add(b); blue_adj[b].add(a)
        else:
            black.add((a,b)); black_adj[a].add(b); black_adj[b].add(a)

    return {'mid':mid,'me':me,'ma':ma,'mi':mi,'blue':blue,'black':black,
            'blue_adj':blue_adj,'black_adj':black_adj,'n':n,'m':m}

# ============================================================================
# TRIPLE-VIEW METRICS
# ============================================================================

def graph_metrics(V, edges, adj, mi, label):
    """Compute standard metrics for a subgraph."""
    E=len(edges)
    if E==0:
        return {'V':V,'E':0,'label':label,'connected':False}

    degs=[len(adj.get(mm,set())) for mm in range(V)]
    active=[mm for mm in range(V) if degs[mm]>0]
    V_active=len(active)

    # Components
    visited=set(); comps=0
    for s in range(V):
        if s in visited or not adj.get(s): continue
        comps+=1; q=deque([s]); visited.add(s)
        while q:
            u=q.popleft()
            for v in adj[u]:
                if v not in visited: visited.add(v); q.append(v)

    # Diameter (sample if large)
    diam=0
    if V_active<=500:
        for s in active[:min(50,len(active))]:
            dist={s:0}; q=deque([s])
            while q:
                u=q.popleft()
                for v in adj.get(u,set()):
                    if v not in dist: dist[v]=dist[u]+1; q.append(v)
            if dist: diam=max(diam,max(dist.values()))

    # Triangles
    tri=0
    for(a,b) in edges:
        tri+=len(adj.get(a,set())&adj.get(b,set()))
    tri//=3

    # Level edges
    lvl=sum(1 for(a,b) in edges if mi[a]['H']==mi[b]['H'])

    # H-increasing vs decreasing
    h_inc=sum(1 for(a,b) in edges if mi[a]['H']!=mi[b]['H'])

    return {
        'V':V, 'V_active':V_active, 'E':E, 'label':label,
        'components':comps,
        'deg_min':min(degs[mm] for mm in active) if active else 0,
        'deg_max':max(degs[mm] for mm in active) if active else 0,
        'deg_avg':sum(degs[mm] for mm in active)/V_active if V_active else 0,
        'diameter':diam,
        'triangles':tri,
        'level_edges':lvl,
        'h_changing':h_inc,
        'density':2*E/(V_active*(V_active-1)) if V_active>1 else 0,
    }

ALL={}

for n in range(3, 9):
    print(f"\n{'='*80}")
    print(f"  n = {n}")
    print(f"{'='*80}")
    t0=time.time()

    D=build(n)
    mid=D['mid']; mi=D['mi']

    # Three views
    combined=graph_metrics(mid, D['me'], D['ma'], mi, 'COMBINED')
    blue_m=graph_metrics(mid, D['blue'], D['blue_adj'], mi, 'BLUE')
    black_m=graph_metrics(mid, D['black'], D['black_adj'], mi, 'BLACK')

    # SC/NS node counts
    sc_count=sum(1 for mm in range(mid) if mi[mm]['sc'])
    ns_count=mid-sc_count

    # Blue subgraph components: should split into SC and NS
    # Count SC-only and NS-only components in blue
    blue_sc_comp=0; blue_ns_comp=0; blue_mixed=0
    visited=set()
    for s in range(mid):
        if s in visited or not D['blue_adj'].get(s): continue
        comp=set(); q=deque([s]); visited.add(s)
        while q:
            u=q.popleft(); comp.add(u)
            for v in D['blue_adj'][u]:
                if v not in visited: visited.add(v); q.append(v)
        has_sc=any(mi[mm]['sc'] for mm in comp)
        has_ns=any(not mi[mm]['sc'] for mm in comp)
        if has_sc and has_ns: blue_mixed+=1
        elif has_sc: blue_sc_comp+=1
        else: blue_ns_comp+=1
    # Isolated nodes (no blue edges)
    blue_isolated=sum(1 for mm in range(mid) if not D['blue_adj'].get(mm))

    print(f"  V={mid} ({sc_count} SC + {ns_count} NS), built in {time.time()-t0:.1f}s")
    print(f"\n  {'METRIC':>20s} {'BLUE':>10s} {'BLACK':>10s} {'COMBINED':>10s}")
    print(f"  {'─'*20} {'─'*10} {'─'*10} {'─'*10}")

    for key in ['E','V_active','components','deg_min','deg_max','density','diameter','triangles','level_edges','h_changing']:
        bv=blue_m.get(key,'?')
        kv=black_m.get(key,'?')
        cv=combined.get(key,'?')
        if isinstance(bv,float): bv=f"{bv:.4f}"
        if isinstance(kv,float): kv=f"{kv:.4f}"
        if isinstance(cv,float): cv=f"{cv:.4f}"
        print(f"  {key:>20s} {str(bv):>10s} {str(kv):>10s} {str(cv):>10s}")

    print(f"\n  Blue components: {blue_sc_comp} SC-only + {blue_ns_comp} NS-only + {blue_mixed} mixed + {blue_isolated} isolated")

    # Degree comparison
    print(f"\n  DEGREE COMPARISON per node type:")
    sc_nodes=[mm for mm in range(mid) if mi[mm]['sc']]
    ns_nodes=[mm for mm in range(mid) if not mi[mm]['sc']]

    if sc_nodes:
        sc_blue_d=[len(D['blue_adj'].get(mm,set())) for mm in sc_nodes]
        sc_black_d=[len(D['black_adj'].get(mm,set())) for mm in sc_nodes]
        sc_total_d=[len(D['ma'].get(mm,set())) for mm in sc_nodes]
        print(f"    SC nodes ({len(sc_nodes)}):")
        print(f"      blue_deg: avg={sum(sc_blue_d)/len(sc_blue_d):.2f}, range=[{min(sc_blue_d)},{max(sc_blue_d)}]")
        print(f"      black_deg: avg={sum(sc_black_d)/len(sc_black_d):.2f}, range=[{min(sc_black_d)},{max(sc_black_d)}]")
        print(f"      total_deg: avg={sum(sc_total_d)/len(sc_total_d):.2f}, range=[{min(sc_total_d)},{max(sc_total_d)}]")
        print(f"      blue/total: {sum(sc_blue_d)/sum(sc_total_d):.4f}" if sum(sc_total_d)>0 else "")

    if ns_nodes:
        ns_blue_d=[len(D['blue_adj'].get(mm,set())) for mm in ns_nodes]
        ns_black_d=[len(D['black_adj'].get(mm,set())) for mm in ns_nodes]
        ns_total_d=[len(D['ma'].get(mm,set())) for mm in ns_nodes]
        print(f"    NS nodes ({len(ns_nodes)}):")
        print(f"      blue_deg: avg={sum(ns_blue_d)/len(ns_blue_d):.2f}, range=[{min(ns_blue_d)},{max(ns_blue_d)}]")
        print(f"      black_deg: avg={sum(ns_black_d)/len(ns_black_d):.2f}, range=[{min(ns_black_d)},{max(ns_black_d)}]")
        print(f"      total_deg: avg={sum(ns_total_d)/len(ns_total_d):.2f}, range=[{min(ns_total_d)},{max(ns_total_d)}]")
        print(f"      blue/total: {sum(ns_blue_d)/sum(ns_total_d):.4f}" if sum(ns_total_d)>0 else "")

    ALL[n]={'combined':combined,'blue':blue_m,'black':black_m,
            'sc':sc_count,'ns':ns_count,
            'blue_sc_comp':blue_sc_comp,'blue_ns_comp':blue_ns_comp,
            'blue_mixed':blue_mixed,'blue_isolated':blue_isolated}

# ============================================================================
# CROSS-n COMPARISON
# ============================================================================

print(f"\n{'='*80}")
print(f"  CROSS-n TRIPLE-VIEW SEQUENCES")
print(f"{'='*80}")

for view in ['BLUE','BLACK','COMBINED']:
    key_map={'BLUE':'blue','BLACK':'black','COMBINED':'combined'}
    vk=key_map[view]
    print(f"\n  {view}:")
    for metric in ['E','components','deg_max','density','diameter','triangles','level_edges']:
        vals=[ALL[n][vk].get(metric,'?') for n in range(3,9)]
        fvals=[]
        for v in vals:
            if isinstance(v,float): fvals.append(f"{v:.4f}")
            else: fvals.append(str(v))
        print(f"    {metric:>15s}: {', '.join(fvals)}")

print(f"\n  BLUE COMPONENT STRUCTURE:")
print(f"    {'n':>3} {'SC_comp':>7} {'NS_comp':>7} {'mixed':>6} {'isolated':>8}")
for n in range(3,9):
    A=ALL[n]
    print(f"    {n:3d} {A['blue_sc_comp']:7d} {A['blue_ns_comp']:7d} "
          f"{A['blue_mixed']:6d} {A['blue_isolated']:8d}")

print(f"\n{'='*80}")
print(f"  DONE")
print(f"{'='*80}")
