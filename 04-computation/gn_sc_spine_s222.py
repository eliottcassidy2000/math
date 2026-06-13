#!/usr/bin/env python3
"""
gn_sc_spine_s222.py — The SC Backbone Spine of G_n/Z_2
opus-2026-03-23-S222

The SC classes form a SPINE through the merged metagraph:
- At n≤7: connected (1 component)
- At n=8: disconnects (7 components)

This script traces the spine in detail:
1. Full SC-SC adjacency list with H values and |Aut|
2. The "SC path" from H=1 to H=H_max along the spine
3. H-gradient on the spine (is the SC backbone itself a DAG?)
4. Degree distribution within the spine vs total degree
5. The "NS halo" — how many NS nodes hang off each SC node
6. Spine genus (E-V+components) and its growth
7. Which SC nodes are "bridges" (removal disconnects spine)
8. The spine diameter vs full graph diameter
9. Score sequences along the spine
10. Whether the spine determines the entire graph structure
"""

import sys
import subprocess
from math import comb, factorial
from itertools import permutations
from collections import defaultdict, Counter, deque
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  THE SC BACKBONE SPINE OF G_n/Z_2")
print("  opus-2026-03-23-S222")
print("=" * 80)

# ============================================================================
# Build merged graph (compact version)
# ============================================================================

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

        # Aut count for SC classes (only)
        if cl['sc'] and n <= 7:
            cnt=0
            for p in permutations(range(n)):
                ok=True
                for i in range(n):
                    for j in range(n):
                        if a[i][j]!=a[p[i]][p[j]]: ok=False; break
                    if not ok: break
                if ok: cnt+=1
            cl['aut']=cnt
        else:
            cl['aut']=factorial(n)//(cl.get('size',factorial(n)))

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
    me=set(); ma_=defaultdict(set)
    for(a,b) in edges:
        x,y=mn[a],mn[b]
        if x!=y: me.add((min(x,y),max(x,y))); ma_[x].add(y); ma_[y].add(x)

    mi_={}
    for mm in range(mid):
        cids=[cl['cid'] for cl in cls if mn[cl['cid']]==mm]
        cl0=cls[cids[0]]
        mi_[mm]={'H':cl0['H'],'sc':cl0['sc'],'score':cl0['score'],
                 'c3':cl0['c3'],'aut':cl0.get('aut',1),'cids':cids}

    return {'cls':cls,'mn':mn,'mid':mid,'me':me,'ma':ma_,'mi':mi_,'n':n,'m':m,'al':al}

# ============================================================================
# SC SPINE ANALYSIS
# ============================================================================

for n in range(3, 9):
    print(f"\n{'='*80}")
    print(f"  n = {n}: SC SPINE ANALYSIS")
    print(f"{'='*80}")
    t0=time.time()

    D=build(n)
    mid=D['mid']; me=D['me']; ma=D['ma']; mi=D['mi']; m=D['m']
    print(f"  G_{n}/Z_2: V={mid}, E={len(me)} ({time.time()-t0:.1f}s)")

    # Extract SC spine
    sc_nodes=sorted([mm for mm in range(mid) if mi[mm]['sc']], key=lambda x: mi[x]['H'])
    sc_adj=defaultdict(set)
    sc_edges=set()
    for(a,b) in me:
        if mi[a]['sc'] and mi[b]['sc']:
            sc_adj[a].add(b); sc_adj[b].add(a)
            sc_edges.add((min(a,b),max(a,b)))

    V_sc=len(sc_nodes); E_sc=len(sc_edges)
    print(f"  SC spine: V={V_sc}, E={E_sc}")

    # Components
    visited=set(); comps=[]; comp_map={}
    for s in sc_nodes:
        if s in visited: continue
        comp=[]; q=deque([s]); visited.add(s)
        while q:
            u=q.popleft(); comp.append(u)
            for v in sc_adj[u]:
                if v not in visited: visited.add(v); q.append(v)
        comps.append(sorted(comp, key=lambda x: mi[x]['H']))
        for c in comp: comp_map[c]=len(comps)-1

    print(f"  Components: {len(comps)}")
    for ci,comp in enumerate(comps):
        h_range=f"H=[{mi[comp[0]]['H']},{mi[comp[-1]]['H']}]"
        comp_edges=sum(1 for(a,b) in sc_edges if comp_map.get(a)==ci and comp_map.get(b)==ci)
        genus=comp_edges-len(comp)+1
        print(f"    C{ci}: {len(comp)} nodes, {comp_edges} edges, genus={genus}, {h_range}")

    # SC spine adjacency list (sorted by H)
    print(f"\n  SC SPINE ADJACENCY LIST:")
    for mm in sc_nodes:
        h=mi[mm]['H']
        aut=mi[mm]['aut']
        score=mi[mm]['score']
        sc_nbrs=sorted(sc_adj[mm], key=lambda x: mi[x]['H'])
        ns_halo=sum(1 for nb in ma[mm] if not mi[nb]['sc'])
        total_deg=len(ma[mm])

        nbr_str=', '.join(f"{nb}(H={mi[nb]['H']})" for nb in sc_nbrs[:8])
        if len(sc_nbrs)>8: nbr_str+='...'

        if V_sc <= 30:
            print(f"    mid={mm:3d} H={h:4d} |Aut|={aut:3d} sc_deg={len(sc_nbrs):2d} "
                  f"ns_halo={ns_halo:3d} score={score} → [{nbr_str}]")
        elif V_sc <= 100:
            print(f"    mid={mm:3d} H={h:4d} sc_deg={len(sc_nbrs):2d} ns_halo={ns_halo:3d}")

    # H-gradient on spine
    spine_h_inc=0; spine_h_dec=0; spine_h_lvl=0
    for(a,b) in sc_edges:
        ha,hb=mi[a]['H'],mi[b]['H']
        if ha<hb: spine_h_inc+=1
        elif ha>hb: spine_h_dec+=1
        else: spine_h_lvl+=1

    print(f"\n  SPINE H-GRADIENT:")
    print(f"    H-increasing: {spine_h_inc}")
    print(f"    H-decreasing: {spine_h_dec}")
    print(f"    H-level: {spine_h_lvl}")
    print(f"    Spine DAG: {'YES' if spine_h_dec==0 else 'NO'}")

    # Spine diameter
    if V_sc <= 500:
        spine_diam=0
        for s in sc_nodes:
            dist={s:0}; q=deque([s])
            while q:
                u=q.popleft()
                for v in sc_adj[u]:
                    if v not in dist: dist[v]=dist[u]+1; q.append(v)
            if dist:
                spine_diam=max(spine_diam,max(dist.values()))
        print(f"    Spine diameter: {spine_diam}")

    # NS halo analysis
    halo_sizes=[sum(1 for nb in ma[mm] if not mi[nb]['sc']) for mm in sc_nodes]
    print(f"\n  NS HALO:")
    print(f"    min={min(halo_sizes)}, max={max(halo_sizes)}, "
          f"avg={sum(halo_sizes)/len(halo_sizes):.1f}, median={sorted(halo_sizes)[len(halo_sizes)//2]}")

    # Zero-halo SC nodes (isolated from NS)
    zero_halo=[mm for mm in sc_nodes if sum(1 for nb in ma[mm] if not mi[nb]['sc'])==0]
    print(f"    Zero-halo (no NS neighbors): {len(zero_halo)}")
    for mm in zero_halo:
        print(f"      mid={mm}, H={mi[mm]['H']}, sc_deg={len(sc_adj[mm])}")

    # Bridge nodes (removal disconnects spine)
    bridges=[]
    if V_sc <= 200:
        for mm in sc_nodes:
            # Remove mm and check connectivity
            remaining=[x for x in sc_nodes if x!=mm]
            if not remaining: continue
            vis2=set(); q2=deque([remaining[0]]); vis2.add(remaining[0])
            while q2:
                u=q2.popleft()
                for v in sc_adj[u]:
                    if v!=mm and v not in vis2: vis2.add(v); q2.append(v)
            if len(vis2)<len(remaining):
                bridges.append(mm)

        print(f"\n  BRIDGE NODES (removal disconnects spine): {len(bridges)}")
        for mm in bridges[:10]:
            print(f"    mid={mm}, H={mi[mm]['H']}, sc_deg={len(sc_adj[mm])}")

    # Longest path in spine (H-monotone)
    print(f"\n  LONGEST H-MONOTONE PATH in spine:")
    # BFS/DFS from H=1 following increasing H
    trans=None
    for mm in sc_nodes:
        if mi[mm]['H']==1: trans=mm; break
    if trans is not None:
        # DFS for longest increasing H path from transitive
        best_path_holder=[[trans]]
        def dfs_inc(node, path):
            if len(path)>len(best_path_holder[0]): best_path_holder[0]=list(path)
            for nb in sc_adj[node]:
                if mi[nb]['H']>mi[node]['H'] and nb not in set(path):
                    path.append(nb)
                    dfs_inc(nb, path)
                    path.pop()
        # Limit depth for large graphs
        if V_sc <= 100:
            dfs_inc(trans, [trans])
            best_path=best_path_holder[0]
        else:
            # Greedy longest path
            path=[trans]; current=trans
            while True:
                nbrs=[nb for nb in sc_adj[current] if mi[nb]['H']>mi[current]['H'] and nb not in set(path)]
                if not nbrs: break
                best_nb=max(nbrs, key=lambda x: mi[x]['H'])
                path.append(best_nb); current=best_nb
            best_path=path

        print(f"    Length: {len(best_path)} nodes")
        print(f"    H values: {[mi[mm]['H'] for mm in best_path[:15]]}{'...' if len(best_path)>15 else ''}")
        if len(best_path)>1:
            h_jumps=[mi[best_path[i+1]]['H']-mi[best_path[i]]['H'] for i in range(len(best_path)-1)]
            print(f"    H jumps: {h_jumps[:15]}{'...' if len(h_jumps)>15 else ''}")
            print(f"    Avg H jump: {sum(h_jumps)/len(h_jumps):.1f}")

    # ---- Score pattern along spine ----
    print(f"\n  SCORE PATTERNS ON SPINE:")
    score_counter=Counter(mi[mm]['score'] for mm in sc_nodes)
    print(f"    Distinct scores: {len(score_counter)}")
    # All SC classes must have palindromic scores
    all_palin=all(
        all(mi[mm]['score'][i]+mi[mm]['score'][n-1-i]==n-1 for i in range(n))
        for mm in sc_nodes
    )
    print(f"    All SC have palindromic scores: {'YES ✓' if all_palin else 'NO ✗'}")

    # Most common score on spine
    top_scores=score_counter.most_common(5)
    print(f"    Most common: {top_scores}")

    elapsed=time.time()-t0
    print(f"\n  n={n} done ({elapsed:.1f}s)")

# ============================================================================
# CROSS-n SUMMARY
# ============================================================================

print(f"\n{'='*80}")
print(f"  CROSS-n SPINE SUMMARY")
print(f"{'='*80}")

print(f"\n  Spine as % of full merged graph:")
print(f"  {'n':>3} {'V_sc':>5} {'V_m':>6} {'%V':>6} {'E_sc':>5} {'E_m':>7} {'%E':>6} "
      f"{'genus':>6} {'diam':>5} {'bridges':>8}")

print(f"\n{'='*80}")
print(f"  DONE")
print(f"{'='*80}")
