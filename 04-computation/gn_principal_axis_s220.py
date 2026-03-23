#!/usr/bin/env python3
"""
gn_principal_axis_s220.py — Principal Axis Analysis of G_n/Z_2
opus-2026-03-23-S220

The transitive tournament class (H=1, always SC) connects to exactly ONE
other SC class via a blue edge. This "principal blue line" is the axis
of symmetry for the merged metagraph.

We line up all classes relative to this axis:
- The transitive class is the "north pole"
- Its unique SC blue neighbor is the "first waypoint"
- Classes on one side vs the other of this line form the bilateral structure

For each n=3..8, compute:
1. The transitive class and its blue SC neighbor
2. BFS from the transitive class along the SC backbone
3. The "bilateral partition" — classes on each side of the principal line
4. Symmetry properties of this partition
5. How the H-gradient flows relative to the axis
"""

import sys
import subprocess
from math import comb, factorial
from itertools import permutations
from collections import defaultdict, Counter, deque
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  PRINCIPAL AXIS ANALYSIS OF G_n/Z_2")
print("  opus-2026-03-23-S220")
print("=" * 80)

# ============================================================================
# Build merged graph (reusable from S218/S219)
# ============================================================================

def build_full_data(n):
    """Build G_n and G_n/Z_2 with full metadata."""
    m = comb(n, 2)
    PAIRS = [(i,j) for i in range(n) for j in range(i+1,n)]

    if n <= 6:
        ALL_PERMS = list(permutations(range(n)))
        def b2a(bits):
            adj=[[0]*n for _ in range(n)]
            for k,(i,j) in enumerate(PAIRS):
                if bits&(1<<k): adj[i][j]=1
                else: adj[j][i]=1
            return adj
        def canon(adj):
            best=None
            for perm in ALL_PERMS:
                form=tuple(adj[perm[i]][perm[j]] for i in range(n) for j in range(n))
                if best is None or form<best: best=form
            return best
        cmap={}
        for bits in range(1<<m):
            adj=b2a(bits); c=canon(adj)
            if c not in cmap: cmap[c]={'bits':bits,'adj':adj,'canon':c,'size':0}
            cmap[c]['size']+=1
        classes=list(cmap.values())
        for i,cl in enumerate(classes): cl['cid']=i
        is_nauty=False
    else:
        res=subprocess.run(['gentourng',str(n)],capture_output=True,text=True)
        lines=[l.strip() for l in (res.stdout or '').split('\n')
               if len(l.strip())==m and all(c in '01' for c in l.strip())]
        def b2a_n(bits):
            adj=[[0]*n for _ in range(n)]
            for k,(i,j) in enumerate(PAIRS):
                if bits&(1<<(m-1-k)): adj[i][j]=1
                else: adj[j][i]=1
            return adj
        classes=[{'bits':int(l,2),'adj':b2a_n(int(l,2)),'cid':i} for i,l in enumerate(lines)]
        is_nauty=True

    def ss(adj): return tuple(sorted(sum(adj[i][j] for j in range(n)) for i in range(n)))
    def c3(adj):
        c=0
        for i in range(n):
            for j in range(i+1,n):
                for k in range(j+1,n):
                    if adj[i][j] and adj[j][k] and adj[k][i]: c+=1
                    if adj[i][k] and adj[k][j] and adj[j][i]: c+=1
        return c
    def Hdp(adj):
        dp={}
        for v in range(n): dp[(1<<v,v)]=1
        for S in range(1,1<<n):
            for v in range(n):
                if not(S&(1<<v)): continue
                val=dp.get((S,v),0)
                if val==0: continue
                for u in range(n):
                    if S&(1<<u): continue
                    if adj[v][u]: key=(S|(1<<u),u); dp[key]=dp.get(key,0)+val
        return sum(dp.get(((1<<n)-1,v),0) for v in range(n))
    def cp(adj):
        scores=[sum(adj[i][j] for j in range(n)) for i in range(n)]
        sg=defaultdict(list)
        for v in range(n): sg[scores[v]].append(v)
        groups=[sg[s] for s in sorted(set(scores))]
        if all(len(g)==1 for g in groups):
            perm=[g[0] for g in groups]
            return tuple(adj[perm[i]][perm[j]] for i in range(n) for j in range(n))
        best=None
        def gen(gs):
            if not gs: yield []; return
            for p in permutations(gs[0]):
                for rest in gen(gs[1:]): yield list(p)+rest
        for perm in gen(groups):
            form=tuple(adj[perm[i]][perm[j]] for i in range(n) for j in range(n))
            if best is None or form<best: best=form
        return best

    hm=defaultdict(list); cm={}
    for cl in classes:
        adj=cl['adj']
        cl['score']=ss(adj); cl['c3']=c3(adj); cl['H']=Hdp(adj)
        cl['own_canon']=cp(adj)
        comp=[[1-adj[i][j] if i!=j else 0 for j in range(n)] for i in range(n)]
        cl['comp_canon']=cp(comp)
        cl['sc']=(cl['own_canon']==cl['comp_canon'])
        hm[(cl['score'],cl['c3'],cl['H'])].append(cl['cid'])
        cm[cl['own_canon']]=cl['cid']
    for cl in classes: cl['comp_cid']=cm.get(cl['comp_canon'],-1)

    def lk(adj):
        s=ss(adj);c=c3(adj);h=Hdp(adj)
        cids=hm.get((s,c,h))
        if not cids: return -1
        if len(cids)==1: return cids[0]
        return cm.get(cp(adj),-1)

    # Build original edges with blue/black classification
    edges=set(); adj_list=defaultdict(set)
    for cl in classes:
        bits=cl['bits']
        for arc in range(m):
            if is_nauty: fb=bits^(1<<(m-1-arc))
            else: fb=bits^(1<<arc)
            if is_nauty:
                fa=[[0]*n for _ in range(n)]
                for k,(i,j) in enumerate(PAIRS):
                    if fb&(1<<(m-1-k)): fa[i][j]=1
                    else: fa[j][i]=1
            else:
                fa=[[0]*n for _ in range(n)]
                for k,(i,j) in enumerate(PAIRS):
                    if fb&(1<<k): fa[i][j]=1
                    else: fa[j][i]=1
            nb=lk(fa)
            if nb>=0 and nb!=cl['cid']:
                edges.add((min(cl['cid'],nb),max(cl['cid'],nb)))
                adj_list[cl['cid']].add(nb)
                adj_list[nb].add(cl['cid'])

    # Classify edges
    edge_color={}
    for(a,b) in edges:
        if classes[a]['sc']==classes[b]['sc']:
            edge_color[(a,b)]='blue'
        else:
            edge_color[(a,b)]='black'

    # Build merged graph
    mn={}; mid=0
    for cl in classes:
        ci=cl['cid']
        if ci in mn: continue
        if cl['sc']: mn[ci]=mid; mid+=1
        else:
            mn[ci]=mid
            comp=cl['comp_cid']
            if comp>=0 and comp!=ci: mn[comp]=mid
            mid+=1
    me=set(); madj=defaultdict(set)
    collapsed=0
    for(a,b) in edges:
        ma,mb=mn[a],mn[b]
        if ma==mb: collapsed+=1
        else:
            me.add((min(ma,mb),max(ma,mb)))
            madj[ma].add(mb); madj[mb].add(ma)

    # Merged edge colors
    me_color={}
    for(a,b) in me:
        # Find original edge representatives
        # An edge in merged graph is blue if connecting same-type nodes
        a_sc=any(classes[ci]['sc'] for ci in range(len(classes)) if mn[ci]==a and classes[ci]['sc'])
        b_sc=any(classes[ci]['sc'] for ci in range(len(classes)) if mn[ci]==b and classes[ci]['sc'])
        # Actually, simpler: check if the merged nodes are SC or NS
        a_is_sc=False; b_is_sc=False
        for cl in classes:
            if mn[cl['cid']]==a and cl['sc']: a_is_sc=True
            if mn[cl['cid']]==b and cl['sc']: b_is_sc=True
        me_color[(min(a,b),max(a,b))]='blue' if a_is_sc==b_is_sc else 'black'

    # Map merged nodes to H values and types
    mn_info={}
    for mi in range(mid):
        cids=[cl['cid'] for cl in classes if mn[cl['cid']]==mi]
        h_vals=set(classes[ci]['H'] for ci in cids)
        is_sc=any(classes[ci]['sc'] for ci in cids)
        mn_info[mi]={
            'H':min(h_vals),'sc':is_sc,'cids':cids,
            'score':classes[cids[0]]['score'],'c3':classes[cids[0]]['c3']
        }

    return {
        'classes':classes,'edges':edges,'adj_list':adj_list,'edge_color':edge_color,
        'mn':mn,'mid':mid,'me':me,'madj':madj,'me_color':me_color,
        'mn_info':mn_info,'collapsed':collapsed,'n':n,'m':m
    }

# ============================================================================
# PRINCIPAL AXIS ANALYSIS
# ============================================================================

for n in range(3, 9):
    print(f"\n{'='*80}")
    print(f"  n = {n}: PRINCIPAL AXIS ANALYSIS")
    print(f"{'='*80}")
    t0=time.time()

    D=build_full_data(n)
    classes=D['classes']; mn=D['mn']; mid=D['mid']
    me=D['me']; madj=D['madj']; me_color=D['me_color']; mn_info=D['mn_info']
    print(f"  Built G_{n}/Z_2: V={mid}, E={len(me)} ({time.time()-t0:.1f}s)")

    # Find the transitive class (H=1, always SC)
    trans_cid=None
    for cl in classes:
        if cl['H']==1:
            trans_cid=cl['cid']; break
    trans_mid=mn[trans_cid]
    print(f"\n  Transitive class: original cid={trans_cid}, merged mid={trans_mid}")
    print(f"    H=1, SC=True, score={classes[trans_cid]['score']}")

    # Find its SC blue neighbors
    sc_blue_nbrs=[]
    black_nbrs=[]
    for nb in madj[trans_mid]:
        edge=(min(trans_mid,nb),max(trans_mid,nb))
        color=me_color.get(edge,'?')
        if color=='blue' and mn_info[nb]['sc']:
            sc_blue_nbrs.append(nb)
        elif color=='black':
            black_nbrs.append(nb)
        elif color=='blue' and not mn_info[nb]['sc']:
            sc_blue_nbrs.append(nb)  # NS-NS blue (trans is SC, so this is SC-SC only if nb is SC)

    # Actually re-classify: transitive is SC. Blue neighbors = other SC nodes.
    # Black neighbors = NS nodes.
    sc_nbrs=[nb for nb in madj[trans_mid] if mn_info[nb]['sc']]
    ns_nbrs=[nb for nb in madj[trans_mid] if not mn_info[nb]['sc']]

    print(f"    SC neighbors (blue edges): {len(sc_nbrs)}")
    for nb in sorted(sc_nbrs, key=lambda x: mn_info[x]['H']):
        print(f"      mid={nb}, H={mn_info[nb]['H']}, score={mn_info[nb]['score']}")
    print(f"    NS neighbors (black edges): {len(ns_nbrs)}")

    # The PRINCIPAL LINE: transitive → first SC blue neighbor
    if len(sc_nbrs)==1:
        axis_partner=sc_nbrs[0]
        print(f"\n  *** PRINCIPAL AXIS: {trans_mid} (H=1) ←→ {axis_partner} (H={mn_info[axis_partner]['H']}) ***")
        print(f"  This is the UNIQUE SC-SC blue edge from the transitive class.")
    elif len(sc_nbrs)==0:
        print(f"\n  *** NO SC BLUE NEIGHBOR — transitive is isolated in SC backbone! ***")
        axis_partner=None
    else:
        print(f"\n  *** MULTIPLE SC BLUE NEIGHBORS — principal axis is NOT unique ***")
        axis_partner=sc_nbrs[0]  # pick lowest H

    # BFS along SC backbone from transitive
    print(f"\n  SC BACKBONE (BFS from transitive):")
    sc_nodes=[mi for mi in range(mid) if mn_info[mi]['sc']]
    sc_adj=defaultdict(set)
    for(a,b) in me:
        if mn_info[a]['sc'] and mn_info[b]['sc']:
            sc_adj[a].add(b); sc_adj[b].add(a)

    visited=set(); queue=deque([trans_mid]); visited.add(trans_mid)
    sc_bfs_order=[]; sc_dist={}; sc_dist[trans_mid]=0
    while queue:
        u=queue.popleft()
        sc_bfs_order.append(u)
        for v in sorted(sc_adj[u], key=lambda x: mn_info[x]['H']):
            if v not in visited:
                visited.add(v); queue.append(v)
                sc_dist[v]=sc_dist[u]+1

    for mi in sc_bfs_order:
        d=sc_dist[mi]
        h=mn_info[mi]['H']
        sc_deg=len(sc_adj[mi])
        ns_deg=sum(1 for nb in madj[mi] if not mn_info[nb]['sc'])
        print(f"    depth={d}: mid={mi}, H={h}, SC_deg={sc_deg}, NS_deg={ns_deg}")

    sc_backbone_edges=sum(1 for(a,b) in me if mn_info[a]['sc'] and mn_info[b]['sc'])
    sc_backbone_components=0
    sc_visited=set()
    for mi in sc_nodes:
        if mi not in sc_visited:
            sc_backbone_components+=1
            q=deque([mi]); sc_visited.add(mi)
            while q:
                u=q.popleft()
                for v in sc_adj[u]:
                    if v not in sc_visited:
                        sc_visited.add(v); q.append(v)

    print(f"\n  SC backbone: {len(sc_nodes)} nodes, {sc_backbone_edges} edges, "
          f"{sc_backbone_components} components")

    # BILATERAL ANALYSIS: partition nodes by which side of the axis they're on
    if axis_partner is not None:
        # Remove the principal axis edge, compute components
        # Side A: connected to transitive (without the axis edge)
        # Side B: connected to axis_partner (without the axis edge)
        # Center: on the axis itself

        # Actually, the "bilateral" structure is:
        # For each non-axis node, compute BFS distance from transitive AND from axis_partner
        # Nodes closer to transitive = "north side"
        # Nodes closer to axis_partner = "south side"
        # Equidistant = "equator"

        def bfs_dist(start, adj, V):
            dist={start:0}; q=deque([start])
            while q:
                u=q.popleft()
                for v in adj[u]:
                    if v not in dist:
                        dist[v]=dist[u]+1; q.append(v)
            return dist

        dist_from_trans=bfs_dist(trans_mid, madj, mid)
        dist_from_partner=bfs_dist(axis_partner, madj, mid)

        north=[]; south=[]; equator=[]; axis_nodes=[trans_mid, axis_partner]
        for mi in range(mid):
            if mi in axis_nodes: continue
            dt=dist_from_trans.get(mi, 999)
            dp=dist_from_partner.get(mi, 999)
            if dt<dp: north.append(mi)
            elif dp<dt: south.append(mi)
            else: equator.append(mi)

        print(f"\n  BILATERAL PARTITION around axis {trans_mid}↔{axis_partner}:")
        print(f"    North (closer to transitive): {len(north)} nodes")
        print(f"    South (closer to axis partner): {len(south)} nodes")
        print(f"    Equator (equidistant): {len(equator)} nodes")
        print(f"    Axis: 2 nodes")
        print(f"    Total: {len(north)+len(south)+len(equator)+2} = {mid}")

        # H distribution on each side
        if north:
            north_H=[mn_info[mi]['H'] for mi in north]
            south_H=[mn_info[mi]['H'] for mi in south]
            equator_H=[mn_info[mi]['H'] for mi in equator]
            print(f"\n    North H: min={min(north_H)}, max={max(north_H)}, "
                  f"avg={sum(north_H)/len(north_H):.1f}")
            if south:
                print(f"    South H: min={min(south_H)}, max={max(south_H)}, "
                      f"avg={sum(south_H)/len(south_H):.1f}")
            if equator:
                print(f"    Equator H: min={min(equator_H)}, max={max(equator_H)}, "
                      f"avg={sum(equator_H)/len(equator_H):.1f}")

        # SC/NS distribution on each side
        north_sc=sum(1 for mi in north if mn_info[mi]['sc'])
        south_sc=sum(1 for mi in south if mn_info[mi]['sc'])
        equator_sc=sum(1 for mi in equator if mn_info[mi]['sc'])
        print(f"\n    North SC: {north_sc}/{len(north)}")
        print(f"    South SC: {south_sc}/{len(south)}")
        print(f"    Equator SC: {equator_sc}/{len(equator)}")

        # Cross-axis edges (connecting north to south)
        cross=0; north_set=set(north); south_set=set(south); eq_set=set(equator)
        for(a,b) in me:
            if (a in north_set and b in south_set) or (a in south_set and b in north_set):
                cross+=1
        print(f"\n    Cross-axis edges (north↔south): {cross}")
        print(f"    Within-north edges: {sum(1 for(a,b) in me if a in north_set and b in north_set)}")
        print(f"    Within-south edges: {sum(1 for(a,b) in me if a in south_set and b in south_set)}")

    # Degree distribution from transitive
    print(f"\n  BFS SHELLS from transitive (merged graph):")
    shells=defaultdict(list)
    for mi,d in dist_from_trans.items():
        shells[d].append(mi)
    for d in sorted(shells):
        nodes=shells[d]
        h_vals=[mn_info[mi]['H'] for mi in nodes]
        sc_count=sum(1 for mi in nodes if mn_info[mi]['sc'])
        print(f"    shell {d}: {len(nodes)} nodes, "
              f"H=[{min(h_vals)},{max(h_vals)}], SC={sc_count}")

    # The H_max class(es) — "south pole"
    h_max=max(mn_info[mi]['H'] for mi in range(mid))
    poles=[mi for mi in range(mid) if mn_info[mi]['H']==h_max]
    print(f"\n  SOUTH POLE (H_max={h_max}): {len(poles)} class(es)")
    for mi in poles:
        d=dist_from_trans.get(mi, -1)
        print(f"    mid={mi}, dist_from_transitive={d}, SC={mn_info[mi]['sc']}")

    elapsed=time.time()-t0
    print(f"\n  n={n} done ({elapsed:.1f}s)")

# ============================================================================
# CROSS-n PATTERNS
# ============================================================================

print(f"\n{'='*80}")
print(f"  CROSS-n PATTERN SUMMARY")
print(f"{'='*80}")

print(f"\n  {'n':>3} {'V_m':>6} {'SC_nbrs_of_trans':>18} {'NS_nbrs_of_trans':>18} "
      f"{'axis_H':>8} {'north':>6} {'south':>6} {'equator':>8}")

# Re-run summary (we'd need to store these, but the output above shows the pattern)
print(f"\n  Look at the output above for the bilateral partition data at each n.")
print(f"  Key patterns to notice:")
print(f"  1. The transitive class always has EXACTLY how many SC blue neighbors?")
print(f"  2. Does the bilateral partition become more symmetric as n grows?")
print(f"  3. How does the equator relate to the H=H_max/2 level?")
print(f"  4. Which side has more SC classes?")

print(f"\n{'='*80}")
print(f"  DONE")
print(f"{'='*80}")
