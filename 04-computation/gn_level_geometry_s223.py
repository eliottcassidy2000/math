#!/usr/bin/env python3
"""
gn_level_geometry_s223.py — Level-by-level geometry of G_n/Z_2
opus-2026-03-23-S223

The user's key insight: the principal blue line (transitive → first SC neighbor)
defines an AXIS. Structures "perpendicular" to this axis appear at each level:

n=4: The NS merged node connects to BOTH SC poles via black edges,
     forming a line perpendicular to the principal blue line.

n=5: The two merged-NS nodes connect via black lines, forming
     perpendicular structure.

For each n, we analyze:
1. The principal blue line (SC backbone path from H=1)
2. At each H-level, what is the "cross-section" perpendicular to the axis?
3. How do self-loops (blueself, blackself) distribute along the axis?
4. What is the "width profile" — how many nodes at each distance from the axis?
5. Level edges as "horizontal" structure within each cross-section
"""

import sys
import subprocess
from math import comb, factorial
from itertools import permutations
from collections import defaultdict, Counter, deque
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  LEVEL-BY-LEVEL GEOMETRY OF G_n/Z_2")
print("  opus-2026-03-23-S223")
print("=" * 80)

# ============================================================================
# Build (compact, from S222)
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
    for cl in cls: cl['comp_cid']=cm2.get(cl['comp_canon'],-1)

    def lk(a):
        s=ss(a);c=c3(a);h=Hf(a)
        cids=hm.get((s,c,h))
        if not cids: return -1
        if len(cids)==1: return cids[0]
        return cm2.get(cp(a),-1)

    # Build edges + detect self-loops
    edges=set(); al=defaultdict(set); self_loop=set()
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
            if nb==cl['cid']:
                self_loop.add(cl['cid'])
            elif nb>=0:
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
        cl0=cls[cids[0]]
        # Self-loop: merged node has self-loop if any original class in it does
        has_sl=any(ci in self_loop for ci in cids)
        mi[mm]={'H':cl0['H'],'sc':cl0['sc'],'score':cl0['score'],
                'c3':cl0['c3'],'self_loop':has_sl,'cids':cids}

    return {'cls':cls,'mn':mn,'mid':mid,'me':me,'ma':ma,'mi':mi,'n':n,'m':m}

# ============================================================================
# LEVEL-BY-LEVEL ANALYSIS
# ============================================================================

for n in range(3, 9):
    print(f"\n{'='*80}")
    print(f"  n = {n}: LEVEL-BY-LEVEL GEOMETRY")
    print(f"{'='*80}")
    t0=time.time()

    D=build(n)
    mid=D['mid']; me=D['me']; ma=D['ma']; mi=D['mi']; m=D['m']
    print(f"  G_{n}/Z_2: V={mid}, E={len(me)} ({time.time()-t0:.1f}s)")

    # Group by H
    h_levels=defaultdict(list)
    for mm in range(mid): h_levels[mi[mm]['H']].append(mm)
    h_sorted=sorted(h_levels.keys())

    # Find transitive
    trans=None
    for mm in range(mid):
        if mi[mm]['H']==1: trans=mm; break

    # BFS distance from transitive
    dist_t={}
    if trans is not None:
        dist_t={trans:0}; q=deque([trans])
        while q:
            u=q.popleft()
            for v in ma[u]:
                if v not in dist_t: dist_t[v]=dist_t[u]+1; q.append(v)

    # ---- CROSS-SECTION AT EACH H-LEVEL ----
    print(f"\n  H-LEVEL CROSS-SECTIONS (the 'perpendicular' structure):")
    print(f"  {'H':>5} {'#nodes':>6} {'#SC':>4} {'#NS':>4} {'#SL':>4} "
          f"{'lvl_E':>5} {'inter_E':>7} {'avg_dist':>8} {'scores':>8}")

    for h in h_sorted:
        nodes=h_levels[h]
        n_sc=sum(1 for mm in nodes if mi[mm]['sc'])
        n_ns=len(nodes)-n_sc
        n_sl=sum(1 for mm in nodes if mi[mm]['self_loop'])

        # Level edges within this H-level
        lvl_e=0
        for(a,b) in me:
            if mi[a]['H']==h and mi[b]['H']==h: lvl_e+=1

        # Inter-level edges (connecting to different H levels)
        inter_e=0
        for mm in nodes:
            for nb in ma[mm]:
                if mi[nb]['H']!=h: inter_e+=1

        # Average distance from transitive
        dists=[dist_t.get(mm,999) for mm in nodes]
        avg_d=sum(dists)/len(dists) if dists else 0

        # Distinct scores
        scores=len(set(mi[mm]['score'] for mm in nodes))

        if len(h_sorted) <= 25 or (len(nodes) >= 2):  # print interesting levels
            print(f"  {h:5d} {len(nodes):6d} {n_sc:4d} {n_ns:4d} {n_sl:4d} "
                  f"{lvl_e:5d} {inter_e:7d} {avg_d:8.1f} {scores:8d}")

    # ---- THE PRINCIPAL LINE vs PERPENDICULAR ----
    print(f"\n  PRINCIPAL LINE ANALYSIS:")

    # Find SC backbone path from transitive (greedy H-increasing)
    sc_adj=defaultdict(set)
    for(a,b) in me:
        if mi[a]['sc'] and mi[b]['sc']:
            sc_adj[a].add(b); sc_adj[b].add(a)

    if trans is not None:
        # Greedy increasing path from transitive along SC backbone
        principal=[trans]; current=trans
        while True:
            nbrs=[nb for nb in sc_adj[current]
                  if mi[nb]['H']>mi[current]['H'] and nb not in set(principal)]
            if not nbrs: break
            # Pick the neighbor with smallest H (stay close to axis)
            best=min(nbrs, key=lambda x: mi[x]['H'])
            principal.append(best); current=best

        print(f"  Principal path (greedy H-increasing along SC backbone):")
        print(f"    Length: {len(principal)} nodes")
        print(f"    H values: {[mi[mm]['H'] for mm in principal]}")
        print(f"    All SC: {all(mi[mm]['sc'] for mm in principal)}")

        # At each node on the principal line, what is the "perpendicular" structure?
        print(f"\n  PERPENDICULAR STRUCTURE at each principal node:")
        for pi, mm in enumerate(principal):
            h=mi[mm]['H']
            # All nodes at this H level
            level_nodes=h_levels[h]
            # Partition: on the principal line vs off it
            on_line=set(principal) & set(level_nodes)
            off_line=[x for x in level_nodes if x not in on_line]

            # Neighbors of mm that are NOT on the principal line and at same H
            perp_nbrs=[nb for nb in ma[mm] if nb in set(level_nodes) and nb!=mm]

            # All connections at this level
            sc_at_level=[x for x in level_nodes if mi[x]['sc']]
            ns_at_level=[x for x in level_nodes if not mi[x]['sc']]

            # Black connections from mm
            black_nbrs=[nb for nb in ma[mm] if not mi[nb]['sc']]  # mm is SC
            blue_nbrs=[nb for nb in ma[mm] if mi[nb]['sc'] and nb!=mm]

            # Self-loop status
            sl='SL' if mi[mm]['self_loop'] else '  '

            if len(principal) <= 15:
                print(f"    [{pi}] mid={mm:3d} H={h:4d} {sl} "
                      f"| level: {len(level_nodes)} ({len(sc_at_level)}SC+{len(ns_at_level)}NS) "
                      f"| perp_nbrs={len(perp_nbrs)} level_E={sum(1 for(a,b) in me if mi[a]['H']==h and mi[b]['H']==h)} "
                      f"| blue→SC:{len(blue_nbrs)} black→NS:{len(black_nbrs)}")

    # ---- SELF-LOOP PATTERN ----
    print(f"\n  SELF-LOOP PATTERN ALONG PRINCIPAL LINE:")
    if trans is not None and len(principal) <= 20:
        sl_pattern=''.join('●' if mi[mm]['self_loop'] else '○' for mm in principal)
        print(f"    {sl_pattern}")
        print(f"    (● = self-loop, ○ = no self-loop)")

    # ---- SYMMETRY OF H-LEVELS AROUND CENTER ----
    h_center=(h_sorted[0]+h_sorted[-1])/2
    print(f"\n  H-LEVEL SYMMETRY around center H={h_center:.1f}:")
    for i in range(min(len(h_sorted)//2, 10)):
        h_lo=h_sorted[i]
        h_hi=h_sorted[-(i+1)]
        n_lo=len(h_levels[h_lo])
        n_hi=len(h_levels[h_hi])
        sc_lo=sum(1 for mm in h_levels[h_lo] if mi[mm]['sc'])
        sc_hi=sum(1 for mm in h_levels[h_hi] if mi[mm]['sc'])
        print(f"    H={h_lo:4d} ({n_lo:3d} nodes, {sc_lo} SC)  ←→  "
              f"H={h_hi:4d} ({n_hi:3d} nodes, {sc_hi} SC)  "
              f"{'MATCH' if n_lo==n_hi else f'diff={n_hi-n_lo}'}")

    # ---- EDGE FLOW: how many edges go UP vs ACROSS ----
    up_edges=0; across_edges=0; level_edges_total=0
    for(a,b) in me:
        ha,hb=mi[a]['H'],mi[b]['H']
        if ha==hb: level_edges_total+=1
        else:
            up_edges+=1

    print(f"\n  EDGE FLOW:")
    print(f"    Vertical (H-changing): {up_edges} ({100*up_edges/len(me):.1f}%)")
    print(f"    Horizontal (H-preserving): {level_edges_total} ({100*level_edges_total/len(me):.1f}%)")

    # ---- WIDEST CROSS-SECTION ----
    widest_h=max(h_sorted, key=lambda h: len(h_levels[h]))
    widest_nodes=h_levels[widest_h]
    print(f"\n  WIDEST LEVEL: H={widest_h} with {len(widest_nodes)} nodes")
    print(f"    SC: {sum(1 for mm in widest_nodes if mi[mm]['sc'])}")
    print(f"    NS: {sum(1 for mm in widest_nodes if not mi[mm]['sc'])}")
    # Internal structure of widest level
    widest_internal=sum(1 for(a,b) in me
                       if mi[a]['H']==widest_h and mi[b]['H']==widest_h)
    print(f"    Internal edges: {widest_internal}")
    if len(widest_nodes)>1:
        max_internal=len(widest_nodes)*(len(widest_nodes)-1)//2
        print(f"    Internal density: {widest_internal}/{max_internal} = "
              f"{widest_internal/max_internal:.3f}" if max_internal>0 else "")

    elapsed=time.time()-t0
    print(f"\n  n={n} done ({elapsed:.1f}s)")

print(f"\n{'='*80}")
print(f"  DONE")
print(f"{'='*80}")
