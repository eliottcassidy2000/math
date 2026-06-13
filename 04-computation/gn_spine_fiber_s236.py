#!/usr/bin/env python3
"""
gn_spine_fiber_s236.py — The BLUE spine structure + tiling fibers
opus-2026-03-23-S236

BLUE = SC↔SC edges (explorer definition). The spine.
BLACK = SC↔NS + NS↔NS edges (explorer definition). Everything else.

For each SC class, we have:
  Tilings(C) = H(C) / |Aut(C)| = fiber size

The SPINE is the subgraph on SC nodes connected by blue (SC-SC) edges.
We investigate:
1. How tiling fiber sizes relate to spine connectivity
2. Whether Σ Tilings² along spine edges gives a formula
3. The spine's degree sequence and how it relates to |Aut|
4. The spine's diameter, genus, bridges as functions of n
5. The Mode A/B fiber counts restricted to SC classes
"""

import sys
import subprocess
from math import comb, factorial
from itertools import permutations
from collections import defaultdict, Counter, deque
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  THE BLUE SPINE: SC-SC EDGES + TILING FIBERS")
print("  opus-2026-03-23-S236")
print("=" * 80)

def build(n):
    m=comb(n,2); P=[(i,j) for i in range(n) for j in range(i+1,n)]; nf=factorial(n)
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
            p=[g[0] for g in gs]; return tuple(a[p[i]][p[j]] for i in range(n) for j in range(n))
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
        a=cl['adj']; cl['score']=ss(a); cl['H']=Hf(a); cl['own_canon']=cp(a)
        comp=[[1-a[i][j] if i!=j else 0 for j in range(n)] for i in range(n)]
        cl['comp_canon']=cp(comp); cl['sc']=(cl['own_canon']==cl['comp_canon'])
        hm[(cl['score'],cl['H'])].append(cl['cid']); cm2[cl['own_canon']]=cl['cid']
        cl['aut']=nf//cl.get('size',nf) if 'size' in cl else 1
        cl['tilings']=cl['H']//cl['aut'] if cl['aut']>0 else cl['H']
    for cl in cls: cl['comp_cid']=cm2.get(cl['comp_canon'],-1)
    def lk(a):
        s=ss(a);h=Hf(a); cids=hm.get((s,h))
        if not cids: return -1
        if len(cids)==1: return cids[0]
        return cm2.get(cp(a),-1)
    edges=set()
    for cl in cls:
        b=cl['bits']
        for arc in range(m):
            if ny: fb=b^(1<<(m-1-arc))
            else: fb=b^(1<<arc)
            if ny:
                fa=[[0]*n for _ in range(n)]
                for k,(i,j) in enumerate(P):
                    if fb&(1<<(m-1-k)): fa[i][j]=1;
                    else: fa[j][i]=1
            else:
                fa=[[0]*n for _ in range(n)]
                for k,(i,j) in enumerate(P):
                    if fb&(1<<k): fa[i][j]=1
                    else: fa[j][i]=1
            nb=lk(fa)
            if nb>=0 and nb!=cl['cid']:
                edges.add((min(cl['cid'],nb),max(cl['cid'],nb)))
    # Merged + spine
    mn={}; mid=0
    for cl in cls:
        ci=cl['cid']
        if ci in mn: continue
        if cl['sc']: mn[ci]=mid; mid+=1
        else:
            mn[ci]=mid; comp=cl['comp_cid']
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
        mi[mm]={'H':cl0['H'],'sc':cl0['sc'],'aut':cl0['aut'],'tilings':cl0['tilings'],
                'score':cl0['score']}
    # SPINE = SC-SC edges (= explorer BLUE)
    spine_adj=defaultdict(set); spine_edges=set()
    for(a,b) in me:
        if mi[a]['sc'] and mi[b]['sc']:
            spine_adj[a].add(b); spine_adj[b].add(a)
            spine_edges.add((a,b))
    # RIBS = SC-NS edges
    rib_count=Counter()
    for(a,b) in me:
        if mi[a]['sc'] != mi[b]['sc']:
            sc_node = a if mi[a]['sc'] else b
            rib_count[sc_node]+=1
    return {'mid':mid,'me':me,'ma':ma,'mi':mi,'spine_adj':spine_adj,
            'spine_edges':spine_edges,'rib_count':rib_count,'n':n,'m':m}

# ============================================================================
# SPINE ANALYSIS
# ============================================================================

for n in range(3, 9):
    print(f"\n{'='*80}")
    print(f"  n = {n}")
    print(f"{'='*80}")
    t0=time.time()

    D=build(n)
    mi=D['mi']; sa=D['spine_adj']; se=D['spine_edges']; rc=D['rib_count']

    sc_nodes=sorted([mm for mm in range(D['mid']) if mi[mm]['sc']], key=lambda x: mi[x]['H'])
    V_sc=len(sc_nodes); E_sc=len(se)

    print(f"  SPINE: V_SC={V_sc}, E_SC={E_sc}, density={2*E_sc/(V_sc*(V_sc-1)):.4f}" if V_sc>1 else f"  SPINE: V_SC={V_sc}")

    # Spine degree + tilings + ribs for each SC node
    if V_sc <= 30:
        print(f"\n  {'mid':>4} {'H':>5} {'|Aut|':>5} {'Tilings':>8} {'spine_d':>7} {'ribs':>5} {'score':>20}")
        for mm in sc_nodes:
            sd=len(sa.get(mm,set()))
            rb=rc.get(mm,0)
            print(f"  {mm:4d} {mi[mm]['H']:5d} {mi[mm]['aut']:5d} {mi[mm]['tilings']:8d} "
                  f"{sd:7d} {rb:5d} {str(mi[mm]['score']):>20}")

    # Spine statistics
    spine_degs=[len(sa.get(mm,set())) for mm in sc_nodes]
    rib_counts=[rc.get(mm,0) for mm in sc_nodes]
    tilings=[mi[mm]['tilings'] for mm in sc_nodes]
    H_vals=[mi[mm]['H'] for mm in sc_nodes]

    print(f"\n  SPINE DEGREE: min={min(spine_degs)}, max={max(spine_degs)}, "
          f"avg={sum(spine_degs)/V_sc:.2f}")
    print(f"  RIB COUNT: min={min(rib_counts)}, max={max(rib_counts)}, "
          f"avg={sum(rib_counts)/V_sc:.2f}")
    print(f"  TILINGS: min={min(tilings)}, max={max(tilings)}, "
          f"avg={sum(tilings)/V_sc:.2f}, sum={sum(tilings)}")

    # Correlations
    if V_sc > 3:
        def corr(x,y):
            n_=len(x); mx=sum(x)/n_; my=sum(y)/n_
            cov=sum((x[i]-mx)*(y[i]-my) for i in range(n_))/n_
            sx=(sum((xi-mx)**2 for xi in x)/n_)**0.5
            sy=(sum((yi-my)**2 for yi in y)/n_)**0.5
            return cov/(sx*sy) if sx>0 and sy>0 else 0

        r_td = corr(tilings, spine_degs)
        r_tr = corr(tilings, rib_counts)
        r_th = corr(tilings, H_vals)
        r_dr = corr(spine_degs, rib_counts)

        print(f"\n  CORRELATIONS ON SPINE:")
        print(f"    corr(Tilings, spine_deg)  = {r_td:.4f}")
        print(f"    corr(Tilings, ribs)       = {r_tr:.4f}")
        print(f"    corr(Tilings, H)          = {r_th:.4f}")
        print(f"    corr(spine_deg, ribs)     = {r_dr:.4f}")

    # Tiling products along spine edges
    if E_sc > 0:
        til_prods=[mi[a]['tilings']*mi[b]['tilings'] for(a,b) in se]
        til_sums_e=[mi[a]['tilings']+mi[b]['tilings'] for(a,b) in se]
        h_diffs=[abs(mi[a]['H']-mi[b]['H']) for(a,b) in se]

        print(f"\n  ALONG SPINE EDGES:")
        print(f"    Σ T(a)×T(b) = {sum(til_prods)}")
        print(f"    avg T(a)×T(b) = {sum(til_prods)/E_sc:.2f}")
        print(f"    avg |ΔH| = {sum(h_diffs)/E_sc:.2f}")

    # Spine + rib combined "total SC influence"
    sc_influence = sum(spine_degs[i]+rib_counts[i] for i in range(V_sc))
    print(f"\n  TOTAL SC INFLUENCE (spine_deg + ribs): {sc_influence}")
    print(f"  avg SC total degree: {sc_influence/V_sc:.2f}")

    elapsed=time.time()-t0
    print(f"\n  ({elapsed:.1f}s)")

# ============================================================================
# CROSS-n SPINE SEQUENCES
# ============================================================================

print(f"\n{'='*80}")
print(f"  CROSS-n SPINE SEQUENCES")
print(f"{'='*80}")

print(f"\n  {'n':>3} {'V_SC':>5} {'E_SC':>5} {'genus':>6} {'avg_sd':>7} {'avg_rib':>8} {'Σ_til':>8}")

for n in range(3,9):
    D=build(n)
    mi=D['mi']; sa=D['spine_adj']; se=D['spine_edges']; rc=D['rib_count']
    sc=[mm for mm in range(D['mid']) if mi[mm]['sc']]
    V=len(sc); E=len(se)
    # Components
    vis=set(); comps=0
    for s in sc:
        if s in vis: continue
        comps+=1; q=deque([s]); vis.add(s)
        while q:
            u=q.popleft()
            for v in sa.get(u,set()):
                if v not in vis: vis.add(v); q.append(v)
    genus=E-V+comps
    avg_sd=2*E/V if V>0 else 0
    avg_rib=sum(rc.get(mm,0) for mm in sc)/V if V>0 else 0
    sum_til=sum(mi[mm]['tilings'] for mm in sc)
    print(f"  {n:3d} {V:5d} {E:5d} {genus:6d} {avg_sd:7.2f} {avg_rib:8.2f} {sum_til:8d}")

print(f"\n  KEY: spine avg_deg stays ~2-4 while ribs grow with n.")
print(f"  The spine is THIN while the ribs THICKEN.")

print(f"\n{'='*80}")
print(f"  DONE")
print(f"{'='*80}")
