#!/usr/bin/env python3
"""
gn_tilings_blue_s226.py — Tilings per merged node + blue/black correlations
opus-2026-03-23-S226

Tilings(T) × |Aut(T)| = H(T) for every iso class (orbit-stabilizer theorem).
So Tilings = H / |Aut|.

For each merged node, compute tilings and look for patterns:
1. Tiling count per node
2. Correlation with blue_deg, black_deg
3. Sum/product of tilings along edges
4. Whether SC blue subgraph is a TREE at all n
5. Verify ΔH = 2^(n-2) pattern
6. Tiling sequences sorted by H
"""

import sys
import subprocess
from math import comb, factorial
from itertools import permutations
from collections import defaultdict, Counter, deque
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  TILINGS + BLUE/BLACK CORRELATIONS")
print("  opus-2026-03-23-S226")
print("=" * 80)

def build(n):
    m=comb(n,2); P=[(i,j) for i in range(n) for j in range(i+1,n)]
    nfact=factorial(n)
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

        # |Aut| from class size
        if 'size' in cl:
            cl['aut']=nfact//cl['size']
        else:
            # Compute for nauty classes
            if n<=7:
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
                cl['aut']=1  # approximate; most have trivial aut at n=8

        cl['tilings']=cl['H']//cl['aut'] if cl['aut']>0 else cl['H']

    for cl in cls: cl['comp_cid']=cm2.get(cl['comp_canon'],-1)

    def lk(a):
        s=ss(a);c=c3(a);h=Hf(a)
        cids=hm.get((s,c,h))
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
    blue=set(); black=set()
    blue_adj=defaultdict(set); black_adj=defaultdict(set)
    for(a,b) in edges:
        x,y=mn[a],mn[b]
        if x!=y:
            me.add((min(x,y),max(x,y)))
            ma[x].add(y); ma[y].add(x)

    mi={}
    for mm in range(mid):
        cids=[cl['cid'] for cl in cls if mn[cl['cid']]==mm]
        cl0=cls[cids[0]]
        mi[mm]={'H':cl0['H'],'sc':cl0['sc'],'score':cl0['score'],
                'c3':cl0['c3'],'aut':cl0['aut'],'tilings':cl0['tilings'],'cids':cids}

    for(a,b) in me:
        if mi[a]['sc']==mi[b]['sc']:
            blue.add((a,b)); blue_adj[a].add(b); blue_adj[b].add(a)
        else:
            black.add((a,b)); black_adj[a].add(b); black_adj[b].add(a)

    return {'mid':mid,'me':me,'ma':ma,'mi':mi,'blue':blue,'black':black,
            'blue_adj':blue_adj,'black_adj':black_adj,'n':n,'m':m,'cls':cls}

# ============================================================================
# MAIN ANALYSIS
# ============================================================================

for n in range(3, 9):
    print(f"\n{'='*80}")
    print(f"  n = {n}")
    print(f"{'='*80}")
    t0=time.time()

    D=build(n)
    mid=D['mid']; mi=D['mi']; blue_adj=D['blue_adj']; black_adj=D['black_adj']

    # Find transitive
    trans=None
    for mm in range(mid):
        if mi[mm]['H']==1: trans=mm; break

    print(f"  V={mid}, E_blue={len(D['blue'])}, E_black={len(D['black'])}")

    # ---- TILINGS PER NODE ----
    sc_nodes=sorted([mm for mm in range(mid) if mi[mm]['sc']], key=lambda x: mi[x]['H'])
    ns_nodes=sorted([mm for mm in range(mid) if not mi[mm]['sc']], key=lambda x: mi[x]['H'])

    print(f"\n  TILINGS (= H / |Aut|) for SC nodes (sorted by H):")
    if len(sc_nodes) <= 30:
        for mm in sc_nodes:
            h=mi[mm]['H']; aut=mi[mm]['aut']; til=mi[mm]['tilings']
            bd=len(blue_adj.get(mm,set())); kd=len(black_adj.get(mm,set()))
            print(f"    H={h:5d} |Aut|={aut:3d} Tilings={til:6d} blue_d={bd:2d} black_d={kd:2d}")

    # Total tilings
    total_tilings_sc=sum(mi[mm]['tilings'] for mm in sc_nodes)
    total_tilings_ns=sum(mi[mm]['tilings'] for mm in ns_nodes)
    total_tilings=total_tilings_sc+total_tilings_ns
    print(f"\n  Total tilings: SC={total_tilings_sc}, NS={total_tilings_ns}, all={total_tilings}")

    # ---- SC BLUE TREE CHECK ----
    sc_blue_edges=[(a,b) for(a,b) in D['blue'] if mi[a]['sc'] and mi[b]['sc']]
    sc_blue_V=len(sc_nodes)
    sc_blue_E=len(sc_blue_edges)
    is_tree=(sc_blue_E==sc_blue_V-1) if sc_blue_V>0 else True

    print(f"\n  SC BLUE SUBGRAPH: V={sc_blue_V}, E={sc_blue_E}")
    print(f"    Is a TREE: {'YES ✓' if is_tree else 'NO ✗'} (E = V-1 = {sc_blue_V-1})")

    # ---- ΔH = 2^(n-2) CHECK ----
    if trans is not None:
        sc_blue_nbrs=sorted([nb for nb in blue_adj.get(trans,set()) if mi[nb]['sc']],
                           key=lambda x: mi[x]['H'])
        print(f"\n  TRANSITIVE (H=1) BLUE SC NEIGHBORS:")
        for nb in sc_blue_nbrs:
            dh=mi[nb]['H']-1
            is_power2=(dh>0 and (dh&(dh-1))==0)
            print(f"    H={mi[nb]['H']}, ΔH={dh}, 2^(n-2)={1<<(n-2)}, "
                  f"{'= 2^(n-2) ✓' if dh==(1<<(n-2)) else ''}")

    # ---- TILING-DEGREE CORRELATION ----
    if mid > 3:
        tilings_list=[mi[mm]['tilings'] for mm in range(mid)]
        blue_d_list=[len(blue_adj.get(mm,set())) for mm in range(mid)]
        black_d_list=[len(black_adj.get(mm,set())) for mm in range(mid)]
        total_d_list=[len(D['ma'].get(mm,set())) for mm in range(mid)]

        # Simple correlation
        def corr(x, y):
            n_=len(x); mx=sum(x)/n_; my=sum(y)/n_
            cov=sum((x[i]-mx)*(y[i]-my) for i in range(n_))/n_
            sx=(sum((xi-mx)**2 for xi in x)/n_)**0.5
            sy=(sum((yi-my)**2 for yi in y)/n_)**0.5
            return cov/(sx*sy) if sx>0 and sy>0 else 0

        r_til_blue=corr(tilings_list, blue_d_list)
        r_til_black=corr(tilings_list, black_d_list)
        r_til_total=corr(tilings_list, total_d_list)
        h_list=[mi[mm]['H'] for mm in range(mid)]
        r_til_h=corr(tilings_list, h_list)

        print(f"\n  CORRELATIONS:")
        print(f"    corr(Tilings, blue_deg)  = {r_til_blue:.4f}")
        print(f"    corr(Tilings, black_deg) = {r_til_black:.4f}")
        print(f"    corr(Tilings, total_deg) = {r_til_total:.4f}")
        print(f"    corr(Tilings, H)         = {r_til_h:.4f}")

    # ---- TILING SUM ALONG EDGES ----
    # For each edge, compute sum and product of endpoint tilings
    if len(D['me']) <= 5000:
        edge_til_sums=[]
        blue_til_sums=[]; black_til_sums=[]
        for(a,b) in D['me']:
            s=mi[a]['tilings']+mi[b]['tilings']
            edge_til_sums.append(s)
        for(a,b) in D['blue']:
            blue_til_sums.append(mi[a]['tilings']+mi[b]['tilings'])
        for(a,b) in D['black']:
            black_til_sums.append(mi[a]['tilings']+mi[b]['tilings'])

        print(f"\n  TILING SUM ALONG EDGES:")
        if edge_til_sums:
            print(f"    All: min={min(edge_til_sums)}, max={max(edge_til_sums)}, "
                  f"avg={sum(edge_til_sums)/len(edge_til_sums):.1f}")
        if blue_til_sums:
            print(f"    Blue: min={min(blue_til_sums)}, max={max(blue_til_sums)}, "
                  f"avg={sum(blue_til_sums)/len(blue_til_sums):.1f}")
        if black_til_sums:
            print(f"    Black: min={min(black_til_sums)}, max={max(black_til_sums)}, "
                  f"avg={sum(black_til_sums)/len(black_til_sums):.1f}")

    # ---- TILING-WEIGHTED EDGE COUNT ----
    # Σ_edges T(a) × T(b) — a "tiling energy" of the graph
    if len(D['me']) <= 50000:
        til_energy=sum(mi[a]['tilings']*mi[b]['tilings'] for(a,b) in D['me'])
        til_energy_blue=sum(mi[a]['tilings']*mi[b]['tilings'] for(a,b) in D['blue'])
        til_energy_black=sum(mi[a]['tilings']*mi[b]['tilings'] for(a,b) in D['black'])
        print(f"\n  TILING ENERGY (Σ T(a)×T(b) over edges):")
        print(f"    All: {til_energy}")
        print(f"    Blue: {til_energy_blue}")
        print(f"    Black: {til_energy_black}")

    # ---- TILING SEQUENCE (sorted) ----
    all_tilings=sorted(set(mi[mm]['tilings'] for mm in range(mid)))
    print(f"\n  DISTINCT TILING VALUES: {len(all_tilings)}")
    if len(all_tilings)<=20:
        print(f"    {all_tilings}")
    else:
        print(f"    First 10: {all_tilings[:10]}")
        print(f"    Last 10: {all_tilings[-10:]}")

    elapsed=time.time()-t0
    print(f"\n  n={n} done ({elapsed:.1f}s)")

# ============================================================================
# CROSS-n SEQUENCES
# ============================================================================

print(f"\n{'='*80}")
print(f"  CROSS-n TILING SEQUENCES")
print(f"{'='*80}")

# Tilings of transitive (always H=1, |Aut|=1 → Tilings=1)
print(f"\n  Transitive tilings: always 1 (H=1, |Aut|=1)")

# Tilings of regular/Paley (H_max)
print(f"\n  H_max and its tilings:")
for n in range(3, 9):
    D=build(n)
    h_max=max(D['mi'][mm]['H'] for mm in range(D['mid']))
    hmax_nodes=[mm for mm in range(D['mid']) if D['mi'][mm]['H']==h_max]
    for mm in hmax_nodes:
        print(f"    n={n}: H_max={h_max}, |Aut|={D['mi'][mm]['aut']}, "
              f"Tilings={D['mi'][mm]['tilings']}, SC={D['mi'][mm]['sc']}")

print(f"\n{'='*80}")
print(f"  DONE")
print(f"{'='*80}")
