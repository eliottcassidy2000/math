#!/usr/bin/env python3
"""
gn_fiber_algebra_s233.py — Fiber algebra of tilings and the metagraph
opus-2026-03-23-S233

FUNDAMENTAL IDENTITY: Σ_C Tilings(C) = 2^{C(n-1,2)} where Tilings(C) = H(C)/|Aut(C)|

The tiling count is the FIBER SIZE of class C in the map:
  {0,1}^{C(n-1,2)} → tournaments → iso classes

For the merged metagraph:
  Total merged tilings = 2^{C(n-1,2)-1} + (1/2) × SC_tilings_sum

We compute:
1. Tiling sums by type (SC, NS, all, merged)
2. H-weighted sums
3. Tiling × degree products
4. Whether tiling products along edges give formulas
5. The "tiling Laplacian" perspective
"""

import sys
import subprocess
from math import comb, factorial, log2
from itertools import permutations
from collections import defaultdict, Counter
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  FIBER ALGEBRA: TILINGS AND THE METAGRAPH")
print("  opus-2026-03-23-S233")
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
        a=cl['adj']; cl['score']=ss(a); cl['H']=Hf(a)
        cl['own_canon']=cp(a)
        comp=[[1-a[i][j] if i!=j else 0 for j in range(n)] for i in range(n)]
        cl['comp_canon']=cp(comp); cl['sc']=(cl['own_canon']==cl['comp_canon'])
        hm[(cl['score'],cl['H'])].append(cl['cid'])
        cm2[cl['own_canon']]=cl['cid']
        if 'size' in cl:
            cl['aut']=nfact//cl['size']
        else:
            cl['aut']=1
        cl['tilings']=cl['H']//cl['aut'] if cl['aut']>0 else cl['H']
    for cl in cls: cl['comp_cid']=cm2.get(cl['comp_canon'],-1)
    def lk(a):
        s=ss(a);h=Hf(a)
        cids=hm.get((s,h))
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
    ba=defaultdict(set); ka=defaultdict(set)
    for(a,b) in edges:
        x,y=mn[a],mn[b]
        if x!=y:
            me.add((min(x,y),max(x,y))); ma[x].add(y); ma[y].add(x)
    mi={}
    for mm in range(mid):
        cids=[cl['cid'] for cl in cls if mn[cl['cid']]==mm]
        cl0=cls[cids[0]]
        mi[mm]={'H':cl0['H'],'sc':cl0['sc'],'aut':cl0['aut'],'tilings':cl0['tilings']}
    for(a,b) in me:
        if mi[a]['sc']==mi[b]['sc']: ba[a].add(b); ba[b].add(a)
        else: ka[a].add(b); ka[b].add(a)
    return {'mid':mid,'me':me,'ma':ma,'mi':mi,'ba':ba,'ka':ka,'n':n,'m':m,'cls':cls}

# ============================================================================
# FIBER ALGEBRA
# ============================================================================

for n in range(3, 9):
    print(f"\n{'='*80}")
    print(f"  n = {n}")
    print(f"{'='*80}")
    t0=time.time()

    D=build(n)
    mid=D['mid']; mi=D['mi']; me=D['me']; ma=D['ma']; ba=D['ba']; ka=D['ka']
    m=D['m']; cls=D['cls']

    # Tiling sums
    total_tilings = sum(mi[mm]['tilings'] for mm in range(mid))
    sc_tilings = sum(mi[mm]['tilings'] for mm in range(mid) if mi[mm]['sc'])
    ns_tilings = total_tilings - sc_tilings

    # Expected: total tilings over ALL classes = 2^{C(n-1,2)}
    expected_all = sum(cl['tilings'] for cl in cls)
    cn1_2 = comb(n-1, 2)
    expected_pow2 = 2**cn1_2

    # For merged: SC classes count once, NS pairs count once each
    # total_merged_tilings = sc_tilings + ns_tilings (NS already merged)
    # But ns_tilings counts each NS PAIR once (both halves have same tiling count)
    # So total_merged_tilings = sc_tilings + ns_tilings

    print(f"  V_merged={mid}, E_merged={len(me)}, m={m}")
    print(f"\n  TILING SUMS:")
    print(f"    Σ_all_classes Tilings(C) = {expected_all}")
    print(f"    2^C(n-1,2) = 2^{cn1_2} = {expected_pow2}")
    print(f"    MATCH: {'YES ✓' if expected_all == expected_pow2 else 'NO ✗'}")

    print(f"\n    Σ_merged Tilings = {total_tilings}")
    print(f"    Σ_SC Tilings = {sc_tilings}")
    print(f"    Σ_NS Tilings = {ns_tilings}")
    print(f"    2^(C(n-1,2)-1) = {expected_pow2//2}")
    print(f"    (2^(C(n-1,2)) + SC_tilings) / 2 = {(expected_pow2 + sc_tilings*2 - sc_tilings)//1}")
    # Actually: merged_total = sc_tilings + ns_tilings
    # All_total = sc_tilings + 2*ns_tilings (each NS pair counted twice)
    # So: ns_tilings = (all_total - sc_tilings) / 2
    ns_from_all = (expected_all - sum(cl['tilings'] for cl in cls if cl['sc'])) // 2
    merged_from_formula = sum(cl['tilings'] for cl in cls if cl['sc']) + ns_from_all
    print(f"    Merged formula: {merged_from_formula}")
    print(f"    = (2^C(n-1,2) + SC_tilings_all) / 2 = ({expected_all} + {sum(cl['tilings'] for cl in cls if cl['sc'])}) / 2 = {(expected_all + sum(cl['tilings'] for cl in cls if cl['sc']))//2}")

    # H-weighted sums
    sum_H = sum(mi[mm]['H'] for mm in range(mid))
    sum_H_sc = sum(mi[mm]['H'] for mm in range(mid) if mi[mm]['sc'])
    sum_H_ns = sum_H - sum_H_sc
    sum_H2 = sum(mi[mm]['H']**2 for mm in range(mid))

    print(f"\n  H SUMS:")
    print(f"    Σ H = {sum_H}")
    print(f"    Σ H (SC) = {sum_H_sc}")
    print(f"    Σ H (NS) = {sum_H_ns}")
    print(f"    Σ H² = {sum_H2}")
    print(f"    avg H = {sum_H/mid:.2f}")

    # Tiling-weighted edge sums
    if len(me) <= 100000:
        sum_TT = sum(mi[a]['tilings']*mi[b]['tilings'] for(a,b) in me)
        sum_TT_blue = sum(mi[a]['tilings']*mi[b]['tilings'] for(a,b) in me if mi[a]['sc']==mi[b]['sc'])
        sum_TT_black = sum_TT - sum_TT_blue

        sum_Tsum = sum(mi[a]['tilings']+mi[b]['tilings'] for(a,b) in me)

        print(f"\n  TILING-WEIGHTED EDGE SUMS:")
        print(f"    Σ T(a)×T(b) = {sum_TT}")
        print(f"    Σ T(a)×T(b) (blue) = {sum_TT_blue}")
        print(f"    Σ T(a)×T(b) (black) = {sum_TT_black}")
        print(f"    Σ (T(a)+T(b)) = {sum_Tsum}")
        print(f"    avg T(a)×T(b) = {sum_TT/len(me):.2f}" if len(me)>0 else "")
        print(f"    avg (T(a)+T(b)) = {sum_Tsum/len(me):.2f}" if len(me)>0 else "")

        # Ratio: Σ T(a)×T(b) / (total_tilings^2 / (2*E))
        # If edges were random: E[T(a)T(b)] ≈ (avg T)^2
        avg_T = total_tilings / mid if mid > 0 else 0
        expected_random = avg_T**2 * len(me)
        print(f"    Random model Σ T×T: {expected_random:.0f}")
        print(f"    Actual/Random ratio: {sum_TT/expected_random:.4f}" if expected_random>0 else "")

    # Degree-tiling correlation in merged
    if mid > 3:
        til_list = [mi[mm]['tilings'] for mm in range(mid)]
        deg_list = [len(ma.get(mm,set())) for mm in range(mid)]
        bdeg_list = [len(ba.get(mm,set())) for mm in range(mid)]

        # Compute Σ deg(v) × T(v)
        sum_dT = sum(deg_list[mm]*til_list[mm] for mm in range(mid))
        sum_bdT = sum(bdeg_list[mm]*til_list[mm] for mm in range(mid))

        print(f"\n  DEGREE-TILING PRODUCTS:")
        print(f"    Σ deg(v)×T(v) = {sum_dT}")
        print(f"    Σ blue_deg(v)×T(v) = {sum_bdT}")
        print(f"    2×Σ T×T over edges = 2×{sum_TT if len(me)<=100000 else '?'}")
        # By handshaking: Σ_v deg(v)×T(v) ≠ Σ_edges T(a)×T(b) in general
        # But: Σ_edges (T(a)+T(b)) = Σ_v deg(v)×T(v)
        print(f"    CHECK: Σ(T(a)+T(b)) = {sum_Tsum}, Σ deg×T = {sum_dT}, "
              f"{'MATCH ✓' if sum_Tsum == sum_dT else 'MISMATCH ✗'}")

    # Key derived quantities
    print(f"\n  KEY RATIOS:")
    print(f"    Σ_merged T / E = {total_tilings/len(me):.4f}" if len(me)>0 else "")
    print(f"    (Σ_merged T)^2 / (mid × E) = {total_tilings**2/(mid*len(me)):.4f}" if len(me)>0 else "")
    print(f"    E × avg_T^2 / Σ T×T = {len(me)*avg_T**2/sum_TT:.4f}" if len(me)<=100000 and sum_TT>0 else "")

    elapsed=time.time()-t0
    print(f"\n  n={n} done ({elapsed:.1f}s)")

# ============================================================================
# CROSS-n SEQUENCES
# ============================================================================

print(f"\n{'='*80}")
print(f"  CROSS-n TILING ALGEBRA SEQUENCES")
print(f"{'='*80}")

print(f"\n  FUNDAMENTAL IDENTITY: Σ_all Tilings = 2^C(n-1,2)")
for n in range(3,9):
    print(f"    n={n}: 2^{comb(n-1,2)} = {2**comb(n-1,2)}")

print(f"\n  THE TILING SUM AS EDGE FORMULA PROXY:")
print(f"  If E ∝ (Σ_merged T)^α for some α, what is α?")

print(f"\n{'='*80}")
print(f"  DONE")
print(f"{'='*80}")
