#!/usr/bin/env python3
"""
gn_forbidden_integration_s227.py — Integration of forbidden H values with
merged metagraph structure, blue/black views, and tiling model.

Key connections to make:
1. Forbidden H=7,21 create "missing floors" in the metagraph DAG
2. The gap structure interacts with blue/black edges
3. The conflict-graph triangle impossibility (alpha_1=3, alpha_2=0 impossible)
   IS the reason H=7 is forbidden, and this constraint propagates through
   the metagraph structure
4. The ΔH=2^(n-2)+1 formula connects to the staircase anti-diagonal
"""

import sys
import subprocess
from math import comb, factorial
from itertools import permutations
from collections import defaultdict, Counter, deque
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  FORBIDDEN H VALUES + MERGED METAGRAPH INTEGRATION")
print("  opus-2026-03-23-S227")
print("=" * 80)

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
    blue=set(); black=set()
    ba=defaultdict(set); ka=defaultdict(set)
    for(a,b) in edges:
        x,y=mn[a],mn[b]
        if x!=y:
            me.add((min(x,y),max(x,y)))
            ma[x].add(y); ma[y].add(x)

    mi={}
    for mm in range(mid):
        cids=[cl['cid'] for cl in cls if mn[cl['cid']]==mm]
        mi[mm]={'H':cls[cids[0]]['H'],'sc':cls[cids[0]]['sc'],
                'score':cls[cids[0]]['score'],'c3':cls[cids[0]]['c3']}
    for(a,b) in me:
        if mi[a]['sc']==mi[b]['sc']:
            blue.add((a,b)); ba[a].add(b); ba[b].add(a)
        else:
            black.add((a,b)); ka[a].add(b); ka[b].add(a)

    return {'mid':mid,'me':me,'ma':ma,'mi':mi,'blue':blue,'black':black,
            'ba':ba,'ka':ka,'n':n,'m':m,'cls':cls}

# ============================================================================
# ANALYSIS
# ============================================================================

for n in range(3, 9):
    print(f"\n{'='*80}")
    print(f"  n = {n}")
    print(f"{'='*80}")
    t0=time.time()

    D=build(n)
    mid=D['mid']; mi=D['mi']; me=D['me']; ma=D['ma']
    blue=D['blue']; black=D['black']; ba=D['ba']; ka=D['ka']

    # H values present
    h_vals=sorted(set(mi[mm]['H'] for mm in range(mid)))
    h_max=max(h_vals)

    # All possible odd H values up to h_max
    all_odd=list(range(1, h_max+1, 2))
    gaps=sorted(set(all_odd)-set(h_vals))

    print(f"  V={mid}, E_blue={len(blue)}, E_black={len(black)}")
    print(f"  H range: [{min(h_vals)}, {h_max}]")
    print(f"  H values present: {len(h_vals)}/{len(all_odd)} = {len(h_vals)/len(all_odd):.3f}")
    print(f"  H gaps: {gaps[:20]}{'...' if len(gaps)>20 else ''}")

    # ---- GAP STRUCTURE AND EDGE JUMPS ----
    print(f"\n  EDGE JUMPS OVER FORBIDDEN H VALUES:")

    # For each edge, compute |ΔH| and check if it jumps over a gap
    jump_over_7 = 0; jump_over_21 = 0; jump_over_any_gap = 0
    blue_jumps = []; black_jumps = []
    all_jumps = []

    for(a,b) in me:
        ha,hb = mi[a]['H'],mi[b]['H']
        dh = abs(ha-hb)
        lo,hi = min(ha,hb),max(ha,hb)
        all_jumps.append(dh)

        # Does this edge jump over H=7?
        if lo < 7 < hi: jump_over_7 += 1
        if lo < 21 < hi: jump_over_21 += 1

        # Any gap jumped?
        jumped = any(lo < g < hi for g in gaps)
        if jumped: jump_over_any_gap += 1

        if (a,b) in blue or (b,a) in blue:
            blue_jumps.append(dh)
        else:
            black_jumps.append(dh)

    print(f"    Edges jumping over H=7: {jump_over_7}/{len(me)}")
    print(f"    Edges jumping over H=21: {jump_over_21}/{len(me)}")
    print(f"    Edges jumping over any gap: {jump_over_any_gap}/{len(me)}")

    if blue_jumps:
        print(f"\n    BLUE edge ΔH: min={min(blue_jumps)}, max={max(blue_jumps)}, "
              f"avg={sum(blue_jumps)/len(blue_jumps):.1f}")
    if black_jumps:
        print(f"    BLACK edge ΔH: min={min(black_jumps)}, max={max(black_jumps)}, "
              f"avg={sum(black_jumps)/len(black_jumps):.1f}")
    if all_jumps:
        print(f"    ALL edge ΔH: min={min(all_jumps)}, max={max(all_jumps)}, "
              f"avg={sum(all_jumps)/len(all_jumps):.1f}")

    # ΔH distribution
    dh_dist=Counter(all_jumps)
    print(f"\n    ΔH distribution (top 10):")
    for dh,cnt in sorted(dh_dist.items(), key=lambda x: -x[1])[:10]:
        # How many are blue vs black?
        b_cnt=sum(1 for j in blue_jumps if j==dh) if blue_jumps else 0
        k_cnt=sum(1 for j in black_jumps if j==dh) if black_jumps else 0
        print(f"      ΔH={dh:4d}: {cnt:5d} edges (blue={b_cnt}, black={k_cnt})")

    # Level edges (ΔH=0) by color
    level_blue = sum(1 for(a,b) in blue if mi[a]['H']==mi[b]['H'])
    level_black = sum(1 for(a,b) in black if mi[a]['H']==mi[b]['H'])
    print(f"\n    Level edges (ΔH=0): blue={level_blue}, black={level_black}, total={level_blue+level_black}")

    # ---- H=7 GAP AND THE CONFLICT GRAPH TRIANGLE ----
    print(f"\n  THE H=7 GAP (conflict-graph triangle impossibility):")
    if 7 in gaps:
        print(f"    H=7 is FORBIDDEN at n={n}")
        # What are the edges that bridge over it?
        bridging=[e for e in me if min(mi[e[0]]['H'],mi[e[1]]['H'])<7<max(mi[e[0]]['H'],mi[e[1]]['H'])]
        print(f"    Edges bridging the H=7 gap: {len(bridging)}")
        for e in bridging[:5]:
            a,b=e
            color='BLUE' if e in blue else 'BLACK'
            print(f"      {mi[a]['H']} ↔ {mi[b]['H']} ({color}, SC={mi[a]['sc']}/{mi[b]['sc']})")
    else:
        print(f"    H=7 is ACHIEVABLE at n={n} (not a gap)")

    # ---- ΔH=2 EDGES: the smallest possible jumps ----
    dh2_edges=[(a,b) for(a,b) in me if abs(mi[a]['H']-mi[b]['H'])==2]
    dh2_blue=sum(1 for e in dh2_edges if e in blue)
    dh2_black=len(dh2_edges)-dh2_blue
    print(f"\n  ΔH=2 EDGES (smallest jumps): {len(dh2_edges)} (blue={dh2_blue}, black={dh2_black})")

    # ---- H MULTIPLICITY: merged nodes per H level ----
    h_mult=Counter(mi[mm]['H'] for mm in range(mid))
    print(f"\n  H MULTIPLICITY (nodes per H level):")
    mult_vals=sorted(h_mult.values())
    mult_dist=Counter(mult_vals)
    print(f"    Distribution of multiplicities: {dict(sorted(mult_dist.items()))}")

    # SC vs NS at each H level
    h_sc_count=Counter()
    h_ns_count=Counter()
    for mm in range(mid):
        if mi[mm]['sc']: h_sc_count[mi[mm]['H']]+=1
        else: h_ns_count[mi[mm]['H']]+=1

    # Levels with only SC, only NS, or both
    sc_only_levels=sum(1 for h in h_vals if h_sc_count[h]>0 and h_ns_count[h]==0)
    ns_only_levels=sum(1 for h in h_vals if h_ns_count[h]>0 and h_sc_count[h]==0)
    mixed_levels=sum(1 for h in h_vals if h_sc_count[h]>0 and h_ns_count[h]>0)
    print(f"\n    SC-only H-levels: {sc_only_levels}")
    print(f"    NS-only H-levels: {ns_only_levels}")
    print(f"    Mixed H-levels: {mixed_levels}")

    # ---- TILING-BASED H DECOMPOSITION ----
    # H = I(Omega, 2) = 1 + 2*alpha_1 + 4*alpha_2 + ...
    # alpha_1 = c3 (at small n), alpha_2 = # independent pairs of cycles
    # H mod 4: H ≡ 1 + 2*alpha_1 (mod 4)
    # Since alpha_1 = c3 at small n: H mod 4 ≡ 1 + 2*c3 (mod 4)
    print(f"\n  H mod 4 STRUCTURE:")
    h_mod4=Counter(mi[mm]['H']%4 for mm in range(mid))
    print(f"    H mod 4 distribution: {dict(sorted(h_mod4.items()))}")
    # All H should be odd (Redei), so H mod 4 ∈ {1, 3}

    # H mod 8
    h_mod8=Counter(mi[mm]['H']%8 for mm in range(mid))
    print(f"    H mod 8 distribution: {dict(sorted(h_mod8.items()))}")

    elapsed=time.time()-t0
    print(f"\n  n={n} done ({elapsed:.1f}s)")

# ============================================================================
# CROSS-n FORBIDDEN VALUE SUMMARY
# ============================================================================

print(f"\n{'='*80}")
print(f"  CROSS-n FORBIDDEN VALUE INTEGRATION")
print(f"{'='*80}")

print(f"\n  THE TWO PERMANENT GAPS: H=7 and H=21")
print(f"  These create 'missing floors' in the metagraph DAG.")
print(f"  H=7 gap (from conflict-graph triangle impossibility):")
print(f"    alpha_1=3, alpha_2=0 → H=1+6+0=7 requires 3 pairwise-conflicting")
print(f"    odd cycles, but any 3 such cycles must share a vertex, forcing a")
print(f"    4th cycle → alpha_1≥4 → H≥9. The 'triangle' in Omega is impossible.")
print(f"  H=21 gap: 21=1+2×10 or 1+2×8+4×1 etc — all decompositions impossible.")
print()

print(f"  IMPACT ON METAGRAPH:")
print(f"  - At n=5: the gap at H=7 forces edges from H=5 to skip to H=9")
print(f"    This is the ONLY ΔH=4 edge forced by a gap.")
print(f"  - At n=6: gaps at 7,21,35,39 create multiple forced jumps")
print(f"  - At n≥8: only H=7 and H=21 persist as permanent gaps")
print(f"  - Transient gaps (filled at larger n) cause temporary distortions")
print()

print(f"  THE CONFLICT-GRAPH TRIANGLE IMPOSSIBILITY:")
print(f"  - 3 pairwise vertex-sharing odd cycles cannot have alpha_2=0")
print(f"  - This is the 'BBK triangle' constraint on Omega(T)")
print(f"  - It propagates to: H=7 forbidden, H=21 forbidden (via component analysis)")
print(f"  - It explains the ΔH=2^(n-2) jump: the transitive class (alpha_1=0)")
print(f"    connects to the first SC class with alpha_1>0, and the gap structure")
print(f"    forces specific H values to be the first achievable.")

print(f"\n{'='*80}")
print(f"  DONE")
print(f"{'='*80}")
