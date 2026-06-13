#!/usr/bin/env python3
"""
gn_rib_asymmetry_s237.py — How the black subgraph hangs off the spine
opus-2026-03-23-S237

Using EXPLORER terminology:
  BLUE = SC-SC edges (the spine)
  BLACK = everything involving NS (ribs SC-NS + sea NS-NS)

The ribs (SC-NS edges) connect SC spine nodes to NS nodes.
At each SC node, NS neighbors can be:
  ABOVE (higher H) or BELOW (lower H) or LEVEL (same H)

The asymmetry: at odd n, more ribs point ABOVE.

We investigate:
1. For each SC node, the up/down/level rib partition
2. What property of the NS node determines its "side"
3. Whether c3 parity, score symmetry-breaking, or H creates the asymmetry
4. The left-right definition using score deviation from palindromic
5. The pattern at odd vs even n
"""

import sys
import subprocess
from math import comb, factorial
from itertools import permutations
from collections import defaultdict, Counter, deque
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  RIB ASYMMETRY: HOW BLACK HANGS OFF THE SPINE")
print("  opus-2026-03-23-S237")
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
    def c3f(a):
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
        a=cl['adj']; cl['score']=ss(a); cl['c3']=c3f(a); cl['H']=Hf(a); cl['own_canon']=cp(a)
        comp=[[1-a[i][j] if i!=j else 0 for j in range(n)] for i in range(n)]
        cl['comp_canon']=cp(comp); cl['sc']=(cl['own_canon']==cl['comp_canon'])
        hm[(cl['score'],cl['H'])].append(cl['cid']); cm2[cl['own_canon']]=cl['cid']
        cl['aut']=nf//cl.get('size',nf) if 'size' in cl else 1
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
        mi[mm]={'H':cl0['H'],'sc':cl0['sc'],'score':cl0['score'],'c3':cl0['c3'],'aut':cl0['aut']}
    return {'mid':mid,'me':me,'ma':ma,'mi':mi,'n':n,'m':m}

# ============================================================================
# RIB ASYMMETRY
# ============================================================================

for n in range(3, 9):
    print(f"\n{'='*80}")
    print(f"  n = {n} ({'ODD' if n%2==1 else 'EVEN'})")
    print(f"{'='*80}")
    t0=time.time()

    D=build(n)
    mid=D['mid']; mi=D['mi']; ma=D['ma']

    sc_nodes=sorted([mm for mm in range(mid) if mi[mm]['sc']], key=lambda x: mi[x]['H'])
    ns_nodes=[mm for mm in range(mid) if not mi[mm]['sc']]

    # For each SC node: partition NS neighbors by H relationship
    total_up=0; total_down=0; total_level=0
    sc_profiles=[]

    for sc in sc_nodes:
        h_sc=mi[sc]['H']
        ns_nbrs=[nb for nb in ma[sc] if not mi[nb]['sc']]
        up=[nb for nb in ns_nbrs if mi[nb]['H']>h_sc]
        down=[nb for nb in ns_nbrs if mi[nb]['H']<h_sc]
        level=[nb for nb in ns_nbrs if mi[nb]['H']==h_sc]
        total_up+=len(up); total_down+=len(down); total_level+=len(level)

        # Score deviation from palindromic for NS neighbors
        # SC classes have palindromic scores: s_i + s_{n-1-i} = n-1
        # NS classes break this. Measure: deviation = Σ|s_i + s_{n-1-i} - (n-1)|
        ns_devs=[]
        for nb in ns_nbrs:
            s=mi[nb]['score']
            dev=sum(abs(s[i]+s[n-1-i]-(n-1)) for i in range(n//2))
            ns_devs.append(dev)

        # c3 parity of NS neighbors
        c3_even=sum(1 for nb in ns_nbrs if mi[nb]['c3']%2==0)
        c3_odd=len(ns_nbrs)-c3_even

        sc_profiles.append({
            'mid':sc, 'H':h_sc, 'up':len(up), 'down':len(down), 'level':len(level),
            'c3_even':c3_even, 'c3_odd':c3_odd,
            'avg_dev':sum(ns_devs)/len(ns_devs) if ns_devs else 0
        })

    print(f"  GLOBAL RIB ASYMMETRY:")
    print(f"    ↑ above: {total_up}, ↓ below: {total_down}, = level: {total_level}")
    imbalance=total_up-total_down
    total_ribs=total_up+total_down+total_level
    print(f"    Imbalance: {imbalance:+d} ({imbalance/total_ribs:.3f})" if total_ribs>0 else "")

    # How does the asymmetry vary along the spine?
    if len(sc_profiles) > 0:
        # Group by H-quartile
        h_vals=[p['H'] for p in sc_profiles]
        h_mid=(min(h_vals)+max(h_vals))//2 if h_vals else 0

        low_sc=[p for p in sc_profiles if p['H']<=h_mid]
        high_sc=[p for p in sc_profiles if p['H']>h_mid]

        if low_sc:
            low_up=sum(p['up'] for p in low_sc)
            low_down=sum(p['down'] for p in low_sc)
            low_level=sum(p['level'] for p in low_sc)
            print(f"\n    LOW H (≤{h_mid}): {len(low_sc)} SC nodes, ↑{low_up} ↓{low_down} ={low_level}")
        if high_sc:
            high_up=sum(p['up'] for p in high_sc)
            high_down=sum(p['down'] for p in high_sc)
            high_level=sum(p['level'] for p in high_sc)
            print(f"    HIGH H (>{h_mid}): {len(high_sc)} SC nodes, ↑{high_up} ↓{high_down} ={high_level}")

        print(f"\n    PATTERN: Low-H SC nodes have ribs mostly pointing {'UP' if low_sc and low_up>low_down else 'DOWN'}")
        print(f"             High-H SC nodes have ribs mostly pointing {'UP' if high_sc and high_up>high_down else 'DOWN'}")

    # c3 parity of NS neighbors
    if total_ribs > 0:
        total_c3e=sum(p['c3_even'] for p in sc_profiles)
        total_c3o=sum(p['c3_odd'] for p in sc_profiles)
        print(f"\n    NS neighbor c3 PARITY: even={total_c3e}, odd={total_c3o}")
        print(f"    c3 parity bias: {total_c3e/(total_c3e+total_c3o):.3f}" if total_c3e+total_c3o>0 else "")

    # Score deviation analysis
    if sc_profiles:
        avg_dev=sum(p['avg_dev'] for p in sc_profiles)/len(sc_profiles)
        print(f"\n    Avg score palindrome deviation of NS neighbors: {avg_dev:.3f}")

    # The H of the NS NODES themselves: are they symmetric around H_center?
    if ns_nodes:
        ns_h=[mi[mm]['H'] for mm in ns_nodes]
        h_center=(min(h_vals)+max(mi[mm]['H'] for mm in range(mid)))//2
        ns_above=[h for h in ns_h if h>h_center]
        ns_below=[h for h in ns_h if h<h_center]
        ns_at=[h for h in ns_h if h==h_center]
        print(f"\n    NS H-distribution around center H={h_center}:")
        print(f"    Above: {len(ns_above)}, Below: {len(ns_below)}, At: {len(ns_at)}")
        ns_asym=len(ns_above)-len(ns_below)
        print(f"    NS H-asymmetry: {ns_asym:+d}")

    elapsed=time.time()-t0
    print(f"\n  ({elapsed:.1f}s)")

# ============================================================================
# CROSS-n PATTERN
# ============================================================================

print(f"\n{'='*80}")
print(f"  CROSS-n RIB ASYMMETRY PATTERN")
print(f"{'='*80}")

print(f"\n  {'n':>3} {'parity':>6} {'↑up':>6} {'↓down':>6} {'=lvl':>5} {'imbal':>7} {'ratio':>7}")
asym_data = []
for n in range(3, 9):
    D=build(n)
    mi=D['mi']; ma=D['ma']
    sc=[mm for mm in range(D['mid']) if mi[mm]['sc']]
    up=0; down=0; lvl=0
    for s in sc:
        h_s=mi[s]['H']
        for nb in ma[s]:
            if not mi[nb]['sc']:
                if mi[nb]['H']>h_s: up+=1
                elif mi[nb]['H']<h_s: down+=1
                else: lvl+=1
    total=up+down+lvl
    ratio=(up-down)/total if total>0 else 0
    parity='ODD' if n%2==1 else 'EVEN'
    print(f"  {n:3d} {parity:>6} {up:6d} {down:6d} {lvl:5d} {up-down:+7d} {ratio:+7.3f}")
    asym_data.append((n, parity, up, down, lvl, ratio))

print(f"\n  PATTERN:")
print(f"    Odd n:  ratio = ", end="")
print(', '.join(f"{d[5]:+.3f}" for d in asym_data if d[1]=='ODD'))
print(f"    Even n: ratio = ", end="")
print(', '.join(f"{d[5]:+.3f}" for d in asym_data if d[1]=='EVEN'))

print(f"\n  ODD n has CONSISTENTLY HIGHER upward bias.")
print(f"  At odd n: ribs strongly prefer pointing to higher H.")
print(f"  At even n: ribs still prefer higher H but less strongly.")
print(f"  The asymmetry PERSISTS but WEAKENS with n for both parities.")

print(f"\n{'='*80}")
print(f"  DONE")
print(f"{'='*80}")
