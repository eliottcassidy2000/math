#!/usr/bin/env python3
"""
gn_lines_parents_s225.py — Lines vs Edges + Parent Class Analysis
opus-2026-03-23-S225

KEY DISTINCTION:
- An EDGE connects two distinct merged nodes (an unordered pair {A,B})
- A LINE is a (tournament, arc-position) transition. Multiple lines can
  connect the same pair of nodes (different arc positions producing same pair)
- A SELF-LINE stays within the same merged node (self-loop)

So we count separately:
  blue_edges, black_edges (distinct pairs connected by blue/black)
  blue_lines, black_lines (total arc-flip transitions, including multiplicity)
  blue_self_lines, black_self_lines (transitions staying within a merged node)

Also: PARENT CLASS analysis — for each iso class at n, which iso class
at n-1 does it "reduce to" when we delete the highest-score vertex?
"""

import sys
import subprocess
from math import comb, factorial
from itertools import permutations
from collections import defaultdict, Counter, deque
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  LINES vs EDGES + PARENT CLASS ANALYSIS")
print("  opus-2026-03-23-S225")
print("=" * 80)

def build_full(n):
    """Build G_n/Z_2 with WEIGHTED edges (line counts)."""
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

    # Count LINES (weighted transitions) not just edges
    edge_weight=Counter()  # (min_cid, max_cid) → count of arc flips producing this
    self_weight=Counter()  # cid → count of arc flips staying in same class
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
                self_weight[cl['cid']]+=1
            elif nb>=0:
                edge_weight[(min(cl['cid'],nb),max(cl['cid'],nb))]+=1

    # Merged graph with weights
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

    mi={}
    for mm in range(mid):
        cids=[cl['cid'] for cl in cls if mn[cl['cid']]==mm]
        mi[mm]={'H':cls[cids[0]]['H'],'sc':cls[cids[0]]['sc'],
                'score':cls[cids[0]]['score'],'c3':cls[cids[0]]['c3'],'cids':cids}

    # Merged edge weights and self-loop weights
    m_edge_weight=Counter()  # (min_mid, max_mid) → total lines
    m_self_weight=Counter()  # mid → total self-lines

    for(a,b),w in edge_weight.items():
        ma,mb=mn[a],mn[b]
        if ma==mb:
            m_self_weight[ma]+=w
        else:
            m_edge_weight[(min(ma,mb),max(ma,mb))]+=w

    for cid,w in self_weight.items():
        m_self_weight[mn[cid]]+=w

    # Classify by blue/black
    blue_edges=set(); black_edges=set()
    blue_lines=0; black_lines=0
    blue_self_lines=0; black_self_lines=0

    for(a,b),w in m_edge_weight.items():
        if mi[a]['sc']==mi[b]['sc']:
            blue_edges.add((a,b)); blue_lines+=w
        else:
            black_edges.add((a,b)); black_lines+=w

    # Self-lines: classify by whether the merged node is SC
    # Self-lines are "blue" because they stay within the same type
    for mm,w in m_self_weight.items():
        blue_self_lines+=w  # Self-loops are within-type → blue

    # Total lines = all (tournament, arc) transitions = 2^m × m
    total_lines = (1<<m)*m if n<=20 else 0
    cross_lines = blue_lines + black_lines
    self_lines_total = sum(m_self_weight.values())

    return {
        'n':n,'m':m,'mid':mid,'mi':mi,'mn':mn,'cls':cls,
        'blue_edges':len(blue_edges),'black_edges':len(black_edges),
        'total_edges':len(blue_edges)+len(black_edges),
        'blue_lines':blue_lines,'black_lines':black_lines,
        'blue_self_lines':blue_self_lines,'black_self_lines':black_self_lines,
        'self_lines_total':self_lines_total,
        'total_lines':total_lines,
        'cross_lines':cross_lines,
        'm_edge_weight':m_edge_weight,
    }

# ============================================================================
# PARENT CLASS ANALYSIS
# ============================================================================

def compute_parents(n, cls, mi, mn):
    """For each class at n, compute its parent class at n-1.
    Parent = iso class of T restricted to vertices {0,...,n-2}
    (delete highest-score vertex, or equivalently, delete vertex n-1
    from the canonical representative)."""
    if n <= 3:
        return {}

    n1 = n - 1
    m1 = comb(n1, 2)
    P1 = [(i,j) for i in range(n1) for j in range(i+1,n1)]

    # Build n-1 class lookup
    if n1 <= 6:
        AP1 = list(permutations(range(n1)))
        def cn1(a):
            best=None
            for p in AP1:
                f=tuple(a[p[i]][p[j]] for i in range(n1) for j in range(n1))
                if best is None or f<best: best=f
            return best
    else:
        def cn1(a):
            scores=[sum(a[i][j] for j in range(n1)) for i in range(n1)]
            sg=defaultdict(list)
            for v in range(n1): sg[scores[v]].append(v)
            gs=[sg[s] for s in sorted(set(scores))]
            if all(len(g)==1 for g in gs):
                p=[g[0] for g in gs]
                return tuple(a[p[i]][p[j]] for i in range(n1) for j in range(n1))
            best=None
            def gen(gs2):
                if not gs2: yield []; return
                for p in permutations(gs2[0]):
                    for rest in gen(gs2[1:]): yield list(p)+rest
            for p in gen(gs):
                f=tuple(a[p[i]][p[j]] for i in range(n1) for j in range(n1))
                if best is None or f<best: best=f
            return best

    # Get parent classes at n-1
    r1=subprocess.run(['gentourng',str(n1)],capture_output=True,text=True)
    ls1=[l.strip() for l in (r1.stdout or '').split('\n') if len(l.strip())==m1 and all(c in '01' for c in l.strip())]
    parent_canons={}
    for pid,line in enumerate(ls1):
        bits=int(line,2)
        a=[[0]*n1 for _ in range(n1)]
        for k,(i,j) in enumerate(P1):
            if bits&(1<<(m1-1-k)): a[i][j]=1
            else: a[j][i]=1
        parent_canons[cn1(a)]=pid

    # For each class at n, delete each vertex, find parent
    parent_map={}  # cid → set of parent pids
    for cl in cls:
        adj=cl['adj']
        parents=set()
        for v in range(n):
            # Delete vertex v
            verts=[i for i in range(n) if i!=v]
            sub=[[adj[verts[i]][verts[j]] for j in range(n1)] for i in range(n1)]
            c=cn1(sub)
            pid=parent_canons.get(c,-1)
            if pid>=0: parents.add(pid)
        parent_map[cl['cid']]=parents

    return parent_map, len(ls1)

# ============================================================================
# MAIN
# ============================================================================

ALL={}

for n in range(3, 9):
    print(f"\n{'='*80}")
    print(f"  n = {n}")
    print(f"{'='*80}")
    t0=time.time()

    D=build_full(n)
    mid=D['mid']; mi=D['mi']; m=D['m']

    print(f"  G_{n}/Z_2: V={mid}, E_total={D['total_edges']}")
    print(f"\n  LINES vs EDGES:")
    print(f"    {'':>20s} {'BLUE':>10s} {'BLACK':>10s} {'TOTAL':>10s}")
    print(f"    {'edges (pairs)':>20s} {D['blue_edges']:10d} {D['black_edges']:10d} {D['total_edges']:10d}")
    print(f"    {'cross-lines':>20s} {D['blue_lines']:10d} {D['black_lines']:10d} {D['cross_lines']:10d}")
    print(f"    {'self-lines':>20s} {D['blue_self_lines']:10d} {D['black_self_lines']:10d} {D['self_lines_total']:10d}")
    print(f"    {'total lines':>20s} {D['blue_lines']+D['blue_self_lines']:10d} {D['black_lines']+D['black_self_lines']:10d} {D['total_lines']:10d}")

    # Average lines per edge
    if D['blue_edges']>0:
        avg_blue_weight=D['blue_lines']/(2*D['blue_edges'])  # bidirectional
    else:
        avg_blue_weight=0
    if D['black_edges']>0:
        avg_black_weight=D['black_lines']/(2*D['black_edges'])
    else:
        avg_black_weight=0
    avg_total_weight=D['cross_lines']/(2*D['total_edges']) if D['total_edges']>0 else 0

    print(f"\n    avg lines/edge: blue={avg_blue_weight:.2f}, black={avg_black_weight:.2f}, total={avg_total_weight:.2f}")
    print(f"    self_lines/total_lines: {D['self_lines_total']/D['total_lines']:.4f}" if D['total_lines']>0 else "")
    print(f"    cross_lines/total_lines: {D['cross_lines']/D['total_lines']:.4f}" if D['total_lines']>0 else "")

    # Weight distribution
    weights=list(D['m_edge_weight'].values())
    if weights:
        wt_dist=Counter(weights)
        print(f"\n    Edge weight distribution: {dict(sorted(wt_dist.items())[:10])}")

    # Parent class analysis
    if n >= 4 and n <= 8:
        print(f"\n  PARENT CLASS ANALYSIS (n={n} → n-1={n-1}):")
        pm, n_parents = compute_parents(n, D['cls'], mi, D['mn'])

        # For each merged node, collect parent set
        merged_parents={}
        for mm in range(mid):
            cids=mi[mm]['cids']
            all_parents=set()
            for ci in cids:
                if ci in pm:
                    all_parents|=pm[ci]
            merged_parents[mm]=all_parents

        # Statistics
        parent_counts=[len(merged_parents[mm]) for mm in range(mid)]
        print(f"    Parent classes at n-1: {n_parents}")
        print(f"    Parents per merged node: min={min(parent_counts)}, max={max(parent_counts)}, "
              f"avg={sum(parent_counts)/len(parent_counts):.2f}")

        # Parent distribution
        p_dist=Counter(parent_counts)
        print(f"    Distribution: {dict(sorted(p_dist.items()))}")

        # How many children does each parent have?
        child_count=Counter()
        for mm in range(mid):
            for p in merged_parents[mm]:
                child_count[p]+=1
        child_vals=list(child_count.values())
        print(f"    Children per parent: min={min(child_vals)}, max={max(child_vals)}, "
              f"avg={sum(child_vals)/len(child_vals):.2f}")

    ALL[n]=D
    print(f"\n  n={n} done ({time.time()-t0:.1f}s)")

# ============================================================================
# SEQUENCES
# ============================================================================

print(f"\n{'='*80}")
print(f"  ALL SEQUENCES")
print(f"{'='*80}")

seq_names=[
    ('blue_edges','Blue edges (merged)'),
    ('black_edges','Black edges (merged)'),
    ('total_edges','Total edges (merged)'),
    ('blue_lines','Blue cross-lines'),
    ('black_lines','Black cross-lines'),
    ('cross_lines','Total cross-lines'),
    ('blue_self_lines','Blue self-lines'),
    ('self_lines_total','Total self-lines'),
    ('total_lines','Total lines (2^m × m)'),
]

for key,label in seq_names:
    vals=[ALL[n][key] for n in range(3,9)]
    print(f"  {label:>25s}: {vals}")

# Derived
print(f"\n  DERIVED RATIOS:")
for n in range(3,9):
    D=ALL[n]
    bl=D['blue_lines']; kl=D['black_lines']; cl=D['cross_lines']
    be=D['blue_edges']; ke=D['black_edges']; te=D['total_edges']
    print(f"  n={n}: blue_lines/blue_edges={bl/(2*be):.2f}" if be>0 else f"  n={n}: no blue edges",
          f"black_lines/black_edges={kl/(2*ke):.2f}" if ke>0 else "no black edges",
          f"self/total={D['self_lines_total']/D['total_lines']:.4f}" if D['total_lines']>0 else "")

print(f"\n{'='*80}")
print(f"  DONE")
print(f"{'='*80}")
