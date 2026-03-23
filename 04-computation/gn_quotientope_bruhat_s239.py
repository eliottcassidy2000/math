#!/usr/bin/env python3
"""
gn_quotientope_bruhat_s239.py — Path 2 (Quotientope) + Path 4 (Bruhat) synthesis
opus-2026-03-23-S239

PATH 2 (QUOTIENTOPE): If the H-gradient poset on G_n/Z_2 is a lattice quotient
of the weak Bruhat order on S_n, then G_n is a quotientope skeleton.
We test: is the H-DAG a LATTICE? Does it have joins and meets?

PATH 4 (BRUHAT ORDER): The Bruhat order on tournaments uses score-increasing
arc reversals. The H-DAG on G_n is related but different (H can increase OR
stay level). We test: how does the H-DAG relate to the score lattice?

KEY QUESTIONS:
1. Is the H-DAG on merged nodes a LATTICE? (unique joins/meets)
2. Does the H-DAG have a RANK FUNCTION? (H itself? or something refined?)
3. What is the MÖBIUS FUNCTION of the H-DAG?
4. Is the H-DAG DISTRIBUTIVE? MODULAR?
5. What is the ZETA POLYNOMIAL? (counts multichains)
"""

import sys
import subprocess
from math import comb, factorial
from itertools import permutations
from collections import defaultdict, Counter, deque
import time

sys.stdout.reconfigure(line_buffering=True)

print("=" * 80)
print("  QUOTIENTOPE + BRUHAT: IS THE H-DAG A LATTICE?")
print("  opus-2026-03-23-S239")
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
        mi[mm]={'H':cl0['H'],'sc':cl0['sc'],'score':cl0['score']}
    return {'mid':mid,'me':me,'ma':ma,'mi':mi,'n':n,'m':m}

# ============================================================================
# PATH 2: LATTICE TESTS ON THE H-DAG
# ============================================================================

for n in range(3, 8):
    print(f"\n{'='*80}")
    print(f"  n = {n}")
    print(f"{'='*80}")
    t0=time.time()

    D=build(n)
    mid=D['mid']; mi=D['mi']; me=D['me']; ma=D['ma']

    # Build the H-DAG: directed edges from lower H to higher H
    # (ignoring level edges)
    dag_adj=defaultdict(set)  # successors (higher H)
    dag_pred=defaultdict(set)  # predecessors (lower H)
    for(a,b) in me:
        ha,hb=mi[a]['H'],mi[b]['H']
        if ha<hb:
            dag_adj[a].add(b); dag_pred[b].add(a)
        elif hb<ha:
            dag_adj[b].add(a); dag_pred[a].add(b)
        # level edges: skip

    # Find min (transitive, H=1) and max (H_max)
    h_min_node=min(range(mid), key=lambda x: mi[x]['H'])
    h_max_nodes=[mm for mm in range(mid) if mi[mm]['H']==max(mi[x]['H'] for x in range(mid))]

    print(f"  V={mid}, E={len(me)}")
    print(f"  H range: [{mi[h_min_node]['H']}, {max(mi[x]['H'] for x in range(mid))}]")

    # TEST 1: Is the DAG a LATTICE?
    # A lattice requires: every pair has a unique join (LUB) and meet (GLB)
    # In the DAG: join(a,b) = lowest common successor
    #             meet(a,b) = highest common predecessor

    # Compute reachability via transitive closure
    # reach_up[v] = set of all nodes reachable by going UP from v
    reach_up={}
    for v in range(mid):
        visited=set(); q=deque([v]); visited.add(v)
        while q:
            u=q.popleft()
            for w in dag_adj[u]:
                if w not in visited: visited.add(w); q.append(w)
        reach_up[v]=visited

    # reach_down[v] = set of all nodes reachable by going DOWN
    reach_down={}
    for v in range(mid):
        visited=set(); q=deque([v]); visited.add(v)
        while q:
            u=q.popleft()
            for w in dag_pred[u]:
                if w not in visited: visited.add(w); q.append(w)
        reach_down[v]=visited

    # Check: does every pair (a,b) have a unique join?
    is_lattice=True
    no_join_count=0; multi_join_count=0
    no_meet_count=0; multi_meet_count=0

    if mid <= 100:  # Only check small graphs
        for a in range(mid):
            for b in range(a+1, mid):
                # Join: smallest element in reach_up[a] ∩ reach_up[b]
                common_up=reach_up[a] & reach_up[b]
                if not common_up:
                    no_join_count+=1; is_lattice=False
                    continue
                # Find minimal elements of common_up
                minimal_up=[x for x in common_up
                           if not any(y in common_up and y!=x for y in reach_down[x]-{x})]
                if len(minimal_up)>1:
                    multi_join_count+=1; is_lattice=False

                # Meet: largest element in reach_down[a] ∩ reach_down[b]
                common_down=reach_down[a] & reach_down[b]
                if not common_down:
                    no_meet_count+=1; is_lattice=False
                    continue
                maximal_down=[x for x in common_down
                             if not any(y in common_down and y!=x for y in reach_up[x]-{x})]
                if len(maximal_down)>1:
                    multi_meet_count+=1; is_lattice=False

    total_pairs=mid*(mid-1)//2
    print(f"\n  LATTICE TEST:")
    print(f"    Total pairs: {total_pairs}")
    print(f"    No join: {no_join_count}, Multi join: {multi_join_count}")
    print(f"    No meet: {no_meet_count}, Multi meet: {multi_meet_count}")
    print(f"    IS LATTICE: {'YES ✓' if is_lattice else 'NO ✗'}")

    if not is_lattice and mid <= 50:
        # Is it at least a SEMILATTICE? (join exists for all, or meet exists for all)
        is_join_sl=(no_join_count==0 and multi_join_count==0)
        is_meet_sl=(no_meet_count==0 and multi_meet_count==0)
        print(f"    Join semilattice: {'YES' if is_join_sl else 'NO'}")
        print(f"    Meet semilattice: {'YES' if is_meet_sl else 'NO'}")

    # TEST 2: RANK FUNCTION
    # Is H a valid rank function? rank(a) < rank(b) whenever a covers b
    # Cover relation: a→b is a cover iff there's no c with a→c→b
    covers=defaultdict(set)
    for a in range(mid):
        for b in dag_adj[a]:
            # Is a→b a cover? No c with a < c < b
            is_cover=True
            for c in dag_adj[a]:
                if c!=b and b in reach_up.get(c,set()):
                    is_cover=False; break
            if is_cover:
                covers[a].add(b)

    # Count covers
    n_covers=sum(len(v) for v in covers.values())
    print(f"\n  COVER RELATION:")
    print(f"    Number of covers: {n_covers}")
    print(f"    Number of edges: {len(me)}")
    print(f"    Covers/edges ratio: {n_covers/len(me):.4f}" if len(me)>0 else "")

    # Check if H increases by CONSTANT amount along covers
    cover_dh=[]
    for a in range(mid):
        for b in covers[a]:
            cover_dh.append(mi[b]['H']-mi[a]['H'])
    if cover_dh:
        dh_dist=Counter(cover_dh)
        print(f"    Cover ΔH distribution: {dict(sorted(dh_dist.items())[:10])}")

    # TEST 3: CHAINS (maximal paths from min to max)
    # Count maximal chains in the H-DAG
    if mid <= 50:
        chain_count=[0]
        def count_chains(v, target_set):
            if v in target_set:
                chain_count[0]+=1; return
            for w in covers[v]:
                count_chains(w, target_set)

        target=set(h_max_nodes)
        count_chains(h_min_node, target)
        print(f"\n  CHAIN COUNT (min→max): {chain_count[0]}")

    # TEST 4: WIDTH via Dilworth (max antichain)
    # An antichain = set of pairwise incomparable elements
    # Already computed as "width" in previous sessions
    h_levels=defaultdict(list)
    for mm in range(mid): h_levels[mi[mm]['H']].append(mm)
    width=max(len(v) for v in h_levels.values())
    print(f"\n  WIDTH (max antichain ≥ max H-level): {width}")

    # TEST 5: MÖBIUS FUNCTION μ(min, max)
    # For a poset P: μ(x,x)=1, Σ_{x≤z≤y} μ(x,z) = 0 for x<y
    if mid <= 100:
        # Compute μ via topological sorting
        topo=[]
        in_deg=Counter()
        for a in range(mid):
            for b in dag_adj[a]: in_deg[b]+=1
        q=deque([v for v in range(mid) if in_deg[v]==0])
        while q:
            u=q.popleft(); topo.append(u)
            for v in dag_adj[u]:
                in_deg[v]-=1
                if in_deg[v]==0: q.append(v)

        # Möbius function from h_min_node
        mu={h_min_node: 1}
        for v in topo:
            if v==h_min_node: continue
            if h_min_node in reach_down[v]:
                # μ(min, v) = -Σ_{min≤w<v} μ(min, w)
                s=0
                for w in reach_down[v]:
                    if w!=v and w in mu:
                        s+=mu[w]
                mu[v]=-s

        print(f"\n  MÖBIUS FUNCTION μ(min, ·):")
        mu_vals=sorted(mu.items(), key=lambda x: mi[x[0]]['H'])
        if len(mu_vals)<=20:
            for v,m in mu_vals:
                print(f"    H={mi[v]['H']:4d}: μ={m}")
        # μ at max
        for mx in h_max_nodes:
            if mx in mu:
                print(f"    μ(min, max) = μ(H=1, H={mi[mx]['H']}) = {mu[mx]}")

    elapsed=time.time()-t0
    print(f"\n  ({elapsed:.1f}s)")

# ============================================================================
# SYNTHESIS
# ============================================================================

print(f"\n{'='*80}")
print(f"  SYNTHESIS: QUOTIENTOPE + BRUHAT")
print(f"{'='*80}")

print(f"""
  IF the H-DAG were a LATTICE: → quotientope theory applies
    → edge count from lattice structure
    → Möbius function gives Euler characteristic
    → Zeta polynomial counts multichains
    → Connects to Pilaud-Santos framework

  IF NOT a lattice: → need weaker structure
    → still a GRADED POSET (H is the rank)
    → Möbius function still defined on intervals
    → but no quotientope skeleton

  THE BRUHAT CONNECTION:
    The weak Bruhat order on S_n has:
    - V = n! permutations
    - Edges: adjacent transpositions
    - Lattice: YES (it's a lattice)
    The score lattice on tournaments:
    - V = tournament score sequences (A000571)
    - Partial order: majorization
    - Lattice: YES (finite distributive)
    G_n/Z_2's H-DAG is a QUOTIENT of the tournament space
    by S_n × Z/2Z, restricted to the H-gradient.
    Whether this quotient preserves lattice structure
    depends on whether H is compatible with the quotient.
""")

print(f"{'='*80}")
print(f"  DONE")
print(f"{'='*80}")
