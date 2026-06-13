"""
tiling_graph_topology.py — Graph topology of tiling connections under various pairings.

The user identified: blue flip-pair graph at n=5 = C4 + 4 pendants (square ring).
This script finds similar constrained structures for ALL natural pairings:
  A. Blue flip-pairs  (tile complement, both grid-sym)
  B. Black flip-pairs (tile complement, both non-grid-sym)
  C. Zeckendorf-cube adjacency (flip one bit, both Zeckendorf)
  D. Transpose-pairs (apply TRANS_MAP, different classes)
  E. Grid-sym Zeckendorf intersections

KEY THEOREMS TO PROVE:
  - Blue count per SC class: ALWAYS ODD (proved: H = |Aut|×#t, all odd)
  - Black count per merged node: ALWAYS EVEN
  - Zeckendorf count per class: varies, but sum = F_{m+2}

oracle-2026-05-10
"""

import sys
from collections import defaultdict, Counter
from itertools import permutations
sys.stdout.reconfigure(encoding='utf-8')

def make_tiles(n):
    return [(x,y) for y in range(1,n-1) for x in range(n,y+1,-1)]
def make_trans_map(n,tiles):
    idx={t:i for i,t in enumerate(tiles)}
    return [idx[(n-y+1,n-x+1)] for x,y in tiles]
def bits_to_adj(bits,n,tiles):
    verts=list(range(n,0,-1))
    A=[[0]*n for _ in range(n)]
    for k in range(n-1): A[k][k+1]=1
    for i,(xl,yl) in enumerate(tiles):
        xi=verts.index(xl); yi=verts.index(yl)
        if bits[i]==0: A[xi][yi]=1
        else: A[yi][xi]=1
    return A
def adj_str(A): return ''.join(''.join(map(str,r)) for r in A)
def canonicalize(A,n,perms):
    best=None
    for p in perms:
        s=adj_str([[A[p[i]][p[j]] for j in range(n)] for i in range(n)])
        if best is None or s<best: best=s
    return best
def full_comp(A,n): return [[1-A[i][j] if i!=j else 0 for j in range(n)] for i in range(n)]
def is_grid_sym(bits,tm): return all(bits[i]==bits[tm[i]] for i in range(len(bits)) if tm[i]!=i)
def is_zeck(bits): return not any(bits[k]==1 and bits[k+1]==1 for k in range(len(bits)-1))
def fib(k):
    if k<=2: return 1
    a,b=1,1
    for _ in range(k-2): a,b=b,a+b
    return b
def H_count(A,n):
    dp={}
    for v in range(n): dp[(1<<v,v)]=1
    for ms in range(2,n+1):
        for mask in range(1<<n):
            if bin(mask).count('1')!=ms: continue
            for v in range(n):
                if not(mask>>v&1): continue
                pm=mask^(1<<v)
                t=sum(dp.get((pm,u),0) for u in range(n) if(pm>>u&1) and A[u][v])
                if t: dp[(mask,v)]=t
    return sum(dp.get(((1<<n)-1,v),0) for v in range(n))

def adj_str_multiline(adj, labels):
    """Print adjacency matrix with labels."""
    n=len(labels)
    rows=[f"  {str(labels[i]):>6}: " + " ".join(f"{adj[i][j]:>3}" for j in range(n)) for i in range(n)]
    return "\n".join(rows)

def describe_graph(adj, labels, title="Graph"):
    """Describe topology of a multigraph with adjacency matrix."""
    n = len(labels)
    degrees = [sum(adj[i][j] for j in range(n)) + 2*adj[i][i] for i in range(n)]
    # Wait: standard convention. Let me just count edges.
    edges = []
    for i in range(n):
        if adj[i][i] > 0:
            for _ in range(adj[i][i]): edges.append((i,i))
        for j in range(i+1,n):
            for _ in range(adj[i][j]): edges.append((i,j))

    deg = [0]*n
    for i,j in edges:
        deg[i] += 1
        if i!=j: deg[j] += 1

    print(f"  {title}:")
    print(f"  Nodes: {n}, Edges: {len(edges)}")
    print(f"  Degree sequence: {sorted(deg)}")
    print(f"  Node degrees: {dict(zip(labels, deg))}")

    # Self-loops
    self_loops = [(i,i) for i,j in edges if i==j]
    if self_loops:
        print(f"  Self-loops at: {[labels[i] for i,_ in self_loops]}")

    # Check connectivity (ignoring self-loops)
    reachable = set()
    if n > 0:
        stack = [0]
        reachable.add(0)
        while stack:
            v = stack.pop()
            for u in range(n):
                if u not in reachable and (adj[v][u]>0 or adj[u][v]>0):
                    reachable.add(u)
                    stack.append(u)
    print(f"  Connected: {len(reachable)==n}")

    # Find bridges/cut vertices would be complex; just note structure
    max_mult = max((adj[i][j] for i in range(n) for j in range(i+1,n)), default=0)
    if max_mult > 1:
        print(f"  Multi-edges (max multiplicity: {max_mult})")

    return edges, deg

# ══════════════════════════════════════════════════════════════════════════
# MAIN ANALYSIS
# ══════════════════════════════════════════════════════════════════════════

for n in range(3, 7):
    print("=" * 70)
    print(f"n={n}")
    print("=" * 70)

    tiles = make_tiles(n)
    m = len(tiles)
    tm = make_trans_map(n, tiles)
    perms = list(permutations(range(n)))

    # Build all tilings
    groups = defaultdict(list)
    comp_groups = {}
    for mask in range(1<<m):
        bits = [(mask>>k)&1 for k in range(m)]
        A = bits_to_adj(bits, n, tiles)
        canon = canonicalize(A, n, perms)
        comp_canon = canonicalize(full_comp(A, n), n, perms)
        gsym = is_grid_sym(bits, tm)
        flip_mask = sum((1-b)<<k for k,b in enumerate(bits))
        zk = is_zeck(bits)
        # trans bits: apply TRANS_MAP
        tb = [0]*m
        for i in range(m): tb[tm[i]] = bits[i]
        trans_mask = sum(b<<k for k,b in enumerate(tb))
        groups[canon].append({'mask':mask,'bits':bits,'A':A,'cc':comp_canon,
                               'gs':gsym,'fm':flip_mask,'zk':zk,'tm':trans_mask})
        comp_groups[canon] = comp_canon

    classes = []
    sig_to_ci = {}
    for ci, (sig, members) in enumerate(sorted(groups.items())):
        A = members[0]['A']
        H = H_count(A, n)
        au = sum(1 for p in perms if all(A[p[i]][p[j]]==A[i][j] for i in range(n) for j in range(n)))
        cc = members[0]['cc']
        sc = (cc == sig)
        Z = sum(1 for t in members if t['zk'])
        gs = sum(1 for t in members if t['gs'])
        classes.append({'ci':ci,'sig':sig,'m':members,'nt':len(members),'H':H,'au':au,'sc':sc,'cc':cc,'Z':Z,'gs':gs})
        sig_to_ci[sig] = ci

    # Update classIndex in members
    mask_to_ci = {}
    for c in classes:
        for t in c['m']:
            t['ci'] = c['ci']
            mask_to_ci[t['mask']] = c['ci']

    # Build merged nodes
    sc_cis = [c['ci'] for c in classes if c['sc']]
    ns_cis = [c['ci'] for c in classes if not c['sc']]
    seen_ns = set()
    ns_pairs = []
    for c in classes:
        if c['sc'] or c['ci'] in seen_ns: continue
        partner_ci = sig_to_ci[c['cc']]
        seen_ns.add(c['ci']); seen_ns.add(partner_ci)
        ns_pairs.append((c['ci'], partner_ci))

    # merged_node_id: SC class -> ci, NS pair -> min(ci_a, ci_b)
    ci_to_merged = {}
    merged_labels = []  # label for each merged node
    for ci in sc_cis:
        ci_to_merged[ci] = ci
        merged_labels.append(f"SC#{ci}")
    for a, b in ns_pairs:
        ci_to_merged[a] = min(a,b)
        ci_to_merged[b] = min(a,b)
        merged_labels.append(f"NS[{a},{b}]")
    merged_ids = sorted(set(ci_to_merged.values()))
    merged_idx = {v: i for i, v in enumerate(merged_ids)}

    print(f"\n  {len(classes)} iso classes, {len(sc_cis)} SC, {len(ns_pairs)} NS-pairs")
    print(f"  m={m}, 2^m={1<<m}, F_{{{m+2}}}={fib(m+2)}")

    # ── A. BLUE FLIP-PAIR GRAPH ──────────────────────────────────────────
    print(f"\n{'─'*60}")
    print("  A. BLUE FLIP-PAIR GRAPH")
    print(f"{'─'*60}")
    print(f"  Blue tilings = grid-symmetric (isGridSym), paired by tile complement")
    print(f"  Total grid-sym: {sum(c['gs'] for c in classes)} = 2^{sum(1 for i in range(m) if tm[i]==i)+(m-sum(1 for i in range(m) if tm[i]==i))//2}")
    print(f"  Blue flip-pairs: {sum(c['gs'] for c in classes)//2}")

    # All SC classes have blue count; verify parity
    blue_per_sc = {c['ci']: c['gs'] for c in classes if c['sc']}
    assert all(v%2==1 for v in blue_per_sc.values() if v>0), "Blue count should be odd for SC!"
    print(f"  Blue counts (SC classes, sorted by ci): {[blue_per_sc[ci] for ci in sorted(blue_per_sc)]}")
    print(f"  All non-zero blue counts odd: {all(v%2==1 for v in blue_per_sc.values() if v>0)} ✓")

    # Build blue adjacency matrix on SC classes
    sc_sorted = sorted(sc_cis)
    sc_idx = {ci: i for i, ci in enumerate(sc_sorted)}
    nsc = len(sc_sorted)
    blue_adj = [[0]*nsc for _ in range(nsc)]

    seen_bp = set()
    for c in classes:
        for t in c['m']:
            if not t['gs']: continue
            fm = t['fm']
            pk = (min(t['mask'], fm), max(t['mask'], fm))
            if pk in seen_bp: continue
            seen_bp.add(pk)
            ci_a = t['ci']
            ci_b = mask_to_ci[fm]
            if ci_a in sc_idx and ci_b in sc_idx:
                ia, ib = sc_idx[ci_a], sc_idx[ci_b]
                blue_adj[ia][ib] += 1
                if ia != ib: blue_adj[ib][ia] += 1

    # Print adjacency
    sc_labels = [f"ci{ci}(gs={blue_per_sc.get(ci,0)})" for ci in sc_sorted]
    print(f"\n  Blue adjacency (SC classes):")
    if nsc <= 12:
        for i in range(nsc):
            row = [blue_adj[i][j] for j in range(nsc)]
            print(f"    {sc_labels[i]:>20}: {row}")

    # Graph properties
    blue_edges, blue_deg = describe_graph(blue_adj, sc_labels, "Blue flip-pair multigraph")

    # Classify structure
    non_zero_nodes = [sc_sorted[i] for i in range(nsc) if blue_per_sc.get(sc_sorted[i],0)>0]
    if nsc <= 12:
        # Find which pairs are adjacent
        print(f"\n  Blue connections (class_ci -> class_ci: count):")
        for ia in range(nsc):
            for ib in range(ia, nsc):
                if blue_adj[ia][ib] > 0:
                    print(f"    ci{sc_sorted[ia]}(H={classes[sc_sorted[ia]]['H']}) -- ci{sc_sorted[ib]}(H={classes[sc_sorted[ib]]['H']}): {blue_adj[ia][ib]} pair(s)")

    # ── B. BLACK FLIP-PAIR GRAPH ─────────────────────────────────────────
    print(f"\n{'─'*60}")
    print("  B. BLACK FLIP-PAIR GRAPH")
    print(f"{'─'*60}")

    total_black = (1<<m) - sum(c['gs'] for c in classes)
    print(f"  Black tilings (non-grid-sym): {total_black}")
    print(f"  Black flip-pairs: {total_black//2}")

    # Black count per merged node
    black_per_merged = defaultdict(int)
    for c in classes:
        mn = ci_to_merged[c['ci']]
        black_per_merged[mn] += c['nt'] - c['gs']

    mn_list = sorted(merged_ids)
    mn_idx = {v: i for i, v in enumerate(mn_list)}
    nmn = len(mn_list)

    print(f"  Black counts per merged node: {[black_per_merged[mn] for mn in mn_list]}")
    print(f"  All even: {all(v%2==0 for v in black_per_merged.values())} ✓")

    # Build black adjacency matrix on merged nodes
    black_adj = [[0]*nmn for _ in range(nmn)]

    seen_blk = set()
    for c in classes:
        for t in c['m']:
            if t['gs']: continue  # skip grid-sym (blue), keep non-grid-sym (black)
            fm = t['fm']
            pk = (min(t['mask'], fm), max(t['mask'], fm))
            if pk in seen_blk: continue
            seen_blk.add(pk)
            mn_a = ci_to_merged[t['ci']]
            mn_b = ci_to_merged[mask_to_ci[fm]]
            ia, ib = mn_idx[mn_a], mn_idx[mn_b]
            black_adj[ia][ib] += 1
            if ia != ib: black_adj[ib][ia] += 1

    mn_labels = []
    for mn in mn_list:
        ci = mn
        if ci in sc_cis:
            mn_labels.append(f"SC#{ci}(bk={black_per_merged[mn]})")
        else:
            a, b = next((a,b) for a,b in ns_pairs if min(a,b)==ci)
            mn_labels.append(f"NS[{a},{b}](bk={black_per_merged[mn]})")

    print(f"\n  Black adjacency (merged nodes with non-zero black):")
    non_zero_mn = [i for i in range(nmn) if black_per_merged[mn_list[i]] > 0]
    if len(non_zero_mn) <= 12:
        sub_labels = [mn_labels[i] for i in non_zero_mn]
        sub_adj = [[black_adj[i][j] for j in non_zero_mn] for i in non_zero_mn]
        for i, row_i in enumerate(non_zero_mn):
            print(f"    {mn_labels[row_i]:>35}: {[black_adj[row_i][j] for j in non_zero_mn]}")

        black_edges, black_deg = describe_graph(sub_adj, sub_labels, "Black flip-pair multigraph")

    # ── C. ZECKENDORF-CUBE ADJACENCY ──────────────────────────────────────
    print(f"\n{'─'*60}")
    print("  C. ZECKENDORF-CUBE ADJACENCY")
    print(f"{'─'*60}")
    print(f"  Fibonacci cube Γ_{m}: edges = pairs of Zeckendorf tilings differing by 1 bit")

    zeck_tilings = [(t, c['ci']) for c in classes for t in c['m'] if t['zk']]
    mask_to_zeck_ci = {t['mask']: c['ci'] for c in classes for t in c['m'] if t['zk']}

    print(f"  Zeckendorf tilings: {len(zeck_tilings)} = F_{{{m+2}}} = {fib(m+2)}")

    # Build Zeckendorf-cube adjacency per class
    z_per_class = defaultdict(int)
    for t, ci in zeck_tilings: z_per_class[ci] += 1

    nz_cis = sorted(ci for ci, cnt in z_per_class.items() if cnt > 0)
    nzc = len(nz_cis)
    nz_idx = {ci: i for i, ci in enumerate(nz_cis)}

    print(f"  Classes with Z-tilings: {nzc}")
    print(f"  Z-counts: {[z_per_class[ci] for ci in nz_cis]}")

    # Fibonacci cube edges
    zeck_adj = [[0]*nzc for _ in range(nzc)]
    zeck_edges_total = 0
    seen_ze = set()
    for t, ci_a in zeck_tilings:
        bits = t['bits']
        # Each neighbor: flip exactly one bit
        for k in range(m):
            new_bits = bits.copy()
            new_bits[k] = 1 - bits[k]
            if is_zeck(new_bits):
                new_mask = sum(b<<j for j,b in enumerate(new_bits))
                pk = (min(t['mask'], new_mask), max(t['mask'], new_mask))
                if pk in seen_ze: continue
                seen_ze.add(pk)
                if new_mask in mask_to_zeck_ci:
                    ci_b = mask_to_zeck_ci[new_mask]
                    if ci_a in nz_idx and ci_b in nz_idx:
                        ia, ib = nz_idx[ci_a], nz_idx[ci_b]
                        zeck_adj[ia][ib] += 1
                        if ia != ib: zeck_adj[ib][ia] += 1
                        zeck_edges_total += 1

    print(f"  Total Fibonacci-cube edges: {zeck_edges_total}")

    # Z-count parities
    z_parities = {ci: z_per_class[ci]%2 for ci in nz_cis}
    print(f"  Z-count parities: {[z_parities[ci] for ci in nz_cis]} (0=even, 1=odd)")
    print(f"  Odd Z-count classes: {sum(1 for v in z_parities.values() if v==1)}")
    print(f"  Even Z-count classes: {sum(1 for v in z_parities.values() if v==0)}")

    # Print Z-cube adjacency
    if nzc <= 12:
        nz_labels = [f"ci{ci}(Z={z_per_class[ci]},{'SC' if classes[ci]['sc'] else 'NS'})" for ci in nz_cis]
        print(f"\n  Zeckendorf-cube adjacency:")
        for i in range(nzc):
            print(f"    {nz_labels[i]:>35}: {[zeck_adj[i][j] for j in range(nzc)]}")

        zeck_deg_list = [sum(zeck_adj[i][j] for j in range(nzc)) + (zeck_adj[i][i] if True else 0) for i in range(nzc)]
        # Actually: degree in Fibonacci cube = number of 1-bits that can be flipped (removing) +
        # 0-bits that can be set to 1 without creating consecutive 1s
        print(f"  Degree sequence in Z-cube class graph: {sorted(sum(zeck_adj[i][j] for j in range(nzc)) for i in range(nzc))}")

    # ── D. GRID-SYM ZECKENDORF INTERSECTION ──────────────────────────────
    print(f"\n{'─'*60}")
    print("  D. GRID-SYM ∩ ZECKENDORF (blue Zeckendorf tilings)")
    print(f"{'─'*60}")

    gz_per_class = defaultdict(int)
    gz_tilings = [(t, c['ci']) for c in classes for t in c['m'] if t['gs'] and t['zk']]
    for t, ci in gz_tilings: gz_per_class[ci] += 1

    total_gz = sum(gz_per_class.values())
    print(f"  Grid-sym ∩ Zeckendorf tilings: {total_gz}")
    print(f"  Per SC class: {[gz_per_class.get(ci,0) for ci in sc_sorted]}")
    print(f"  Sum: {sum(gz_per_class.get(ci,0) for ci in sc_sorted)}")

    # Parity of gz counts
    gz_counts = [gz_per_class.get(ci,0) for ci in sc_sorted if gz_per_class.get(ci,0)>0]
    print(f"  Non-zero gz counts (all odd?): {gz_counts} → all odd: {all(v%2==1 for v in gz_counts)}")

    # Flip-pairs within gz tilings
    gz_adj = [[0]*nsc for _ in range(nsc)]
    seen_gz = set()
    for t, ci_a in gz_tilings:
        fm = t['fm']
        pk = (min(t['mask'], fm), max(t['mask'], fm))
        if pk in seen_gz: continue
        seen_gz.add(pk)
        ci_b = mask_to_ci[fm]
        # Is T_flip also gz?
        t_flip = next((x for x in classes[ci_b]['m'] if x['mask']==fm), None)
        if t_flip and t_flip['gs'] and t_flip['zk']:
            if ci_a in sc_idx and ci_b in sc_idx:
                ia, ib = sc_idx[ci_a], sc_idx[ci_b]
                gz_adj[ia][ib] += 1
                if ia != ib: gz_adj[ib][ia] += 1

    gz_pairs = sum(gz_adj[i][j] for i in range(nsc) for j in range(i,nsc))
    print(f"  Blue-Zeckendorf flip-pairs (both gz): {gz_pairs}")

    print()

# ══════════════════════════════════════════════════════════════════════════
# PATTERN SUMMARY
# ══════════════════════════════════════════════════════════════════════════
print("=" * 70)
print("PATTERN SUMMARY")
print("=" * 70)

print("""
UNIVERSAL CONSTRAINTS (proved theorems):

1. Blue tiling count per SC class: ALWAYS ODD
   Proof: H(C) = |Aut(C)| × #tilings(C), H always odd → #tilings odd.
   Grid-sym tilings = fixed points of TRANS_MAP on SC class.
   #grid-sym = #tilings - 2×(internal trans-pairs) = odd - even = ODD.
   [Since internal trans-pairs = (#tilings - #grid-sym)/2 must be integer,
    and #tilings is odd, #grid-sym must be odd too.]

2. Black tiling count per MERGED node: ALWAYS EVEN
   SC class: #black = #tilings - #grid-sym = odd - odd = EVEN.
   NS merged pair: #black = #t_A + #t_B (both black, no grid-sym)
                           = odd + odd = EVEN.

3. Sum of all blue tilings = 2^{free_bits} = 2^{SC + non-self-pairs}
   where free_bits = #self-paired-tiles + #non-self-pairs.

4. Sum of all black tilings = 2^m - 2^{free_bits}

TOPOLOGY EVOLUTION OF BLUE FLIP-PAIR GRAPH:

n=3: m=1, 2 blue tilings, 1 pair. 2 SC classes with 1 blue each.
     Graph: P_2 (single edge between the two SC classes).

n=4: m=3, 4 blue tilings, 2 pairs. 2 SC classes with blue={1,3}.
     Graph: 1 cross-edge + 1 self-loop (tadpole at degree-3 class).

n=5: m=6, 16 blue tilings, 8 pairs. 8 SC classes with blue={1,1,1,1,3,3,3,3}.
     Graph: C4 (square ring on 4-degree-3 nodes) + 4 pendant nodes (1-degree).
     → The square ring connects the four 3-classes cyclically.
     → Each 3-class has 1 pendant to a 1-class.

n=6: m=10, 64 blue tilings, 32 pairs. 12 SC classes with blue sorted={1,1,3,3,5,5,7,7,7,7,9,9}.

ZECKENDORF OBSERVATION:
  Z-count per class: no universal parity (can be 0,1,2,3,...)
  But: classes with Z=0 tend to have the highest H values
  And: sum = F_{m+2} = Fibonacci number
  The Fibonacci cube Γ_m IS the topology of the Zeckendorf tilings!
""")
