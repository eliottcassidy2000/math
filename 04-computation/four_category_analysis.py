"""
four_category_analysis.py — Systematic analysis of the four tiling categories.

Categories per merged node:
  PURE BLUE:  SC class, blue>0, black=0  (all tilings grid-sym)
  MIXED:      SC class, blue>0, black>0  (some grid-sym, some not)
  PURE BLACK: NS merged pair, blue=0, black>0 (no grid-sym)
  ZERO:       blue=0, black=0 (impossible for any node)

KEY SUMS at n=5:
  pure_blue [1,1,3]=5,  mixed_blue [1,1,3,3,3]=11  → total_blue=16=2^4
  mixed_black [4,6,8,8,10]=36, pure_black [2,10]=12 → total_black=48=64-16

For each mixed node: a (black, blue) PAIR where black+blue=#tilings.

oracle-2026-05-10
"""

import sys
from collections import defaultdict, Counter
from itertools import permutations
sys.stdout.reconfigure(encoding='utf-8')

def make_tiles(n): return [(x,y) for y in range(1,n-1) for x in range(n,y+1,-1)]
def make_trans_map(n,tiles):
    idx={t:i for i,t in enumerate(tiles)}
    return [idx[(n-y+1,n-x+1)] for x,y in tiles]
def bits_to_adj(bits,n,tiles):
    verts=list(range(n,0,-1)); A=[[0]*n for _ in range(n)]
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
def aut_size(A,n,perms):
    return sum(1 for p in perms if all(A[p[i]][p[j]]==A[i][j] for i in range(n) for j in range(n)))
def fib(k):
    if k<=2: return 1
    a,b=1,1
    for _ in range(k-2): a,b=b,a+b
    return b

ALL_DATA = {}

for n in range(3, 7):
    tiles=make_tiles(n); m=len(tiles); tm=make_trans_map(n,tiles); perms=list(permutations(range(n)))
    groups=defaultdict(list)
    for mask in range(1<<m):
        bits=[(mask>>k)&1 for k in range(m)]
        A=bits_to_adj(bits,n,tiles)
        canon=canonicalize(A,n,perms); cc=canonicalize(full_comp(A,n),n,perms)
        gs=is_grid_sym(bits,tm); fm=sum((1-b)<<k for k,b in enumerate(bits)); zk=is_zeck(bits)
        tb=[0]*m
        for i in range(m): tb[tm[i]]=bits[i]
        tmask=sum(b<<k for k,b in enumerate(tb))
        groups[canon].append({'mask':mask,'bits':bits,'A':A,'cc':cc,'gs':gs,'fm':fm,'zk':zk,'tm':tmask})

    classes=[]; sig_to_ci={}
    for ci,(sig,members) in enumerate(sorted(groups.items())):
        A=members[0]['A']; H=H_count(A,n)
        au=aut_size(A,n,perms); sc=(members[0]['cc']==sig)
        Z=sum(1 for t in members if t['zk']); gs=sum(1 for t in members if t['gs'])
        classes.append({'ci':ci,'sig':sig,'m':members,'nt':len(members),'H':H,'au':au,'sc':sc,'cc':members[0]['cc'],'Z':Z,'gs':gs})
        sig_to_ci[sig]=ci
    for c in classes:
        for t in c['m']: t['ci']=c['ci']
    mask_to_ci={t['mask']:c['ci'] for c in classes for t in c['m']}

    # NS pairs
    sc_cis=[c['ci'] for c in classes if c['sc']]
    seen_ns=set(); ns_pairs=[]
    for c in classes:
        if c['sc'] or c['ci'] in seen_ns: continue
        p=next(x for x in classes if x['sig']==c['cc']); seen_ns.add(c['ci']); seen_ns.add(p['ci'])
        ns_pairs.append((c['ci'],p['ci']))
    ci_to_merged={ci:ci for ci in sc_cis}
    for a,b in ns_pairs: ci_to_merged[a]=min(a,b); ci_to_merged[b]=min(a,b)

    # Blue/black per class
    sc_data=[]
    for c in classes:
        if not c['sc']: continue
        blue=c['gs']; black=c['nt']-c['gs']
        cat='pure_blue' if black==0 else 'mixed'
        sc_data.append({'ci':c['ci'],'H':c['H'],'au':c['au'],'nt':c['nt'],'blue':blue,'black':black,'cat':cat,'Z':c['Z']})
    ns_data=[]
    for a,b in ns_pairs:
        ca=next(x for x in classes if x['ci']==a)
        cb=next(x for x in classes if x['ci']==b)
        nt_merged=ca['nt']+cb['nt']
        ns_data.append({'a':a,'b':b,'nt':nt_merged,'black':nt_merged,'H_a':ca['H'],'H_b':cb['H'],'au_a':ca['au'],'au_b':cb['au']})

    # Blue graph adjacency for SC classes
    sc_sorted=sorted(sc_cis); sc_idx={ci:i for i,ci in enumerate(sc_sorted)}; nsc=len(sc_sorted)
    blue_adj=[[0]*nsc for _ in range(nsc)]; seen_bp=set()
    for c in classes:
        for t in c['m']:
            if not t['gs']: continue
            pk=(min(t['mask'],t['fm']),max(t['mask'],t['fm']))
            if pk in seen_bp: continue; seen_bp.add(pk)
            ci_a=t['ci']; ci_b=mask_to_ci[t['fm']]
            if ci_a in sc_idx and ci_b in sc_idx:
                ia,ib=sc_idx[ci_a],sc_idx[ci_b]
                blue_adj[ia][ib]+=1
                if ia!=ib: blue_adj[ib][ia]+=1

    # Blue graph role (ring/pendant/isolated/self-loop)
    blue_deg=[sum(blue_adj[i][j] for j in range(nsc)) for i in range(nsc)]
    self_loop=[blue_adj[i][i]>0 for i in range(nsc)]
    for d in sc_data:
        idx=sc_idx[d['ci']]
        d['blue_deg']=d['blue']  # degree = blue count
        d['self_loop']=self_loop[idx]
        # Find connections
        d['connections']=[sc_sorted[j] for j in range(nsc) if blue_adj[idx][j]>0 and j!=idx]
        d['self_loop_count']=blue_adj[idx][idx]

    # Collect free_bits
    self_paired=sum(1 for i,ti in enumerate(tm) if ti==i)
    non_self=(m-self_paired)//2
    free_bits=self_paired+non_self

    ALL_DATA[n]={'classes':classes,'sc_data':sc_data,'ns_data':ns_data,'m':m,'free_bits':free_bits,
                 'sc_sorted':sc_sorted,'blue_adj':blue_adj,'sc_idx':sc_idx}

# ══════════════════════════════════════════════════════════════════════════
# MAIN OUTPUT
# ══════════════════════════════════════════════════════════════════════════

for n in range(3, 7):
    d=ALL_DATA[n]; sc=d['sc_data']; ns=d['ns_data']
    m=d['m']; fb=d['free_bits']
    total_t=1<<m; total_blue=1<<fb; total_black=total_t-total_blue

    pb=[x for x in sc if x['cat']=='pure_blue']
    mx=[x for x in sc if x['cat']=='mixed']

    print("="*70)
    print(f"n={n}  m={m}  2^m={total_t}  free_bits={fb}  blue={total_blue}  black={total_black}")
    print("="*70)

    # Four categories
    pb_vals=sorted(x['blue'] for x in pb)
    mx_blue=sorted(x['blue'] for x in mx)
    mx_black=sorted(x['black'] for x in mx)
    ns_black=sorted(x['black'] for x in ns)

    print(f"\nFOUR CATEGORIES:")
    print(f"  Pure blue:  {pb_vals}  sum={sum(pb_vals)}")
    print(f"  Mixed blue: {mx_blue}  sum={sum(mx_blue)}")
    print(f"  Mixed blk:  {mx_black}  sum={sum(mx_black)}")
    print(f"  Pure blk:   {ns_black}  sum={sum(ns_black)}")
    print(f"  Check: {sum(pb_vals)}+{sum(mx_blue)}={total_blue}✓  {sum(mx_black)}+{sum(ns_black)}={total_black}{'✓' if sum(mx_black)+sum(ns_black)==total_black else '✗'}")

    # Mixed pairs (black, blue) sorted by black then blue
    mx_pairs=sorted([(x['black'],x['blue'],x['H'],x['au'],x['nt'],x['Z'],x['cat'],x['ci']) for x in mx])
    print(f"\nMIXED PAIRS (black, blue) [total {len(mx_pairs)}]:")
    print(f"  {'(black,blue)':>12} | {'#t':>4} | {'H':>5} | {'|Aut|':>5} | {'Z':>3} | {'role'}")
    for bk,bl,H,au,nt,Z,cat,ci in mx_pairs:
        role='self-loop' if d['sc_data'][next(i for i,x in enumerate(d['sc_data']) if x['ci']==ci)]['self_loop'] else ''
        print(f"  ({bk:>4},{bl:>3}) | {nt:>4} | {H:>5} | {au:>5} | {Z:>3} | {role}")

    # Pure blue details
    print(f"\nPURE BLUE CLASSES:")
    for x in sorted(pb, key=lambda x:x['blue']):
        print(f"  ci{x['ci']}: blue={x['blue']}, H={x['H']}, |Aut|={x['au']}, Z={x['Z']}")
        conns=x['connections']
        if conns:
            conn_details=[(c, next(y['H'] for y in d['sc_data'] if y['ci']==c)) for c in conns]
            print(f"    → blue-connects to: {[(f'ci{c}(H={H})',) for c,H in conn_details]}")

    # Pure black details
    print(f"\nPURE BLACK (NS pairs):")
    for x in sorted(ns, key=lambda x:x['black']):
        print(f"  [{x['a']},{x['b']}]: black={x['black']}, H_a={x['H_a']}(au={x['au_a']}), H_b={x['H_b']}(au={x['au_b']})")

    # Constraint analysis
    print(f"\nCONSTRAINT ANALYSIS:")
    print(f"  Blue total={total_blue}: sum of {len(sc)} odd numbers")
    print(f"  = {sum(pb_vals)} (pure) + {sum(mx_blue)} (mixed)")
    print(f"  = {len(pb)} values summing to {sum(pb_vals)} + {len(mx)} values summing to {sum(mx_blue)}")

    print(f"  Black total={total_black}: sum of {len(sc)+len(ns)} even numbers")
    print(f"  = {sum(mx_black)} (mixed) + {sum(ns_black)} (pure)")
    print(f"  = {len(mx)} positive even values + {len(ns)} positive even values + {sum(1 for x in sc if x['black']==0)} zeros")

    # H-value distribution by category
    print(f"\nH-DISTRIBUTION BY CATEGORY:")
    H_pb=Counter(x['H'] for x in pb)
    H_mx=Counter(x['H'] for x in mx)
    H_ns_a=Counter(x['H_a'] for x in ns)
    H_ns_b=Counter(x['H_b'] for x in ns)
    all_H=sorted(set(list(H_pb.keys())+list(H_mx.keys())+list(H_ns_a.keys())))
    print(f"  H val | pure_blue | mixed | NS_class")
    for H in all_H:
        print(f"  {H:>5} | {H_pb.get(H,0):>9} | {H_mx.get(H,0):>5} | {H_ns_a.get(H,0)+H_ns_b.get(H,0):>8}")

    # Mixed pairs sorted by H
    print(f"\nMIXED PAIRS SORTED BY H THEN |AUT|:")
    for bk,bl,H,au,nt,Z,cat,ci in sorted(mx_pairs, key=lambda x:(x[2],x[3])):
        ratio=bk/bl if bl>0 else float('inf')
        print(f"  H={H:>4}, |Aut|={au}: (bk={bk:>4}, bl={bl:>2}) #t={nt:>4}, bk/bl={ratio:.2f}, bk+bl={bk+bl:>4}")

    print()

# ══════════════════════════════════════════════════════════════════════════
# CROSS-N PATTERNS
# ══════════════════════════════════════════════════════════════════════════
print("="*70)
print("CROSS-n PATTERNS")
print("="*70)

print("\n1. PURE BLUE CLASSES across n:")
print(f"   {'n':>3} | {'pure_blue vals':>30} | {'sum':>6} | {'H values'}")
for n in range(3,7):
    d=ALL_DATA[n]; sc=d['sc_data']
    pb=[x for x in sc if x['cat']=='pure_blue']
    vals=sorted(x['blue'] for x in pb)
    H_vals=sorted(x['H'] for x in pb)
    print(f"   {n:>3} | {str(vals):>30} | {sum(vals):>6} | {H_vals}")

print("\n2. PURE BLACK (NS pairs) across n:")
print(f"   {'n':>3} | {'pure_black vals (merged)':>40} | {'count':>6}")
for n in range(3,7):
    d=ALL_DATA[n]; ns=d['ns_data']
    vals=sorted(x['black'] for x in ns)
    print(f"   {n:>3} | {str(vals):>40} | {len(vals):>6}")

print("\n3. MIXED PAIR COUNTS and SUMS across n:")
print(f"   {'n':>3} | {'#mixed':>7} | {'Σ_blue':>8} | {'Σ_black':>9} | {'tot_blue':>9} | {'tot_black'}")
for n in range(3,7):
    d=ALL_DATA[n]; sc=d['sc_data']
    mx=[x for x in sc if x['cat']=='mixed']
    ns=d['ns_data']; m=d['m']; fb=d['free_bits']
    tot_b=1<<fb; tot_bk=(1<<m)-tot_b
    sb=sum(x['blue'] for x in mx); sbk=sum(x['black'] for x in mx)
    print(f"   {n:>3} | {len(mx):>7} | {sb:>8} | {sbk:>9} | {tot_b:>9} | {tot_bk}")

print("\n4. THE (black, blue) PAIR STRUCTURE:")
print("   At each n, sorted mixed pairs show progression:")
for n in range(3,7):
    d=ALL_DATA[n]; sc=d['sc_data']
    mx=[x for x in sc if x['cat']=='mixed']
    pairs=sorted([(x['black'],x['blue']) for x in mx])
    bl_only=sorted(set(x['blue'] for x in mx))
    print(f"   n={n}: blue values {bl_only}, pairs (bk,bl): {pairs}")

print("\n5. PURE BLUE CONDITION: when is a SC class pure blue?")
print("   (black=0 iff ALL #tilings are grid-symmetric)")
for n in range(3,7):
    d=ALL_DATA[n]; sc=d['sc_data']
    pb=[x for x in sc if x['cat']=='pure_blue']
    print(f"   n={n}: {len(pb)} pure blue classes")
    for x in sorted(pb,key=lambda y:y['nt']):
        reason=[]
        if x['nt']==1: reason.append('#t=1')
        if x['nt']==x['blue']: reason.append('all_grid_sym')
        if x['au']>1: reason.append(f'|Aut|={x["au"]}')
        print(f"     ci{x["ci"]}: H={x["H"]}, |Aut|={x["au"]}, #t={x["nt"]}, blue={x["blue"]}, Z={x["Z"]} [{",".join(reason)}]")

print("\n6. RATIO bk/bl FOR MIXED PAIRS (measures 'how black' each mixed node is):")
for n in range(3,7):
    d=ALL_DATA[n]; sc=d['sc_data']
    mx=sorted([x for x in sc if x['cat']=='mixed'],key=lambda y:(y['black'],y['blue']))
    ratios=[round(x['black']/x['blue'],2) for x in mx]
    print(f"   n={n}: {ratios}")

print("\n7. MIXED PAIR STRUCTURE: (black, blue) pairs with #tilings:")
print("   Key identity: black + blue = #tilings = H/|Aut| (always odd)")
for n in range(3,7):
    d=ALL_DATA[n]; sc=d['sc_data']
    mx=sorted([x for x in sc if x['cat']=='mixed'],key=lambda y:(y['nt'],y['blue']))
    print(f"\n   n={n}:")
    for x in mx:
        print(f"   H={x['H']:>4}, |Aut|={x['au']}, #t={x['nt']:>4}: black={x['black']:>4}+blue={x['blue']:>2}={x['nt']:>4}, Z={x['Z']}, self_loop={x['self_loop']}")

print("\n8. PURE BLUE VALUES as DIVISORS of 2^{free_bits}:")
print("   The pure_blue sum must divide into the total in a specific way.")
for n in range(3,7):
    d=ALL_DATA[n]; sc=d['sc_data']
    pb=[x for x in sc if x['cat']=='pure_blue']
    fb=d['free_bits']; tot_b=1<<fb
    pb_sum=sum(x['blue'] for x in pb)
    mx_sum=tot_b-pb_sum
    print(f"   n={n}: 2^{fb}={tot_b} = {pb_sum}(pure) + {mx_sum}(mixed)")

print("\n9. COUNTS SUMMARY TABLE:")
print(f"   {'n':>2} | {'SC':>4} | {'NS-pr':>6} | {'#PureB':>7} | {'#Mixed':>7} | {'#PureBk':>8} | {'PureB_sum':>9} | {'PureBk_sum':>11}")
for n in range(3,7):
    d=ALL_DATA[n]; sc=d['sc_data']; ns=d['ns_data']
    pb=[x for x in sc if x['cat']=='pure_blue']
    mx=[x for x in sc if x['cat']=='mixed']
    print(f"   {n:>2} | {len(sc):>4} | {len(ns):>6} | {len(pb):>7} | {len(mx):>7} | {len(ns):>8} | {sum(x['blue'] for x in pb):>9} | {sum(x['black'] for x in ns):>11}")

print("\n10. ZECKENDORF vs BLUE COUNT CORRELATION:")
print("   Does higher blue count → more Zeckendorf tilings?")
for n in range(3,7):
    d=ALL_DATA[n]; sc=d['sc_data']
    by_blue=defaultdict(list)
    for x in sc: by_blue[x['blue']].append(x['Z'])
    print(f"   n={n}: blue→Z distribution: {dict(sorted({k:sorted(v) for k,v in by_blue.items()}.items()))}")

print("\n11. PURE BLUE SELF-COMPLEMENT STRUCTURE:")
print("   For pure blue class C (all tilings grid-sym): what is C's flip-partner?")
print("   The flip-partner of a pure-blue class is always MIXED or PURE-BLUE.")
for n in range(3,7):
    d=ALL_DATA[n]; sc=d['sc_data']
    pb=[x for x in sc if x['cat']=='pure_blue']
    for x in pb:
        conns=x['connections']
        conn_cats=[next(y['cat'] for y in sc if y['ci']==c) for c in conns]
        print(f"   n={n} ci{x['ci']}(H={x['H']},blue={x['blue']}): connects to {[(c,next(y['cat'] for y in sc if y['ci']==c)) for c in conns]}")
