"""
zeckendorf_deep2.py — Follow-up exploration based on deep_exploration.py results.

Key findings to pursue:
  1. Two-tile H formula via window intersection geometry
  2. H(chained) = 1 + 2^{r1-1} + 2^{r2-1} - 2^{r1+r2-2} — prove/verify
  3. q=3 Jacobsthal connection: a(m,q=3) = I(P_m,2)
  4. G_n(2)/A000568 = 7/4 at n=4 and n=5
  5. H-distribution GF over Z_n has symmetry?
  6. NEW: two-tile H depends only on span and overlap, not individual positions?

oracle-2026-05-09
"""

import sys
from math import comb, factorial, gcd, isqrt
from fractions import Fraction
from collections import Counter, defaultdict
sys.stdout.reconfigure(encoding='utf-8')

def fib(n):
    if n<=0: return 0
    if n==1: return 1
    a,b=1,1
    for _ in range(n-2): a,b=b,a+b
    return b

def h_val(r): return 1+(1<<(r-1))

def tiles_for_n(n):
    return [(b,b+r) for r in range(2,n) for b in range(n-r)]

def tiling_to_adj(bits, n):
    adj=[[0]*n for _ in range(n)]
    for k in range(n-1): adj[k][k+1]=1
    tiles=tiles_for_n(n)
    for k,(b,a) in enumerate(tiles):
        if (bits>>k)&1: adj[a][b]=1
        else: adj[b][a]=1
    return adj

def compute_H(adj,n):
    dp={}
    for v in range(n): dp[(1<<v,v)]=1
    for ms in range(2,n+1):
        for mask in range(1<<n):
            if bin(mask).count('1')!=ms: continue
            for v in range(n):
                if not(mask&(1<<v)): continue
                pm=mask^(1<<v)
                t=sum(dp.get((pm,u),0) for u in range(n) if(pm&(1<<u))and adj[u][v])
                if t: dp[(mask,v)]=t
    return sum(dp.get(((1<<n)-1,v),0) for v in range(n))

def two_tile_H(n, i, j):
    """Exact H for tournament with exactly tiles i and j active."""
    tiles=tiles_for_n(n)
    bits=(1<<i)|(1<<j)
    return compute_H(tiling_to_adj(bits,n),n)


# ══════════════════════════════════════════════════════════════════════════════
# THREAD 1: Two-tile H formula — depends on span and overlap only?
# ══════════════════════════════════════════════════════════════════════════════
print("="*70)
print("THREAD 1: Two-tile H formula — characterizing the correction")
print("="*70)

print("""
For two tiles (b1,a1) and (b2,a2), define:
  r1 = a1-b1, r2 = a2-b2 (ranges)
  b_min = min(b1,b2), a_max = max(a1,a2)  (total span = a_max - b_min)
  overlap = max(0, min(a1,a2) - max(b1,b2))  (window intersection size)

HYPOTHESIS: H(two tiles) = f(r1, r2, overlap) depends only on ranges and overlap.
""")

# Collect all two-tile H values organized by (r1, r2, overlap)
results = {}
for n in [5, 6, 7]:
    tiles = tiles_for_n(n)
    for i in range(len(tiles)):
        for j in range(i+1, len(tiles)):
            b1,a1 = tiles[i]
            b2,a2 = tiles[j]
            r1,r2 = a1-b1, a2-b2
            if r1>r2: r1,r2,b1,a1,b2,a2 = r2,r1,b2,a2,b1,a1  # ensure r1<=r2
            overlap = max(0, min(a1,a2)-max(b1,b2))
            key = (r1,r2,overlap)
            H = two_tile_H(n, tiles.index((b1 if a1-b1==r1 else b2), (a1 if a1-b1==r1 else a2)),
                               tiles.index((b2 if a1-b1==r1 else b1), (a2 if a1-b1==r1 else a1)))
            if key not in results:
                results[key] = set()
            results[key].add(H)

# Actually let me redo this cleanly
results2 = defaultdict(set)
for n in [5, 6, 7]:
    tiles = tiles_for_n(n)
    for i in range(len(tiles)):
        for j in range(i+1, len(tiles)):
            b1,a1 = tiles[i]
            b2,a2 = tiles[j]
            r1,r2 = a1-b1, a2-b2
            overlap = max(0, min(a1,a2)-max(b1,b2))
            if r1>r2: r1,r2 = r2,r1
            bits = (1<<i)|(1<<j)
            H = compute_H(tiling_to_adj(bits,n),n)
            h_prod = h_val(r1)*h_val(r2)
            results2[(r1,r2,overlap)].add((H, h_prod))

print(f"  {'(r1,r2,overlap)':>20} | {'H values seen':>30} | {'H_prod':>8} | {'H unique?'}")
print("-"*75)
for key in sorted(results2.keys()):
    r1,r2,ov = key
    vals = results2[key]
    H_vals = sorted(set(v[0] for v in vals))
    H_prod = vals.__iter__().__next__()[1]
    unique = len(H_vals)==1
    note = "UNIQUE" if unique else f"{len(H_vals)} values"
    if r1+r2 <= 7:
        print(f"  ({r1},{r2},{ov:>2}) {'':>8} | {str(H_vals):>30} | {H_prod:>8} | {note}")

print("""
KEY OBSERVATION:
For each (r1, r2, overlap), is H unique? If YES, then H = f(r1, r2, overlap).
""")

# Find cases where H is NOT unique for same (r1,r2,overlap)
non_unique = {k:v for k,v in results2.items() if len(set(x[0] for x in v))>1}
if non_unique:
    print("  Cases where H is NOT unique for same (r1,r2,overlap):")
    for key, vals in sorted(non_unique.items()):
        r1,r2,ov = key
        H_vals = sorted(set(v[0] for v in vals))
        print(f"    (r1={r1},r2={r2},overlap={ov}): H ∈ {H_vals}")
    print()
    print("  → H depends on MORE than just (r1, r2, overlap)!")
    print("  → What additional geometric parameter distinguishes these cases?")
else:
    print("  All (r1,r2,overlap) keys give UNIQUE H values! H = f(r1,r2,overlap). ✓")

# If non-unique, classify by more refined geometry
print()
print("  Refined analysis: classify by (r1, r2, ov, b_min_shift) where")
print("  b_min_shift = max(b1,b2) - min(b1,b2) = relative starting position shift")

results3 = defaultdict(set)
for n in [5, 6, 7]:
    tiles = tiles_for_n(n)
    for i in range(len(tiles)):
        for j in range(i+1, len(tiles)):
            b1,a1 = tiles[i]
            b2,a2 = tiles[j]
            r1,r2 = a1-b1, a2-b2
            overlap = max(0, min(a1,a2)-max(b1,b2))
            shift = abs(b2-b1)  # relative shift
            if r1>r2: r1,r2 = r2,r1
            bits=(1<<i)|(1<<j)
            H = compute_H(tiling_to_adj(bits,n),n)
            results3[(r1,r2,overlap,shift)].add(H)

non_unique3 = {k:v for k,v in results3.items() if len(v)>1}
if non_unique3:
    print("  Still non-unique after adding shift:")
    for key,vals in sorted(non_unique3.items()):
        print(f"    {key}: H ∈ {sorted(vals)}")
else:
    print("  With (r1,r2,overlap,shift): ALL keys give unique H! H = f(r1,r2,ov,shift). ✓")

# Find the actual formula
print()
print("  Two-tile H formula table:")
print(f"  {'(r1,r2,ov,shift)':>25} | {'H':>5} | {'formula?'}")
formula_checks = {}
for key in sorted(results3.keys()):
    r1,r2,ov,shift = key
    if r1+r2>7: continue
    H_vals = sorted(results3[key])
    H = H_vals[0]
    # Try formula: H = 1 + 2^{r1-1} + 2^{r2-1} - 2^{r1+r2-ov-2}?
    span = r1+r2-ov  # effective span considering overlap
    formula1 = 1 + h_val(r1) + h_val(r2) - h_val(r1+r2-ov)
    formula2 = h_val(span)
    formula3 = 1 + 2**(r1-1) + 2**(r2-1) - 2**(r1+r2-ov-2+1)
    ok1 = '✓' if H==formula1 else f'?({formula1})'
    ok2 = '✓' if H==formula2 else f'?({formula2})'
    ok3 = '✓' if H==formula3 else f'?({formula3})'
    print(f"  ({r1},{r2},{ov},{shift}) | H={H:>4} | "
          f"span={r1+r2-ov} h(span)={ok2} f1={ok1}")

# ══════════════════════════════════════════════════════════════════════════════
# THREAD 2: Prove H(two tiles) = h(span) formula (h=1+2^{r-1})
# ══════════════════════════════════════════════════════════════════════════════
print()
print("="*70)
print("THREAD 2: Is H(two tiles) = h(effective_span) ?")
print("="*70)

print("""
HYPOTHESIS: For two tiles of ranges r1, r2 with overlap d (in vertex count):
  effective_span = r1 + r2 - d
  H(two tiles) = h(effective_span) = 1 + 2^{effective_span - 1}

Note: d=0 (disjoint) → span = r1+r2, H_prod = h(r1)*h(r2)
      So h(r1+r2) ≠ h(r1)*h(r2) in general.
      h(r1)*h(r2) = (1+2^{r1-1})(1+2^{r2-1}) = 1+2^{r1-1}+2^{r2-1}+2^{r1+r2-2}
      h(r1+r2) = 1+2^{r1+r2-1}

For disjoint tiles: the hypothesis H=h(span) would give 1+2^{r1+r2-1},
but the PRODUCT FORMULA gives h(r1)*h(r2) = 1+2^{r1-1}+2^{r2-1}+2^{r1+r2-2}.
These differ! (e.g., r1=r2=2: h(4)=9 vs h(2)*h(2)=9 — equal! But r1=r2=3: h(6)=33 vs h(3)*h(3)=25 ≠ 33.)

So for disjoint tiles, H = h(r1)*h(r2) (product formula), NOT h(r1+r2).
For overlapping tiles, let's check h(r1+r2-overlap):
""")

# Systematic check
print("  Checking H = h(r1+r2-overlap) for all two-tile pairs:")
formula_ok = True
for n in [5,6,7]:
    tiles=tiles_for_n(n)
    for i in range(len(tiles)):
        for j in range(i+1,len(tiles)):
            b1,a1=tiles[i]; b2,a2=tiles[j]
            r1,r2=a1-b1,a2-b2
            ov=max(0,min(a1,a2)-max(b1,b2))
            bits=(1<<i)|(1<<j)
            H=compute_H(tiling_to_adj(bits,n),n)
            span=r1+r2-ov
            pred=h_val(span)
            if H!=pred and r1+r2<=8:
                if ov>0:  # only check overlapping
                    print(f"    n={n} tiles {tiles[i]},{tiles[j]}: r1={r1},r2={r2},ov={ov} "
                          f"span={span} H={H} h(span)={pred} {'✓' if H==pred else '✗'}")

# Check for disjoint tiles separately
print()
print("  Checking H = h(r1)*h(r2) (product formula) for disjoint tiles:")
disjoint_ok = True
for n in [5,6,7]:
    tiles=tiles_for_n(n)
    for i in range(len(tiles)):
        for j in range(i+1,len(tiles)):
            b1,a1=tiles[i]; b2,a2=tiles[j]
            r1,r2=a1-b1,a2-b2
            ov=max(0,min(a1,a2)-max(b1,b2))
            if ov>0: continue  # only disjoint
            bits=(1<<i)|(1<<j)
            H=compute_H(tiling_to_adj(bits,n),n)
            pred=h_val(r1)*h_val(r2)
            if H!=pred:
                print(f"    FAILS: n={n} tiles {tiles[i]},{tiles[j]}: r1={r1},r2={r2} H={H} prod={pred}")
                disjoint_ok=False

if disjoint_ok:
    print("    ALL DISJOINT: H = h(r1)*h(r2) ✓")


# ══════════════════════════════════════════════════════════════════════════════
# THREAD 3: The UNIFIED two-tile formula
# ══════════════════════════════════════════════════════════════════════════════
print()
print("="*70)
print("THREAD 3: Unified two-tile formula — OCF interpretation")
print("="*70)

print("""
From the data, for two tiles of ranges r1, r2 with overlap d vertices:
  H_approx = h(r1)*h(r2) [product of single-tile values]
  H_true = h(r1)*h(r2) if d=0 [disjoint: product formula exact]
  H_true = ? if d>0 [overlap: formula below]

Let's find the formula. The key quantity is:
  α_1 = (H-1)/2 [since H is always odd for tournaments]
  α_1_approx = (H_approx-1)/2 = α_1(r1) + α_1(r2) + α_1(r1)*α_1(r2)
             = 2^{r1-2} + 2^{r2-2} + 2^{r1+r2-4}

For overlapping tiles, the correction is to α_1 and α_2 separately.
""")

for n in [6,7]:
    tiles = tiles_for_n(n)
    print(f"\n  n={n}: Two-tile α values:")
    seen = set()
    for i in range(len(tiles)):
        for j in range(i+1,len(tiles)):
            b1,a1=tiles[i]; b2,a2=tiles[j]
            r1,r2=a1-b1,a2-b2
            ov=max(0,min(a1,a2)-max(b1,b2))
            if r1>r2: r1,r2=r2,r1
            key=(r1,r2,ov)
            if key in seen or r1+r2>8: continue
            seen.add(key)
            bits=(1<<i)|(1<<j)
            H=compute_H(tiling_to_adj(bits,n),n)
            # H = 1 + 2*α_1 + 4*α_2 (no higher for 2 tiles)
            # We know α_2 from: H = 1 + 2*(α_1(r1)+α_1(r2)+Δα_1) + 4*(α_1(r1)*α_1(r2)-Δα_2)
            # Actually: need to find α_1_true and α_2_true
            # For now: H = 1 + 2α_1 + 4α_2. But H is the total.
            # We know H_approx = (1+2α_1^a)(1+2α_1^b) = 1+2(α_1^a+α_1^b)+4α_1^aα_1^b
            # where α_1^a = 2^{r1-2}, α_1^b = 2^{r2-2}
            a1a = 2**(r1-2) if r1>=2 else 0
            a1b = 2**(r2-2) if r2>=2 else 0
            H_prod = h_val(r1)*h_val(r2)
            # True H: what are α_1_true and α_2_true?
            # 1 + 2*α_1_true + 4*α_2_true = H
            # We don't know individual values, just their combination.
            # But: α_1_true <= α_1^a + α_1^b (can gain or lose cycles)
            # And: α_2_true <= α_1_true*(α_1_true-1)/2 (independence condition)
            delta_H = H_prod - H
            delta_a1 = (H_prod-H)//4 if ov==0 else None
            if ov>0:
                # H = 1 + 2*a1_total + 4*a2_total
                # a1_total ≤ a1^a + a1^b (some cycles may be shared/lost)
                # For overlap: some cycles from r1 and r2 share vertices → lose α_2 contribution
                pass
            print(f"    (r1={r1},r2={r2},ov={ov}): H={H}, H_prod={H_prod}, ΔH={H_prod-H}, "
                  f"H-1={H-1}, (H-1)/2={(H-1)//2}, span={r1+r2-ov}")

# Key test: for overlapping tiles, does H = h(span)?
print()
print("  Testing H = h(span) for ALL overlapping pairs:")
all_match = True
for n in [5,6,7]:
    tiles=tiles_for_n(n)
    for i in range(len(tiles)):
        for j in range(i+1,len(tiles)):
            b1,a1=tiles[i]; b2,a2=tiles[j]
            r1,r2=a1-b1,a2-b2
            ov=max(0,min(a1,a2)-max(b1,b2))
            if ov==0: continue
            bits=(1<<i)|(1<<j)
            H=compute_H(tiling_to_adj(bits,n),n)
            span=r1+r2-ov
            pred=h_val(span)
            if H!=pred:
                if r1+r2<=9:
                    print(f"  FAILS at n={n}: ({tiles[i]},{tiles[j]}) "
                          f"r1={r1},r2={r2},ov={ov},span={span}: H={H}≠h({span})={pred}")
                all_match=False
if all_match:
    print("  ALL OVERLAPPING: H = h(span) = h(r1+r2-overlap) ✓")


# ══════════════════════════════════════════════════════════════════════════════
# THREAD 4: Jacobsthal connection — q=3 count = I(P_m, 2)
# ══════════════════════════════════════════════════════════════════════════════
print()
print("="*70)
print("THREAD 4: Jacobsthal-Zeckendorf-OCF triangle")
print("="*70)

print("""
From Thread H: count of q-state 'no-two-adjacent-nonzero' strings of length m:
  a(m,q) = a(m-1,q) + (q-1)*a(m-2,q)
  For q=3: a(m,3) = a(m-1,3) + 2*a(m-2,3) → Jacobsthal recurrence
  a(m,3) = (2^{m+2} - (-1)^m) / 3 = I(P_m, 2)

AND I(P_m, 2) = H(T) where Omega(T) = path graph P_m.

QUESTION: What tournament T has Omega(T) = P_m?

For Omega = path, we need m odd cycles C_1,...,C_m where consecutive
pairs (C_i, C_{i+1}) conflict but non-consecutive don't.

At n=5 (m=C(4,2)=6 tiles, conflict graph Omega has up to C(alpha_1, 2) edges):
  For a single range-5 tile at n=6 (range r=5):
    α_1 = 2^{5-2} = 8 cycles
    I(Omega_{range5}, 2) = H(range-5 tile) = 1+2^4 = 17 (α_2=0 for single tile? or not)

Actually, for a single backward arc of range r, the conflict graph Omega_r
has α_1(r) = 2^{r-2} vertices. What's the edge structure?
Omega_r = WHAT GRAPH?
""")

# Compute the OCF conflict graph for single tiles
def single_tile_conflict_graph(n, tile_idx):
    """Find odd cycles and their conflicts for a single tile."""
    tiles = tiles_for_n(n)
    b, a = tiles[tile_idx]
    r = a - b
    bits = 1 << tile_idx
    adj = tiling_to_adj(bits, n)
    H = compute_H(adj, n)
    print(f"  Tile ({b},{a}) range {r} at n={n}: H={H}, h_val={h_val(r)}, α_1={2**(r-2) if r>=2 else 0}")
    return H

for n in [5, 6, 7]:
    tiles = tiles_for_n(n)
    print(f"\n  Single-tile H values at n={n}:")
    for i, (b,a) in enumerate(tiles):
        r = a-b
        bits = 1<<i
        adj = tiling_to_adj(bits, n)
        H = compute_H(adj, n)
        h = h_val(r)
        # alpha_1 = 2^{r-2}, alpha_2 = 0? Let's verify via H = 1+2*alpha_1
        alpha_1_from_H = (H-1)//2
        alpha_1_formula = 2**(r-2) if r>=2 else 0
        ok = '✓' if H==h and alpha_1_from_H==alpha_1_formula else '✗'
        print(f"    ({b},{a}) r={r}: H={H}=h({r})={h} α_1={alpha_1_from_H}=2^{{{r-2}}}={alpha_1_formula} {ok}")

print("""
For a single tile of range r: H = 1+2^{r-1} = 1+2*2^{r-2}.
So α_1 = 2^{r-2} (the number of independent odd cycles).
And α_2 = 0: all cycles from a single tile CONFLICT pairwise (Omega = complete graph K_{α_1}).

I(K_{α_1}, 2) = 1 + α_1*2 = 1 + 2*2^{r-2} = 1+2^{r-1} ✓

So Omega_r = K_{2^{r-2}} (complete graph on 2^{r-2} vertices).
This makes sense: all cycles from a single arc use the SAME arc, so they all conflict.
""")

# The Jacobsthal connection:
# I(P_m, 2) = (2^{m+2} - (-1)^m)/3
# If we have m SUCCESSIVE single-range-2 tiles forming a PATH of conflicts:
# Omega = P_m where each vertex represents one 3-cycle from a range-2 tile
# (range-2 tiles have exactly 1 cycle each = K_1 = single vertex in Omega)
# P_m has m vertices and m-1 edges.
# I(P_m, 2) = (2^{m+2}-(-1)^m)/3 ← Jacobsthal

print()
print("  Jacobsthal at x=2 interpretation:")
print("  I(P_m, 2) = #{independent sets in path P_m, weighted by 2^size}")
for m in range(1, 10):
    ip = (2**(m+2)-(-1)**m)//3
    print(f"    m={m}: I(P_m,2) = {ip}", end="")
    # What tournament would have Omega = P_m?
    # P_m path: m vertices, edges between consecutive.
    # If m = α_1 of a tournament with linear cycle structure...
    print()

print("""
KEY INSIGHT:
  If we could construct a tournament T with Omega(T) = P_m (path on m vertices),
  then H(T) = I(P_m, 2) = Jacobsthal(m).

  A tournament with m backward tiles of range 2, CHAINED:
    Tile_1: (0,2)r=2, Tile_2: (2,4)r=2, ..., Tile_m: (2m-2, 2m)r=2
    All chained: endpoint of tile k = start of tile k+1.

  Does such a tournament have Omega = P_m?

  Each range-2 tile contributes 1 cycle (C_k = cycle using tile k's arc).
  Consecutive chained tiles share a vertex → conflict.
  Non-consecutive chained tiles don't share vertices → no conflict.
  → Omega = PATH P_m! ✓

  H = I(P_m, 2) = Jacobsthal(m). ✓
""")

# Verify: m chained range-2 tiles give H = Jacobsthal(m)
def jacobsthal(m):
    return (2**(m+2)-(-1)**m)//3

print("  Verification: m chained range-2 tiles give H = Jacobsthal(m)?")
for n in range(3, 10):
    tiles = tiles_for_n(n)
    # Find chained range-2 tiles: (0,2),(2,4),(4,6),...
    chained = []
    for k in range(n-2):
        tile = (2*k, 2*k+2)  # range-2 tiles at positions 0,2,4,...
        if tile in tiles:
            chained.append(tiles.index(tile))
    m_chain = len(chained)
    if m_chain == 0: continue
    # Activate all chained tiles
    bits = sum(1<<k for k in chained)
    adj = tiling_to_adj(bits, n)
    H = compute_H(adj, n)
    J = jacobsthal(m_chain)
    ok = '✓' if H==J else '✗'
    print(f"  n={n}: {m_chain} chained tiles → H={H}, J({m_chain})={J} {ok}")


# ══════════════════════════════════════════════════════════════════════════════
# THREAD 5: The G_n(2)/A000568 = 7/4 mystery
# ══════════════════════════════════════════════════════════════════════════════
print()
print("="*70)
print("THREAD 5: G_n(2)/A000568(n) = 7/4 at n=4 AND n=5")
print("="*70)

from math import factorial
from collections import Counter

def cycle_types_with_counts(n):
    def partitions(n, max_val=None):
        if max_val is None: max_val=n
        if n==0: yield (); return
        for k in range(min(n,max_val),0,-1):
            for rest in partitions(n-k,k): yield (k,)+rest
    for p in partitions(n):
        c=Counter(p)
        cnt=factorial(n)
        for l,freq in c.items(): cnt//=(l**freq*factorial(freq))
        yield p,cnt

def B_of_partition(p):
    parts=list(p)
    return sum((l-1)//2 for l in parts)+sum(gcd(parts[i],parts[j]) for i in range(len(parts)) for j in range(i+1,len(parts)))

def G_n(n,q):
    total=0
    for p,cnt in cycle_types_with_counts(n): total+=cnt*q**B_of_partition(p)
    return total//factorial(n)

A000568={1:1,2:1,3:2,4:4,5:12,6:56,7:456,8:6880}

print()
print("  G_n(2) = (1/n!) Σ_σ 2^{B(σ)} where sum is over ALL σ")
print("  A000568(n) = (1/n!) Σ_{σ: all cycles odd} 2^{B(σ)}")
print()
print("  The difference D_n = G_n(2) - A000568(n) counts:")
print("  D_n = (1/n!) Σ_{σ: some even cycle} 2^{B(σ)}")

for n in range(2, 9):
    g2 = G_n(n,2)
    a = A000568.get(n,0)
    D = g2 - a
    r = Fraction(g2, a) if a else None
    print(f"  n={n}: G_n(2)={g2}, A000568={a}, D={D}, ratio={r}")

print("""
OBSERVATION: ratio=7/4 at n=4 and n=5. Let's understand WHY.

At n=4: G_4(2)=7, A000568(4)=4. Ratio = 7/4.
At n=5: G_5(2)=21, A000568(5)=12. Ratio = 21/12 = 7/4.

This means: G_5(2)/G_4(2) = 21/7 = 3 = A000568(5)/A000568(4) = 12/4 = 3.

So both G_n(2) and A000568(n) multiply by 3 from n=4 to n=5!
""")

print("  Consecutive ratios:")
g_prev = G_n(1,2)
a_prev = A000568[1]
for n in range(2, 9):
    g = G_n(n,2)
    a = A000568.get(n,0)
    rg = Fraction(g, g_prev)
    ra = Fraction(a, a_prev)
    print(f"  n={n}: G_n/G_{{n-1}}={rg}, A000568(n)/A000568(n-1)={ra} {'SAME' if rg==ra else ''}")
    g_prev, a_prev = g, a


# ══════════════════════════════════════════════════════════════════════════════
# THREAD 6: H-distribution GF over Z_n — generating function structure
# ══════════════════════════════════════════════════════════════════════════════
print()
print("="*70)
print("THREAD 6: GF Σ_{T∈Z_n} x^{H(T)} — coefficient structure")
print("="*70)

def H_dist_Zn(n):
    from itertools import product
    tiles=tiles_for_n(n)
    m=len(tiles)
    dist=Counter()
    def gen(pos, last):
        if pos==m:
            yield ()
            return
        for bit in (0,1):
            if bit==1 and last==1: continue
            for rest in gen(pos+1,bit):
                yield (bit,)+rest
    for bits_tuple in gen(0,0):
        bits_int=sum(b<<k for k,b in enumerate(bits_tuple))
        adj=tiling_to_adj(bits_int,n)
        H=compute_H(adj,n)
        dist[H]+=1
    return dist

for n in [4,5,6]:
    dist=H_dist_Zn(n)
    sorted_items=sorted(dist.items())
    H_vals=[h for h,c in sorted_items]
    coeffs=[c for h,c in sorted_items]
    total=sum(coeffs)
    sumH=sum(h*c for h,c in sorted_items)
    maxH=max(H_vals)
    print(f"\n  n={n}: GF = Σ c_H * x^H, |Z_n|={total}:")
    print(f"  H values: {H_vals}")
    print(f"  Coefficients: {coeffs}")
    print(f"  Sum H = {sumH} = {total} * {Fraction(sumH,total)}")
    print(f"  Max H = {maxH}, Min H = {min(H_vals)}")

    # Check: do coefficients at H and H_MAX+2-H match? (Palindrome test)
    print(f"  Palindrome test (c_H = c_{{maxH+2-H}}):")
    is_palindrome = True
    for h,c in sorted_items:
        complement = maxH + 2 - h
        c_comp = dist.get(complement, 0)
        if c != c_comp:
            is_palindrome = False
            print(f"    c_{h}={c} vs c_{complement}={c_comp} ✗")
    if is_palindrome:
        print(f"    PALINDROME SYMMETRIC! ✓ c_H = c_{{{maxH}+2-H}}")

    # Check: H values are all ≡ 1 (mod 2) [all odd]
    all_odd = all(h%2==1 for h in H_vals)
    print(f"  All H odd: {all_odd}")

    # Check: H ≡ 1 (mod 4) for all?
    mod4 = set(h%4 for h in H_vals)
    print(f"  H mod 4: {mod4}")

    # Check: H values form arithmetic/geometric progression?
    diffs = [H_vals[i+1]-H_vals[i] for i in range(len(H_vals)-1)]
    print(f"  Gaps between consecutive H: {diffs}")


# ══════════════════════════════════════════════════════════════════════════════
# THREAD 7: New result — H-distribution over Z_n has palindrome symmetry?
# ══════════════════════════════════════════════════════════════════════════════
print()
print("="*70)
print("THREAD 7: Complement involution on Z_n — H symmetry")
print("="*70)

print("""
The 'Zeckendorf complement' of a tiling z = (z_0,...,z_{m-1}) is
z* = (1-z_{m-1},...,1-z_0) — reverse and flip.

Is z* Zeckendorf whenever z is? (Z_n closed under this involution?)
What is H(T_{z*}) vs H(T_z)?
""")

for n in [5, 6]:
    tiles = tiles_for_n(n)
    m = len(tiles)
    pairs = []
    seen = set()
    def gen_zeck(pos, last):
        if pos == m:
            yield ()
            return
        for bit in (0,1):
            if bit==1 and last==1: continue
            for rest in gen_zeck(pos+1,bit):
                yield (bit,)+rest
    for z in gen_zeck(0,0):
        # Compute H(z)
        bits_int=sum(b<<k for k,b in enumerate(z))
        adj=tiling_to_adj(bits_int,n)
        H=compute_H(adj,n)
        # Reverse-and-flip
        zstar=tuple(1-b for b in reversed(z))
        # Is zstar Zeckendorf?
        valid_star=not any(zstar[k]==1 and zstar[k+1]==1 for k in range(m-1))
        if valid_star:
            bits_star=sum(b<<k for k,b in enumerate(zstar))
            adj_star=tiling_to_adj(bits_star,n)
            Hstar=compute_H(adj_star,n)
            pairs.append((z,H,zstar,Hstar))

    closed=len(pairs)==fib(m+2)
    print(f"\n  n={n}: Z_n closed under reverse-flip? {closed} ({len(pairs)}/{fib(m+2)})")
    if pairs:
        h_vs_hstar=[(p[1],p[3]) for p in pairs]
        same=sum(1 for h,hs in h_vs_hstar if h==hs)
        print(f"  H(z) = H(z*) always? {same}/{fib(m+2)} cases")
        # Distribution of H(z*) vs H(z)
        diff_dist=Counter(hs-h for h,hs in h_vs_hstar)
        print(f"  H(z*)-H(z) distribution: {dict(sorted(diff_dist.items()))}")


# ══════════════════════════════════════════════════════════════════════════════
# THREAD 8: The KEY FORMULA — H for ALL Zeckendorf tilings via pair structure
# ══════════════════════════════════════════════════════════════════════════════
print()
print("="*70)
print("THREAD 8: H formula for ALL Zeckendorf tilings at n=5")
print("="*70)

print("""
At n=5 with tiles (0,2),(1,3),(2,4),(0,3),(1,4),(0,4) [ranges 2,2,2,3,3,4]:

For any Zeckendorf tiling, the ACTIVE tiles have windows that may overlap.
Using the formulas:
  - Single tile range r: H = h(r) = 1+2^{r-1}
  - Two disjoint tiles: H = h(r1)*h(r2)
  - Two overlapping tiles with overlap d: H = h(r1+r2-d)

Can we compute H for ALL tilings from this?

For 3+ tiles, we need higher-order corrections...
""")

# Compute and verify using the two-tile formula
def H_from_window_geometry(active_tiles):
    """Try to compute H from window geometry (2-tile formula extended)."""
    if len(active_tiles) == 0: return 1
    if len(active_tiles) == 1:
        b,a = active_tiles[0]
        return h_val(a-b)
    if len(active_tiles) == 2:
        (b1,a1),(b2,a2) = active_tiles
        ov = max(0,min(a1,a2)-max(b1,b2))
        r1,r2 = a1-b1,a2-b2
        if ov == 0:
            return h_val(r1)*h_val(r2)  # product formula
        else:
            return h_val(r1+r2-ov)  # overlap formula
    return None  # no formula yet for 3+

# Test at n=5
n=5
tiles=tiles_for_n(n)
m=len(tiles)
print(f"  Testing 2-tile formula predictions at n={n}:")
for z in gen_zeck(0,0) if False else []:
    pass

total_ok=0; total_tested=0
for bits_tuple in list(__import__('itertools').product([0,1], repeat=m)):
    if any(bits_tuple[k]==1 and bits_tuple[k+1]==1 for k in range(m-1)):
        continue  # not Zeckendorf
    active=[tiles[k] for k,b in enumerate(bits_tuple) if b]
    if len(active)!=2: continue
    total_tested+=1
    bits_int=sum(b<<k for k,b in enumerate(bits_tuple))
    adj=tiling_to_adj(bits_int,n)
    H=compute_H(adj,n)
    pred=H_from_window_geometry(active)
    if H==pred: total_ok+=1

print(f"  Two-tile Zeckendorf tilings: {total_ok}/{total_tested} match formula")
print()

# What about THREE-tile tilings? Do they satisfy H = h(span) or product formula?
print("  Three-tile Zeckendorf tilings at n=5:")
for bits_tuple in __import__('itertools').product([0,1], repeat=m):
    if any(bits_tuple[k]==1 and bits_tuple[k+1]==1 for k in range(m-1)):
        continue
    active=[tiles[k] for k,b in enumerate(bits_tuple) if b]
    if len(active)!=3: continue
    bits_int=sum(b<<k for k,b in enumerate(bits_tuple))
    adj=tiling_to_adj(bits_int,n)
    H=compute_H(adj,n)
    ranges=sorted([a-b for b,a in active])
    # Compute span and overlaps
    all_b=min(b for b,a in active)
    all_a=max(a for b,a in active)
    span=all_a-all_b
    h_span=h_val(span)
    print(f"  active={active} ranges={ranges} H={H} span={span} h(span)={h_span} match={'✓' if H==h_span else '✗'}")
