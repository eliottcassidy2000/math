"""
Connections between the summand graph, metagraph structure, and recursive formulas.
"""
import sys, os
sys.path.insert(0, '/home/ubuntu/math/03-artifacts/code')
from tournament_lib import all_tournaments, hamiltonian_path_count
from itertools import permutations
from collections import defaultdict, Counter
from math import factorial, comb

def H(T): return hamiltonian_path_count(T)

# ─────────────────────────────────────────────
# PART 1: H-SPECTRUM AND SUMMAND STRUCTURE
# ─────────────────────────────────────────────
print("=== H-SPECTRUM AT EACH n ===")
spectra = {}
for n in range(2, 8):
    h_vals = sorted(set(H(T) for T in all_tournaments(n)))
    spectra[n] = h_vals
    forb = [x for x in [7,21] if x <= h_vals[-1] and x not in h_vals]
    print(f"n={n}: spectrum={h_vals}")
    print(f"  count={len(h_vals)}, range=[{h_vals[0]},{h_vals[-1]}], forbidden in range: {forb}")

# What new H values appear at each step?
print("\n=== NEW H-VALUES INTRODUCED AT EACH n ===")
prev = set(spectra[2])
for n in range(3, 8):
    curr = set(spectra[n])
    new_vals = sorted(curr - prev)
    print(f"n={n}: NEW values = {new_vals}")
    # Which of these are NOT products of smaller H-values?
    smaller_prods = set()
    for a in prev:
        for b in prev:
            smaller_prods.add(a*b)
    non_product_new = [v for v in new_vals if v not in smaller_prods]
    print(f"       Not products of earlier H-values: {non_product_new}")
    prev = curr

# ─────────────────────────────────────────────
# PART 2: SUMMAND PAIRS FROM DELETION-CONTRACTION
# ─────────────────────────────────────────────
print("\n=== DELETION-CONTRACTION SUMMAND PAIRS ===")
# H(T) = H(T\e) + H(T/e). Here T\e is T with arc (u,v) removed (0 in both directions)
# and T/e is T contracted at (u,v): merge u,v into new vertex w.

def ham_paths_digraph(adj):
    """Count Hamiltonian paths in a digraph (adj = list of lists, 0/1 entries)."""
    n = len(adj)
    count = 0
    for perm in permutations(range(n)):
        valid = all(adj[perm[i]][perm[i+1]] == 1 for i in range(n-1))
        if valid:
            count += 1
    return count

def delete_arc(T, u, v):
    n = len(T)
    T2 = [row[:] for row in T]
    T2[u][v] = 0  # remove the arc; no arc in either direction
    return T2

def contract_arc(T, u, v):
    """Contract u→v: new vertex w inherits IN-edges from u, OUT-edges from v."""
    n = len(T)
    others = [i for i in range(n) if i != u and i != v]
    m = n - 1
    new_idx = {x: i+1 for i, x in enumerate(others)}
    w = 0
    T2 = [[0]*m for _ in range(m)]
    for x in others:
        xi = new_idx[x]
        if T[x][u]: T2[xi][w] = 1      # x beats u → x beats w
        if T[v][x]: T2[w][xi] = 1      # v beats x → w beats x
        for y in others:
            if y != x: T2[xi][new_idx[y]] = T[x][y]
    return T2

print("Verifying and cataloguing DC pairs:")
for n in range(3, 7):
    pair_counter = Counter()
    for T in all_tournaments(n):
        HT = H(T)
        # Test one representative arc per tournament
        for u in range(n):
            for v in range(n):
                if T[u][v]:
                    T_del = delete_arc(T, u, v)
                    T_con = contract_arc(T, u, v)
                    H_del = ham_paths_digraph(T_del)
                    H_con = ham_paths_digraph(T_con)
                    assert H_del + H_con == HT, f"DC violated! {H_del}+{H_con}≠{HT}"
                    pair_counter[(min(H_del,H_con), max(H_del,H_con))] += 1
                    break
            else: continue
            break
    print(f"n={n}: Top DC summand pairs (a,b) with a+b=H(T):")
    for (a,b), cnt in sorted(pair_counter.items(), key=lambda x:-x[1])[:8]:
        print(f"  {a}+{b}={a+b}  (occurs {cnt}x)")

# ─────────────────────────────────────────────
# PART 3: METAGRAPH ΔH DISTRIBUTION
# ─────────────────────────────────────────────
print("\n=== METAGRAPH ΔH DISTRIBUTION (arc flip changes) ===")

def flip(T, u, v):
    T2 = [row[:] for row in T]
    T2[u][v] = 0; T2[v][u] = 1
    return T2

for n in range(3, 8):
    n_verts = n
    dH_vals = Counter()
    total_flips = 0
    for T in all_tournaments(n):
        HT = H(T)
        for u in range(n):
            for v in range(n):
                if T[u][v]:
                    T2 = flip(T, u, v)
                    dH_vals[H(T2) - HT] += 1
                    total_flips += 1
    # Each flip is double-counted (u→v and v→u are inverses)
    pos = {k: v//2 for k,v in dH_vals.items() if k > 0}
    zero = dH_vals[0] // 2
    print(f"n={n}: ΔH distribution (positive changes): {dict(sorted(pos.items()))}")
    print(f"  Level edges (ΔH=0): {zero}")
    
    # Check: are all |ΔH| values EVEN?
    all_even = all(k % 2 == 0 for k in dH_vals.keys())
    print(f"  All |ΔH| even: {all_even}  (max |ΔH|={max(abs(k) for k in dH_vals)})")

# ─────────────────────────────────────────────
# PART 4: THE ΣH RECURSION
# ─────────────────────────────────────────────
print("\n=== ΣH OVER ALL TILINGS — RECURSION ANALYSIS ===")
# From THM-292: ΣH(n) = Σ_{succession-free σ} 2^{C(n-1,2)-n+1+bp(σ)}
# A000255 recursion: a(n) = (n+1)*a(n-1) - (n-2)*a(n-2)

def succession_free_bp(n):
    """Count succession-free perms of [0..n-1] by bp (base-path descent-by-1 count)."""
    bp_count = Counter()
    for perm in permutations(range(n)):
        # Succession: perm[i+1] = perm[i]+1
        if any(perm[i+1] == perm[i]+1 for i in range(n-1)):
            continue
        # bp: perm[i] = perm[i+1]+1 (descent by 1 = base-path arc)
        bp = sum(1 for i in range(n-1) if perm[i] == perm[i+1]+1)
        bp_count[bp] += 1
    return bp_count

print("Succession-free perm distribution by bp:")
for n in range(2, 9):
    bpd = succession_free_bp(n)
    total = sum(bpd.values())
    sum_H = sum(cnt * (2 ** (comb(n-1,2)-n+1+bp)) for bp, cnt in bpd.items())
    print(f"n={n}: A000255={total}, bp_dist={dict(sorted(bpd.items()))}, ΣH={sum_H}")

# Verify A000255 recurrence: a(n) = (n+1)*a(n-1) - (n-2)*a(n-2)
print("\nA000255 recurrence verification:")
a = [0, 1, 3, 11, 53, 309, 2119, 16687]  # A000255
for n in range(4, 8):
    lhs = a[n]
    rhs = (n+1)*a[n-1] - (n-2)*a[n-2]
    print(f"n={n}: A000255={lhs}, recurrence gives {rhs} {'✓' if lhs==rhs else '✗'}")

# The ΣH in terms of A000255_eval(n, x=2):
print("\nΣH as weighted A000255 sum:")
for n in range(2, 9):
    bpd = succession_free_bp(n)
    m = comb(n-1,2)
    # ΣH = 2^{m-n+1} * Σ_bp N(n,bp)*2^bp = 2^{m-n+1} * A255_eval(n,2)
    a255_eval_at_2 = sum(cnt * (2**bp) for bp, cnt in bpd.items())
    sigma_H = (2**(m-n+1)) * a255_eval_at_2
    print(f"n={n}: 2^{{m-n+1={m-n+1}}} * A255_eval(n,2)={a255_eval_at_2} = ΣH={sigma_H}")

# ─────────────────────────────────────────────
# PART 5: THE SUMMAND GRAPH LAYERS AND METAGRAPH LEVELS
# ─────────────────────────────────────────────
print("\n=== SUMMAND GRAPH LAYERS vs H-SPECTRUM GROWTH ===")

def min_depth(n, memo={}):
    if n in memo: return memo[n]
    if n <= 2: memo[n]=0; return 0
    pairs = [(a, n-a) for a in range(1, n//2 + (1 if n%2==1 else 0))]
    if not pairs: memo[n]=float('inf'); return float('inf')
    d = 1 + min(max(min_depth(a), min_depth(b)) for a,b in pairs)
    memo[n] = d; return d

print("H-values by summand-graph depth:")
for n in range(3, 8):
    by_depth = defaultdict(list)
    for h in spectra[n]:
        by_depth[min_depth(h)].append(h)
    print(f"n={n}: spectrum by depth:")
    for d in sorted(by_depth):
        print(f"  depth {d}: {by_depth[d]}")

# ─────────────────────────────────────────────
# PART 6: THE {1,4,6} MODULE IN METAGRAPH TERMS
# ─────────────────────────────────────────────
print("\n=== THE {1,4,6} MODULE IN METAGRAPH TERMS ===")
print("H=1: transitive tournament (unique)")
print("H=4: even, never achieved (Rédei — H always odd)")
print("H=6: even, never achieved\n")
print("Reinterpretation: in the METAGRAPH, which H-values play the role of {1,4,6}?")
print("= which H-values are UNREACHABLE by arc-flips from the {2,3}-generated H-values?")
print()

# In the metagraph: the 'root' H-values are those achievable at n=3: {1,3}
# Starting from {1,3} and applying ΔH steps, which H-values are reachable?
# ΔH values form a set of EVEN integers

for n in range(3, 8):
    reachable = set(spectra[3])  # Start from n=3 values
    added = True
    while added:
        added = False
        for h in list(reachable):
            for delta in range(-200, 201, 2):  # Even steps only
                h_new = h + delta
                if h_new > 0 and h_new % 2 == 1 and h_new in spectra[n]:
                    if h_new not in reachable:
                        reachable.add(h_new)
                        added = True
    missing = [h for h in spectra[n] if h not in reachable]
    print(f"n={n}: Starting from {{1,3}}, reachable by even steps: {sorted(reachable)[:10]}...")
    print(f"       H-values NOT reachable from {{1,3}}: {missing}")

# ─────────────────────────────────────────────
# PART 7: THE RECURSIVE FORMULA FOR H_MAX(n)
# ─────────────────────────────────────────────
print("\n=== RECURSIVE STRUCTURE OF H_MAX(n) ===")
h_max = {n: max(spectra[n]) for n in spectra}
print("H_max values:", h_max)
print()
for n in range(4, 8):
    hm_n = h_max[n]
    hm_prev = h_max[n-1]
    ratio = hm_n / hm_prev
    # Check if H_max(n) = H_max(n-1) * (n) / some_constant
    print(f"n={n}: H_max={hm_n}, prev={hm_prev}, ratio={ratio:.4f}")
    # Look for H_max(n) as a summand-graph combination of H_max(n-1)
    for a in range(1, hm_n):
        b = hm_n - a
        if b > 0 and b < hm_n:
            # Is a reachable at n-1 or smaller?
            for prev_n in [n-2, n-1]:
                if prev_n in spectra and a in spectra[prev_n]:
                    print(f"  H_max({n})={hm_n} = {a} (from n={prev_n}) + {b}")
                    break
            break

# ─────────────────────────────────────────────
# PART 8: SUMMAND GRAPH DEPTH = TOURNAMENT HIERARCHY LEVEL
# ─────────────────────────────────────────────
print("\n=== DEPTH CORRESPONDENCE: SUMMAND GRAPH ↔ TOURNAMENT LEVEL ===")
print("Key observation: H=1 (depth 0), H=3 (depth 1), H=5 (depth 2)")
print("H=9,11,13,15 (depth 3-4), matching metagraph 'level' from transitive\n")

for n in range(3, 8):
    h_to_depth = {h: min_depth(h) for h in spectra[n]}
    # How many distinct depths appear?
    depths = sorted(set(h_to_depth.values()))
    print(f"n={n}: depths of achievable H-values = {depths}")
    for d in depths:
        vals = [h for h,dd in h_to_depth.items() if dd==d]
        print(f"  Depth {d}: {sorted(vals)}")

# ─────────────────────────────────────────────
# PART 9: ΣH GENERATING FUNCTION AND SUMMAND RECURRENCE
# ─────────────────────────────────────────────
print("\n=== ΣH GENERATING FUNCTION RECURRENCE ===")
# ΣH(n) satisfies: what recurrence?
sigma_H_vals = []
for n in range(2, 10):
    bpd = succession_free_bp(n) if n <= 9 else None
    if bpd is not None:
        m = comb(n-1,2)
        sh = sum(cnt * (2**(m-n+1+bp)) for bp,cnt in bpd.items())
        sigma_H_vals.append((n, sh))
        print(f"ΣH({n}) = {sh}")

# Look for linear recurrence
vals = [v for n,v in sigma_H_vals]
print("\nRatios ΣH(n)/ΣH(n-1):")
for i in range(1, len(vals)):
    print(f"  ΣH({i+2})/ΣH({i+1}) = {vals[i]}/{vals[i-1]} = {vals[i]/vals[i-1]:.4f}")

# Try recurrence ΣH(n) = c1*n*ΣH(n-1) + c2*ΣH(n-2)?
print("\nLooking for recurrence ΣH(n) = a(n)*ΣH(n-1) + b(n)*ΣH(n-2):")
for i in range(2, len(vals)):
    n = i+2
    # c1 * vals[i-1] + c2 * vals[i-2] = vals[i]
    # Two unknowns, try simple forms
    for a in range(1, 20):
        for b in range(-20, 20):
            if a * vals[i-1] + b * vals[i-2] == vals[i]:
                print(f"  ΣH({n}) = {a}*ΣH({n-1}) + {b}*ΣH({n-2}) ✓")

