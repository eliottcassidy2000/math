"""
Doubling/Blowup study for tournaments - final version.
Uses networkx for isomorphism, DP for HP counting.
"""

import itertools
from collections import defaultdict, Counter
import math
import networkx as nx

# ============================================================
# Canonical form (for small n<=7 only for enumeration)
# ============================================================

def canonical_form(n, arcs):
    arc_set = frozenset(arcs)
    best = None
    for perm in itertools.permutations(range(n)):
        relabeled = tuple(sorted((perm[u], perm[v]) for u,v in arc_set))
        if best is None or relabeled < best:
            best = relabeled
    return best

def is_SC(n, arcs):
    """Check if tournament is self-complementary using networkx."""
    G = nx.DiGraph()
    G.add_nodes_from(range(n))
    G.add_edges_from(arcs)
    Gop = nx.DiGraph()
    Gop.add_nodes_from(range(n))
    Gop.add_edges_from((v,u) for u,v in arcs)
    return nx.is_isomorphic(G, Gop)

def score_seq(n, arcs):
    out = [0]*n
    for u,v in arcs:
        out[u] += 1
    return tuple(sorted(out))

# ============================================================
# HP counting via DP
# ============================================================

def count_hp(n, arcs):
    arc_set = set(arcs)
    dp = [[0]*n for _ in range(1<<n)]
    for v in range(n):
        dp[1<<v][v] = 1
    for mask in range(1, 1<<n):
        for v in range(n):
            if not (mask & (1<<v)) or dp[mask][v] == 0:
                continue
            for w in range(n):
                if not (mask & (1<<w)) and (v,w) in arc_set:
                    dp[mask|(1<<w)][w] += dp[mask][v]
    full = (1<<n)-1
    return sum(dp[full][v] for v in range(n))

# ============================================================
# Blowups
# ============================================================

def lex_blowup(n, arcs):
    new_arcs = set()
    for v in range(n):
        new_arcs.add((2*v, 2*v+1))
    for u,v in arcs:
        for i in range(2):
            for j in range(2):
                new_arcs.add((2*u+i, 2*v+j))
    return 2*n, frozenset(new_arcs)

def sc_blowup(n, arcs):
    new_arcs = set()
    for v in range(n):
        new_arcs.add((2*v, 2*v+1))
    for u,v in arcs:
        new_arcs.add((2*u,   2*v))    # u0->v0 (lane)
        new_arcs.add((2*u+1, 2*v+1))  # u1->v1 (lane)
        new_arcs.add((2*v,   2*u+1))  # v0->u1 (cross)
        new_arcs.add((2*v+1, 2*u))    # v1->u0 (cross)
    return 2*n, frozenset(new_arcs)

# ============================================================
# Enumerate iso classes for n
# ============================================================

def all_iso_classes(n):
    pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
    seen = set()
    reps = []
    for choices in itertools.product([0,1], repeat=len(pairs)):
        arcs = frozenset((i,j) if c==0 else (j,i) for c,(i,j) in zip(choices,pairs))
        can = canonical_form(n, arcs)
        if can not in seen:
            seen.add(can)
            reps.append(arcs)
    return reps

# ============================================================
# Circulant
# ============================================================

def circulant(n, S):
    return frozenset((i, (i+s)%n) for i in range(n) for s in S)

def is_circulant(n, arcs):
    arc_set = frozenset(arcs)
    return frozenset(((u+1)%n, (v+1)%n) for u,v in arc_set) == arc_set

def connection_set(n, arcs):
    if not is_circulant(n, arcs):
        return None
    return tuple(sorted((v - 0) % n for u,v in arcs if u == 0))

# ============================================================
# PART 1: Main table for n=3,4,5
# ============================================================
print("="*90)
print("PART 1: H values for all iso classes at n=3,4,5")
print("="*90)

all_results = {}

for n in [3, 4, 5]:
    classes = all_iso_classes(n)
    print(f"\n--- n={n}: {len(classes)} iso classes (blowups -> n={2*n}) ---")
    hdr = f"{'Cls':<5} {'Score seq':<18} {'H(T)':<7} {'H(Lex)':<8} {'H(SC)':<9} {'T SC?':<7} {'Lex SC?':<8} {'SC SC?':<8}"
    print(hdr)
    print("-"*72)

    all_results[n] = []

    for idx, arcs in enumerate(classes):
        ss = score_seq(n, arcs)
        H_T = count_hp(n, arcs)

        Nl, lex_arcs = lex_blowup(n, arcs)
        H_lex = count_hp(Nl, lex_arcs)

        Ns, sc_arcs = sc_blowup(n, arcs)
        H_sc = count_hp(Ns, sc_arcs)

        t_sc   = is_SC(n, arcs)
        lex_sc = is_SC(Nl, lex_arcs)
        sc_sc  = is_SC(Ns, sc_arcs)

        print(f"T{idx:<4} {str(ss):<18} {H_T:<7} {H_lex:<8} {H_sc:<9} {str(t_sc):<7} {str(lex_sc):<8} {str(sc_sc)}")

        all_results[n].append({
            'idx': idx, 'arcs': arcs, 'ss': ss,
            'H_T': H_T, 'H_lex': H_lex, 'H_sc': H_sc,
            't_sc': t_sc, 'lex_sc': lex_sc, 'sc_sc': sc_sc
        })

# ============================================================
# PART 2: SC blowup of SC tournaments
# ============================================================
print("\n" + "="*90)
print("PART 2: SC blowup of SC tournaments")
print("="*90)
for n in [3,4,5]:
    sc_ts = [r for r in all_results[n] if r['t_sc']]
    print(f"\nn={n}: {len(sc_ts)} SC tournament(s)")
    for r in sc_ts:
        print(f"  T{r['idx']} ss={r['ss']}: H(T)={r['H_T']}, H(T_SC)={r['H_sc']}, T_SC is SC: {r['sc_sc']}")

# ============================================================
# PART 3: Circulant analysis
# ============================================================
print("\n" + "="*90)
print("PART 3: Circulant tournaments")
print("="*90)

# C3
c3 = circulant(3, {1})
Nc3, c3sc = sc_blowup(3, c3)
H_c3 = count_hp(3, c3)
H_c3sc = count_hp(6, c3sc)
Nc3l, c3lex = lex_blowup(3, c3)
H_c3lex = count_hp(6, c3lex)
print(f"\nC3 = circ(3,{{1}}): H={H_c3}, SC={is_SC(3,c3)}")
print(f"  C3_Lex (n=6): H={H_c3lex}, circ={is_circulant(6,c3lex)}, SC={is_SC(6,c3lex)}")
print(f"  C3_SC  (n=6): H={H_c3sc},  circ={is_circulant(6,c3sc)}, SC={is_SC(6,c3sc)}")
if is_circulant(6, c3sc):
    print(f"    conn set: {connection_set(6,c3sc)}")

# Paley(5)
c5 = circulant(5, {1,2})
Nc5, c5sc = sc_blowup(5, c5)
H_c5 = count_hp(5, c5)
H_c5sc = count_hp(10, c5sc)
print(f"\nPaley(5) = circ(5,{{1,2}}): H={H_c5}, SC={is_SC(5,c5)}")
print(f"  Paley(5)_SC (n=10): H={H_c5sc}, circ={is_circulant(10,c5sc)}, SC={is_SC(10,c5sc)}")
if is_circulant(10, c5sc):
    print(f"    conn set: {connection_set(10,c5sc)}")

# All circulant n=5 tournaments
print("\n  All circulant tournaments on n=5:")
for S in [{1,2},{1,3},{2,3},{1,4},{2,4},{3,4}]:
    ct = circulant(5, S)
    if is_circulant(5, ct):  # verify
        H = count_hp(5, ct)
        Nct, ctsc = sc_blowup(5, ct)
        Hsc = count_hp(10, ctsc)
        circ_sc = is_circulant(10, ctsc)
        isc = is_SC(10, ctsc)
        cs = connection_set(10, ctsc) if circ_sc else 'N/A'
        print(f"  circ(5,{sorted(S)}): H={H} -> SC blowup: H={Hsc}, circ={circ_sc}, SC={isc}, conn={cs}")

# Paley(7)
c7 = circulant(7, {1,2,4})
Nc7, c7sc = sc_blowup(7, c7)
H_c7 = count_hp(7, c7)
H_c7sc = count_hp(14, c7sc)
print(f"\nPaley(7) = circ(7,{{1,2,4}}): H={H_c7}, SC={is_SC(7,c7)}")
print(f"  Paley(7)_SC (n=14): H={H_c7sc}, circ={is_circulant(14,c7sc)}, SC={is_SC(14,c7sc)}")
if is_circulant(14, c7sc):
    cs = connection_set(14, c7sc)
    print(f"    conn set: {cs}")

# ============================================================
# PART 4: Tower computation
# ============================================================
print("\n" + "="*90)
print("PART 4: Tower computation")
print("="*90)

# C3 tower
print("\nC3 tower (SC blowup):")
H0 = count_hp(3, c3)
N1, t1 = sc_blowup(3, c3)
H1 = count_hp(N1, t1)
N2, t2 = sc_blowup(N1, t1)
H2 = count_hp(N2, t2)
print(f"  n=3:  H = {H0}")
print(f"  n=6:  H = {H1}")
print(f"  n=12: H = {H2}")
print(f"  Ratios: {H1}/{H0} = {H1/H0:.4f}, {H2}/{H1} = {H2/H1:.4f}")
print(f"  H1/H0^2 = {H1}/{H0**2} = {H1/H0**2:.6f}")
print(f"  H2/H1^2 = {H2}/{H1**2} = {H2/H1**2:.6f}")

# Transitive n=3 tower
print("\nTransitive-3 tower (SC blowup):")
trans3 = frozenset([(0,1),(0,2),(1,2)])
Ht0 = count_hp(3, trans3)
Nt1, tt1 = sc_blowup(3, trans3)
Ht1 = count_hp(Nt1, tt1)
Nt2, tt2 = sc_blowup(Nt1, tt1)
Ht2 = count_hp(Nt2, tt2)
print(f"  n=3:  H = {Ht0}")
print(f"  n=6:  H = {Ht1}")
print(f"  n=12: H = {Ht2}")
print(f"  Ratios: {Ht1}/{Ht0} = {Ht1/Ht0:.4f}, {Ht2}/{Ht1} = {Ht2/Ht1:.4f}")

# Paley(5) tower
print("\nPaley(5) tower (SC blowup):")
Hp0 = count_hp(5, c5)
Np1, pt1 = sc_blowup(5, c5)
Hp1 = count_hp(Np1, pt1)
Np2, pt2 = sc_blowup(Np1, pt1)
Hp2 = count_hp(Np2, pt2)
print(f"  n=5:  H = {Hp0}")
print(f"  n=10: H = {Hp1}")
print(f"  n=20: H = {Hp2}")
print(f"  Ratios: {Hp1}/{Hp0} = {Hp1/Hp0:.4f}, {Hp2}/{Hp1} = {Hp2/Hp1:.4f}")
print(f"  H1/H0^2 = {Hp1/Hp0**2:.6f}")
print(f"  H2/H1^2 = {Hp2/Hp1**2:.6f}")

# ============================================================
# PART 5: Coarse path structure (n=3 and n=4 cases only)
# ============================================================
print("\n" + "="*90)
print("PART 5: Coarse path structure in SC blowup")
print("="*90)

def coarse_analysis(n, base_arcs, sc_arcs, label, max_n=10):
    N = 2*n
    if N > max_n:
        # Use DP to enumerate paths and extract coarse pattern
        arc_set = set(sc_arcs)
        base_arc_set = set(base_arcs)

        # DP tracking coarse sequence -- too complex, skip
        print(f"\n{label}: N={N} too large for brute force coarse analysis")
        return None

    arc_set = set(sc_arcs)
    base_arc_set = set(base_arcs)

    coarse_counts = defaultdict(int)
    for perm in itertools.permutations(range(N)):
        if all((perm[i], perm[i+1]) in arc_set for i in range(N-1)):
            coarse = tuple(v//2 for v in perm)
            coarse_counts[coarse] += 1

    total = sum(coarse_counts.values())

    t_hps = set()
    for perm in itertools.permutations(range(n)):
        if all((perm[i],perm[i+1]) in base_arc_set for i in range(n-1)):
            t_hps.add(perm)

    print(f"\n{label}: total HP = {total}, |HP(T)| = {len(t_hps)}")
    cnt_dist = Counter(coarse_counts.values())
    print(f"  Distinct coarse seqs: {len(coarse_counts)}, count distribution: {dict(sorted(cnt_dist.items()))}")

    hp_coarse = {c: coarse_counts[c] for c in coarse_counts if c in t_hps}
    non_hp_coarse = {c: coarse_counts[c] for c in coarse_counts if c not in t_hps}
    print(f"  Coarse seqs = HP in T: {len(hp_coarse)}, total paths: {sum(hp_coarse.values())}")
    print(f"  Coarse seqs != HP in T: {len(non_hp_coarse)}, total paths: {sum(non_hp_coarse.values())}")
    print(f"  All coarse patterns:")
    for coarse, cnt in sorted(coarse_counts.items()):
        tag = '[HP in T]' if coarse in t_hps else ''
        print(f"    {coarse}: {cnt} {tag}")

    return coarse_counts

# C3 (n=3 -> N=6)
coarse_analysis(3, c3, c3sc, "C3_SC (n=6)")

# Transitive-3 (n=3 -> N=6)
_, trans3sc = sc_blowup(3, trans3)
coarse_analysis(3, trans3, trans3sc, "Trans3_SC (n=6)")

# All n=4 classes (N=8, brute force feasible)
print("\n--- n=4 cases ---")
for r in all_results[4]:
    _, sc4 = sc_blowup(4, r['arcs'])
    coarse_analysis(4, r['arcs'], sc4, f"T{r['idx']}(ss={r['ss']})_SC", max_n=8)

# ============================================================
# PART 6: Formula search
# ============================================================
print("\n" + "="*90)
print("PART 6: Formula search")
print("="*90)

print("\nAll data:")
print(f"{'n':<4} {'T':<5} {'ss':<22} {'H_T':<8} {'H_Lex':<9} {'H_SC':<9} {'H_SC/H_T^2':<14} {'H_Lex/H_T^2':<14}")
for n in [3,4,5]:
    for r in all_results[n]:
        r1 = r['H_sc'] / r['H_T']**2 if r['H_T'] > 0 else float('nan')
        r2 = r['H_lex'] / r['H_T']**2 if r['H_T'] > 0 else float('nan')
        print(f"{n:<4} T{r['idx']:<4} {str(r['ss']):<22} {r['H_T']:<8} {r['H_lex']:<9} {r['H_sc']:<9} {r1:<14.4f} {r2:<14.4f}")

# Key check: is H(T_SC) the same for all tournaments with same H(T)?
print("\n--- H_SC grouped by H_T ---")
ht_to_hsc = defaultdict(list)
for n in [3,4,5]:
    for r in all_results[n]:
        ht_to_hsc[(n,r['H_T'])].append((r['H_sc'], r['ss'], r['idx']))

for key in sorted(ht_to_hsc.keys()):
    n, H_T = key
    entries = ht_to_hsc[key]
    hsc_vals = [e[0] for e in entries]
    print(f"  n={n}, H_T={H_T}: H_SC values = {hsc_vals} -- {'SAME' if len(set(hsc_vals))==1 else 'DIFFERENT'}")

# Key check: is H(Lex) the same for all tournaments with same H(T)?
print("\n--- H_Lex grouped by H_T ---")
ht_to_hlex = defaultdict(list)
for n in [3,4,5]:
    for r in all_results[n]:
        ht_to_hlex[(n,r['H_T'])].append((r['H_lex'], r['ss'], r['idx']))

for key in sorted(ht_to_hlex.keys()):
    n, H_T = key
    entries = ht_to_hlex[key]
    hlex_vals = [e[0] for e in entries]
    print(f"  n={n}, H_T={H_T}: H_Lex values = {hlex_vals} -- {'SAME' if len(set(hlex_vals))==1 else 'DIFFERENT'}")

print("\n--- H_SC as function of (n, score_seq) ---")
for n in [3,4,5]:
    by_ss = defaultdict(list)
    for r in all_results[n]:
        by_ss[r['ss']].append((r['H_T'], r['H_sc'], r['idx']))
    for ss, vals in sorted(by_ss.items()):
        print(f"  n={n}, ss={ss}: {vals}")

print("\n--- Looking for H_Lex = (2n-1)!! * H_T? ---")
for n in [3,4,5]:
    dbl_fact = math.prod(range(1, 2*n, 2))  # (2n-1)!! = 1*3*5*...*(2n-1)
    print(f"  n={n}: (2n-1)!! = (2*{n}-1)!! = {dbl_fact}")
    for r in all_results[n]:
        if r['H_T'] > 0:
            print(f"    T{r['idx']}: H_Lex={r['H_lex']}, (2n-1)!!*H_T = {dbl_fact*r['H_T']}, ratio = {r['H_lex']/(dbl_fact*r['H_T']):.4f}")

print("\nDone!")
