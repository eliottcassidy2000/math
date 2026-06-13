"""
Doubling/Blowup study for tournaments.
Fixed canonical form using sorted-tuple comparison.
"""

import itertools
from collections import defaultdict, Counter
import math

# ============================================================
# Canonical form: min over relabelings using sorted-tuple
# ============================================================

def canonical_form(n, arcs):
    """Canonical form: min over all relabelings (lex on sorted arc list)."""
    arc_set = frozenset(arcs)
    best = None
    for perm in itertools.permutations(range(n)):
        relabeled = tuple(sorted((perm[u], perm[v]) for u,v in arc_set))
        if best is None or relabeled < best:
            best = relabeled
    return best

def is_isomorphic(n, a1, a2):
    return canonical_form(n, a1) == canonical_form(n, a2)

def is_SC(n, arcs):
    rev = frozenset((v,u) for u,v in arcs)
    return canonical_form(n, arcs) == canonical_form(n, rev)

# ============================================================
# Score sequence
# ============================================================

def score_seq(n, arcs):
    out = [0]*n
    for u,v in arcs:
        out[u] += 1
    return tuple(sorted(out))

# ============================================================
# Hamiltonian path counting via DP
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
# Blowup constructions
# ============================================================

def lex_blowup(n, arcs):
    """T[K2]: u->v gives all 4 arcs ui->vj; v0->v1 internal."""
    new_arcs = set()
    for v in range(n):
        new_arcs.add((2*v, 2*v+1))
    for u,v in arcs:
        for i in range(2):
            for j in range(2):
                new_arcs.add((2*u+i, 2*v+j))
    return 2*n, frozenset(new_arcs)

def sc_blowup(n, arcs):
    """T_SC: u->v gives lane arcs u0->v0, u1->v1 AND cross arcs v0->u1, v1->u0."""
    new_arcs = set()
    for v in range(n):
        new_arcs.add((2*v, 2*v+1))
    for u,v in arcs:
        new_arcs.add((2*u,   2*v))      # u0->v0
        new_arcs.add((2*u+1, 2*v+1))    # u1->v1
        new_arcs.add((2*v,   2*u+1))    # v0->u1
        new_arcs.add((2*v+1, 2*u))      # v1->u0
    return 2*n, frozenset(new_arcs)

# ============================================================
# Enumerate iso classes
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
# Circulant utilities
# ============================================================

def circulant(n, S):
    return frozenset((i, (i+s)%n) for i in range(n) for s in S)

def is_circulant(n, arcs):
    arc_set = frozenset(arcs)
    rotated = frozenset(((u+1)%n, (v+1)%n) for u,v in arc_set)
    return rotated == arc_set

def connection_set(n, arcs):
    if not is_circulant(n, arcs):
        return None
    return tuple(sorted((v - 0) % n for u,v in arcs if u == 0))

# ============================================================
# PART 1: Main table for n=3,4,5
# ============================================================
print("="*80)
print("PART 1: H values for all iso classes at n=3,4,5")
print("="*80)

all_results = {}

for n in [3, 4, 5]:
    print(f"\n{'='*50}")
    print(f"n={n} ({2**( n*(n-1)//2)} tournaments, blowups have n={2*n})")
    hdr = f"{'Cls':<5} {'Score seq':<16} {'H(T)':<7} {'H(Lex)':<8} {'H(SC)':<8} {'T SC?':<7} {'Lex SC?':<8} {'SC SC?':<8} {'H(SC)/H(T)^2':<14}"
    print(hdr)
    print("-"*85)

    classes = all_iso_classes(n)
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

        ratio = H_sc / H_T**2 if H_T > 0 else float('nan')

        print(f"T{idx:<4} {str(ss):<16} {H_T:<7} {H_lex:<8} {H_sc:<8} {str(t_sc):<7} {str(lex_sc):<8} {str(sc_sc):<8} {ratio:.4f}")

        all_results[n].append({
            'idx': idx, 'arcs': arcs, 'ss': ss,
            'H_T': H_T, 'H_lex': H_lex, 'H_sc': H_sc,
            't_sc': t_sc, 'lex_sc': lex_sc, 'sc_sc': sc_sc
        })

# ============================================================
# PART 2: Iso-class consolidation check
# ============================================================
print("\n" + "="*80)
print("PART 2: Number of iso classes and H values summary")
print("="*80)
for n in [3,4,5]:
    classes = all_results[n]
    print(f"\nn={n}: {len(classes)} iso classes")
    # Group by score seq
    by_ss = defaultdict(list)
    for r in classes:
        by_ss[r['ss']].append(r)
    for ss, rs in sorted(by_ss.items()):
        hs = [r['H_T'] for r in rs]
        print(f"  ss={ss}: {len(rs)} class(es), H(T) values: {hs}")

# ============================================================
# PART 3: Circulant analysis
# ============================================================
print("\n" + "="*80)
print("PART 3: Circulant tournaments")
print("="*80)

# C3 (the unique cyclic 3-tournament)
c3 = circulant(3, {1})
Nc3, c3sc = sc_blowup(3, c3)
H_c3 = count_hp(3, c3)
H_c3sc = count_hp(6, c3sc)
print(f"\nC3 = circ(3,{{1}}): H={H_c3}, SC={is_SC(3,c3)}, score={score_seq(3,c3)}")
print(f"C3_SC (n=6): H={H_c3sc}, is_circulant={is_circulant(6,c3sc)}, SC={is_SC(6,c3sc)}")
if is_circulant(6, c3sc):
    cs = connection_set(6, c3sc)
    print(f"  C3_SC connection set: {cs}")
print(f"  score: {score_seq(6, c3sc)}")

# Lex blowup of C3
Nlc3, c3lex = lex_blowup(3, c3)
H_c3lex = count_hp(6, c3lex)
print(f"C3_Lex (n=6): H={H_c3lex}, is_circulant={is_circulant(6,c3lex)}, SC={is_SC(6,c3lex)}")

# Paley(5)
c5 = circulant(5, {1,2})
Nc5, c5sc = sc_blowup(5, c5)
H_c5 = count_hp(5, c5)
H_c5sc = count_hp(10, c5sc)
print(f"\nPaley(5) = circ(5,{{1,2}}): H={H_c5}, SC={is_SC(5,c5)}, score={score_seq(5,c5)}")
print(f"Paley(5)_SC (n=10): H={H_c5sc}, is_circulant={is_circulant(10,c5sc)}, SC={is_SC(10,c5sc)}")
if is_circulant(10, c5sc):
    cs = connection_set(10, c5sc)
    print(f"  Paley(5)_SC connection set: {cs}")

# C5 with different connection set
for S in [{1,2},{1,3},{2,3}]:
    ct = circulant(5, S)
    if is_circulant(5, ct):
        H = count_hp(5, ct)
        Nct, ctsc = sc_blowup(5, ct)
        Hsc = count_hp(10, ctsc)
        print(f"circ(5,{sorted(S)}): H={H}, SC={is_SC(5,ct)} -> SC blowup: H={Hsc}, circ={is_circulant(10,ctsc)}, SC={is_SC(10,ctsc)}")

# Paley(7)
c7 = circulant(7, {1,2,4})
Nc7, c7sc = sc_blowup(7, c7)
H_c7 = count_hp(7, c7)
H_c7sc = count_hp(14, c7sc)
print(f"\nPaley(7) = circ(7,{{1,2,4}}): H={H_c7}, SC={is_SC(7,c7)}, score={score_seq(7,c7)}")
print(f"Paley(7)_SC (n=14): H={H_c7sc}, is_circulant={is_circulant(14,c7sc)}, SC={is_SC(14,c7sc)}")
if is_circulant(14, c7sc):
    cs = connection_set(14, c7sc)
    print(f"  Paley(7)_SC connection set: {cs}")

# ============================================================
# PART 4: Tower computation
# ============================================================
print("\n" + "="*80)
print("PART 4: Tower C3 -> C3_SC -> (C3_SC)_SC")
print("="*80)

# Level 0
H0 = count_hp(3, c3)

# Level 1
H1 = count_hp(6, c3sc)

# Level 2
N12, c3sc2 = sc_blowup(6, c3sc)
H2 = count_hp(12, c3sc2)

print(f"Level 0 (n=3):  H = {H0}")
print(f"Level 1 (n=6):  H = {H1}")
print(f"Level 2 (n=12): H = {H2}")
print(f"Ratios: H1/H0 = {H1/H0:.4f}, H2/H1 = {H2/H1:.4f}")
print(f"H1/(H0^2) = {H1/H0**2:.4f}")
print(f"H2/(H1^2) = {H2/H1**2:.4f}")
print(f"H1/H0^2 exact: {H1}/{H0**2} = {H1}/ {H0**2}")
# n grows as 3,6,12,... = 3*2^k
# H: 3, ?, ?
print(f"Sequence: {H0}, {H1}, {H2}")

# Tower for transitive n=3
trans3 = frozenset([(0,1),(0,2),(1,2)])
H_t0 = count_hp(3, trans3)
N, trans3sc = sc_blowup(3, trans3)
H_t1 = count_hp(6, trans3sc)
N2, trans3sc2 = sc_blowup(6, trans3sc)
H_t2 = count_hp(12, trans3sc2)
print(f"\nTower for Transitive-3:")
print(f"Level 0 (n=3):  H = {H_t0}")
print(f"Level 1 (n=6):  H = {H_t1}")
print(f"Level 2 (n=12): H = {H_t2}")
print(f"Ratios: H1/H0 = {H_t1/H_t0:.4f}, H2/H1 = {H_t2/H_t1:.4f}")

# ============================================================
# PART 5: Coarse path analysis
# ============================================================
print("\n" + "="*80)
print("PART 5: Coarse path analysis")
print("="*80)

def coarse_analysis(n, arcs, sc_arcs, label):
    N = 2*n
    arc_set = set(sc_arcs)
    base_arcs = set(arcs)

    coarse_counts = defaultdict(int)

    for perm in itertools.permutations(range(N)):
        if all((perm[i], perm[i+1]) in arc_set for i in range(N-1)):
            coarse = tuple(v//2 for v in perm)
            coarse_counts[coarse] += 1

    total = sum(coarse_counts.values())

    # Hamiltonian paths in T
    t_hps = set()
    for perm in itertools.permutations(range(n)):
        if all((perm[i],perm[i+1]) in base_arcs for i in range(n-1)):
            t_hps.add(perm)

    print(f"\n{label}: total HP = {total}, |HP(T)| = {len(t_hps)}")
    print(f"  Distinct coarse seqs: {len(coarse_counts)}")

    cnt_dist = Counter(coarse_counts.values())
    print(f"  Count distribution: {dict(sorted(cnt_dist.items()))}")

    # How many coarse seqs correspond to HPs in T?
    hp_coarse = sum(1 for c in coarse_counts if c in t_hps)
    non_hp_coarse = sum(1 for c in coarse_counts if c not in t_hps)
    print(f"  Coarse seqs that ARE HPs in T: {hp_coarse}")
    print(f"  Coarse seqs that are NOT HPs in T: {non_hp_coarse}")

    if hp_coarse > 0:
        hp_total = sum(coarse_counts[c] for c in coarse_counts if c in t_hps)
        non_hp_total = sum(coarse_counts[c] for c in coarse_counts if c not in t_hps)
        print(f"  HP paths via 'T-coarse' patterns: {hp_total}")
        print(f"  HP paths via 'non-T-coarse' patterns: {non_hp_total}")

    print(f"  All coarse patterns:")
    for coarse, cnt in sorted(coarse_counts.items()):
        is_hp = '[HP in T]' if coarse in t_hps else ''
        print(f"    {coarse}: {cnt} paths {is_hp}")

    return coarse_counts

# C3 case (n=6, feasible)
cp_c3 = coarse_analysis(3, c3, c3sc, "C3_SC (n=6)")

# Transitive-3 case
cp_trans3 = coarse_analysis(3, trans3, trans3sc, "Trans3_SC (n=6)")

# All n=4 cases (n=8 blowup, 8 vertices, 8! paths to check is fine with DP)
# But brute force 8! = 40320 is also fine
for r in all_results[4]:
    _, sc4 = sc_blowup(4, r['arcs'])
    coarse_analysis(4, r['arcs'], sc4, f"T{r['idx']}(n=4,ss={r['ss']})_SC (n=8)")

# ============================================================
# PART 6: Formula search
# ============================================================
print("\n" + "="*80)
print("PART 6: Formula search")
print("="*80)

print("\nAll data (n, class, H_T, H_Lex, H_SC):")
for n in [3,4,5]:
    for r in all_results[n]:
        print(f"  n={n} T{r['idx']} ss={r['ss']}: H_T={r['H_T']}, H_Lex={r['H_lex']}, H_SC={r['H_sc']}")

# H_Lex pattern
print("\nLex blowup analysis:")
print("  For transitive (H=1): H_Lex should relate to (2n)!/(2n) type")
# For transitive n-tournament: H=1 (unique HP), the lex blowup is transitive 2n-tournament, H=1
# For cyclic C3: H=3, H_Lex=45

print("\nH_SC / H_T ratios:")
for n in [3,4,5]:
    for r in all_results[n]:
        if r['H_T'] > 0:
            ratio = r['H_sc'] / r['H_T']
            print(f"  n={n} T{r['idx']} ss={r['ss']}: {r['H_sc']}/{r['H_T']} = {ratio:.4f}")

# Is H_SC constant within iso class of same score seq but different H?
print("\nWithin n=4, different classes with same score seq:")
by_ss4 = defaultdict(list)
for r in all_results[4]:
    by_ss4[r['ss']].append(r)
for ss, rs in sorted(by_ss4.items()):
    if len(rs) > 1:
        print(f"  ss={ss}: {[(r['H_T'], r['H_lex'], r['H_sc']) for r in rs]}")

# ============================================================
# PART 7: SC blowup of SC tournaments
# ============================================================
print("\n" + "="*80)
print("PART 7: SC blowup of SC tournaments")
print("="*80)
for n in [3,4,5]:
    sc_ts = [r for r in all_results[n] if r['t_sc']]
    print(f"\nn={n}: SC tournaments: {len(sc_ts)}")
    for r in sc_ts:
        print(f"  T{r['idx']} ss={r['ss']}: H(T)={r['H_T']}, H(T_SC)={r['H_sc']}, T_SC is SC: {r['sc_sc']}")

print("\nDone!")
