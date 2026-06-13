"""
Doubling/Blowup study for tournaments.
Uses nauty (via pynauty or networkx) or hash-based canonical form.
Fast canonical form using trace of powers of adjacency matrix.
"""

import itertools
from collections import defaultdict, Counter
import subprocess
import sys

# ============================================================
# Fast canonical form via nauty command line
# ============================================================

def tournament_to_dreadnaut(n, arcs):
    """Convert tournament to dreadnaut input for nauty."""
    lines = []
    lines.append(f"n={n} g")
    for v in range(n):
        nbrs = [str(w) for w in range(n) if (v,w) in arcs]
        if nbrs:
            lines.append(f"{v}: {' '.join(nbrs)};")
        else:
            lines.append(f"{v};")
    lines.append(".")
    lines.append("cx")
    lines.append("##")
    return "\n".join(lines)

def canonical_via_nauty(n, arcs):
    """Use nauty/dreadnaut to get canonical form."""
    arc_set = frozenset(arcs)
    inp = tournament_to_dreadnaut(n, arc_set)
    try:
        result = subprocess.run(
            ['dreadnaut'],
            input=inp, capture_output=True, text=True, timeout=10
        )
        return result.stdout.strip()
    except Exception:
        return None

# ============================================================
# Canonical form: use sorted adjacency rows (good for small n)
# For correctness we need actual isomorphism; use fast hash.
# ============================================================

def adj_matrix(n, arcs):
    A = [[0]*n for _ in range(n)]
    for u,v in arcs:
        A[u][v] = 1
    return A

def score_seq(n, arcs):
    A = adj_matrix(n, arcs)
    return tuple(sorted(sum(A[i]) for i in range(n)))

def canonical_form_small(n, arcs):
    """Canonical form by min over all relabelings. Only for n<=8."""
    arc_set = frozenset(arcs)
    best = None
    for perm in itertools.permutations(range(n)):
        relabeled = frozenset((perm[u], perm[v]) for u,v in arc_set)
        if best is None or relabeled < best:
            best = relabeled
    return best

def canonical_form(n, arcs):
    if n <= 8:
        return canonical_form_small(n, arcs)
    # For larger n, use nauty if available
    c = canonical_via_nauty(n, arcs)
    if c is not None:
        return c
    # Fallback: use score-based invariant (not perfect but fast)
    return frozenset(arcs)

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
# SC check (isomorphism with reverse)
# ============================================================

def reverse_arcs(arcs):
    return frozenset((v,u) for u,v in arcs)

def is_SC(n, arcs):
    rev = reverse_arcs(arcs)
    return canonical_form(n, arcs) == canonical_form(n, rev)

# ============================================================
# Blowup constructions
# ============================================================

def lex_blowup(n, arcs):
    """T[K2]: u->v gives all 4 arcs ui->vj; v0->v1 internal."""
    arc_set = set(arcs)
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
        new_arcs.add((2*v,   2*u+1))    # v0->u1 (cross, T^op)
        new_arcs.add((2*v+1, 2*u))      # v1->u0 (cross, T^op)
    return 2*n, frozenset(new_arcs)

# ============================================================
# Enumerate all iso classes for n
# ============================================================

def all_tournaments_isoclass(n):
    pairs = [(i,j) for i in range(n) for j in range(i+1,n)]
    seen = set()
    reps = []
    for choices in itertools.product([0,1], repeat=len(pairs)):
        arcs = set()
        for k,(i,j) in enumerate(pairs):
            arcs.add((i,j) if choices[k]==0 else (j,i))
        arcs = frozenset(arcs)
        can = canonical_form(n, arcs)
        if can not in seen:
            seen.add(can)
            reps.append(arcs)
    return reps

# ============================================================
# PART 1: H values for n=3,4,5
# ============================================================
print("="*70)
print("PART 1: H values for all iso classes at n=3,4,5")
print("="*70)

all_results = {}

for n in [3, 4, 5]:
    print(f"\n{'='*40}")
    print(f"n={n} (blowups have n={2*n} vertices)")
    print(f"{'Class':<6} {'Score seq':<20} {'H(T)':<8} {'H(Lex)':<8} {'H(SC)':<8} {'T SC?':<7} {'Lex SC?':<8} {'TscSC?':<8} {'H(SC)/H(T)^2':<14}")
    print("-"*90)

    classes = all_tournaments_isoclass(n)
    all_results[n] = []

    for idx, arcs in enumerate(classes):
        ss = score_seq(n, arcs)
        H_T = count_hp(n, arcs)

        Nl, lex_arcs = lex_blowup(n, arcs)
        H_lex = count_hp(Nl, lex_arcs)

        Ns, sc_arcs = sc_blowup(n, arcs)
        H_sc = count_hp(Ns, sc_arcs)

        t_sc  = is_SC(n, arcs)
        lex_sc = is_SC(Nl, lex_arcs)
        sc_sc  = is_SC(Ns, sc_arcs)

        ratio = H_sc / H_T**2 if H_T > 0 else float('nan')

        print(f"T{idx:<5} {str(ss):<20} {H_T:<8} {H_lex:<8} {H_sc:<8} {str(t_sc):<7} {str(lex_sc):<8} {str(sc_sc):<8} {ratio:<14.4f}")

        all_results[n].append({
            'idx': idx, 'arcs': arcs, 'ss': ss,
            'H_T': H_T, 'H_lex': H_lex, 'H_sc': H_sc,
            't_sc': t_sc, 'lex_sc': lex_sc, 'sc_sc': sc_sc
        })

# ============================================================
# PART 2: Circulant tournaments
# ============================================================
print("\n" + "="*70)
print("PART 2: Circulant tournaments")
print("="*70)

def circulant(n, S):
    return frozenset((i, (i+s)%n) for i in range(n) for s in S)

def is_circulant(n, arcs):
    rotated = frozenset(((u+1)%n, (v+1)%n) for u,v in arcs)
    return rotated == frozenset(arcs)

def find_connection_set(n, arcs):
    """If circulant, find S."""
    if not is_circulant(n, arcs):
        return None
    S = {(v-0)%n for u,v in arcs if u==0}
    return tuple(sorted(S))

# C3
c3 = circulant(3, {1})
N, c3sc = sc_blowup(3, c3)
print(f"\nC3 = circ(3,{{1}}): H={count_hp(3,c3)}, SC={is_SC(3,c3)}")
print(f"C3_SC (n=6): H={count_hp(6,c3sc)}, is_circulant={is_circulant(6,c3sc)}, SC={is_SC(6,c3sc)}")
if is_circulant(6, c3sc):
    print(f"  connection set: {find_connection_set(6,c3sc)}")
print(f"  score seq: {score_seq(6,c3sc)}")

# C5(1,2) = Paley(5)
c5 = circulant(5, {1,2})
N5, c5sc = sc_blowup(5, c5)
print(f"\nC5({{1,2}}) = Paley(5): H={count_hp(5,c5)}, SC={is_SC(5,c5)}")
print(f"C5_SC (n=10): H={count_hp(10,c5sc)}, is_circulant={is_circulant(10,c5sc)}, SC={is_SC(10,c5sc)}")
if is_circulant(10, c5sc):
    print(f"  connection set: {find_connection_set(10,c5sc)}")
print(f"  score seq: {score_seq(10,c5sc)}")

# C7(1,2,4) = Paley(7)
c7 = circulant(7, {1,2,4})
N7, c7sc = sc_blowup(7, c7)
H_c7 = count_hp(7, c7)
H_c7sc = count_hp(14, c7sc)
print(f"\nC7({{1,2,4}}) = Paley(7): H={H_c7}, SC={is_SC(7,c7)}")
print(f"C7_SC (n=14): H={H_c7sc}, is_circulant={is_circulant(14,c7sc)}, SC={is_SC(14,c7sc)}")
if is_circulant(14, c7sc):
    print(f"  connection set: {find_connection_set(14,c7sc)}")

# ============================================================
# PART 3: Tower computation C3 -> C3_SC -> (C3_SC)_SC
# ============================================================
print("\n" + "="*70)
print("PART 3: Tower computation")
print("="*70)

H0 = count_hp(3, c3)
H1 = count_hp(6, c3sc)
N12, c3sc2 = sc_blowup(6, c3sc)
H2 = count_hp(12, c3sc2)

print(f"C3 (n=3): H = {H0}")
print(f"C3_SC (n=6): H = {H1}")
print(f"(C3_SC)_SC (n=12): H = {H2}")
print(f"Ratios: H1/H0 = {H1/H0:.4f}, H2/H1 = {H2/H1:.4f}")
print(f"H1/H0^2 = {H1/H0**2:.4f}")
print(f"H2/H1^2 = {H2/H1**2:.4f}")

# ============================================================
# PART 4: Coarse path structure in SC blowup
# ============================================================
print("\n" + "="*70)
print("PART 4: Coarse path structure in SC blowup (C3 case)")
print("="*70)

def coarse_path_analysis(n, arcs, sc_arcs, label=""):
    N = 2*n
    arc_set = set(sc_arcs)
    base_arcs = set(arcs)

    coarse_counts = defaultdict(int)

    # DP-based enumeration with path tracking
    # dp[mask][v] = list of partial paths or just count
    # Too expensive for paths, use brute force for n<=5
    if N > 10:
        print(f"  Skipping {label}, N={N} too large for brute force")
        return

    for perm in itertools.permutations(range(N)):
        if all((perm[i], perm[i+1]) in arc_set for i in range(N-1)):
            # Extract coarse pattern (sequence of pairs)
            coarse = tuple(v//2 for v in perm)
            coarse_counts[coarse] += 1

    total = sum(coarse_counts.values())
    print(f"\n{label}: total HP = {total}")

    # Group by coarse sequence
    unique_coarses = list(coarse_counts.keys())
    print(f"Distinct coarse sequences: {len(unique_coarses)}")

    # Check which coarse sequences are themselves HPs in T
    t_hps = set()
    for perm in itertools.permutations(range(n)):
        if all((perm[i],perm[i+1]) in base_arcs for i in range(n-1)):
            t_hps.add(perm)
    print(f"HPs in T: {len(t_hps)}")

    # Distinct counts per coarse pattern
    cnt_dist = Counter(coarse_counts.values())
    print(f"Distribution of (coarse->count): {dict(sorted(cnt_dist.items()))}")

    print(f"All coarse patterns and their HP counts:")
    for coarse, cnt in sorted(coarse_counts.items()):
        is_hp_in_T = coarse in t_hps
        print(f"  {coarse}: {cnt} paths {'[HP in T]' if is_hp_in_T else ''}")

    return coarse_counts

cp3 = coarse_path_analysis(3, c3, c3sc, "C3_SC (n=6)")

# Also do the transitive tournament n=3
trans3 = frozenset({(0,1),(0,2),(1,2)})
_, trans3sc = sc_blowup(3, trans3)
cp_trans3 = coarse_path_analysis(3, trans3, trans3sc, "Trans3_SC (n=6)")

# n=4 cases
print("\n" + "="*40)
for r in all_results[4]:
    _, sc4 = sc_blowup(4, r['arcs'])
    coarse_path_analysis(4, r['arcs'], sc4, f"T{r['idx']}(n=4)_SC (n=8)")

# ============================================================
# PART 5: Formula search - deeper
# ============================================================
print("\n" + "="*70)
print("PART 5: Formula analysis")
print("="*70)

print("\nAll data points (n, H_T, H_SC, H_Lex):")
for n in [3,4,5]:
    for r in all_results[n]:
        # H_SC / H_T, and H_SC / n!
        import math
        nfact = math.factorial(n)
        print(f"  n={n}, T{r['idx']}, ss={r['ss']}: H_T={r['H_T']}, H_Lex={r['H_lex']}, H_SC={r['H_sc']}")
        # Check if H_SC = H_T * (2*H_T - something)?
        # Check if H_SC = (2n-1)!! * something?

# Key pattern: for regular tournaments, H_T = n! / n (for cyclic) or similar
# Let's check H_SC vs 2^k * H_T
print("\nChecking H_SC / (2^k * H_T) for various k:")
for n in [3,4,5]:
    for r in all_results[n]:
        for k in range(0, 8):
            if r['H_T'] > 0 and r['H_sc'] % (2**k * r['H_T']) == 0:
                q = r['H_sc'] // (2**k * r['H_T'])
                print(f"  n={n} T{r['idx']}: H_SC = 2^{k} * H_T * {q}  (H_T={r['H_T']}, H_SC={r['H_sc']})")

print("\nChecking H_SC / H_T^2:")
for n in [3,4,5]:
    for r in all_results[n]:
        if r['H_T'] > 0:
            ratio = r['H_sc'] / r['H_T']**2
            print(f"  n={n} T{r['idx']}: H_SC/H_T^2 = {r['H_sc']}/{r['H_T']**2} = {ratio:.6f}")

print("\nChecking H_Lex / H_T^2 and H_Lex / H_T:")
for n in [3,4,5]:
    for r in all_results[n]:
        if r['H_T'] > 0:
            r1 = r['H_lex'] / r['H_T']**2
            r2 = r['H_lex'] / r['H_T']
            print(f"  n={n} T{r['idx']}: H_Lex={r['H_lex']}, H_Lex/H_T^2={r1:.4f}, H_Lex/H_T={r2:.4f}")

print("\nDone!")
