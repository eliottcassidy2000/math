#!/usr/bin/env python3
"""
GLMY-Walsh Bridge Investigation
================================
Investigates connections between GLMY path homology (Betti numbers)
and the Walsh-Fourier framework for tournaments.

For n=5 (all 1024) and n=6 (all 32768):
- Compute GLMY Betti numbers
- Compute H(T), t3(T), t5(T), score sequence, signed HP count S(T)
- Test correlations and formulas linking β_1, β_3 to Walsh invariants
- Verify complement duality β(T) = β(T^op)
- Compute Euler characteristic and its relation to H, t3
"""

import numpy as np
from itertools import permutations, combinations
from collections import defaultdict, Counter
import sys

# ============================================================
# CORE: Path Homology (from path_homology_v2.py)
# ============================================================

def enumerate_allowed_paths(A, n, p):
    if p < 0:
        return []
    if p == 0:
        return [(v,) for v in range(n)]
    adj = [[] for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if A[i][j] == 1:
                adj[i].append(j)
    paths = []
    stack = []
    for start in range(n):
        stack.append(([start], 1 << start))
        while stack:
            path, visited = stack.pop()
            if len(path) == p + 1:
                paths.append(tuple(path))
                continue
            v = path[-1]
            for u in adj[v]:
                if not (visited & (1 << u)):
                    stack.append((path + [u], visited | (1 << u)))
    return paths

def boundary_coeffs(path):
    p = len(path) - 1
    result = []
    for i in range(p + 1):
        face = path[:i] + path[i+1:]
        result.append(((-1)**i, face))
    return result

def build_full_boundary_matrix(allowed_p, allowed_pm1):
    if not allowed_p or not allowed_pm1:
        return np.zeros((max(len(allowed_pm1), 0), max(len(allowed_p), 0)))
    idx_pm1 = {path: i for i, path in enumerate(allowed_pm1)}
    M = np.zeros((len(allowed_pm1), len(allowed_p)))
    for j, path in enumerate(allowed_p):
        for sign, face in boundary_coeffs(path):
            if face in idx_pm1:
                M[idx_pm1[face], j] += sign
    return M

def compute_omega_basis(A, n, p, allowed_p, allowed_pm1):
    dim_Ap = len(allowed_p)
    if dim_Ap == 0:
        return np.zeros((0, 0))
    if p == 0:
        return np.eye(dim_Ap)
    allowed_pm1_set = set(allowed_pm1)
    non_allowed_faces = {}
    na_count = 0
    for j, path in enumerate(allowed_p):
        for sign, face in boundary_coeffs(path):
            if len(set(face)) == len(face) and face not in allowed_pm1_set:
                if face not in non_allowed_faces:
                    non_allowed_faces[face] = na_count
                    na_count += 1
    if na_count == 0:
        return np.eye(dim_Ap)
    P = np.zeros((na_count, dim_Ap))
    for j, path in enumerate(allowed_p):
        for sign, face in boundary_coeffs(path):
            if face in non_allowed_faces:
                P[non_allowed_faces[face], j] += sign
    U, S, Vt = np.linalg.svd(P, full_matrices=True)
    rank = sum(s > 1e-10 for s in S)
    null_space = Vt[rank:].T
    if null_space.shape[1] == 0:
        return np.zeros((dim_Ap, 0))
    return null_space

def path_betti_numbers(A, n, max_dim=None):
    if max_dim is None:
        max_dim = n - 1
    allowed = {}
    for p in range(-1, max_dim + 2):
        if p < 0:
            allowed[p] = []
        else:
            allowed[p] = enumerate_allowed_paths(A, n, p)
    omega = {}
    for p in range(max_dim + 2):
        omega[p] = compute_omega_basis(A, n, p, allowed[p], allowed[p-1])
    betti = []
    for p in range(max_dim + 1):
        dim_omega_p = omega[p].shape[1] if omega[p].ndim == 2 else 0
        if dim_omega_p == 0:
            betti.append(0)
            continue
        bd_p = build_full_boundary_matrix(allowed[p], allowed[p-1])
        bd_p_omega = bd_p @ omega[p]
        if bd_p_omega.shape[0] > 0 and bd_p_omega.shape[1] > 0:
            S_p = np.linalg.svd(bd_p_omega, compute_uv=False)
            rank_p = sum(s > 1e-8 for s in S_p)
        else:
            rank_p = 0
        ker_dim = dim_omega_p - rank_p
        dim_omega_p1 = omega[p+1].shape[1] if omega[p+1].ndim == 2 else 0
        if dim_omega_p1 > 0:
            bd_p1 = build_full_boundary_matrix(allowed[p+1], allowed[p])
            bd_p1_omega = bd_p1 @ omega[p+1]
            S_p1 = np.linalg.svd(bd_p1_omega, compute_uv=False)
            im_dim = sum(s > 1e-8 for s in S_p1)
        else:
            im_dim = 0
        beta_p = ker_dim - im_dim
        betti.append(max(0, beta_p))
    return betti

# ============================================================
# Tournament Invariants
# ============================================================

def all_tournaments(n):
    edges = [(i, j) for i in range(n) for j in range(i+1, n)]
    m = len(edges)
    for mask in range(1 << m):
        A = [[0]*n for _ in range(n)]
        for idx, (i, j) in enumerate(edges):
            if (mask >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
        yield A, mask

def complement_tournament(A, n):
    """T^op: reverse all arcs."""
    B = [[0]*n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            if i != j:
                B[i][j] = 1 - A[i][j]
    return B

def score_sequence(A, n):
    scores = sorted([sum(A[i]) for i in range(n)])
    return tuple(scores)

def count_3cycles(A, n):
    t3 = 0
    for i in range(n):
        for j in range(i+1, n):
            for k in range(j+1, n):
                if (A[i][j] and A[j][k] and A[k][i]) or \
                   (A[j][i] and A[k][j] and A[i][k]):
                    t3 += 1
    return t3

def count_5cycles(A, n):
    """Count directed 5-cycles (up to rotation)."""
    count = 0
    for combo in combinations(range(n), 5):
        for perm in permutations(combo):
            is_cycle = True
            for idx in range(5):
                if A[perm[idx]][perm[(idx+1) % 5]] != 1:
                    is_cycle = False
                    break
            if is_cycle:
                count += 1
    return count // 5  # each cycle counted 5 times (rotations)

def ham_path_count(A, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)) or (mask, v) not in dp:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if A[v][u]:
                    key = (mask | (1 << u), u)
                    dp[key] = dp.get(key, 0) + dp[(mask, v)]
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

def signed_hp_count(A, n):
    """S(T) = sum_P prod_{i} B[P_i][P_{i+1}] where B = 2A - J."""
    B = [[2*A[i][j] - 1 if i != j else 0 for j in range(n)] for i in range(n)]
    S = 0
    for perm in permutations(range(n)):
        prod = 1
        for i in range(n - 1):
            prod *= B[perm[i]][perm[i+1]]
        S += prod
    return S

def isomorphism_type(A, n):
    """Return canonical score sequence as a rough type identifier."""
    return score_sequence(A, n)

# ============================================================
# Main Analysis
# ============================================================

def analyze_n(n, max_betti_dim=None):
    if max_betti_dim is None:
        max_betti_dim = n - 1

    print(f"\n{'='*70}")
    print(f"  ANALYSIS: n = {n}  (all {1 << (n*(n-1)//2)} tournaments)")
    print(f"{'='*70}")

    records = []
    total = 1 << (n*(n-1)//2)

    for idx, (A, mask) in enumerate(all_tournaments(n)):
        if (idx + 1) % 1000 == 0:
            print(f"  ... processed {idx+1}/{total}", file=sys.stderr)

        betti = path_betti_numbers(A, n, max_dim=max_betti_dim)
        H = ham_path_count(A, n)
        t3 = count_3cycles(A, n)
        score = score_sequence(A, n)

        # Signed HP count (expensive at n>=7, fine for n<=6)
        S = signed_hp_count(A, n) if n <= 6 else None

        # 5-cycles only for small n
        t5 = count_5cycles(A, n) if n <= 6 else None

        # Complement
        Ac = complement_tournament(A, n)

        # Euler characteristic
        chi = sum((-1)**m * betti[m] for m in range(len(betti)))

        records.append({
            'A': A, 'mask': mask, 'betti': tuple(betti),
            'H': H, 't3': t3, 't5': t5, 'score': score,
            'S': S, 'chi': chi, 'Ac': Ac
        })

    return records


def report_section(title):
    print(f"\n{'─'*70}")
    print(f"  {title}")
    print(f"{'─'*70}")


def main():
    print("=" * 70)
    print("  GLMY PATH HOMOLOGY ↔ WALSH-FOURIER BRIDGE")
    print("  Investigating connections between Betti numbers and Walsh invariants")
    print("=" * 70)

    # ── n = 5 ──
    records5 = analyze_n(5, max_betti_dim=4)

    report_section("1. BETTI DISTRIBUTION (n=5)")
    betti_dist = Counter(r['betti'] for r in records5)
    for bt in sorted(betti_dist.keys()):
        print(f"  β = {list(bt)}: {betti_dist[bt]} tournaments")

    report_section("2. DETAILED TABLE BY ISOMORPHISM TYPE (n=5)")
    # Group by (score, t3) as rough type
    type_groups = defaultdict(list)
    for r in records5:
        key = (r['score'], r['t3'])
        type_groups[key].append(r)

    print(f"  {'Score':<20} {'t3':>4} {'t5':>4} {'H':>4} {'S':>6} {'β':>20} {'χ':>4}  count")
    print(f"  {'─'*20} {'─'*4} {'─'*4} {'─'*4} {'─'*6} {'─'*20} {'─'*4}  ─────")
    for key in sorted(type_groups.keys()):
        reps = type_groups[key]
        r0 = reps[0]
        # Check if all in group have same betti
        bettis = set(r['betti'] for r in reps)
        Hs = set(r['H'] for r in reps)
        Ss = set(r['S'] for r in reps)
        t5s = set(r['t5'] for r in reps)
        for bt in sorted(bettis):
            cnt = sum(1 for r in reps if r['betti'] == bt)
            rep = next(r for r in reps if r['betti'] == bt)
            print(f"  {str(r0['score']):<20} {r0['t3']:>4} {rep['t5']:>4} {rep['H']:>4} {rep['S']:>6} {str(list(bt)):<20} {rep['chi']:>4}  ×{cnt}")

    report_section("3. β_1 vs (t3, H, score, S) CORRELATION (n=5)")
    # At n=5, β_1 ∈ {0, 1}
    b1_vals = [r['betti'][1] for r in records5]
    t3_vals = [r['t3'] for r in records5]
    H_vals = [r['H'] for r in records5]
    S_vals = [r['S'] for r in records5]

    # Group by β_1
    for b1 in sorted(set(b1_vals)):
        subset = [r for r in records5 if r['betti'][1] == b1]
        t3s = [r['t3'] for r in subset]
        Hs = [r['H'] for r in subset]
        Ss = [r['S'] for r in subset]
        scores = Counter(r['score'] for r in subset)
        print(f"\n  β_1 = {b1}: {len(subset)} tournaments")
        print(f"    t3 range: {min(t3s)} to {max(t3s)}, values: {sorted(set(t3s))}")
        print(f"    H range:  {min(Hs)} to {max(Hs)}, values: {sorted(set(Hs))}")
        print(f"    S range:  {min(Ss)} to {max(Ss)}, values: {sorted(set(Ss))}")
        print(f"    scores: {dict(scores)}")

    # Test: is β_1 determined by score sequence alone?
    score_to_b1 = defaultdict(set)
    for r in records5:
        score_to_b1[r['score']].add(r['betti'][1])
    ambiguous = {s: vals for s, vals in score_to_b1.items() if len(vals) > 1}
    print(f"\n  Score sequence determines β_1? {'YES' if not ambiguous else 'NO'}")
    if ambiguous:
        print(f"    Ambiguous scores: {ambiguous}")

    # Test: is β_1 determined by t3?
    t3_to_b1 = defaultdict(set)
    for r in records5:
        t3_to_b1[r['t3']].add(r['betti'][1])
    ambiguous_t3 = {t: vals for t, vals in t3_to_b1.items() if len(vals) > 1}
    print(f"  t3 determines β_1? {'YES' if not ambiguous_t3 else 'NO'}")
    if ambiguous_t3:
        print(f"    Ambiguous t3 values: {ambiguous_t3}")

    # Test: is β_1 determined by (score, t3)?
    st_to_b1 = defaultdict(set)
    for r in records5:
        st_to_b1[(r['score'], r['t3'])].add(r['betti'][1])
    ambiguous_st = {k: v for k, v in st_to_b1.items() if len(v) > 1}
    print(f"  (score, t3) determines β_1? {'YES' if not ambiguous_st else 'NO'}")

    # Test: is β_1 determined by S?
    S_to_b1 = defaultdict(set)
    for r in records5:
        S_to_b1[r['S']].add(r['betti'][1])
    ambiguous_S = {k: v for k, v in S_to_b1.items() if len(v) > 1}
    print(f"  S (signed HP) determines β_1? {'YES' if not ambiguous_S else 'NO'}")

    # Try simple formulas
    report_section("4. TESTING FORMULAS FOR β_1 (n=5)")

    # β_1 = 1 iff t3 == 0?
    formula1 = all((r['betti'][1] == 1) == (r['t3'] == 0) for r in records5)
    print(f"  β_1 = [t3 == 0]?  {formula1}")

    # β_1 = 1 iff score = (1,1,2,3,3)?  (transitive has score (0,1,2,3,4))
    trans_score = (0, 1, 2, 3, 4)
    formula2 = all((r['betti'][1] == 1) == (r['score'] == trans_score) for r in records5)
    print(f"  β_1 = [score == (0,1,2,3,4)]?  {formula2}")

    # β_1 = 1 iff H == max(H)?
    max_H = max(r['H'] for r in records5)
    formula3 = all((r['betti'][1] == 1) == (r['H'] == max_H) for r in records5)
    print(f"  β_1 = [H == {max_H}]?  {formula3}")

    # Euler characteristic
    report_section("5. EULER CHARACTERISTIC (n=5)")
    print(f"  chi = Σ (-1)^m β_m")
    chi_vals = Counter(r['chi'] for r in records5)
    for c in sorted(chi_vals.keys()):
        print(f"    χ = {c}: {chi_vals[c]} tournaments")

    # chi vs t3
    chi_t3 = defaultdict(set)
    for r in records5:
        chi_t3[r['t3']].add(r['chi'])
    print(f"\n  χ by t3:")
    for t3 in sorted(chi_t3.keys()):
        print(f"    t3={t3}: χ ∈ {chi_t3[t3]}")

    # chi = 1 - β_1 check
    chi_formula = all(r['chi'] == 1 - r['betti'][1] for r in records5)
    print(f"\n  χ = 1 - β_1?  {chi_formula}  (expected: β_0=1, β_2=β_3=β_4=0)")

    # ── Complement Duality ──
    report_section("6. COMPLEMENT DUALITY β(T) = β(T^op) (n=5)")
    duality_holds = True
    duality_failures = 0
    for r in records5:
        betti_comp = path_betti_numbers(r['Ac'], 5, max_dim=4)
        if tuple(betti_comp) != r['betti']:
            duality_holds = False
            duality_failures += 1
    print(f"  β(T) = β(T^op) for ALL n=5 tournaments?  {duality_holds}")
    if not duality_holds:
        print(f"    Failures: {duality_failures}/{len(records5)}")

    # ── n = 6 ──
    print(f"\n\n{'='*70}")
    print(f"  NOW ANALYZING n = 6  (32768 tournaments, may take a few minutes)")
    print(f"{'='*70}")

    records6 = analyze_n(6, max_betti_dim=5)

    report_section("7. BETTI DISTRIBUTION (n=6)")
    betti_dist6 = Counter(r['betti'] for r in records6)
    for bt in sorted(betti_dist6.keys()):
        print(f"  β = {list(bt)}: {betti_dist6[bt]} tournaments")

    report_section("8. β_1 ANALYSIS (n=6)")
    for b1 in sorted(set(r['betti'][1] for r in records6)):
        subset = [r for r in records6 if r['betti'][1] == b1]
        t3s = sorted(set(r['t3'] for r in subset))
        Hs = sorted(set(r['H'] for r in subset))
        Ss = sorted(set(r['S'] for r in subset))
        scores = Counter(r['score'] for r in subset)
        print(f"\n  β_1 = {b1}: {len(subset)} tournaments")
        print(f"    t3 values: {t3s}")
        print(f"    H values: {Hs}")
        print(f"    S values: {Ss}")

    # Determinacy tests at n=6
    score_to_b1_6 = defaultdict(set)
    t3_to_b1_6 = defaultdict(set)
    st_to_b1_6 = defaultdict(set)
    S_to_b1_6 = defaultdict(set)
    for r in records6:
        score_to_b1_6[r['score']].add(r['betti'][1])
        t3_to_b1_6[r['t3']].add(r['betti'][1])
        st_to_b1_6[(r['score'], r['t3'])].add(r['betti'][1])
        S_to_b1_6[r['S']].add(r['betti'][1])

    print(f"\n  Score determines β_1? {'YES' if all(len(v)==1 for v in score_to_b1_6.values()) else 'NO'}")
    print(f"  t3 determines β_1? {'YES' if all(len(v)==1 for v in t3_to_b1_6.values()) else 'NO'}")
    print(f"  (score, t3) determines β_1? {'YES' if all(len(v)==1 for v in st_to_b1_6.values()) else 'NO'}")
    print(f"  S determines β_1? {'YES' if all(len(v)==1 for v in S_to_b1_6.values()) else 'NO'}")

    report_section("9. β_3 ANALYSIS (n=6)")
    for b3_val in sorted(set(r['betti'][3] for r in records6)):
        subset = [r for r in records6 if r['betti'][3] == b3_val]
        print(f"\n  β_3 = {b3_val}: {len(subset)} tournaments")
        if len(subset) <= 50:
            t3s = sorted(set(r['t3'] for r in subset))
            Hs = sorted(set(r['H'] for r in subset))
            Ss = sorted(set(r['S'] for r in subset))
            scores = Counter(r['score'] for r in subset)
            print(f"    t3 values: {t3s}")
            print(f"    H values: {Hs}")
            print(f"    S values: {Ss}")
            print(f"    scores: {dict(scores)}")

    # β_3 Walsh characterization
    S_to_b3 = defaultdict(set)
    t3_to_b3 = defaultdict(set)
    for r in records6:
        S_to_b3[r['S']].add(r['betti'][3])
        t3_to_b3[r['t3']].add(r['betti'][3])
    print(f"\n  S determines β_3? {'YES' if all(len(v)==1 for v in S_to_b3.values()) else 'NO'}")
    print(f"  t3 determines β_3? {'YES' if all(len(v)==1 for v in t3_to_b3.values()) else 'NO'}")

    # β_1 and β_3 mutual exclusivity
    report_section("10. β_1 AND β_3 MUTUAL EXCLUSIVITY (n=6)")
    both_nonzero = sum(1 for r in records6 if r['betti'][1] > 0 and r['betti'][3] > 0)
    print(f"  Tournaments with β_1 > 0 AND β_3 > 0: {both_nonzero}")
    print(f"  Mutual exclusivity holds? {'YES' if both_nonzero == 0 else 'NO'}")

    report_section("11. EULER CHARACTERISTIC (n=6)")
    chi_vals6 = Counter(r['chi'] for r in records6)
    for c in sorted(chi_vals6.keys()):
        print(f"    χ = {c}: {chi_vals6[c]} tournaments")

    # chi = 1 - β_1 + β_3 ?  (β_0=1, β_2=0)
    chi_formula6 = all(r['chi'] == sum((-1)**m * r['betti'][m] for m in range(len(r['betti']))) for r in records6)
    print(f"\n  χ = Σ(-1)^m β_m? {chi_formula6} (tautology check)")

    # Check relationship to H, t3
    chi_to_invariants = defaultdict(lambda: {'H': set(), 't3': set(), 'S': set()})
    for r in records6:
        chi_to_invariants[r['chi']]['H'].add(r['H'])
        chi_to_invariants[r['chi']]['t3'].add(r['t3'])
        chi_to_invariants[r['chi']]['S'].add(r['S'])
    print(f"\n  χ vs invariants:")
    for c in sorted(chi_to_invariants.keys()):
        d = chi_to_invariants[c]
        print(f"    χ={c}: H ∈ {sorted(d['H'])}, t3 ∈ {sorted(d['t3'])}")

    # ── Complement Duality n=6 ──
    report_section("12. COMPLEMENT DUALITY β(T) = β(T^op) (n=6)")
    # For efficiency: compute complement mask and look up
    edges6 = [(i, j) for i in range(6) for j in range(i+1, 6)]
    m6 = len(edges6)
    full_mask = (1 << m6) - 1

    # Build mask -> record lookup
    mask_to_betti = {}
    for r in records6:
        mask_to_betti[r['mask']] = r['betti']

    duality_holds_6 = True
    duality_failures_6 = 0
    for r in records6:
        comp_mask = full_mask ^ r['mask']
        if comp_mask in mask_to_betti:
            if mask_to_betti[comp_mask] != r['betti']:
                duality_holds_6 = False
                duality_failures_6 += 1
    print(f"  β(T) = β(T^op) for ALL n=6 tournaments?  {duality_holds_6}")
    if not duality_holds_6:
        print(f"    Failures: {duality_failures_6}/{len(records6)}")

    # ── Formula search ──
    report_section("13. FORMULA SEARCH: β_1 = f(score, t3, t5, H, S)")

    # n=5: try β_1 = some function
    print("\n  n=5:")
    # Check: β_1 = 1 iff t3 = 0
    check_t3_zero_5 = all((r['betti'][1] == 1) == (r['t3'] == 0) for r in records5)
    print(f"    β_1 = [t3==0]?  {check_t3_zero_5}")

    # Check: β_1 is related to score being strictly monotone (transitive)
    check_transitive_5 = all((r['betti'][1] == 1) == (len(set(r['score'])) == 5) for r in records5)
    print(f"    β_1 = [all scores distinct]?  {check_transitive_5}")

    # n=6: same checks
    print("\n  n=6:")
    check_t3_zero_6 = all((r['betti'][1] == 1) == (r['t3'] == 0) for r in records6)
    print(f"    β_1 = [t3==0]?  {check_t3_zero_6}")

    # Actually at n=5 and n=6, t3=0 iff transitive (acyclic). Check.
    t3_zero_scores_5 = set(r['score'] for r in records5 if r['t3'] == 0)
    t3_zero_scores_6 = set(r['score'] for r in records6 if r['t3'] == 0)
    print(f"\n    n=5 t3=0 scores: {t3_zero_scores_5}")
    print(f"    n=6 t3=0 scores: {t3_zero_scores_6}")

    # For β_1 at n=6 (if β_1 can be >1):
    b1_vals_6 = sorted(set(r['betti'][1] for r in records6))
    print(f"\n    n=6 β_1 values: {b1_vals_6}")

    # Look for β_1 formula at n=6 involving t3
    if max(b1_vals_6) > 1:
        # More complex formula needed
        for b1 in b1_vals_6:
            if b1 > 0:
                subset = [r for r in records6 if r['betti'][1] == b1]
                t3s = sorted(set(r['t3'] for r in subset))
                print(f"    β_1={b1}: t3 ∈ {t3s}")

    # ── Walsh spectral characterization ──
    report_section("14. WALSH SPECTRAL CHARACTERIZATION ATTEMPT")
    print("  Testing if S (signed HP count) separates Betti classes...")

    # n=5
    S_betti_5 = defaultdict(Counter)
    for r in records5:
        S_betti_5[r['S']][r['betti']] += 1
    print(f"\n  n=5: S → β distribution:")
    for S_val in sorted(S_betti_5.keys()):
        print(f"    S={S_val:>4}: {dict(S_betti_5[S_val])}")

    # n=6
    S_betti_6 = defaultdict(Counter)
    for r in records6:
        S_betti_6[r['S']][r['betti']] += 1
    print(f"\n  n=6: S → β distribution:")
    for S_val in sorted(S_betti_6.keys()):
        print(f"    S={S_val:>4}: {dict(S_betti_6[S_val])}")

    # ── t5 analysis ──
    report_section("15. 5-CYCLE (t5) ANALYSIS")
    print("\n  n=5:")
    t5_betti_5 = defaultdict(Counter)
    for r in records5:
        t5_betti_5[r['t5']][r['betti']] += 1
    for t5 in sorted(t5_betti_5.keys()):
        print(f"    t5={t5}: {dict(t5_betti_5[t5])}")

    print("\n  n=6:")
    t5_betti_6 = defaultdict(Counter)
    for r in records6:
        t5_betti_6[r['t5']][r['betti']] += 1
    for t5 in sorted(t5_betti_6.keys()):
        print(f"    t5={t5}: {dict(t5_betti_6[t5])}")

    # ── Combined formula attempt ──
    report_section("16. COMBINED FORMULA SEARCH")

    # At n=6, test: β_1 = C(n,2) - n - t3_related?
    # Actually let's think about it: β_1 counts independent directed 1-cycles
    # not bounding any 2-chain. Since β_2 = 0, every 2-cycle bounds.
    # β_1 = dim(ker ∂_1 on Ω_1) - dim(im ∂_2 on Ω_2)

    # Simple numerical correlation
    if len(records6) > 0:
        b1_arr = np.array([r['betti'][1] for r in records6])
        t3_arr = np.array([r['t3'] for r in records6])
        H_arr = np.array([r['H'] for r in records6])
        S_arr = np.array([r['S'] for r in records6])

        if b1_arr.std() > 0:
            corr_t3 = np.corrcoef(b1_arr, t3_arr)[0, 1]
            corr_H = np.corrcoef(b1_arr, H_arr)[0, 1]
            corr_S = np.corrcoef(b1_arr, S_arr)[0, 1]
            print(f"\n  n=6 correlations with β_1:")
            print(f"    corr(β_1, t3) = {corr_t3:.6f}")
            print(f"    corr(β_1, H)  = {corr_H:.6f}")
            print(f"    corr(β_1, S)  = {corr_S:.6f}")

    # Same for β_3 at n=6
    b3_arr = np.array([r['betti'][3] for r in records6])
    if b3_arr.std() > 0:
        corr_t3_b3 = np.corrcoef(b3_arr, t3_arr)[0, 1]
        corr_H_b3 = np.corrcoef(b3_arr, H_arr)[0, 1]
        corr_S_b3 = np.corrcoef(b3_arr, S_arr)[0, 1]
        t5_arr = np.array([r['t5'] for r in records6])
        corr_t5_b3 = np.corrcoef(b3_arr, t5_arr)[0, 1]
        print(f"\n  n=6 correlations with β_3:")
        print(f"    corr(β_3, t3) = {corr_t3_b3:.6f}")
        print(f"    corr(β_3, H)  = {corr_H_b3:.6f}")
        print(f"    corr(β_3, S)  = {corr_S_b3:.6f}")
        print(f"    corr(β_3, t5) = {corr_t5_b3:.6f}")

    # ── Summary ──
    report_section("SUMMARY OF FINDINGS")
    print("""
  Key questions answered:
  1. Does score sequence determine β_1? (checked n=5, n=6)
  2. Does t3 determine β_1? (checked n=5, n=6)
  3. Does signed HP count S determine β? (checked n=5, n=6)
  4. β_1 and β_3 mutual exclusivity (checked n=6)
  5. Complement duality β(T) = β(T^op) (checked n=5, n=6)
  6. Euler characteristic vs tournament invariants
  7. Correlations between Betti numbers and Walsh invariants
    """)

    print("\nDone.")


if __name__ == '__main__':
    main()
