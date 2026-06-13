"""
sigma_spectral_detail.py -- kind-pasteur-2026-03-13-S61

The sigma graph spectrum resolves ALL simplex profile ambiguities.
This script extracts the EXACT spectral features that differ.

Key questions:
1. What eigenvalues differ between lo and hi?
2. Is it a single eigenvalue that shifts?
3. Can we identify a simpler invariant (trace of sigma^k for some k)?
4. Is it the "labeled sigma profile" = sigma degree sequence per vertex?
5. Connection to the tournament's skew structure
"""

import numpy as np
from itertools import combinations
from collections import Counter, defaultdict

def bits_to_adj(bits, n):
    A = np.zeros((n, n), dtype=int)
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def lambda_graph(A, n):
    L = np.zeros((n, n), dtype=int)
    for u in range(n):
        for v in range(u+1, n):
            for w in range(n):
                if w == u or w == v:
                    continue
                if (A[u][v] and A[v][w] and A[w][u]) or (A[v][u] and A[u][w] and A[w][v]):
                    L[u][v] += 1
                    L[v][u] += 1
    return L

def count_directed_k_cycles(A, n, k):
    Ak = np.linalg.matrix_power(A, k)
    return int(np.trace(Ak)) // k

def get_pair_data(A, n):
    L = lambda_graph(A, n)
    A2 = A @ A
    pairs = []
    for u in range(n):
        for v in range(u+1, n):
            sig = n - 2 - int(A2[u][v]) - int(A2[v][u])
            lam = int(L[u][v])
            delta = n - 2 - sig - lam
            pairs.append((u, v, sig, lam, delta))
    return pairs

n = 7
total_bits = n * (n-1) // 2

print("=" * 60)
print("SIGMA SPECTRAL DETAIL ANALYSIS AT n=7")
print("=" * 60)

np.random.seed(42)

# Collect data
profile_groups = defaultdict(list)

for trial in range(10000):
    bits = np.random.randint(0, 1 << total_bits)
    A = bits_to_adj(bits, n)
    pair_data = get_pair_data(A, n)
    c7 = count_directed_k_cycles(A, n, 7)
    profile = tuple(sorted((s, l, d) for u, v, s, l, d in pair_data))

    # Build sigma matrix
    sig_mat = np.zeros((n, n), dtype=int)
    lam_mat = np.zeros((n, n), dtype=int)
    del_mat = np.zeros((n, n), dtype=int)
    for u, v, s, l, d in pair_data:
        sig_mat[u][v] = sig_mat[v][u] = s
        lam_mat[u][v] = lam_mat[v][u] = l
        del_mat[u][v] = del_mat[v][u] = d

    # Sigma degree sequence (sum of sigma weights per vertex)
    sig_deg = tuple(sorted(int(sig_mat[i].sum()) for i in range(n)))
    # Lambda degree sequence
    lam_deg = tuple(sorted(int(lam_mat[i].sum()) for i in range(n)))
    # Delta degree sequence
    del_deg = tuple(sorted(int(del_mat[i].sum()) for i in range(n)))

    profile_groups[profile].append({
        'bits': bits,
        'c7': c7,
        'sig_mat': sig_mat.copy(),
        'lam_mat': lam_mat.copy(),
        'del_mat': del_mat.copy(),
        'sig_deg': sig_deg,
        'lam_deg': lam_deg,
        'del_deg': del_deg,
        'A': A.copy(),
    })

# Find ambiguous profiles
ambig = {p: entries for p, entries in profile_groups.items()
         if len(set(e['c7'] for e in entries)) > 1}

print(f"\nAmbiguous profiles: {len(ambig)}")

for prof_idx, (profile, entries) in enumerate(sorted(ambig.items(),
                                                      key=lambda x: -len(x[1]))):
    c7_vals = sorted(set(e['c7'] for e in entries))
    if len(c7_vals) != 2:
        continue
    c7_lo, c7_hi = c7_vals

    lo = [e for e in entries if e['c7'] == c7_lo]
    hi = [e for e in entries if e['c7'] == c7_hi]

    print(f"\n{'='*50}")
    print(f"Profile {prof_idx+1}: c7={c7_lo} vs {c7_hi} (dc7={c7_hi-c7_lo})")

    # 1. Sigma degree sequence
    lo_sdeg = set(e['sig_deg'] for e in lo)
    hi_sdeg = set(e['sig_deg'] for e in hi)
    print(f"  Sigma deg seq: lo={lo_sdeg}, hi={hi_sdeg}")
    if lo_sdeg != hi_sdeg:
        print(f"    -> SIGMA DEGREE SEQUENCE DISTINGUISHES!")
    else:
        print(f"    -> Sigma deg seq same (need deeper invariant)")

    # 2. Lambda degree sequence
    lo_ldeg = set(e['lam_deg'] for e in lo)
    hi_ldeg = set(e['lam_deg'] for e in hi)
    print(f"  Lambda deg seq: lo={lo_ldeg}, hi={hi_ldeg}")

    # 3. Delta degree sequence
    lo_ddeg = set(e['del_deg'] for e in lo)
    hi_ddeg = set(e['del_deg'] for e in hi)
    print(f"  Delta deg seq: lo={lo_ddeg}, hi={hi_ddeg}")

    # 4. Trace of sigma^k for k=2,3,4
    for label, group in [("lo", lo[:5]), ("hi", hi[:5])]:
        e = group[0]
        sig = e['sig_mat'].astype(float)
        for k in [2, 3, 4]:
            tr = int(round(np.trace(np.linalg.matrix_power(sig, k))))
            print(f"    {label}: tr(sigma^{k}) = {tr}", end="")
        print()

    # 5. Check: do ALL lo have same sigma degree sequence? Same for hi?
    lo_unique_sdeg = set(e['sig_deg'] for e in lo)
    hi_unique_sdeg = set(e['sig_deg'] for e in hi)
    print(f"  Sigma deg seq varieties: lo has {len(lo_unique_sdeg)}, hi has {len(hi_unique_sdeg)}")

    # 6. The full sigma matrix spectrum (detailed eigenvalues)
    e_lo = lo[0]
    e_hi = hi[0]
    sig_eigs_lo = sorted(np.linalg.eigvalsh(e_lo['sig_mat'].astype(float)))
    sig_eigs_hi = sorted(np.linalg.eigvalsh(e_hi['sig_mat'].astype(float)))
    print(f"  Sigma eigenvalues:")
    print(f"    lo: [{', '.join(f'{x:.4f}' for x in sig_eigs_lo)}]")
    print(f"    hi: [{', '.join(f'{x:.4f}' for x in sig_eigs_hi)}]")
    diffs = [abs(a - b) for a, b in zip(sig_eigs_lo, sig_eigs_hi)]
    print(f"    diffs: [{', '.join(f'{x:.4f}' for x in diffs)}]")

    # 7. Score sequences (from the original tournament)
    lo_scores = set(tuple(sorted(int(sum(e['A'][i])) for i in range(n))) for e in lo)
    hi_scores = set(tuple(sorted(int(sum(e['A'][i])) for i in range(n))) for e in hi)
    print(f"  Tournament scores: lo={lo_scores}, hi={hi_scores}")

# Now the BIG question: does sigma degree sequence ALONE determine c7
# (within fixed simplex profile)?
print(f"\n\n{'='*60}")
print("SIGMA DEGREE SEQUENCE AS HIDDEN DIMENSION")
print("=" * 60)

# For ALL profiles (not just ambiguous), check if simplex profile + sigma deg seq
# determines c7
total_combos = 0
ambig_combos = 0
for profile, entries in profile_groups.items():
    if len(entries) < 2:
        continue
    sdeg_c7 = defaultdict(set)
    for e in entries:
        sdeg_c7[e['sig_deg']].add(e['c7'])
    for sdeg, c7s in sdeg_c7.items():
        total_combos += 1
        if len(c7s) > 1:
            ambig_combos += 1

print(f"  (profile, sigma_deg_seq) pairs: {total_combos}")
print(f"  Ambiguous for c7: {ambig_combos}")
if ambig_combos == 0:
    print(f"  -> SIMPLEX PROFILE + SIGMA DEGREE SEQUENCE DETERMINES c7!")

# Even simpler: does sigma degree sequence ALONE determine c7?
print(f"\n--- Sigma Degree Sequence Alone ---")
sdeg_c7_global = defaultdict(set)
for entries in profile_groups.values():
    for e in entries:
        sdeg_c7_global[e['sig_deg']].add(e['c7'])
ambig_sdeg = sum(1 for c7s in sdeg_c7_global.values() if len(c7s) > 1)
print(f"  Distinct sigma deg seqs: {len(sdeg_c7_global)}")
print(f"  Ambiguous for c7: {ambig_sdeg}")

# What about lambda degree sequence alone?
print(f"\n--- Lambda Degree Sequence Alone ---")
ldeg_c7_global = defaultdict(set)
for entries in profile_groups.values():
    for e in entries:
        ldeg_c7_global[e['lam_deg']].add(e['c7'])
ambig_ldeg = sum(1 for c7s in ldeg_c7_global.values() if len(c7s) > 1)
print(f"  Distinct lambda deg seqs: {len(ldeg_c7_global)}")
print(f"  Ambiguous for c7: {ambig_ldeg}")

# What about (sigma_deg, lambda_deg) together?
print(f"\n--- (Sigma, Lambda) Degree Sequences Together ---")
sl_deg_c7 = defaultdict(set)
for entries in profile_groups.values():
    for e in entries:
        sl_deg_c7[(e['sig_deg'], e['lam_deg'])].add(e['c7'])
ambig_sl = sum(1 for c7s in sl_deg_c7.values() if len(c7s) > 1)
print(f"  Distinct (sig_deg, lam_deg) pairs: {len(sl_deg_c7)}")
print(f"  Ambiguous for c7: {ambig_sl}")

# What about the tournament score sequence + c3?
print(f"\n--- Score Sequence + c3 ---")
score_c3_c7 = defaultdict(set)
c3_data = {}
for entries in profile_groups.values():
    for e in entries:
        scores = tuple(sorted(int(sum(e['A'][i])) for i in range(n)))
        c3 = count_directed_k_cycles(e['A'], n, 3)
        score_c3_c7[(scores, c3)].add(e['c7'])
ambig_sc = sum(1 for c7s in score_c3_c7.values() if len(c7s) > 1)
print(f"  Distinct (scores, c3): {len(score_c3_c7)}")
print(f"  Ambiguous for c7: {ambig_sc}")

# Score sequence + c3 + c5?
print(f"\n--- Score Sequence + c3 + c5 ---")
sc35_c7 = defaultdict(set)
for entries in profile_groups.values():
    for e in entries:
        scores = tuple(sorted(int(sum(e['A'][i])) for i in range(n)))
        c3 = count_directed_k_cycles(e['A'], n, 3)
        c5 = count_directed_k_cycles(e['A'], n, 5)
        sc35_c7[(scores, c3, c5)].add(e['c7'])
ambig_sc35 = sum(1 for c7s in sc35_c7.values() if len(c7s) > 1)
print(f"  Distinct (scores, c3, c5): {len(sc35_c7)}")
print(f"  Ambiguous for c7: {ambig_sc35}")

print("\n\nDone.")
