#!/usr/bin/env python3
"""
beta3_seesaw_approach.py — opus-2026-03-09-S53

ALGEBRAIC APPROACH via adjacent-odd seesaw:
  beta_1 * beta_3 = 0 for all tournaments

This means: if beta_3 > 0, then beta_1 = 0 (T is strongly connected).

Also from S45/S46: beta_{2k-1} * beta_{2k+1} = 0 generally.

Key question: can we use the seesaw together with the LES to get beta_3 ≤ 1?

The LES approach alone gives: beta_3(T) ≤ beta_3(T\v) + 1
(from H_3(T,T\v) ≤ 1 and rank-nullity)

The seesaw gives: beta_3 > 0 → beta_1 = 0 (strongly connected)

Combining: if beta_3(T) = 2, then:
1. T is strongly connected (from seesaw)
2. ALL deletions T\v have beta_3(T\v) ≥ 1 (from LES)
3. Each T\v with beta_3=1 is also strongly connected (from seesaw at n-1)

Can we get a contradiction from (1) + (2) + (3)?

IDEA: Strong connectivity constrains the score sequence.
A strongly connected tournament has min score ≥ 1 and max score ≤ n-2.
(Vertices with score 0 or n-1 make T not strongly connected.)

Also: if beta_3(T\v) ≥ 1 for all v, then ALL T\v are strongly connected.
This means: for every v, T\v has no vertex with score 0 or n-2 (=5 at n=7).

What does this mean for T's score sequence?
If v has score d in T, and v beats vertices S⊂V\{v} (|S|=d):
- Each w∈S has score d_w in T, and d_w-1 in T\v (they lost one arc to v)
  Wait no: in T\v, the score of w is (score in T) - T[w][v].
  If v→w (w∈S, w lost to v), then T[w][v]=0, so score_w(T\v) = d_w
  If w→v (w∉S, w beat v), then T[w][v]=1, so score_w(T\v) = d_w - 1

So in T\v:
- vertices that v beats: keep their score
- vertices that beat v: lose 1 from their score

For T\v to be strongly connected (min score ≥ 1 in T\v):
- vertices beating v: d_w - 1 ≥ 1, so d_w ≥ 2
- vertices beaten by v: d_w ≥ 1 (they keep score ≥ 1)
- Also max score ≤ n-3 (= 4 at n=7) in T\v:
  - vertices beating v: d_w - 1 ≤ 4, so d_w ≤ 5
  - vertices beaten by v: d_w ≤ 4

Wait, strong connectivity doesn't mean min score ≥ 1. A tournament
can be strongly connected with a vertex of score 1. It just can't have
score 0 (sink — no outgoing edges means it can't reach others).

Actually: score 0 = receives from everyone = sink. Not in any outpath.
NOT strongly connected only if vertex has in-degree n-1 or out-degree n-1?
No, strong connectivity is more subtle.

Let me just compute: for which score sequences is beta_1 > 0 at n=6?
(beta_1 > 0 means NOT strongly connected)
"""

import numpy as np
from itertools import combinations, permutations
from collections import Counter, defaultdict

def tournament_from_bits(n, bits):
    T = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if bits & (1 << idx):
                T[i][j] = 1
            else:
                T[j][i] = 1
            idx += 1
    return T

def is_strongly_connected(T, n):
    """Check if tournament T on n vertices is strongly connected."""
    # BFS from vertex 0
    visited = {0}
    queue = [0]
    while queue:
        v = queue.pop(0)
        for w in range(n):
            if w not in visited and T[v][w]:
                visited.add(w)
                queue.append(w)
    if len(visited) < n:
        return False
    # Also check reverse reachability
    visited = {0}
    queue = [0]
    while queue:
        v = queue.pop(0)
        for w in range(n):
            if w not in visited and T[w][v]:
                visited.add(w)
                queue.append(w)
    return len(visited) == n

def is_allowed_path(T, path):
    for i in range(len(path)-1):
        if not T[path[i]][path[i+1]]:
            return False
    return len(path) == len(set(path))

def compute_beta3(T, n):
    tol = 1e-8
    all_paths = {}
    for p in range(0, 6):
        paths = []
        if p + 1 <= n:
            for verts in combinations(range(n), p+1):
                for perm in permutations(verts):
                    if is_allowed_path(T, perm):
                        paths.append(perm)
        all_paths[p] = paths

    def build_omega(pd, mp):
        omega = {}
        for p in range(0, mp + 1):
            a_p = pd[p]
            if not a_p:
                omega[p] = np.zeros((0, 0))
                continue
            a_pm1_set = set(pd[p-1]) if p > 0 else set()
            na = {}
            for sigma in a_p:
                for i in range(1, len(sigma)-1):
                    face = sigma[:i] + sigma[i+1:]
                    if p > 0 and face not in a_pm1_set:
                        if face not in na:
                            na[face] = len(na)
            if not na:
                omega[p] = np.eye(len(a_p))
            else:
                mat = np.zeros((len(na), len(a_p)))
                for j, sigma in enumerate(a_p):
                    for i in range(1, len(sigma)-1):
                        face = sigma[:i] + sigma[i+1:]
                        if face in na:
                            mat[na[face], j] += (-1)**i
                U, S, Vt = np.linalg.svd(mat, full_matrices=True)
                rank = int(np.sum(S > tol))
                if len(a_p) - rank == 0:
                    omega[p] = np.zeros((len(a_p), 0))
                else:
                    omega[p] = Vt[rank:].T
        return omega

    omega = build_omega(all_paths, 4)
    boundary = {}
    for p in range(1, 5):
        a_p = all_paths[p]
        a_pm1 = all_paths[p-1]
        if not a_p or not a_pm1:
            boundary[p] = np.zeros((0, 0))
            continue
        idx = {path: i for i, path in enumerate(a_pm1)}
        mat = np.zeros((len(a_pm1), len(a_p)))
        for j, sigma in enumerate(a_p):
            for i in range(len(sigma)):
                face = sigma[:i] + sigma[i+1:]
                if face in idx:
                    mat[idx[face], j] += (-1)**i
        boundary[p] = mat

    ranks = {}
    dims = {}
    for p in range(5):
        dims[p] = omega[p].shape[1] if omega[p].ndim == 2 else 0
    for p in range(1, 5):
        Om_p = omega[p]
        Om_pm1 = omega[p-1]
        if Om_p.ndim < 2 or Om_p.shape[1] == 0 or Om_pm1.ndim < 2 or Om_pm1.shape[1] == 0:
            ranks[p] = 0
            continue
        dp = Om_pm1.T @ boundary[p] @ Om_p
        S = np.linalg.svd(dp, compute_uv=False)
        ranks[p] = int(np.sum(S > tol))

    ker3 = dims[3] - ranks.get(3, 0)
    im4 = ranks.get(4, 0)
    return ker3 - im4

def main():
    print("=" * 70)
    print("SEESAW + LES APPROACH: structural constraints on beta_3=2")
    print("=" * 70)

    # At n=6 exhaustive: if beta_3=1, then T is strongly connected
    n = 6
    num_arcs = n*(n-1)//2
    n_total = 1 << num_arcs

    b3_sc = Counter()  # (beta_3, strongly_connected)
    b3_1_deletion_sc = Counter()  # for beta_3=1 tournaments: #deletions strongly connected

    print(f"\n  Exhaustive n={n} analysis...")
    for bits in range(n_total):
        if bits % 5000 == 0:
            print(f"    ... {bits}/{n_total}")
        T = tournament_from_bits(n, bits)
        sc = is_strongly_connected(T, n)
        try:
            b3 = compute_beta3(T, n)
        except:
            continue
        b3_sc[(b3, sc)] += 1

        if b3 == 1:
            # Check how many deletions are strongly connected
            sc_del_count = 0
            for v in range(n):
                remaining = [i for i in range(n) if i != v]
                T_sub = [[T[remaining[a]][remaining[b]] for b in range(n-1)] for a in range(n-1)]
                if is_strongly_connected(T_sub, n-1):
                    sc_del_count += 1
            b3_1_deletion_sc[sc_del_count] += 1

    print(f"\n  (beta_3, strongly_connected) distribution at n=6:")
    for key, cnt in sorted(b3_sc.items()):
        print(f"    beta_3={key[0]}, SC={key[1]}: {cnt}")

    print(f"\n  For beta_3=1 tournaments: #deletions that are strongly connected:")
    for key, cnt in sorted(b3_1_deletion_sc.items()):
        print(f"    {key}/6 SC: {cnt}")

    # Key question: at n=6, beta_3=1 → T strongly connected →
    # all deletions also strongly connected?
    # If so, then beta_3=2 at n=7 would require:
    # all 7 deletions have beta_3=1 → all 7 deletions strongly connected
    # → T must be "doubly strongly connected"

    # Now check the Euler characteristic constraint
    print(f"\n{'='*70}")
    print("EULER CHARACTERISTIC ANALYSIS")
    print("=" * 70)

    # For tournaments: chi = 1 - beta_1 + beta_2 - beta_3 + ...
    # Since beta_2 = 0: chi = 1 - beta_1 - beta_3 + beta_4 - ...
    # From seesaw: beta_1 * beta_3 = 0.
    #
    # If beta_3 = 2: chi = 1 - 0 - 2 + beta_4 - ... = -1 + beta_4 - ...
    # At n=7: chi ∈ {0, 1} for generic tournaments, but Paley has chi=7.
    #
    # Actually, let's check: what values does chi take at n=7?
    n = 7
    num_arcs = n*(n-1)//2
    n_total = 1 << num_arcs
    rng = np.random.RandomState(42)

    print(f"\n  Sampling chi at n={n}...")
    chi_dist = Counter()
    for trial in range(300):
        bits = rng.randint(0, n_total)
        T = tournament_from_bits(n, bits)

        # Need full Betti for chi
        tol = 1e-8
        all_paths = {}
        for p in range(0, 7):
            paths = []
            if p + 1 <= n:
                for verts in combinations(range(n), p+1):
                    for perm in permutations(verts):
                        if is_allowed_path(T, perm):
                            paths.append(perm)
            all_paths[p] = paths

        def build_omega(pd, mp):
            omega = {}
            for p in range(0, mp + 1):
                a_p = pd.get(p, [])
                if not a_p:
                    omega[p] = np.zeros((0, 0))
                    continue
                a_pm1_set = set(pd.get(p-1, [])) if p > 0 else set()
                na = {}
                for sigma in a_p:
                    for i in range(1, len(sigma)-1):
                        face = sigma[:i] + sigma[i+1:]
                        if p > 0 and face not in a_pm1_set:
                            if face not in na:
                                na[face] = len(na)
                if not na:
                    omega[p] = np.eye(len(a_p))
                else:
                    mat = np.zeros((len(na), len(a_p)))
                    for j, sigma in enumerate(a_p):
                        for i in range(1, len(sigma)-1):
                            face = sigma[:i] + sigma[i+1:]
                            if face in na:
                                mat[na[face], j] += (-1)**i
                    try:
                        U, S, Vt = np.linalg.svd(mat, full_matrices=True)
                    except:
                        omega[p] = np.eye(len(a_p))
                        continue
                    rank = int(np.sum(S > tol))
                    if len(a_p) - rank == 0:
                        omega[p] = np.zeros((len(a_p), 0))
                    else:
                        omega[p] = Vt[rank:].T
            return omega

        omega = build_omega(all_paths, 6)
        boundary = {}
        for p in range(1, 7):
            a_p = all_paths.get(p, [])
            a_pm1 = all_paths.get(p-1, [])
            if not a_p or not a_pm1:
                boundary[p] = np.zeros((0, 0))
                continue
            idx_map = {path: i for i, path in enumerate(a_pm1)}
            mat = np.zeros((len(a_pm1), len(a_p)))
            for j, sigma in enumerate(a_p):
                for i in range(len(sigma)):
                    face = sigma[:i] + sigma[i+1:]
                    if face in idx_map:
                        mat[idx_map[face], j] += (-1)**i
            boundary[p] = mat

        ranks = {}
        dims = {}
        for p in range(7):
            dims[p] = omega[p].shape[1] if omega[p].ndim == 2 else 0
        for p in range(1, 7):
            Om_p = omega.get(p, np.zeros((0,0)))
            Om_pm1 = omega.get(p-1, np.zeros((0,0)))
            if Om_p.ndim < 2 or Om_p.shape[1] == 0 or Om_pm1.ndim < 2 or Om_pm1.shape[1] == 0:
                ranks[p] = 0
                continue
            dp = Om_pm1.T @ boundary[p] @ Om_p
            try:
                S = np.linalg.svd(dp, compute_uv=False)
                ranks[p] = int(np.sum(S > tol))
            except:
                ranks[p] = 0

        betti = []
        for p in range(7):
            ker = dims[p] - ranks.get(p, 0)
            im_next = ranks.get(p+1, 0)
            betti.append(ker - im_next)

        chi = sum((-1)**p * betti[p] for p in range(len(betti)))
        chi_dist[chi] += 1

    print(f"  Chi distribution at n=7: {dict(sorted(chi_dist.items()))}")

if __name__ == '__main__':
    main()
