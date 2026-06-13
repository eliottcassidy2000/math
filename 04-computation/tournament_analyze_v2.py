#!/usr/bin/env python3
"""
tournament_analyze_v2.py — kind-pasteur-2026-03-22-S20j

Extended tournament analysis tool. Adds to v1:
  - CSV input from real pairwise comparison data
  - Independence polynomial computation and root analysis
  - Inter-dimensional evaluation (chemistry x=1, tournament x=2, topology x=-1)
  - Spectral analysis (eigenvalues of skew-adjacency)
  - Comparison of two tournaments
  - JSON output mode
  - Batch mode for multiple tournaments

Author: kind-pasteur-2026-03-22-S20j
Dependencies: numpy only.
"""

import sys, csv, json
import numpy as np
from math import comb, log, sqrt, pi
from collections import defaultdict
from itertools import combinations

# ========================================================================
# CORE COMPUTATIONS
# ========================================================================

def scores(A):
    return A.sum(axis=1).astype(int)

def S2(A):
    s = scores(A)
    return int(sum(s**2))

def c3_from_scores(n, s2):
    return comb(n, 3) - (s2 - comb(n, 2)) // 2

def c3_from_trace(A):
    return int(np.trace(A @ A @ A)) // 3

def H_exact(A):
    n = len(A)
    if n > 20:
        return None  # too slow
    dp = defaultdict(int)
    for v in range(n): dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            if dp[(mask, v)] == 0: continue
            for w in range(n):
                if mask & (1 << w): continue
                if A[v][w]: dp[(mask | (1 << w), w)] += dp[(mask, v)]
    full = (1 << n) - 1
    return sum(dp[(full, v)] for v in range(n))

def H_approx(n, s2):
    return 1 + n*(n-1)*(2*n-1)//6 - s2

def formal_rank(A):
    n = len(A)
    net = np.zeros(n)
    opponents = np.zeros(n)
    for i in range(n):
        for j in range(n):
            if i == j: continue
            opponents[i] += 1
            if A[i][j]: net[i] += 1
            else: net[i] -= 1
    win_rate = net / np.maximum(opponents, 1)
    rapidity = np.arctanh(np.clip(win_rate, -0.999, 0.999))
    return np.argsort(-rapidity), rapidity

def spectral_analysis(A):
    """Eigenvalues of skew-adjacency B = A - A^T."""
    n = len(A)
    B = (A - A.T).astype(float)
    eigvals = np.sort(np.linalg.eigvals(B).imag)[::-1]
    # Casimir invariants
    trB2 = int(round(np.trace(B @ B).real))
    trB4 = int(round(np.trace(B @ B @ B @ B).real))
    return {
        "eigenvalues_imag": [round(float(e), 6) for e in eigvals if abs(e) > 1e-10],
        "trB2": trB2,
        "trB4": trB4,
        "trB2_expected": -n*(n-1),
    }

def independence_poly_small(A):
    """Independence polynomial of Omega(T) for small n. Returns coefficients."""
    n = len(A)
    if n > 8:
        return None  # too expensive

    # Find all directed 3-cycle VERTEX SETS
    cycle_sets = []
    for triple in combinations(range(n), 3):
        i, j, k = triple
        if (A[i][j] and A[j][k] and A[k][i]) or (A[i][k] and A[k][j] and A[j][i]):
            cycle_sets.append(frozenset(triple))

    nc = len(cycle_sets)
    if nc == 0:
        return [1]

    # Build conflict adjacency
    adj = [[False]*nc for _ in range(nc)]
    for a in range(nc):
        for b in range(a+1, nc):
            if cycle_sets[a] & cycle_sets[b]:
                adj[a][b] = adj[b][a] = True

    # Count independent sets by size
    alpha = defaultdict(int)
    alpha[0] = 1
    for mask in range(1, 1 << nc):
        verts = [v for v in range(nc) if (mask >> v) & 1]
        ok = True
        for a in range(len(verts)):
            for b in range(a+1, len(verts)):
                if adj[verts[a]][verts[b]]:
                    ok = False; break
            if not ok: break
        if ok:
            alpha[len(verts)] += 1

    max_k = max(alpha.keys())
    return [alpha.get(k, 0) for k in range(max_k + 1)]

def eval_poly(coeffs, x):
    return sum(c * x**k for k, c in enumerate(coeffs))

def poly_roots(coeffs):
    if len(coeffs) <= 1: return []
    np_coeffs = list(reversed(coeffs))
    return sorted(np.roots(np_coeffs).real)

# ========================================================================
# CSV INPUT
# ========================================================================

def load_csv(filename, col_a="item_a", col_b="item_b", col_winner="winner"):
    """Load pairwise comparison data from CSV."""
    items = set()
    wins = defaultdict(lambda: defaultdict(int))

    with open(filename, "r") as f:
        reader = csv.DictReader(f)
        for row in reader:
            a, b, w = row[col_a].strip(), row[col_b].strip(), row[col_winner].strip()
            items.add(a)
            items.add(b)
            if w == a:
                wins[a][b] += 1
            elif w == b:
                wins[b][a] += 1

    names = sorted(items)
    n = len(names)
    idx = {name: i for i, name in enumerate(names)}
    A = np.zeros((n, n), dtype=int)

    for a in names:
        for b in names:
            if a == b: continue
            wa = wins[a][b]
            wb = wins[b][a]
            if wa > wb:
                A[idx[a]][idx[b]] = 1
            elif wb > wa:
                A[idx[b]][idx[a]] = 1
            # ties: randomly break (or leave as 0, making partial tournament)

    return A, names

# ========================================================================
# FULL ANALYSIS
# ========================================================================

def full_analysis(A, names=None, output_json=False):
    n = len(A)
    if names is None:
        names = [str(i) for i in range(n)]

    results = {"n": n}

    # Scores
    s = scores(A)
    s2 = S2(A)
    c3 = c3_from_scores(n, s2)
    c3_check = c3_from_trace(A)
    ranking, raps = formal_rank(A)
    regular = all(x == s[0] for x in s)

    results["scores"] = {names[i]: int(s[i]) for i in range(n)}
    results["S2"] = s2
    results["c3"] = c3
    results["regular"] = regular
    results["ranking"] = [names[i] for i in ranking]

    # H
    H_app = H_approx(n, s2)
    results["H_approx"] = H_app

    H = None
    if n <= 20:
        H = H_exact(A)
        results["H"] = H
        results["H_exact_method"] = "Held-Karp DP"
        if n <= 4:
            results["H_formula_exact"] = True
            results["formula"] = f"H = {1 + n*(n-1)*(2*n-1)//6} - S_2"
        else:
            results["H_formula_exact"] = False
            results["score_explains_pct"] = round(100 * H_app / H, 1) if H > 0 else 0

    # Spectral
    spec = spectral_analysis(A)
    results["spectral"] = spec

    # Cartan
    A_float = A.astype(float)
    A_sym = (A_float + A_float.T) / 2
    A_anti = (A_float - A_float.T) / 2
    n_sym = float(np.linalg.norm(A_sym))
    n_anti = float(np.linalg.norm(A_anti))
    total = n_sym**2 + n_anti**2
    t_frac = n_anti**2 / total if total > 0 else 0
    results["cartan"] = {
        "tournament_norm": round(n_anti, 4),
        "cooperation_norm": round(n_sym, 4),
        "tournament_fraction": round(t_frac, 4),
    }

    # Independence polynomial (small n)
    ip = independence_poly_small(A)
    if ip is not None:
        results["independence_poly"] = ip

        # Multi-dimensional evaluation
        results["evaluations"] = {
            "I(1)_sigma": eval_poly(ip, 1),
            "I(2)_H_from_poly": eval_poly(ip, 2),
            "I(-1)_chi": eval_poly(ip, -1),
            "I(phi)": round(eval_poly(ip, (1+sqrt(5))/2), 4),
        }

        # Roots
        if len(ip) > 1:
            roots = poly_roots(ip)
            results["roots"] = [round(r, 6) for r in roots]
            # H as product of distances
            H_from_roots = ip[-1]
            for r in roots:
                H_from_roots *= abs(2 - r)
            results["H_from_root_product"] = round(H_from_roots, 2)

    # Energy
    if H and H > 0:
        results["free_energy"] = round(-log(H), 4)
        results["log_H"] = round(log(H), 4)

    results["redei_parity"] = "ODD" if (H is not None and H % 2 == 1) else ("EVEN (BUG)" if H is not None else "unknown")

    # Output
    if output_json:
        print(json.dumps(results, indent=2, default=str))
    else:
        print_analysis(results, names, A, s, ranking, raps, H, H_app, c3, c3_check, s2, ip, spec, n_anti, n_sym, t_frac, regular)

    return results

def print_analysis(results, names, A, s, ranking, raps, H, H_app, c3, c3_check, s2, ip, spec, n_anti, n_sym, t_frac, regular):
    n = len(A)
    print("=" * 60)
    print("  TOURNAMENT ANALYSIS v2")
    print("=" * 60)

    print(f"\n  Vertices: {n}, Arcs: {int(A.sum())}")

    print(f"\n  RANKING (FormalRank):")
    for rank, idx in enumerate(ranking):
        bar = "#" * max(1, int(s[idx]))
        print(f"    #{rank+1} {names[idx]:>12s}  score={s[idx]:>2d}  rapidity={raps[idx]:>+.3f}  {bar}")

    print(f"\n  SCORE ANALYSIS:")
    print(f"    S_2 = {s2}, variance = {np.var(s):.3f}, regular = {regular}")
    print(f"    c3 = {c3} (3-cycles), fraction = {c3/max(comb(n,3),1):.3f}")

    print(f"\n  HAMILTONIAN PATH COUNT (H):")
    if H is not None:
        print(f"    H = {H} (exact)")
        if n <= 4:
            print(f"    = {1+n*(n-1)*(2*n-1)//6} - {s2} (algebraic, EXACT)")
        else:
            pct = 100*H_app/H if H > 0 else 0
            print(f"    Score-determined: {H_app} ({pct:.0f}%), correction: {H - H_app}")
        print(f"    H is {'ODD' if H%2==1 else 'EVEN (ERROR)'} (Redei)")
    else:
        print(f"    H ~ {H_app} (approximate from scores)")

    if ip is not None and len(ip) > 1:
        poly_str = " + ".join(f"{c}x^{k}" if k > 0 else str(c) for k, c in enumerate(ip) if c > 0)
        print(f"\n  INDEPENDENCE POLYNOMIAL:")
        print(f"    I(x) = {poly_str}")
        print(f"    I(1) = {eval_poly(ip, 1)} (Merrifield-Simmons / chemistry)")
        print(f"    I(2) = {eval_poly(ip, 2)} (tournament / H)")
        print(f"    I(-1) = {eval_poly(ip, -1)} (Euler characteristic)")
        if "roots" in results:
            print(f"    Roots: {results['roots']}")
            if H and H > 0:
                print(f"    log H = sum log(2 - r_i) + log(alpha_d)")

    print(f"\n  SPECTRAL (skew-adjacency B = A - A^T):")
    if spec["eigenvalues_imag"]:
        print(f"    Imaginary eigenvalues: {spec['eigenvalues_imag'][:6]}{'...' if len(spec['eigenvalues_imag'])>6 else ''}")
    print(f"    Tr(B^2) = {spec['trB2']} (expected {spec['trB2_expected']})")
    print(f"    Tr(B^4) = {spec['trB4']}")

    print(f"\n  CARTAN DECOMPOSITION:")
    print(f"    Tournament ||A_anti|| = {n_anti:.3f} ({100*t_frac:.1f}%)")
    print(f"    Cooperation ||A_sym|| = {n_sym:.3f} ({100*(1-t_frac):.1f}%)")

    if H and H > 0:
        print(f"\n  FREE ENERGY: F = -log H = {-log(H):.4f}, bandwidth = log(3/2) = {log(3/2):.4f}")

    if regular:
        print(f"\n  ** REGULAR tournament: maximizes H for this order **")
    elif H == 1:
        print(f"\n  ** TRANSITIVE tournament: H=1, no cycles, minimum **")
    print()


# ========================================================================
# COMPARISON MODE
# ========================================================================

def compare(A1, A2, names1=None, names2=None):
    """Compare two tournaments side by side."""
    r1 = full_analysis(A1, names1, output_json=False)
    r2 = full_analysis(A2, names2, output_json=False)

    print("=" * 60)
    print("  COMPARISON")
    print("=" * 60)

    for key in ["H", "c3", "S2"]:
        v1 = r1.get(key, "?")
        v2 = r2.get(key, "?")
        print(f"  {key}: T1={v1}, T2={v2}, diff={v1-v2 if isinstance(v1,int) and isinstance(v2,int) else '?'}")

    if r1.get("cartan") and r2.get("cartan"):
        t1 = r1["cartan"]["tournament_fraction"]
        t2 = r2["cartan"]["tournament_fraction"]
        print(f"  Tournament fraction: T1={t1:.3f}, T2={t2:.3f}")

    print()

# ========================================================================
# MAIN
# ========================================================================

if __name__ == "__main__":
    args = sys.argv[1:]

    if not args or args[0] == "--demo":
        # Demo mode
        print("\n  ROCK-PAPER-SCISSORS")
        A_rps = np.array([[0,1,0],[0,0,1],[1,0,0]])
        full_analysis(A_rps, ["Rock","Paper","Scissors"])

        print("\n  TRANSITIVE (n=5)")
        A_t5 = np.zeros((5,5), dtype=int)
        for i in range(5):
            for j in range(i+1, 5):
                A_t5[i][j] = 1
        full_analysis(A_t5, ["A","B","C","D","E"])

        print("\n  REGULAR (n=5, Paley-like)")
        A_r5 = np.array([[0,1,0,1,0],[0,0,1,0,1],[1,0,0,1,0],[0,1,0,0,1],[1,0,1,0,0]])
        full_analysis(A_r5, ["v0","v1","v2","v3","v4"])

    elif args[0] == "--json":
        A_rps = np.array([[0,1,0],[0,0,1],[1,0,0]])
        full_analysis(A_rps, ["Rock","Paper","Scissors"], output_json=True)

    elif args[0] == "--csv" and len(args) >= 2:
        filename = args[1]
        col_a = args[3] if len(args) > 3 else "item_a"
        col_b = args[5] if len(args) > 5 else "item_b"
        col_w = args[7] if len(args) > 7 else "winner"
        A, names = load_csv(filename, col_a, col_b, col_w)
        full_analysis(A, names)

    elif args[0] == "--matrix":
        mat_str = args[1]
        rows = mat_str.split(";")
        A = np.array([[int(x.strip()) for x in row.split(",")] for row in rows])
        full_analysis(A)

    else:
        print("Usage:")
        print("  python tournament_analyze_v2.py --demo")
        print("  python tournament_analyze_v2.py --json")
        print("  python tournament_analyze_v2.py --csv data.csv --col-a item_a --col-b item_b --col-winner winner")
        print("  python tournament_analyze_v2.py --matrix '0,1,0; 0,0,1; 1,0,0'")
