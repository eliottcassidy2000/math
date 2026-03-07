#!/usr/bin/env python3
"""
Special Evaluations of the Bivariate Tournament Eulerian Polynomial
====================================================================
Phi_T(x, y) = sum_P x^{fwd(P)} y^{bwd(P)}
            = sum_k a_k(T) x^k y^{n-1-k}

where a_k = #{permutations P of V with exactly k forward arcs in T}.

Investigate special evaluations and compare with I(Omega(T), z).

Author: opus-2026-03-07
"""

import sys, os, random
from itertools import permutations, combinations
from fractions import Fraction

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '03-artifacts', 'code'))
from tournament_lib import (
    all_tournaments, random_tournament, find_odd_cycles,
    conflict_graph, independence_poly_at, hamiltonian_path_count,
    independence_poly_at_fast,
)


# ---------------------------------------------------------------------------
# Core: compute a_k coefficients and Phi evaluations
# ---------------------------------------------------------------------------

def eulerian_coefficients(T):
    """Compute a_k(T) for k = 0, ..., n-1.
    a_k = number of vertex permutations with exactly k forward arcs."""
    n = len(T)
    counts = [0] * n
    for perm in permutations(range(n)):
        fwd = sum(1 for i in range(n - 1) if T[perm[i]][perm[i + 1]])
        counts[fwd] += 1
    return counts


def phi_eval(a_k, x, y):
    """Evaluate Phi_T(x, y) = sum_k a_k * x^k * y^{n-1-k}."""
    n = len(a_k)  # n-1+1 = n coefficients for degree n-1
    total = 0
    for k in range(n):
        total += a_k[k] * (x ** k) * (y ** (n - 1 - k))
    return total


def phi_eval_frac(a_k, x, y):
    """Evaluate with Fraction arithmetic for exact rational results."""
    n = len(a_k)
    total = Fraction(0)
    x, y = Fraction(x), Fraction(y)
    for k in range(n):
        total += Fraction(a_k[k]) * (x ** k) * (y ** (n - 1 - k))
    return total


# ---------------------------------------------------------------------------
# Independence polynomial coefficients
# ---------------------------------------------------------------------------

def indep_poly_coefficients(cycles):
    """Compute independence polynomial coefficients [alpha_0, alpha_1, ..., alpha_d]
    for the conflict graph Omega(T)."""
    m = len(cycles)
    if m == 0:
        return [1]

    vsets = [frozenset(c) for c in cycles]
    nbr = [0] * m
    for a in range(m):
        for b in range(a + 1, m):
            if vsets[a] & vsets[b]:
                nbr[a] |= 1 << b
                nbr[b] |= 1 << a

    max_size = len(set(v for c in cycles for v in c)) // 3 + 1
    counts = [0] * (max_size + 1)

    for mask in range(1 << m):
        ok = True
        seen = 0
        temp = mask
        sz = 0
        while temp:
            v = (temp & -temp).bit_length() - 1
            if nbr[v] & seen:
                ok = False
                break
            seen |= 1 << v
            temp &= temp - 1
            sz += 1
        if ok:
            if sz < len(counts):
                counts[sz] += 1
            else:
                counts.extend([0] * (sz - len(counts) + 1))
                counts[sz] = 1

    while len(counts) > 1 and counts[-1] == 0:
        counts.pop()
    return counts


def indep_poly_eval(alphas, z):
    """Evaluate I(Omega, z) = sum_k alpha_k * z^k."""
    return sum(a * z**k for k, a in enumerate(alphas))


# ---------------------------------------------------------------------------
# Analysis
# ---------------------------------------------------------------------------

def analyze_tournament(T, label=""):
    n = len(T)
    a_k = eulerian_coefficients(T)
    H = hamiltonian_path_count(T)
    cycles = find_odd_cycles(T)
    alphas = indep_poly_coefficients(cycles)

    # Verify basics
    import math
    assert a_k[n - 1] == H, f"a_{{n-1}} should be H(T): {a_k[n-1]} != {H}"
    assert sum(a_k) == math.factorial(n), f"sum(a_k)={sum(a_k)} != {n}!"

    # Compute I(Omega, z) for various z
    I_at = {}
    for z in [0, 1, 2, 3, -1, Fraction(1, 2)]:
        I_at[z] = indep_poly_eval(alphas, z)

    # Verify OCF: H(T) = I(Omega, 2)
    assert H == I_at[2], f"OCF fails: H={H}, I(Omega,2)={I_at[2]}"

    # Compute special Phi evaluations
    evals = {}

    # 1. Phi(x, 0) = a_{n-1} * x^{n-1} = H * x^{n-1}
    # At x=1: Phi(1, 0) = H
    evals['Phi(1,0)'] = phi_eval(a_k, 1, 0)
    assert evals['Phi(1,0)'] == H

    # 2. Phi(1, 1) = n!
    evals['Phi(1,1)'] = phi_eval(a_k, 1, 1)
    assert evals['Phi(1,1)'] == math.factorial(n)

    # 3. Phi(x, 1) = E_T(x) = sum_k a_k x^k  (the tournament Eulerian poly)
    # 4. Phi(1, y) = sum_k a_k y^{n-1-k}

    # 5. Phi(x, -x) = sum_k a_k x^k (-x)^{n-1-k} = (-1)^{n-1} x^{n-1} sum_k a_k (-1)^k
    alt_sum = sum(a_k[k] * ((-1) ** k) for k in range(n))
    evals['Phi(1,-1)'] = phi_eval(a_k, 1, -1)
    assert evals['Phi(1,-1)'] == alt_sum * ((-1) ** (n - 1))
    evals['alt_sum'] = alt_sum  # sum_k a_k (-1)^k

    # 6. Phi(x, x-1) on the line y = x - 1.
    #    When x-y = 1, the OCF decomposition simplifies: (x-y)^{n-1-f} = 1^{n-1-f} = 1
    #    So Phi(x, x-1) = A_n(x,x-1) + sum_I 2^parts(I) A_{f+1}(x,x-1) I(T)
    # At x = 1: Phi(1, 0) = H(T)  (already known)
    # At x = 3/2: Phi(3/2, 1/2)
    evals['Phi(3/2,1/2)'] = phi_eval_frac(a_k, Fraction(3, 2), Fraction(1, 2))
    evals['Phi(2,1)'] = phi_eval(a_k, 2, 1)  # = E_T(2)
    evals['Phi(0,-1)'] = phi_eval(a_k, 0, -1)  # only k=0 term survives = a_0 * (-1)^{n-1}

    # 7. Phi(2, -1): x - y = 3
    evals['Phi(2,-1)'] = phi_eval(a_k, 2, -1)

    # 8. Phi(x, 1/x) / x^{n-1} = sum_k a_k x^{k-(n-1-k)} / x^{n-1} ...
    # Actually: Phi(x, 1/x) = sum_k a_k x^k (1/x)^{n-1-k} = sum_k a_k x^{2k-(n-1)}
    # So Phi(x, 1/x) * x^{n-1} = sum_k a_k x^{2k}  (the even part)
    # At x=1: Phi(1, 1) / 1 = n!
    # At x = i (imaginary): gives alternating sum... skip imaginary for now

    # Palindromic at x=1: just n!

    # More evaluations along y = x - 1 line (where x - y = 1):
    for xv in [0, 1, 2, 3, -1, Fraction(1, 2)]:
        key = f'Phi({xv},{xv}-1)=Phi({xv},{xv-1})'
        evals[key] = phi_eval_frac(a_k, Fraction(xv), Fraction(xv - 1))

    # Phi at points targeting I(Omega, z) for z = 0, 1, 2, 3
    # We know Phi(1, 0) = H = I(Omega, 2).
    # Can we find (x, y) such that Phi(x, y) = I(Omega, 1)?
    # I(Omega, 1) = total number of independent sets
    # I(Omega, 0) = 1
    # I(Omega, 3) = sum alpha_k * 3^k

    # Additional structured evaluations
    evals['Phi(-1,1)'] = phi_eval(a_k, -1, 1)
    evals['Phi(0,1)'] = phi_eval(a_k, 0, 1)   # = a_0
    evals['Phi(1,-1)_v2'] = phi_eval(a_k, 1, -1)
    evals['Phi(-1,0)'] = phi_eval(a_k, -1, 0)  # = a_{n-1} * (-1)^{n-1}
    evals['Phi(2,0)'] = phi_eval(a_k, 2, 0)    # = a_{n-1} * 2^{n-1} = H * 2^{n-1}
    evals['Phi(3,1)'] = phi_eval(a_k, 3, 1)    # = E_T(3)
    evals['Phi(1,2)'] = phi_eval(a_k, 1, 2)    # sum a_k 2^{n-1-k}

    # Ratio-based: Phi(r, 1) / Phi(1, 0)  for various r
    # Phi(r, 1) = E_T(r) = sum a_k r^k
    E_at = {}
    for r in [0, 1, 2, 3, -1, Fraction(1, 2)]:
        E_at[r] = phi_eval_frac(a_k, Fraction(r), Fraction(1))

    return {
        'label': label,
        'n': n,
        'a_k': a_k,
        'H': H,
        'cycles': len(cycles),
        'alphas': alphas,
        'I_at': I_at,
        'evals': evals,
        'E_at': E_at,
        'alt_sum': alt_sum,
    }


def print_report(results):
    """Print analysis for a list of tournament results."""
    if not results:
        return

    n = results[0]['n']
    print(f"\n{'='*80}")
    print(f"BIVARIATE TOURNAMENT EULERIAN POLYNOMIAL: SPECIAL EVALUATIONS (n={n})")
    print(f"{'='*80}")
    print(f"Number of tournaments analyzed: {len(results)}")

    # Print each tournament's data
    for r in results:
        print(f"\n--- {r['label']} ---")
        print(f"  a_k = {r['a_k']}")
        print(f"  H(T) = {r['H']},  #odd_cycles = {r['cycles']}")
        print(f"  I(Omega) coeffs (alphas) = {r['alphas']}")
        print(f"  I(Omega, 0)=1, I(Omega, 1)={r['I_at'][1]}, I(Omega, 2)={r['I_at'][2]}, I(Omega, 3)={r['I_at'][3]}")

    # Tabulate key evaluations
    print(f"\n{'='*80}")
    print("KEY EVALUATIONS TABLE")
    print(f"{'='*80}")

    # Header
    key_evals = [
        'Phi(1,0)', 'Phi(1,1)', 'Phi(1,-1)', 'Phi(-1,1)', 'Phi(0,1)',
        'Phi(2,1)', 'Phi(1,2)', 'Phi(2,-1)', 'Phi(2,0)', 'Phi(3,1)',
        'Phi(3/2,1/2)', 'Phi(0,-1)', 'Phi(-1,0)',
    ]

    for ev in key_evals:
        vals = [r['evals'].get(ev, '?') for r in results]
        # Check if all same
        all_same = len(set(str(v) for v in vals)) == 1
        marker = " [CONST]" if all_same else ""
        vals_str = [str(v) for v in vals[:min(8, len(vals))]]
        print(f"  {ev:20s} = {', '.join(vals_str)}{marker}")

    # y = x - 1 line
    print(f"\n{'='*80}")
    print("EVALUATIONS ON THE LINE y = x - 1 (where x - y = 1)")
    print(f"{'='*80}")
    for xv in [0, 1, 2, 3, -1, Fraction(1, 2)]:
        key = f'Phi({xv},{xv}-1)=Phi({xv},{xv-1})'
        vals = [r['evals'].get(key, '?') for r in results]
        all_same = len(set(str(v) for v in vals)) == 1
        marker = " [CONST]" if all_same else ""
        vals_str = [str(v) for v in vals[:min(8, len(vals))]]
        print(f"  x={xv}: Phi({xv},{xv-1}) = {', '.join(vals_str)}{marker}")

    # E_T(r) = Phi(r, 1) = sum a_k r^k
    print(f"\n{'='*80}")
    print("E_T(r) = Phi(r, 1) — Tournament Eulerian polynomial at various r")
    print(f"{'='*80}")
    for r in [0, 1, 2, 3, -1, Fraction(1, 2)]:
        vals = [res['E_at'].get(r, '?') for res in results]
        all_same = len(set(str(v) for v in vals)) == 1
        marker = " [CONST]" if all_same else ""
        vals_str = [str(v) for v in vals[:min(8, len(vals))]]
        print(f"  E_T({r}) = {', '.join(vals_str)}{marker}")

    # Now the big question: does any Phi evaluation match I(Omega, z)?
    print(f"\n{'='*80}")
    print("MATCHING Phi_T(x, y) TO I(Omega(T), z)")
    print(f"{'='*80}")

    print("\n  Known: Phi(1, 0) = H(T) = I(Omega, 2)  [the OCF identity]")

    # For each evaluation, check if it equals I(Omega, z) for z in {0,1,2,3,-1,1/2}
    all_eval_keys = list(key_evals)
    for xv in [0, 1, 2, 3, -1, Fraction(1, 2)]:
        key = f'Phi({xv},{xv}-1)=Phi({xv},{xv-1})'
        all_eval_keys.append(key)

    i_targets = [0, 1, 2, 3, -1, Fraction(1, 2)]

    for ev in all_eval_keys:
        for z in i_targets:
            # Check if Phi eval == I(Omega, z) for ALL tournaments
            match = True
            for r in results:
                phi_val = r['evals'].get(ev, None)
                if phi_val is None:
                    match = False
                    break
                i_val = r['I_at'].get(z, None)
                if i_val is None:
                    match = False
                    break
                if phi_val != i_val:
                    match = False
                    break
            if match:
                print(f"  *** MATCH: {ev} = I(Omega, {z}) for ALL tournaments! ***")

    # Also check ratio relationships: Phi(...) / I(Omega, z) constant?
    print(f"\n{'='*80}")
    print("RATIO ANALYSIS: Phi(x,y) / I(Omega, z)")
    print(f"{'='*80}")

    for ev in all_eval_keys:
        for z in [1, 3]:
            ratios = []
            valid = True
            for r in results:
                phi_val = r['evals'].get(ev, None)
                i_val = r['I_at'].get(z, None)
                if phi_val is None or i_val is None or i_val == 0:
                    valid = False
                    break
                ratios.append(Fraction(phi_val) / Fraction(i_val))
            if valid and len(set(ratios)) == 1:
                print(f"  CONSTANT RATIO: {ev} / I(Omega, {z}) = {ratios[0]} for all tournaments")

    # Alternating sum analysis
    print(f"\n{'='*80}")
    print("ALTERNATING SUM: sum_k a_k(-1)^k")
    print(f"{'='*80}")
    for r in results:
        print(f"  {r['label']}: alt_sum = {r['alt_sum']}, H = {r['H']}, I(Omega,1) = {r['I_at'][1]}")

    # Check: is alt_sum always odd? always 1?
    alt_sums = [r['alt_sum'] for r in results]
    print(f"\n  All alt_sums: {alt_sums}")
    if all(a == alt_sums[0] for a in alt_sums):
        print(f"  -> CONSTANT = {alt_sums[0]}")
    if all(a % 2 == 1 for a in alt_sums):
        print(f"  -> Always odd")

    # a_0 analysis (permutations with 0 forward arcs = number of Hamiltonian paths in T^op)
    print(f"\n{'='*80}")
    print("a_0 (backward Hamiltonian paths) vs a_{n-1} (forward = H(T))")
    print(f"{'='*80}")
    for r in results:
        print(f"  {r['label']}: a_0 = {r['a_k'][0]}, a_{{n-1}} = {r['a_k'][n-1]}, sum = {r['a_k'][0] + r['a_k'][n-1]}")

    # Phi on diagonal x = y
    print(f"\n{'='*80}")
    print("Phi(x, x) = x^{n-1} * n!  (since sum a_k = n!)")
    print(f"{'='*80}")
    import math
    nf = math.factorial(n)
    for r in results:
        val = phi_eval(r['a_k'], 2, 2)
        expected = 2**(n-1) * nf
        print(f"  {r['label']}: Phi(2,2) = {val}, expected 2^{n-1}*{n}! = {expected}, match = {val == expected}")

    # Phi(1+t, 1-t) — symmetric perturbation around (1,1)
    # Phi(1+t, 1-t) = sum_k a_k (1+t)^k (1-t)^{n-1-k}
    # At t=0: n!.  Derivative at t=0?
    print(f"\n{'='*80}")
    print("Phi(1+t, 1-t) for small t  (x-y = 2t)")
    print(f"{'='*80}")
    for t_val in [Fraction(1, 2), 1, Fraction(1, 4)]:
        x_v = 1 + t_val
        y_v = 1 - t_val
        key = f"Phi({x_v},{y_v})"
        vals = []
        for r in results:
            v = phi_eval_frac(r['a_k'], Fraction(x_v), Fraction(y_v))
            vals.append(v)
        vals_str = [str(v) for v in vals[:min(8, len(vals))]]
        all_same = len(set(str(v) for v in vals)) == 1
        marker = " [CONST]" if all_same else ""
        print(f"  t={t_val}: {key} = {', '.join(vals_str)}{marker}")

    # KEY: derivative d/dt Phi(1+t, 1-t)|_{t=0}
    # = sum_k a_k [k (1+t)^{k-1}(1-t)^{n-1-k} - (n-1-k)(1+t)^k(1-t)^{n-2-k}]|_{t=0}
    # = sum_k a_k [k - (n-1-k)] = sum_k a_k (2k - n + 1)
    print(f"\n  d/dt Phi(1+t,1-t)|_{{t=0}} = sum_k a_k(2k-n+1):")
    for r in results:
        deriv = sum(r['a_k'][k] * (2*k - n + 1) for k in range(n))
        print(f"    {r['label']}: {deriv}")

    # Search for LINEAR relationships: Phi(x,y) = c1 * I(Omega, z1) + c0
    print(f"\n{'='*80}")
    print("SEARCHING FOR AFFINE MATCH: Phi(x,y) = c * I(Omega, z) + d")
    print(f"{'='*80}")

    # For each (Phi eval, z), check if linear relationship exists across all tournaments
    # Need at least 3 tournaments for this to be meaningful
    if len(results) >= 3:
        for ev in all_eval_keys:
            for z in [1, 3]:
                # Collect (I_val, Phi_val) pairs
                pairs = []
                valid = True
                for r in results:
                    phi_val = r['evals'].get(ev, None)
                    i_val = r['I_at'].get(z, None)
                    if phi_val is None or i_val is None:
                        valid = False
                        break
                    pairs.append((Fraction(i_val), Fraction(phi_val)))
                if not valid or len(pairs) < 2:
                    continue

                # Check if all points are collinear
                # (x1, y1), (x2, y2): slope = (y2-y1)/(x2-x1) if x2 != x1
                x1, y1 = pairs[0]
                x2, y2 = pairs[1]
                if x1 == x2:
                    if all(p[0] == x1 for p in pairs):
                        # All same I value - degenerate
                        continue
                    else:
                        continue

                slope = (y2 - y1) / (x2 - x1)
                intercept = y1 - slope * x1
                all_collinear = True
                for xi, yi in pairs[2:]:
                    if yi != slope * xi + intercept:
                        all_collinear = False
                        break
                if all_collinear and slope != 0:
                    print(f"  AFFINE: {ev} = {slope} * I(Omega, {z}) + {intercept}")

    # Final: look for ANY (a, b) such that Phi(a, b) = I(Omega, z)
    # Try a systematic grid — use floats for speed, confirm with Fraction
    print(f"\n{'='*80}")
    print("SYSTEMATIC GRID SEARCH: Phi(a, b) = I(Omega, z)?")
    print(f"{'='*80}")
    grid_half = [a / 2.0 for a in range(-4, 9)]  # -2, -1.5, ..., 4
    found_any = False
    for z in [0, 1, 2, 3, -1]:
        for xa in grid_half:
            for yb in grid_half:
                match = True
                for r in results:
                    phi_v = phi_eval(r['a_k'], xa, yb)
                    i_v = r['I_at'][z]
                    if abs(phi_v - i_v) > 0.01:
                        match = False
                        break
                if match:
                    # Confirm with exact arithmetic
                    exact_match = True
                    xa_f, yb_f = Fraction(int(xa * 2), 2), Fraction(int(yb * 2), 2)
                    for r in results:
                        phi_v = phi_eval_frac(r['a_k'], xa_f, yb_f)
                        if phi_v != Fraction(r['I_at'][z]):
                            exact_match = False
                            break
                    if not exact_match:
                        continue
                    # Skip known
                    if int(xa*2) == 2 and int(yb*2) == 0 and z == 2:
                        continue
                    print(f"  *** FOUND: Phi({xa_f}, {yb_f}) = I(Omega, {z}) for ALL {len(results)} tournaments ***")
                    found_any = True
    if not found_any:
        print("  No new matches found on the grid (besides Phi(1,0) = I(Omega,2)).")

    # Check: Phi(a, b) always an integer on integer grid? Find constants.
    print(f"\n{'='*80}")
    print("INTEGRALITY & CONSTANTS on integer grid")
    print(f"{'='*80}")
    for xa in range(-2, 5):
        for yb in range(-2, 5):
            val_set = set(phi_eval(r['a_k'], xa, yb) for r in results)
            if len(val_set) == 1:
                print(f"  Phi({xa},{yb}) = {val_set.pop()} [CONST across all T]")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    import math

    # ---- n = 5: exhaustive ----
    print("Computing n=5 (exhaustive: all 2^10 = 1024 tournaments)...")
    results_5 = []
    seen_a = set()  # deduplicate by a_k profile
    for idx, T in enumerate(all_tournaments(5)):
        a_k = eulerian_coefficients(T)
        key = tuple(a_k)
        if key not in seen_a:
            seen_a.add(key)
            r = analyze_tournament(T, label=f"n5-type{len(seen_a)}")
            results_5.append(r)

    print(f"  Found {len(results_5)} distinct a_k profiles out of 1024 tournaments")
    print_report(results_5)

    # ---- n = 7: 10 random ----
    print(f"\n\nComputing n=7 (10 random tournaments)...")
    rng = random.Random(42)
    results_7 = []
    for idx in range(10):
        T = random_tournament(7, rng=rng)
        r = analyze_tournament(T, label=f"n7-rand{idx}")
        results_7.append(r)

    print_report(results_7)

    # ---- Summary of findings ----
    print(f"\n\n{'='*80}")
    print("SUMMARY OF FINDINGS")
    print(f"{'='*80}")

    print("""
1. Phi(1, 0) = H(T) = I(Omega, 2)   [the OCF — always holds]
2. Phi(1, 1) = n!                     [trivial — all perms contribute 1]
3. Phi(x, x) = x^{n-1} * n!          [trivial — factors out]
4. Phi(0, 1) = a_0 = H(T^{op})       [backward Hamiltonian paths]
5. Alt sum = sum_k a_k(-1)^k          [see values above]

KEY QUESTION RESULTS:
- Grid search for Phi(a, b) = I(Omega, z) beyond z=2 is reported above.
- Affine relationships Phi(...) = c * I(Omega, z) + d are reported above.
""")

    # Extra: for n=5, print the FULL table of all distinct types
    print(f"\n{'='*80}")
    print("FULL n=5 TYPE TABLE")
    print(f"{'='*80}")
    print(f"{'Type':>10s} | {'a_k':>30s} | {'H':>5s} | {'alphas':>20s} | I(O,1) | I(O,3) | alt_sum")
    print("-" * 110)
    for r in results_5:
        print(f"{r['label']:>10s} | {str(r['a_k']):>30s} | {r['H']:>5d} | {str(r['alphas']):>20s} | {r['I_at'][1]:>6} | {r['I_at'][3]:>6} | {r['alt_sum']:>6}")


if __name__ == "__main__":
    main()
