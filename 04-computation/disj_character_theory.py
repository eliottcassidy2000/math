#!/usr/bin/env python3
"""
disj_character_theory.py -- kind-pasteur-2026-03-13-S60

CHARACTER-THEORETIC EXPLANATION of disj(k1,k2) = linear(c_k).

For circulant tournament on Z_p with connection set S:
  σ(d) = +1 if d ∈ S, -1 if p-d ∈ S.

A k-subset V = {v_1,...,v_k} supports a directed k-cycle iff
the subtournament on V has a Hamiltonian cycle.

For 3-cycles: {a,b,c} is active iff σ(b-a)*σ(c-b)*σ(a-c) = -1.

Define for each k-subset V: ind(V) = 1 if V supports a k-cycle, 0 otherwise.
Then c_k = sum_V ind(V).

disj(k1,k2) = sum_{V1∩V2=∅, |V1|=k1, |V2|=k2} ind(V1)*ind(V2)

KEY FORMULA: ind(V) depends on σ values, which can be decomposed:
  σ(d) = sum_t a_t * ω^(td)   (Fourier expansion of σ on Z_p)

where ω = e^{2πi/p} and a_t are the Fourier coefficients of σ.
Since σ(d) = -σ(p-d), σ is "odd": a_t = -a_{-t} (= purely imaginary part).

The eigenvalue D_t = Im(sum_{s∈S} ω^{st}) satisfies:
  σ(d) = (2/p) * sum_{t=1}^{m} D_t * sin(2π td/p)

So σ is a linear function of the D_t values.

For 3-cycles, ind(V) involves PRODUCT of σ values (degree 3 in σ).
For c_3, the sum of ind(V) over all 3-subsets gives a degree-3 polynomial in D_t.
By the power-sum / moment connection, this becomes a function of S_2, S_4, S_6, etc.

For disj(3,3), the "disjoint pair product" is degree 6 in σ (3+3).
If it were unconstrained, this would be a function of S_2,...,S_12 (up to degree 6 moments).

BUT: the constraint V1∩V2=∅ introduces correlations that REDUCE the effective degree.
Specifically, when V1 and V2 are disjoint, their combined degree-6 product factors
into LOWER-degree sums over the complementary vertex set.

Let me test: what degree moments determine disj(k1,k2)?
Hypothesis: disj(k1,k2) depends on moments up to degree max(k1,k2), NOT k1+k2.
"""

import cmath
import numpy as np
from itertools import combinations
from collections import defaultdict
from fractions import Fraction


def build_adj(p, S):
    S_set = set(S)
    A = [[0]*p for _ in range(p)]
    for i in range(p):
        for s in S_set:
            A[i][(i+s) % p] = 1
    return A


def count_ham_cycles(A, verts):
    k = len(verts)
    if k == 3:
        a, b, c = verts
        return (A[a][b]*A[b][c]*A[c][a]) + (A[a][c]*A[c][b]*A[b][a])
    start = 0
    dp = {(1 << start, start): 1}
    for mask in range(1, 1 << k):
        if not (mask & (1 << start)):
            continue
        for v in range(k):
            if not (mask & (1 << v)):
                continue
            key = (mask, v)
            if key not in dp or dp[key] == 0:
                continue
            cnt = dp[key]
            for w in range(k):
                if mask & (1 << w):
                    continue
                if A[verts[v]][verts[w]]:
                    nkey = (mask | (1 << w), w)
                    dp[nkey] = dp.get(nkey, 0) + cnt
    full = (1 << k) - 1
    total = 0
    for v in range(k):
        if v == start:
            continue
        key = (full, v)
        if key in dp and dp[key] > 0:
            if A[verts[v]][verts[start]]:
                total += dp[key]
    return total


def compute_moments(p, S):
    m = (p - 1) // 2
    omega = cmath.exp(2j * cmath.pi / p)
    D_vals = []
    for t in range(1, m + 1):
        lam = sum(omega ** (s * t) for s in S)
        D_vals.append(lam.imag)
    moments = {}
    for k in range(2, p + 1):
        moments[k] = sum(d**k for d in D_vals)
    return moments, D_vals


def test_moment_degree(p):
    """For each disj(k1,k2), find the MINIMUM set of moments needed."""
    m = (p - 1) // 2
    N = 1 << m

    print(f"\n{'='*70}")
    print(f"MOMENT DEGREE ANALYSIS at p={p}, m={m}, N={N}")
    print(f"{'='*70}")

    data = []
    for bits in range(N):
        S = []
        for j in range(m):
            if bits & (1 << j):
                S.append(j + 1)
            else:
                S.append(p - (j + 1))

        A = build_adj(p, S)
        moments, D_vals = compute_moments(p, S)

        # Enumerate cycles
        cycles_by_k = {}
        for k in range(3, p + 1, 2):
            cyc_list = []
            for subset in combinations(range(p), k):
                nc = count_ham_cycles(A, list(subset))
                for _ in range(nc):
                    cyc_list.append(frozenset(subset))
            if cyc_list:
                cycles_by_k[k] = cyc_list

        c_k = {k: len(v) for k, v in cycles_by_k.items()}

        # Disjoint pairs
        disj = {}
        for k1 in sorted(cycles_by_k):
            for k2 in sorted(cycles_by_k):
                if k1 > k2 or k1 + k2 > p:
                    continue
                cyc1 = cycles_by_k[k1]
                cyc2 = cycles_by_k[k2]
                count = 0
                if k1 == k2:
                    for i in range(len(cyc1)):
                        for j in range(i + 1, len(cyc1)):
                            if not (cyc1[i] & cyc1[j]):
                                count += 1
                else:
                    for c1 in cyc1:
                        for c2 in cyc2:
                            if not (c1 & c2):
                                count += 1
                disj[(k1, k2)] = count

        data.append({
            'bits': bits, 'c_k': c_k, 'disj': disj,
            'moments': moments, 'D_vals': D_vals
        })

    # For each disj(k1,k2), test which SINGLE moment S_{2j} determines it
    active_disj_types = [(k1, k2) for k1, k2 in sorted(set(
        k for d in data for k in d['disj'].keys()
    )) if any(d['disj'].get((k1, k2), 0) > 0 for d in data)
        for k1, k2 in [(k1, k2)]]  # This doesn't work, fix below

    active_disj_types = sorted(set(
        key for d in data for key in d['disj'].keys()
        if d['disj'][key] > 0
    ))

    all_moment_orders = list(range(2, p + 1))
    even_moment_orders = [k for k in all_moment_orders if k % 2 == 0]

    print(f"\n  Active disjoint pair types: {active_disj_types}")
    print(f"\n  Testing moment sufficiency:")

    for k1, k2 in active_disj_types:
        y = np.array([d['disj'].get((k1, k2), 0) for d in data], dtype=float)
        if len(set(y)) <= 1:
            print(f"  disj({k1},{k2}) = CONSTANT = {int(y[0])}")
            continue

        # Test: disj = f(S_4), f(S_4, S_6), f(S_4, S_6, S_8), ...
        for n_mom in range(1, len(even_moment_orders)):
            mom_indices = even_moment_orders[:n_mom]
            cols = [np.array([d['moments'][mi] for d in data], dtype=float)
                    for mi in mom_indices]
            X = np.column_stack(cols + [np.ones(N)])
            coeffs, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
            pred = X @ coeffs
            max_err = np.max(np.abs(y - pred))

            mom_str = ','.join(f'S{k}' for k in mom_indices)
            exact = max_err < 0.5
            if exact:
                coeff_str = ' + '.join(f'{c:.8f}*S{k}' for c, k in
                                       zip(coeffs[:-1], mom_indices))
                coeff_str += f' + {coeffs[-1]:.2f}'
                print(f"  disj({k1},{k2}) = f({mom_str}): {coeff_str}  "
                      f"err={max_err:.6f} EXACT")
                break
        else:
            print(f"  disj({k1},{k2}): NOT linear in any moment subset tested!")

    # KEY TEST: Is the moment degree for disj(k1,k2) = max(k1,k2) or k1+k2?
    print(f"\n  MOMENT DEGREE SUMMARY:")
    print(f"  {'type':>10} {'max(k)':>7} {'k1+k2':>7} {'min moments needed':>20}")

    for k1, k2 in active_disj_types:
        y = np.array([d['disj'].get((k1, k2), 0) for d in data], dtype=float)
        if len(set(y)) <= 1:
            print(f"  ({k1},{k2}):>10 {'N/A':>7} {'N/A':>7} {'CONSTANT':>20}")
            continue

        min_mom = 'not found'
        for n_mom in range(1, len(even_moment_orders)):
            mom_indices = even_moment_orders[:n_mom]
            cols = [np.array([d['moments'][mi] for d in data], dtype=float)
                    for mi in mom_indices]
            X = np.column_stack(cols + [np.ones(N)])
            coeffs, _, _, _ = np.linalg.lstsq(X, y, rcond=None)
            pred = X @ coeffs
            max_err = np.max(np.abs(y - pred))
            if max_err < 0.5:
                min_mom = f'S4..S{mom_indices[-1]}'
                break

        print(f"  ({k1},{k2}){' ':>5} {max(k1,k2):>7} {k1+k2:>7} {min_mom:>20}")

    # Also test: what is the CHARACTER EXPANSION of disj(3,3)?
    # disj(3,3) = sum_{disj V1,V2} ind(V1)*ind(V2)
    # ind(V) = (1 - prod_V)/2 where prod_V = product of σ values around the 3-cycle
    # So disj = (1/4) sum_{disj} (1-p1)(1-p2)
    #         = (1/4) [D - 2P + Q]
    # D = constant (number of disjoint 3-subset pairs)
    # P = sum_{disj(V1,V2)} prod(V1)
    # Q = sum_{disj(V1,V2)} prod(V1)*prod(V2)

    print(f"\n  CHARACTER DECOMPOSITION of disj(3,3):")

    # Compute prod(V) = σ(b-a)*σ(c-b)*σ(a-c mod p) for each 3-subset
    for d in data[:4]:  # First 4 orientations
        S = d['bits']
        S_set = set()
        for j in range(m):
            if S & (1 << j):
                S_set.add(j + 1)
            else:
                S_set.add(p - (j + 1))

        # σ function
        sigma = {}
        for j in range(1, p):
            sigma[j] = 1 if j in S_set else -1

        # Compute D, P, Q for 3-subsets
        triples = list(combinations(range(p), 3))
        prods = {}
        for V in triples:
            a, b, c = V
            d1 = (b - a) % p
            d2 = (c - b) % p
            d3 = (a - c) % p
            prods[frozenset(V)] = sigma[d1] * sigma[d2] * sigma[d3]

        # Count disjoint pairs
        D_count = 0
        P_sum = 0
        Q_sum = 0
        for i in range(len(triples)):
            V1 = frozenset(triples[i])
            for j in range(i + 1, len(triples)):
                V2 = frozenset(triples[j])
                if not (V1 & V2):
                    D_count += 1
                    P_sum += prods[V1] + prods[V2]
                    Q_sum += prods[V1] * prods[V2]

        disj_33 = (D_count - P_sum + Q_sum) // 4
        disj_actual = d['disj'].get((3, 3), 0)

        c3 = sum(1 for V, p_val in prods.items() if p_val == -1)

        # Total prod sum (over ALL 3-subsets)
        total_prod = sum(prods.values())
        # c3 = (C(p,3) - total_prod) / 2

        print(f"  bits={d['bits']:>2}: D={D_count}, P={P_sum}, Q={Q_sum}, "
              f"(D-P+Q)/4={disj_33}, actual={disj_actual}, "
              f"c3={c3}, total_prod={total_prod}")

    return data


# Run for p=7 and p=11
for p_val in [7, 11]:
    test_moment_degree(p_val)

print("\nDONE.")
