#!/usr/bin/env python3
"""
K_trace_formula.py -- Express K = c5 - 2*ov1 - 2*ov2 in terms of
matrix traces and show why it's constant for regular tournaments.

Key insight: tr(A^5)/5 counts directed 5-cycles. We need matrix
expressions for ov1 and ov2 in terms of A.

ov2 = sum_{i->j} C(lambda_{ij}, 2) where lambda_{ij} = (A^2)_{ji}
    = (1/2) sum_{i->j} [(A^2)_{ji}^2 - (A^2)_{ji}]

ov1 can be derived from: sum_v C(c3(v), 2) = ov1 + 2*ov2

We seek: K = f(tr(A^k), other trace quantities) that simplifies for regular A.

Author: kind-pasteur-2026-03-12-S60
"""

import numpy as np
from itertools import combinations
from collections import defaultdict
import random


def random_regular_tournament(n):
    assert n % 2 == 1
    m = (n - 1) // 2
    A = np.zeros((n, n), dtype=int)
    for i in range(n):
        for j in range(i+1, n):
            if random.random() < 0.5:
                A[i][j] = 1
            else:
                A[j][i] = 1
    for _ in range(10000):
        scores = A.sum(axis=1)
        if scores.max() == m and scores.min() == m:
            return A
        high = int(scores.argmax())
        low = int(scores.argmin())
        if high != low and A[high][low]:
            A[high][low] = 0
            A[low][high] = 1
        else:
            for j in range(n):
                if j != high and A[high][j]:
                    if scores[j] < m:
                        A[high][j] = 0
                        A[j][high] = 1
                        break
    scores = A.sum(axis=1)
    if all(s == m for s in scores):
        return A
    return None


def compute_all_quantities(A, n):
    """Compute all relevant trace and combinatorial quantities."""
    A2 = A @ A
    A3 = A @ A2
    A4 = A @ A3
    A5 = A @ A4

    # Basic traces
    trA3 = int(np.trace(A3))
    trA5 = int(np.trace(A5))
    c3 = trA3 // 3
    c5 = trA5 // 5

    # Score-related
    scores = A.sum(axis=1)
    S1 = int(scores.sum())
    S2 = int((scores**2).sum())
    S3 = int((scores**3).sum())
    S4 = int((scores**4).sum())

    # ov2 via lambda values
    # lambda_{ij} = (A^2)_{ji} for edge i->j
    ov2 = 0
    sum_lambda = 0
    sum_lambda_sq = 0
    for i in range(n):
        for j in range(n):
            if A[i][j]:
                lam = int(A2[j][i])
                sum_lambda += lam
                sum_lambda_sq += lam * lam
                ov2 += lam * (lam - 1) // 2

    # ov1 via vertex 3-cycle counts
    # c3(v) for each vertex
    c3v = np.zeros(n, dtype=int)
    for v in range(n):
        for a in range(n):
            if a == v:
                continue
            for b in range(a+1, n):
                if b == v:
                    continue
                if (A[v][a] and A[a][b] and A[b][v]) or (A[v][b] and A[b][a] and A[a][v]):
                    c3v[v] += 1

    sum_c3v = int(c3v.sum())
    # Each 3-cycle counted at 3 vertices: sum c3v = 3*c3
    assert sum_c3v == 3 * c3, f"sum_c3v={sum_c3v} != 3*c3={3*c3}"

    sum_c3v_choose2 = int(sum(cv * (cv - 1) // 2 for cv in c3v))
    # sum C(c3v, 2) = ov1 + 2*ov2
    ov1 = sum_c3v_choose2 - 2 * ov2

    K = c5 - 2 * ov1 - 2 * ov2

    # Additional trace quantities
    # Butterfly sum: sum_{i->j} (A^2)_{ji}^2
    butterfly = sum_lambda_sq  # same thing

    # sum_{i->j} (A^2)_{ji} = tr(A^3) = 3*c3
    assert sum_lambda == 3 * c3, f"sum_lambda={sum_lambda} != 3*c3={3*c3}"

    # Hadamard products (element-wise)
    # A ⊙ A^2^T : entry (i,j) = A[i][j] * A2[j][i] = indicator(i->j) * lambda_{ij}
    AoA2T = A * A2.T
    # sum of AoA2T = sum_lambda = 3*c3 (verified above)

    # A ⊙ (A2.T)^2 : entry (i,j) = A[i][j] * (A2[j][i])^2
    AoA2T_sq = A * (A2.T ** 2)
    sum_AoA2T_sq = int(AoA2T_sq.sum())  # = butterfly

    # tr(A * diag(A2.T)... hmm various combinations
    # Let's try: what combination of traces gives K?

    # Known: ov2 = (butterfly - 3*c3) / 2
    # Known: ov1 = sum_c3v_choose2 - 2*ov2

    # For REGULAR tournaments: c3(v) should be constant = 3*c3/n
    c3v_constant = all(cv == c3v[0] for cv in c3v)

    # Compute: sum C(c3v, 2) for regular
    if c3v_constant:
        c3v_val = int(c3v[0])
        sum_c3v_c2_formula = n * c3v_val * (c3v_val - 1) // 2
    else:
        c3v_val = None
        sum_c3v_c2_formula = None

    return {
        'n': n, 'c3': c3, 'c5': c5, 'K': K,
        'ov1': ov1, 'ov2': ov2,
        'trA3': trA3, 'trA5': trA5,
        'S2': S2, 'S3': S3, 'S4': S4,
        'sum_lambda': sum_lambda, 'butterfly': butterfly,
        'c3v': c3v, 'c3v_constant': c3v_constant, 'c3v_val': c3v_val,
        'sum_c3v_c2': sum_c3v_choose2,
    }


def main():
    print("=" * 70)
    print("K IN TERMS OF MATRIX TRACES")
    print("=" * 70)

    # ====== PART 1: Trace formulas for regular n=5 ======
    print("\n" + "=" * 60)
    print("PART 1: All regular n=5 tournaments")
    print("=" * 60)

    n = 5
    m = 2
    for bits in range(1 << 10):
        A = np.zeros((n, n), dtype=int)
        idx = 0
        for i in range(n):
            for j in range(i+1, n):
                if (bits >> idx) & 1:
                    A[i][j] = 1
                else:
                    A[j][i] = 1
                idx += 1
        scores = A.sum(axis=1)
        if not all(s == m for s in scores):
            continue

        r = compute_all_quantities(A, n)
        print(f"  bits={bits:>4}: c3={r['c3']}, c5={r['c5']}, "
              f"ov1={r['ov1']}, ov2={r['ov2']}, K={r['K']}, "
              f"butterfly={r['butterfly']}, c3v_const={r['c3v_constant']}")

    # ====== PART 2: Regular n=7 trace analysis ======
    print("\n" + "=" * 60)
    print("PART 2: Regular n=7 trace analysis")
    print("=" * 60)

    n = 7
    m = 3
    random.seed(42)

    K_expected = -3 * n * (n**2 - 1) * (n**2 - 9) // 320

    results_by_class = defaultdict(list)
    for trial in range(500):
        A = random_regular_tournament(n)
        if A is None:
            continue
        r = compute_all_quantities(A, n)
        key = (r['c5'], r['ov1'], r['ov2'])
        if len(results_by_class[key]) < 3:
            results_by_class[key].append(r)

    print(f"  K_expected = {K_expected}")
    print(f"\n  Classes found:")
    for key in sorted(results_by_class):
        r = results_by_class[key][0]
        print(f"\n    c5={r['c5']}, ov1={r['ov1']}, ov2={r['ov2']}, K={r['K']}")
        print(f"    butterfly={r['butterfly']}, c3v_const={r['c3v_constant']}, c3v_val={r['c3v_val']}")
        print(f"    Trace quantities: trA5={r['trA5']}")

    # ====== PART 3: Is c3(v) always constant for regular? ======
    print("\n" + "=" * 60)
    print("PART 3: Is c3(v) constant for regular tournaments?")
    print("=" * 60)

    # For regular tournament: is c3(v) = 3c3/n for ALL v?
    c3_total = n * (n**2 - 1) // 24  # = 14 at n=7
    c3v_expected = 3 * c3_total // n  # = 6 at n=7

    print(f"  n=7: c3={c3_total}, 3c3/n = {3*c3_total/n}")

    # Check: is 3c3/n an integer?
    # 3c3/n = 3*n(n^2-1)/24/n = (n^2-1)/8
    print(f"  (n^2-1)/8 = {(n**2-1)/8}")

    c3v_all_constant = True
    random.seed(123)
    for trial in range(100):
        A = random_regular_tournament(n)
        if A is None:
            continue
        r = compute_all_quantities(A, n)
        if not r['c3v_constant']:
            c3v_all_constant = False
            print(f"    Trial {trial}: c3(v) NOT constant! c3v = {r['c3v']}")
            break

    if c3v_all_constant:
        print(f"  c3(v) = {c3v_expected} for ALL regular n=7 tournaments (100 samples)")
        print(f"  This is a well-known result: in regular tournaments, every vertex")
        print(f"  is in the same number of 3-cycles = (n^2-1)/8")
    else:
        print(f"  c3(v) is NOT always constant!")

    # ====== PART 4: The algebraic identity ======
    print("\n" + "=" * 60)
    print("PART 4: Expressing K algebraically")
    print("=" * 60)

    # For regular tournaments with c3(v) = (n^2-1)/8 = L:
    # sum C(c3v, 2) = n * C(L, 2) = n * L*(L-1)/2
    # ov1 + 2*ov2 = n * L*(L-1)/2
    # ov2 = (butterfly - 3c3) / 2
    # ov1 = n*L*(L-1)/2 - 2*ov2 = n*L*(L-1)/2 - butterfly + 3c3
    # K = c5 - 2*ov1 - 2*ov2
    #   = c5 - 2*(n*L*(L-1)/2 - butterfly + 3c3) - 2*(butterfly - 3c3)/2
    #   = c5 - n*L*(L-1) + 2*butterfly - 6c3 - butterfly + 3c3
    #   = c5 + butterfly - n*L*(L-1) - 3c3
    #   = tr(A^5)/5 + butterfly - n*L*(L-1) - c3

    # Where butterfly = sum_{i->j} (A^2)_{ji}^2
    # = sum_{i,j} A[i][j] * (A^2)_{ji}^2

    # Wait, let me redo this carefully:
    # ov2 = (1/2) sum_{i->j} [lambda^2 - lambda] = (butterfly - sum_lambda)/2 = (butterfly - 3c3)/2
    # ov1 = sum_v C(c3v,2) - 2*ov2 = n*L(L-1)/2 - (butterfly - 3c3)
    # 2*ov1 = n*L(L-1) - 2*(butterfly - 3c3) = n*L(L-1) - 2*butterfly + 6c3
    # 2*ov2 = butterfly - 3c3
    # K = c5 - 2*ov1 - 2*ov2
    #   = c5 - [n*L(L-1) - 2*butterfly + 6c3] - [butterfly - 3c3]
    #   = c5 - n*L(L-1) + 2*butterfly - 6c3 - butterfly + 3c3
    #   = c5 + butterfly - n*L(L-1) - 3c3

    L = (n**2 - 1) // 8
    c3_val = n * (n**2 - 1) // 24
    constant_part = -n * L * (L - 1) - 3 * c3_val

    print(f"\n  For n=7: L = c3(v) = {L}, c3 = {c3_val}")
    print(f"  K = c5 + butterfly - n*L*(L-1) - 3*c3")
    print(f"  K = c5 + butterfly + ({constant_part})")
    print(f"  K_expected = {K_expected}")
    print(f"  So: c5 + butterfly = {K_expected - constant_part}")

    # Verify
    print(f"\n  Verification:")
    random.seed(77)
    for trial in range(10):
        A = random_regular_tournament(n)
        if A is None:
            continue
        r = compute_all_quantities(A, n)
        c5_plus_butterfly = r['c5'] + r['butterfly']
        K_from_formula = c5_plus_butterfly + constant_part
        print(f"    c5={r['c5']}, butterfly={r['butterfly']}, "
              f"c5+butterfly={c5_plus_butterfly}, K_formula={K_from_formula}, K_actual={r['K']}, "
              f"{'OK' if K_from_formula == r['K'] else 'FAIL'}")

    # ====== PART 5: Why is c5 + butterfly constant? ======
    print("\n" + "=" * 60)
    print("PART 5: Why is c5 + butterfly constant for regular tournaments?")
    print("=" * 60)

    # K = c5 + butterfly + constant
    # So we need: c5 + butterfly = constant for regular tournaments.
    #
    # c5 = tr(A^5)/5
    # butterfly = sum_{i,j} A[i][j] * (A^2[j][i])^2
    #           = sum_{i,j} A_{ij} * [sum_k A_{jk} A_{ki}]^2
    #           = sum_{i,j,k,l} A_{ij} A_{jk} A_{ki} A_{jl} A_{li}
    #
    # This is a sum over configurations:
    #   i->j->k->i and i->j->l->i (two 3-cycles through edge i->j)
    # which is exactly what ov2 counts (with the right normalization).
    #
    # But can we express butterfly in terms of traces?
    # butterfly = sum_{ij} A_{ij} [(A^2)_{ji}]^2
    # = sum_{ij} A_{ij} (A^T A^T)_{ij}^2... hmm
    # = tr(A * diag... )... no
    #
    # Actually: (A^2)_{ji} = (A^T)^2)_{ij} = (A^{T2})_{ij}
    # Wait no: (A^2)_{ji} = sum_k A_{jk}A_{ki} = sum_k (A^T)_{kj}(A^T)_{ik}
    # Hmm.

    # Let B = A^2. Then butterfly = sum_{ij} A_{ij} * B_{ji}^2.
    # This can be written as: sum_{ij} A_{ij} * (B ⊙ B)_{ji}
    # = tr(A * (B ⊙ B)^T) = tr(A * (B^T ⊙ B^T))
    # where ⊙ is Hadamard (element-wise) product.

    # So: butterfly = tr(A * (A^{2T} ⊙ A^{2T}))

    # Verify this matrix formula:
    random.seed(55)
    for trial in range(5):
        A = random_regular_tournament(n)
        if A is None:
            continue
        A2 = A @ A
        B = A2.T * A2.T  # Hadamard square of A^2 transposed
        butterfly_matrix = int(np.trace(A @ B))
        r = compute_all_quantities(A, n)
        print(f"  Trial {trial}: butterfly_direct={r['butterfly']}, butterfly_matrix={butterfly_matrix}, "
              f"{'OK' if butterfly_matrix == r['butterfly'] else 'FAIL'}")

    # ====== PART 6: c5 + butterfly as a single trace ======
    print("\n" + "=" * 60)
    print("PART 6: Can c5 + butterfly be expressed as a single trace formula?")
    print("=" * 60)

    # We have: c5 + butterfly constant for regular tournaments
    # c5 = tr(A^5)/5
    # butterfly = tr(A * (A^{2T} ⊙ A^{2T}))
    #
    # For regular A: A*J = m*J and J*A = m*J where J is all-ones.
    # Also A + A^T = J - I.
    # So A^T = J - I - A.
    # And A^2 = A*A, (A^2)^T = (A*A)^T = A^T * A^T = (J-I-A)^2 = ...
    #
    # This substitution should allow expressing everything in terms of A and J.

    # Let's try: for regular tournament, what are the eigenvalues of A?
    # A has eigenvalues m (for eigenvector J) and eigenvalues of the form
    # -1/2 + i*sqrt(n)/2 (for the other n-1 eigenvalues if DRT).
    # For general regular tournament, eigenvalues have real part -1/2.

    # For regular: A*1 = m*1 where 1 is all-ones vector.
    # A^T*1 = m*1 (since in-degree = m too).
    # So (A+A^T)*1 = 2m*1 = (n-1)*1. Check: A+A^T = J-I, (J-I)*1 = (n-1)*1. Yes.

    # The key identity for regular tournaments:
    # A^2 = A*A. (A^2)*1 = A*(A*1) = A*(m*1) = m*(A*1) = m^2*1.
    # So A^2 also has all-ones as eigenvector with eigenvalue m^2.

    # tr(A^5) = sum of 5th powers of eigenvalues
    # For regular: one eigenvalue is m, others are lambda_1,...,lambda_{n-1} with sum -m.
    # tr(A^5) = m^5 + sum lambda_i^5

    # For the butterfly: it involves Hadamard products, which don't simplify with eigenvalues.

    # Let me try a DIFFERENT approach: express everything in terms of A, J, I.

    # For regular: A*J = J*A = m*J. A + A^T = J - I.

    # tr(A^5) = tr(A^5)
    # We need a direct algebraic proof that tr(A^5)/5 + butterfly = const.

    # Let's check if butterfly = tr(something simple):
    random.seed(42)
    A = random_regular_tournament(n)
    r = compute_all_quantities(A, n)
    A2 = A @ A
    A3 = A @ A2
    A4 = A @ A3
    A5 = A @ A4
    AT = A.T
    AT2 = AT @ AT
    ATA = AT @ A
    AAT = A @ AT

    # Try various combinations:
    candidates = {
        'tr(A^3)': int(np.trace(A3)),
        'tr(A^5)': int(np.trace(A5)),
        'tr(A^2 * A^T)': int(np.trace(A2 @ AT)),
        'tr(A * A^T * A)': int(np.trace(A @ AT @ A)),
        'tr(A^2 * A^T^2)': int(np.trace(A2 @ AT2)),
        'tr((A*A^T)^2)': int(np.trace(AAT @ AAT)),
        'tr(A * A^T * A * A^T)': int(np.trace(A @ AT @ A @ AT)),
        'tr(A * A^T * A^T * A)': int(np.trace(A @ AT @ AT @ A)),
        'tr(A^3 * A^T)': int(np.trace(A3 @ AT)),
        'tr(A * A^T * A^2)': int(np.trace(A @ AT @ A2)),
        'tr(A^2 * A * A^T)': int(np.trace(A2 @ A @ AT)),
        'tr(A^2 * A^T * A)': int(np.trace(A2 @ AT @ A)),
        'sum(A ⊙ A2^T ⊙ A2^T)': int((A * A2.T * A2.T).sum()),
        'butterfly': r['butterfly'],
        'c5': r['c5'],
        'c5+butterfly': r['c5'] + r['butterfly'],
    }

    print(f"  For one regular n=7 tournament:")
    for name, val in candidates.items():
        print(f"    {name} = {val}")

    # Now check which are constant across regular tournaments:
    print(f"\n  Checking constancy across 50 regular n=7 tournaments:")
    values = defaultdict(set)
    random.seed(42)
    for trial in range(500):
        A = random_regular_tournament(n)
        if A is None:
            continue
        A2 = A @ A
        A3 = A @ A2
        A5 = A @ A3 @ A
        AT = A.T
        AT2 = AT @ AT
        AAT = A @ AT

        vals = {
            'tr(A^5)': int(np.trace(A5)),
            'tr(A^2*A^T^2)': int(np.trace(A2 @ AT2)),
            'tr((AAT)^2)': int(np.trace(AAT @ AAT)),
            'tr(A^3*A^T)': int(np.trace(A3 @ AT)),
            'butterfly': int((A * A2.T * A2.T).sum()),
            'c5+butterfly': int(np.trace(A5))//5 + int((A * A2.T * A2.T).sum()),
        }
        for name, val in vals.items():
            values[name].add(val)

    for name in values:
        v = values[name]
        status = "CONSTANT" if len(v) == 1 else f"VARIES ({len(v)} values)"
        print(f"    {name}: {sorted(v)} -- {status}")

    # ====== PART 7: Which trace expressions are constant? ======
    print("\n" + "=" * 60)
    print("PART 7: Systematic search for constant trace expressions")
    print("=" * 60)

    # For regular tournaments, which trace monomials are constant?
    # Level 1: tr(A^k) for k=1,...,5
    # Level 2: tr(products of A and A^T of total length 4 or 5)

    random.seed(42)
    tournaments = []
    for _ in range(500):
        A = random_regular_tournament(n)
        if A is not None:
            tournaments.append(A)
        if len(tournaments) >= 30:
            break

    print(f"  Testing {len(tournaments)} regular n=7 tournaments")

    def trace_expr(A, name):
        AT = A.T
        A2 = A @ A
        A3 = A @ A2
        if name == 'trA': return int(np.trace(A))
        if name == 'trA2': return int(np.trace(A2))
        if name == 'trA3': return int(np.trace(A3))
        if name == 'trA4': return int(np.trace(A3 @ A))
        if name == 'trA5': return int(np.trace(A3 @ A2))
        if name == 'trAAT': return int(np.trace(A @ AT))
        if name == 'trA2AT': return int(np.trace(A2 @ AT))
        if name == 'trAATA': return int(np.trace(A @ AT @ A))
        if name == 'trA2AT2': return int(np.trace(A2 @ AT @ AT))
        if name == 'trAAT2': return int(np.trace((A @ AT) @ (A @ AT)))
        if name == 'trA3AT': return int(np.trace(A3 @ AT))
        if name == 'trA2ATA': return int(np.trace(A2 @ AT @ A))
        if name == 'trAA2T': return int(np.trace(A @ A2.T))
        if name == 'butterfly': return int((A * (A2.T)**2).sum())
        if name == 'trA3AT2': return int(np.trace(A3 @ AT @ AT))
        if name == 'trA2ATA2': return int(np.trace(A2 @ AT @ A2))
        if name == 'trA2ATAT': return int(np.trace(A2 @ AT @ AT))
        return None

    expr_names = ['trA', 'trA2', 'trA3', 'trA4', 'trA5',
                  'trAAT', 'trA2AT', 'trAATA', 'trA2AT2', 'trAAT2',
                  'trA3AT', 'trA2ATA', 'trAA2T', 'butterfly',
                  'trA3AT2', 'trA2ATAT']

    expr_values = {name: set() for name in expr_names}
    expr_first = {}

    for A in tournaments:
        for name in expr_names:
            val = trace_expr(A, name)
            if val is not None:
                expr_values[name].add(val)
                if name not in expr_first:
                    expr_first[name] = val

    print(f"\n  Constant expressions (across all regular n=7):")
    for name in expr_names:
        v = expr_values[name]
        if len(v) == 1:
            print(f"    {name} = {list(v)[0]} -- CONSTANT")

    print(f"\n  Variable expressions:")
    for name in expr_names:
        v = expr_values[name]
        if len(v) > 1:
            print(f"    {name} = {sorted(v)} -- {len(v)} values")

    # Now check: is c5 + butterfly = trA5/5 + butterfly constant?
    combo_vals = set()
    for A in tournaments:
        A2 = A @ A
        A5 = A @ A2 @ A2
        c5 = int(np.trace(A5)) // 5
        butterfly = int((A * (A2.T)**2).sum())
        combo_vals.add(c5 + butterfly)

    print(f"\n  c5 + butterfly = {sorted(combo_vals)} -- {'CONSTANT' if len(combo_vals)==1 else 'VARIES'}")

    # Try: is tr(A^5) + 5*butterfly constant?
    combo2 = set()
    for A in tournaments:
        A2 = A @ A
        A5 = A @ A2 @ A2
        val = int(np.trace(A5)) + 5 * int((A * (A2.T)**2).sum())
        combo2.add(val)
    print(f"  tr(A^5) + 5*butterfly = {sorted(combo2)} -- {'CONSTANT' if len(combo2)==1 else 'VARIES'}")

    # ====== PART 8: Expanding butterfly in terms of standard traces ======
    print("\n" + "=" * 60)
    print("PART 8: Express butterfly as linear combination of traces")
    print("=" * 60)

    # butterfly = sum_{i,j,k,l} A_{ij}*A_{jk}*A_{ki}*A_{jl}*A_{li}
    # This is a degree-5 monomial in A entries.
    # Can it be a linear combination of tr(M1*M2*...) for products of A and A^T?

    # Build a system: for each tournament, compute traces and butterfly.
    # Find coefficients c_i such that butterfly = sum c_i * trace_i

    # Variable traces (level 5): trA5, and various mixed products
    # Let's identify all distinct trace monomials of degree 5 in {A, A^T}

    # Actually, butterfly involves 5 A-entries but NOT A^T. So it should be
    # expressible as tr() of products of A only... but which?

    # butterfly = sum_{i,j,k,l} A_{ij}A_{jk}A_{ki}A_{jl}A_{li}
    # Fix indices: sum over i,j,k,l.
    # = sum_j [sum_i A_{ij}A_{ki}A_{li}] * [sum_{} A_{jk}A_{jl}]... hmm that doesn't factor nicely.

    # Let's think of it graph-theoretically. butterfly counts:
    # Start at i, go to j (edge i->j)
    # From j, take two paths back to i: j->k->i and j->l->i
    # This is a "butterfly" rooted at i with center j, wings k and l (possibly k=l).

    # When k != l: this is two distinct 3-cycles through edge i->j.
    # When k = l: this is a single 3-cycle squared, counting lambda_{ij}.

    # butterfly = sum_{i,j: i->j} lambda_{ij}^2
    # = sum_{i,j: i->j} [lambda_{ij}*(lambda_{ij}-1) + lambda_{ij}]
    # = 2*ov2 + 3*c3

    # Verify!
    random.seed(42)
    for trial in range(5):
        A = random_regular_tournament(n)
        if A is None:
            continue
        r = compute_all_quantities(A, n)
        butterfly_from_ov2 = 2 * r['ov2'] + 3 * r['c3']
        print(f"  Trial {trial}: butterfly={r['butterfly']}, 2*ov2+3*c3={butterfly_from_ov2}, "
              f"{'OK' if r['butterfly'] == butterfly_from_ov2 else 'FAIL'}")

    # So: butterfly = 2*ov2 + 3*c3 (always, by definition)
    # Therefore: K = c5 + butterfly - n*L*(L-1) - 3*c3
    #            = c5 + 2*ov2 + 3*c3 - n*L*(L-1) - 3*c3
    #            = c5 + 2*ov2 - n*L*(L-1)

    print(f"\n  SIMPLIFIED: K = c5 + 2*ov2 - n*L*(L-1)")
    print(f"  For n=7: K = c5 + 2*ov2 - 7*6*5 = c5 + 2*ov2 - 210")

    random.seed(42)
    for trial in range(10):
        A = random_regular_tournament(n)
        if A is None:
            continue
        r = compute_all_quantities(A, n)
        K_simplified = r['c5'] + 2*r['ov2'] - n*L*(L-1)
        print(f"    c5={r['c5']}, ov2={r['ov2']}, c5+2*ov2={r['c5']+2*r['ov2']}, "
              f"K_simplified={K_simplified}, K_actual={r['K']}, "
              f"{'OK' if K_simplified == r['K'] else 'FAIL'}")

    # Wait: K = c5 - 2*ov1 - 2*ov2
    # And: K = c5 + 2*ov2 - nL(L-1)  (from butterfly substitution)
    # So: -2*ov1 - 2*ov2 = 2*ov2 - nL(L-1)
    # => -2*ov1 = 4*ov2 - nL(L-1)
    # => ov1 = nL(L-1)/2 - 2*ov2
    # And: ov1 + 2*ov2 = nL(L-1)/2 = sum C(c3v, 2) = n*C(L,2)
    # Check: nL(L-1)/2 = n*C(L,2). Yes!
    # So this is consistent with our definition.

    # So: K = c5 - 2*ov1 - 2*ov2 = c5 + 2*ov2 - n*L(L-1)
    # Both forms are equivalent via ov1 = n*C(L,2) - 2*ov2.
    #
    # But the SECOND form is more illuminating because:
    # K = c5 + 2*ov2 - const
    # So K constant iff c5 + 2*ov2 = const.
    #
    # In matrix terms:
    # c5 = tr(A^5)/5
    # ov2 = (butterfly - 3c3)/2
    # c5 + 2*ov2 = tr(A^5)/5 + butterfly - 3c3
    # And butterfly = tr(A*(A^{2T})^{⊙2})

    # So we need: tr(A^5)/5 + tr(A*(A^{2T})^{⊙2}) = const for regular A.
    # That is: tr(A^5) + 5*tr(A*(A^{2T})^{⊙2}) = 5*const + 15*c3

    print(f"\n  So the core identity is: c5 + 2*ov2 = constant for regular tournaments")
    print(f"  where ov2 = sum_{{i->j}} C(lambda_{{ij}}, 2) and lambda_{{ij}} = (A^2)_{{ji}}")

    # ====== PART 9: c5 + 2*ov2 in terms of traces ======
    print("\n" + "=" * 60)
    print("PART 9: Verify c5 + 2*ov2 = constant")
    print("=" * 60)

    for nn in [5, 7, 9]:
        mm = (nn - 1) // 2
        LL = (nn**2 - 1) // 8
        K_exp = -3 * nn * (nn**2 - 1) * (nn**2 - 9) // 320
        target = K_exp + nn * LL * (LL - 1)
        print(f"\n  n={nn}: L={LL}, n*L*(L-1)={nn*LL*(LL-1)}, K={K_exp}")
        print(f"  c5 + 2*ov2 should be {target}")

        random.seed(42)
        combo_vals = set()
        for _ in range(500 if nn <= 7 else 200):
            A = random_regular_tournament(nn)
            if A is None:
                continue
            A2 = A @ A
            c5 = int(np.trace(A @ A2 @ A2)) // 5
            ov2 = 0
            for i in range(nn):
                for j in range(nn):
                    if A[i][j]:
                        lam = int(A2[j][i])
                        ov2 += lam * (lam - 1) // 2
            combo_vals.add(c5 + 2 * ov2)

        print(f"  c5 + 2*ov2 = {sorted(combo_vals)} -- {'CONSTANT = '+str(list(combo_vals)[0]) if len(combo_vals)==1 else 'VARIES!'}")

    # Final: express in terms of tr(A^5) and another trace
    # c5 = tr(A^5)/5
    # 2*ov2 = butterfly - 3*c3 = sum_{ij} A_{ij}*(A2_{ji})^2 - 3*c3
    # = sum_{ij} A_{ij}*(A2_{ji})^2 - tr(A^3)
    #
    # So c5 + 2*ov2 = tr(A^5)/5 + sum A⊙(A^{2T})^2 - tr(A^3)
    #
    # For regular, this equals:
    # n*L*(L-1) + K = [specific number depending only on n]

    # Can we prove this algebraically?
    # Using A + A^T = J - I for tournaments...
    # A^T = J - I - A
    # A^2 = A*A
    # (A^2)^T = (J-I-A)^2 = J^2 - 2J + I + JA + AJ - A - AJ - JA + A^2
    #         = nJ - 2J + I - A + A^2 (using JA = AJ = mJ for regular)
    #         = (n-2)J + I - A + A^2
    #         Hmm wait: J^2 = nJ, JA = mJ, AJ = mJ.
    #         (J-I-A)^2 = J^2 - J*I - J*A - I*J + I^2 + I*A - A*J + A*I + A^2
    #         = nJ - J - mJ - J + I + A - mJ + A + A^2
    #         = (n-1-2m)J + I + 2A + A^2
    #         = 0*J + I + 2A + A^2 (since n=2m+1)
    #         = I + 2A + A^2
    # So: A^{2T} = I + 2A + A^2

    print(f"\n  For regular tournaments: A^{{2T}} = I + 2A + A^2")
    print(f"  Verification:")
    random.seed(42)
    A = random_regular_tournament(n)
    A2 = A @ A
    AT = A.T
    AT2 = AT @ AT
    check = np.eye(n, dtype=int) + 2*A + A2
    print(f"  A^{{2T}} == I + 2A + A^2: {np.array_equal(AT2, check)}")

    print("\nDONE.")


if __name__ == '__main__':
    main()
