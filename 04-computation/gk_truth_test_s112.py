#!/usr/bin/env python3
"""
gk_truth_test_s112.py — Which g_k is correct? Transfer matrix (degree k) vs opus (degree 3)?
kind-pasteur-2026-03-15-S112

Test: compute W(n)/n! - 1 = CV^2 both directly and via each g_k family.
The correct one matches. The other is wrong.
"""

from fractions import Fraction
from math import factorial, comb
from itertools import permutations

def compute_W_direct(n):
    """W(n) = sum over all permutations sigma of 2^{adj1(sigma)} * NUD(sigma)
    where adj1 = number of unit ascents, NUD = no unit descent indicator.
    Actually: W(n)/n! = E[prod_{j=0}^{n-2} (1+Z_j)]."""
    nfact = factorial(n)
    total = Fraction(0)
    for perm in permutations(range(n)):
        # Check NUD: no position where sigma(j+1) = sigma(j) - 1
        nud = True
        adj1 = 0
        for j in range(n-1):
            if perm[j+1] == perm[j] - 1:
                nud = False
                break
            if perm[j+1] == perm[j] + 1:
                adj1 += 1
        if nud:
            total += Fraction(2**adj1)
    return total  # This is W(n)

def compute_W_from_product(n):
    """W(n)/n! = E[prod (1+Z_j)]"""
    nfact = factorial(n)
    total = Fraction(0)
    for perm in permutations(range(n)):
        prod = Fraction(1)
        skip = False
        for j in range(n-1):
            x = 1 if perm[j+1] == perm[j] + 1 else 0
            y = 1 if perm[j+1] == perm[j] - 1 else 0
            z = x - y
            prod *= (1 + z)
            if prod == 0:
                skip = True
                break
        if not skip:
            total += prod
    return total  # = W(n)

def falling_factorial(n, k):
    r = Fraction(1)
    for i in range(k):
        r *= (n - i)
    return r

# Transfer matrix g_k computation
def transfer_matrix_gk_all(n, k_max=None):
    """Compute g_k(n-2k) for all k via transfer matrix."""
    num_edges = n - 2
    if k_max is None:
        k_max = (n - 1) // 2

    # State: [A, B, C] as polynomials in x
    max_deg = k_max
    state = [[Fraction(1)] + [Fraction(0)]*max_deg,
             [Fraction(0)]*(max_deg+1),
             [Fraction(0)]*(max_deg+1)]

    def poly_add(a, b):
        result = [Fraction(0)] * max(len(a), len(b))
        for i in range(len(a)):
            result[i] += a[i]
        for i in range(len(b)):
            result[i] += b[i]
        return result[:max_deg+1]

    def poly_shift(a):
        return [Fraction(0)] + a[:max_deg]

    def poly_scale(a, c):
        return [v * c for v in a]

    for step in range(num_edges):
        A, B, C = state
        new_A = poly_add(A, C)
        new_B = poly_add(poly_scale(poly_shift(A), Fraction(2)), poly_shift(C))
        new_C = list(B)
        state = [new_A, new_B, new_C]

    total = poly_add(poly_add(state[0], state[1]), state[2])

    # g_k = total[k] / 2
    result = {}
    for k in range(1, min(max_deg+1, len(total))):
        if total[k] != 0:
            result[k] = total[k] / 2
    return result

# Opus g_k polynomials
def opus_gk(k, m):
    """Opus g_k (claimed degree 3 for k >= 3)."""
    if k == 1: return Fraction(m)
    if k == 2: return Fraction(m * m)
    if k == 3: return Fraction(2*m**3 + m, 3)
    if k == 4: return Fraction(10*m**3 - 33*m**2 + 50*m - 24, 3)
    if k == 5: return Fraction(388*m**3 - 2040*m**2 + 3431*m - 1776, 3)
    if k == 6: return Fraction(69660*m**3 - 380445*m**2 + 653748*m - 342960, 3)
    if k == 7: return Fraction(19826270*m**3 - 109486152*m**2 + 189674605*m - 100014720, 3)
    if k == 8: return Fraction(7309726742*m**3 - 40641958545*m**2 + 70757788486*m - 37425556680, 3)
    if k == 9: return Fraction(3262687720240*m**3 - 18232387983408*m**2 + 31858349908595*m - 16888649645424, 3)
    return None

print("="*70)
print("TRUTH TEST: Transfer Matrix vs Opus g_k")
print("="*70)

for n in range(3, 11):  # up to n=10 (10! = 3.6M, borderline)
    if n > 8:
        # Skip direct computation for large n, use product formula
        print(f"\nn={n}: (skipping direct W computation, too large)")
        continue

    W = compute_W_direct(n)
    nfact = factorial(n)
    cv2_true = W / nfact - 1

    print(f"\nn={n}: W={W}, W/n!={W/nfact}, CV^2={float(cv2_true):.10f}")

    # Transfer matrix sum
    tm_gk = transfer_matrix_gk_all(n)
    cv2_tm = Fraction(0)
    for k, gk in tm_gk.items():
        m = n - 2*k
        if m < 0:
            continue
        term = 2 * gk / falling_factorial(n, 2*k)
        cv2_tm += term

    # Opus sum
    cv2_opus = Fraction(0)
    for k in range(1, n//2 + 1):
        m = n - 2*k
        if m < 0:
            continue
        gk = opus_gk(k, m)
        if gk is None:
            cv2_opus = None
            break
        term = 2 * gk / falling_factorial(n, 2*k)
        cv2_opus += term

    tm_match = "OK" if cv2_tm == cv2_true else f"FAIL (diff={cv2_tm - cv2_true})"
    opus_match = "OK" if cv2_opus is not None and cv2_opus == cv2_true else \
                 f"FAIL (diff={cv2_opus - cv2_true})" if cv2_opus is not None else "incomplete"

    print(f"  Transfer matrix: CV^2 = {float(cv2_tm):.10f} {tm_match}")
    opus_val = f"{float(cv2_opus):.10f}" if cv2_opus is not None else "?"
    print(f"  Opus cubic:      CV^2 = {opus_val} {opus_match}")

    # Show individual terms
    print(f"  Transfer terms:")
    for k in sorted(tm_gk.keys()):
        m = n - 2*k
        if m < 0: continue
        print(f"    k={k}: g_{k}({m}) = {tm_gk[k]}")

    if cv2_opus is not None:
        print(f"  Opus terms:")
        for k in range(1, n//2 + 1):
            m = n - 2*k
            if m < 0: continue
            gk = opus_gk(k, m)
            if gk is not None:
                print(f"    k={k}: g_{k}({m}) = {gk}")

# Now use the nud_weight.c program or W recurrence for larger n
# W(n) satisfies W(n) = n*W(n-1) + 2*(n-1)*W(n-2)
# With W(1) = 1, W(2) = 2

print("\n" + "="*70)
print("LARGER n VIA RECURRENCE: W(n) = n*W(n-1) + 2*(n-1)*W(n-2)")
print("="*70)

W_cache = {1: Fraction(1), 2: Fraction(2)}
def W_rec(n):
    if n in W_cache:
        return W_cache[n]
    result = n * W_rec(n-1) + 2*(n-1) * W_rec(n-2)
    W_cache[n] = result
    return result

for n in range(3, 22):
    W = W_rec(n)
    nfact = factorial(n)
    cv2_true = Fraction(W, nfact) - 1

    # Transfer matrix
    tm_gk = transfer_matrix_gk_all(n)
    cv2_tm = Fraction(0)
    for k, gk in tm_gk.items():
        m = n - 2*k
        if m < 0: continue
        cv2_tm += 2 * gk / falling_factorial(n, 2*k)

    # Opus
    cv2_opus = Fraction(0)
    complete = True
    for k in range(1, n//2 + 1):
        m = n - 2*k
        if m < 0: continue
        gk = opus_gk(k, m)
        if gk is None:
            complete = False
            break
        cv2_opus += 2 * gk / falling_factorial(n, 2*k)

    tm_ok = cv2_tm == cv2_true
    opus_ok = complete and cv2_opus == cv2_true

    tm_str = "OK" if tm_ok else f"FAIL (err={float(cv2_tm-cv2_true):.2e})"
    opus_str = ("OK" if opus_ok else f"FAIL (err={float(cv2_opus-cv2_true):.2e})") if complete else "incomplete (k too high)"

    print(f"n={n:2d}: CV^2={float(cv2_true):.8f}  TM:{tm_str:20s}  Opus:{opus_str}")

print("\nDone!")
