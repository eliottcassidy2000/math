#!/usr/bin/env python3
"""
Q(x) DIVISIBILITY — THE x(x-2) FACTOR AND H ≥ |Pf| PROOF PATH
opus-2026-03-14-S68

From x_deformation_pfaffian.py we found:
  Q(x) = I(CG,x)² - det(I+xA) is always divisible by x.

For the single-cycle case: Q(x) = -x(x-2)(x+1) — x=2 is a root!

KEY QUESTION: Is Q(x) always divisible by x(x-2)?
If yes, then at x=2: Q(2) = 0 implies H² = det(I+2A) = Pf(S)² for ALL
tournaments, not just single-cycle ones. That would mean H = |Pf| always.

But we KNOW H > |Pf| for many tournaments. So x(x-2) does NOT always divide Q(x).

NEW QUESTION: What IS the x=2 substitution structure?
  Q(2) = I(CG,2)² - det(I+2A) = H² - Pf² = 8q with q ≥ 0.
  Q(x) = x · R(x) for some R(x) (since Q(0) = 0 always).
  So Q(2) = 2 · R(2) = 8q, meaning R(2) = 4q.
  R(2) ≡ 0 (mod 4)???

Let's investigate R(x) = Q(x)/x systematically.
"""

import numpy as np
from itertools import combinations, permutations
from sympy import symbols, Matrix, det, factor, expand, Poly, Symbol, div, Rational

x = Symbol('x')

def adj_matrix_bits(bits, n):
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> idx) & 1: A[i][j] = 1
            else: A[j][i] = 1
            idx += 1
    return A

def count_hp(A_np, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)): continue
            if (mask, v) not in dp: continue
            for u in range(n):
                if mask & (1 << u): continue
                if A_np[v][u]:
                    key = (mask | (1 << u), u)
                    dp[key] = dp.get(key, 0) + dp[(mask, v)]
    return sum(dp.get(((1 << n) - 1, v), 0) for v in range(n))

def find_odd_cycles(A, n):
    cycles = []
    for length in range(3, n+1, 2):
        for verts in combinations(range(n), length):
            for perm in permutations(verts):
                if all(A[perm[i]][perm[(i+1) % length]] for i in range(length)):
                    mi = perm.index(min(perm))
                    canon = perm[mi:] + perm[:mi]
                    cycles.append(canon)
    return list(set(cycles))

def build_cg(cycles):
    nc = len(cycles)
    adj = [[False]*nc for _ in range(nc)]
    for i in range(nc):
        for j in range(i+1, nc):
            if set(cycles[i]) & set(cycles[j]):
                adj[i][j] = adj[j][i] = True
    return adj

def indep_poly(adj, nc):
    alpha = [0]*(nc+1)
    for mask in range(1 << nc):
        verts = [i for i in range(nc) if mask & (1 << i)]
        ok = True
        for a in range(len(verts)):
            for b in range(a+1, len(verts)):
                if adj[verts[a]][verts[b]]:
                    ok = False; break
            if not ok: break
        if ok:
            alpha[len(verts)] += 1
    return alpha

print("=" * 78)
print("  Q(x) DIVISIBILITY ANALYSIS — THE PATH TO H ≥ |Pf|")
print("=" * 78)

# ============================================================================
# PART 1: Q(x)/x = R(x), what is R(2)?
# ============================================================================

print("""
╔══════════════════════════════════════════════════════════════════════════════╗
║  PART 1: R(x) = Q(x)/x — THE REDUCED GAP POLYNOMIAL                       ║
╚══════════════════════════════════════════════════════════════════════════════╝

Q(x) = I(CG,x)² - det(I+xA) = x · R(x)
At x=2: Q(2) = 2·R(2) = H² - Pf²
So R(2) = (H² - Pf²)/2 = (H-|Pf|)(H+|Pf|)/2

Since H and |Pf| are both odd: H-|Pf| and H+|Pf| are both even.
So (H-|Pf|)(H+|Pf|)/2 is even·even/2 = even.
Actually R(2) = (H²-Pf²)/2. We know H²-Pf² ≡ 0 (mod 8), so R(2) ≡ 0 (mod 4).
""")

for n in [3, 4, 5]:
    m = n*(n-1)//2
    print(f"\nn={n}:")
    seen = {}

    for bits in range(1 << m):
        A_list = adj_matrix_bits(bits, n)
        A_np = np.array(A_list)
        H = count_hp(A_np, n)

        cycles = find_odd_cycles(A_np, n)
        if len(cycles) > 12: continue

        cg = build_cg(cycles)
        alpha = indep_poly(cg, len(cycles))

        Hx = sum(alpha[k] * x**k for k in range(len(alpha)))
        A_sym = Matrix(A_list)
        det_x = expand(det(Matrix.eye(n) + x * A_sym))
        Q = expand(Hx**2 - det_x)

        key = str(expand(Q))
        if key in seen: continue
        seen[key] = True

        if Q == 0:
            print(f"  H={H}: Q(x)=0 identically → R(x)=0")
            continue

        # Divide Q by x
        R = expand(Q / x)  # symbolic division
        # Verify it's a polynomial
        R_poly = Poly(R, x)

        R_at_2 = int(R.subs(x, 2))
        Q_at_2 = int(Q.subs(x, 2))

        # det at x=2
        det_at_2 = int(det_x.subs(x, 2))
        pf_sq = det_at_2
        pf = int(round(abs(pf_sq)**0.5))

        print(f"  H={H}: Q(x) = {factor(Q)}")
        print(f"    R(x) = Q(x)/x = {factor(R)}")
        print(f"    R(2) = {R_at_2}, Q(2) = {Q_at_2} = H²-Pf² = {H}²-{pf}² = {H**2-pf**2}")
        print(f"    R(2) mod 4 = {R_at_2 % 4}")

        # Check: is R(x) divisible by (x-2)?
        R_at_2_check = R.subs(x, 2)
        if R_at_2_check == 0:
            # R(x) = (x-2) · S(x)
            S = expand(R / (x - 2))
            print(f"    R(2)=0! So Q(x) = x(x-2)·S(x) with S(x) = {factor(S)}")
            print(f"    S(2) = {S.subs(x, 2)}")
        else:
            # Factor out what we can from R(2)
            # R(2) = (H²-Pf²)/2, should be ≡ 0 (mod 4)
            print(f"    R(2) ≠ 0 → Q(x) NOT divisible by x(x-2)")

# ============================================================================
# PART 2: THE x=2 SUBSTITUTION — WHY IS Q(2) ≥ 0?
# ============================================================================

print("""
╔══════════════════════════════════════════════════════════════════════════════╗
║  PART 2: WHY Q(2) ≥ 0? — THE ALGEBRAIC MECHANISM                           ║
╚══════════════════════════════════════════════════════════════════════════════╝

We have Q(x) = I(CG,x)² - det(I+xA).

Method: write det(I+xA) in terms of I(CG,x) directly.

At x=2: I(CG,2) = H, and det(I+2A) = Pf(S)². We need H² ≥ Pf².

Approach: Use the matrix-tree theorem generalization.
det(I+xA) involves ALL substructures of A, while I(CG,x) only involves
independent sets of odd cycles. The "extra" terms in det should not
overwhelm the independence polynomial.

Let's look at det(I+xA) as a polynomial in x and compare its coefficients
with those of I(CG,x)².
""")

print("Coefficient-by-coefficient comparison: I(CG,x)² vs det(I+xA)")
print()

for n in [3, 4, 5]:
    m = n*(n-1)//2
    print(f"n={n}:")
    seen_pairs = {}

    for bits in range(1 << m):
        A_list = adj_matrix_bits(bits, n)
        A_np = np.array(A_list)
        H = count_hp(A_np, n)

        cycles = find_odd_cycles(A_np, n)
        if len(cycles) > 12: continue

        cg = build_cg(cycles)
        alpha = indep_poly(cg, len(cycles))

        Hx = sum(alpha[k] * x**k for k in range(len(alpha)))
        A_sym = Matrix(A_list)
        det_x = expand(det(Matrix.eye(n) + x * A_sym))

        Hx2 = expand(Hx**2)

        # Get coefficients
        p_H2 = Poly(Hx2, x) if Hx2.has(x) else Poly(Hx2, x)
        p_det = Poly(det_x, x) if det_x.has(x) else Poly(det_x, x)

        h2_coeffs = dict(p_H2.as_dict()) if Hx2 != 1 else {(0,): 1}
        det_coeffs = dict(p_det.as_dict()) if det_x != 1 else {(0,): 1}

        key = (tuple(sorted(h2_coeffs.items())), tuple(sorted(det_coeffs.items())))
        if key in seen_pairs: continue
        seen_pairs[key] = True

        # Build aligned coefficient list
        max_deg = max(
            max((k[0] for k in h2_coeffs), default=0),
            max((k[0] for k in det_coeffs), default=0)
        )

        print(f"  H={H}, α={alpha[:max(k for k,a in enumerate(alpha) if a)+1] if any(alpha) else [1]}:")
        print(f"    {'deg':>4}  {'[x^d] I²':>10}  {'[x^d] det':>10}  {'diff':>8}")
        for d in range(max_deg + 1):
            c_h2 = h2_coeffs.get((d,), 0)
            c_det = det_coeffs.get((d,), 0)
            print(f"    {d:4d}  {int(c_h2):10d}  {int(c_det):10d}  {int(c_h2 - c_det):8d}")
        print()
    print()

# ============================================================================
# PART 3: THE SIGN PATTERN — WHEN IS Q(x) > 0 vs < 0?
# ============================================================================

print("""
╔══════════════════════════════════════════════════════════════════════════════╗
║  PART 3: SIGN PATTERN OF Q(x) — IS x=2 ALWAYS IN THE POSITIVE REGION?      ║
╚══════════════════════════════════════════════════════════════════════════════╝

For Q(x) not identically zero:
- Q(0) = 0 always
- Q(2) ≥ 0 always (this is H ≥ |Pf|)
- Q(x) can be negative for some x > 0

WHERE are the sign changes? Is x=2 always between roots (or at a root)?
""")

from sympy import solve, real_roots, nroots

for n in [3, 4, 5]:
    m = n*(n-1)//2
    print(f"\nn={n}:")
    seen = {}

    for bits in range(1 << m):
        A_list = adj_matrix_bits(bits, n)
        A_np = np.array(A_list)
        H = count_hp(A_np, n)

        cycles = find_odd_cycles(A_np, n)
        if len(cycles) > 12: continue

        cg = build_cg(cycles)
        alpha = indep_poly(cg, len(cycles))
        Hx = sum(alpha[k] * x**k for k in range(len(alpha)))
        A_sym = Matrix(A_list)
        det_x = expand(det(Matrix.eye(n) + x * A_sym))
        Q = expand(Hx**2 - det_x)

        key = str(Q)
        if key in seen or Q == 0: continue
        seen[key] = True

        # Find real roots of Q(x)
        try:
            roots = nroots(Poly(Q, x), n=6)
            real_r = sorted([complex(r).real for r in roots if abs(complex(r).imag) < 1e-6])
        except:
            real_r = []

        roots_str = ", ".join(f"{r:.4f}" for r in real_r)
        print(f"  H={H}: Q(x) = {factor(Q)}")
        print(f"    Real roots: {roots_str}")
        print(f"    Q(2) = {int(Q.subs(x, 2))}")

        # Check sign in [0, 2] interval
        signs = []
        for xv in [0.5, 1.0, 1.5, 2.0, 2.5, 3.0]:
            val = float(Q.subs(x, xv))
            signs.append(f"Q({xv})={'+'if val>0 else '-' if val<0 else '0'}{abs(val):.1f}")
        print(f"    Signs: {', '.join(signs)}")
        print()

# ============================================================================
# PART 4: THE CRITICAL INSIGHT — det(I+xA) AS A PERMANENT?
# ============================================================================

print("""
╔══════════════════════════════════════════════════════════════════════════════╗
║  PART 4: SEARCHING FOR THE BRIDGE — det vs perm vs I²                      ║
╚══════════════════════════════════════════════════════════════════════════════╝

At x=2:
  det(I+2A) = Pf(S)²   (THM-174)
  I(CG,2)² = H²        (OCF squared)
  H² ≥ Pf²             (observed, need proof)

Could there be a matrix M(x) such that:
  perm(M) = I(CG,x)²  and  det(M) = det(I+xA)?

Then |det| ≤ perm (Hadamard-type inequality) would give H² ≥ |Pf|² = det.

Let's check: for the "2-or-0" matrix M = I+2A at x=2:
  M_{ij} ∈ {0, 1, 2}
  det(M) = Pf(S)²
  perm(M) = ???
""")

def permanent(M):
    """Compute permanent of matrix M (Ryser formula for small matrices)."""
    n = len(M)
    if n == 0: return 1

    # Ryser's formula
    total = 0
    for S in range(1, 1 << n):
        # S encodes subset of columns
        col_sums = [0] * n
        bits_S = bin(S).count('1')
        for j in range(n):
            if S & (1 << j):
                for i in range(n):
                    col_sums[i] += M[i][j]
        prod = 1
        for i in range(n):
            prod *= col_sums[i]
        total += ((-1)**(n - bits_S)) * prod

    return total * ((-1)**n)

print("Comparing det(I+2A) vs perm(I+2A) vs H² for small tournaments:")
print(f"{'n':>3}  {'H':>5}  {'det(I+2A)':>10}  {'perm(I+2A)':>12}  {'H²':>8}  {'perm=H²?':>10}")
print("-" * 55)

for n in [3, 4, 5]:
    m = n*(n-1)//2
    seen_H = set()
    for bits in range(1 << m):
        b = [(bits >> i) & 1 for i in range(m)]
        A_np = np.zeros((n,n), dtype=int)
        idx = 0
        for i in range(n):
            for j in range(i+1, n):
                if b[idx]: A_np[i][j] = 1
                else: A_np[j][i] = 1
                idx += 1

        H = count_hp(A_np, n)
        if H in seen_H: continue
        seen_H.add(H)

        M = np.eye(n, dtype=int) + 2 * A_np
        d = int(round(np.linalg.det(M.astype(float))))
        p = permanent(M.tolist())

        print(f"{n:3d}  {H:5d}  {d:10d}  {p:12d}  {H**2:8d}  {'✓' if p == H**2 else '✗':>10}")

# ============================================================================
# PART 5: THE PERMANENT CONNECTION — DEEPER
# ============================================================================

print("""
╔══════════════════════════════════════════════════════════════════════════════╗
║  PART 5: perm(I+2A) vs H² — IS THE PERMANENT EQUAL TO H²?                  ║
╚══════════════════════════════════════════════════════════════════════════════╝
""")

# Exhaustive check for n=3,4,5
for n in [3, 4, 5, 6]:
    m = n*(n-1)//2
    total = 0
    matches = 0
    max_test = min(1 << m, 2**15)

    for bits in range(max_test):
        b = [(bits >> i) & 1 for i in range(m)]
        A_np = np.zeros((n,n), dtype=int)
        idx = 0
        for i in range(n):
            for j in range(i+1, n):
                if b[idx]: A_np[i][j] = 1
                else: A_np[j][i] = 1
                idx += 1

        H = count_hp(A_np, n)
        M = np.eye(n, dtype=int) + 2 * A_np
        p = permanent(M.tolist())

        total += 1
        if p == H**2:
            matches += 1

    print(f"  n={n}: {matches}/{total} have perm(I+2A) = H² ({'ALL' if matches == total else f'{matches}/{total}'})")

print("""
╔══════════════════════════════════════════════════════════════════════════════╗
║  PART 6: IF perm = H², THEN THE PROOF IS:                                   ║
╚══════════════════════════════════════════════════════════════════════════════╝

If perm(I+2A) = H² for all tournaments, then:

  H² = perm(I+2A)         (conjectured identity)
  Pf² = det(I+2A)         (THM-174)
  perm(M) ≥ |det(M)|      (for nonneg real M: trivially perm ≥ |det| by
                            definition since perm uses + while det uses ±)

Wait: perm ≥ |det| is NOT true in general!
For M with nonneg entries: det can be negative while perm is always ≥ 0.
But |det| ≤ perm is also not always true.

HOWEVER: for our M = I+2A with entries in {0, 1, 2}:
  det(M) = Pf² ≥ 0 (always nonneg! a perfect square!)
  perm(M) = H² ≥ 0 (always nonneg!)

The question is: does perm(M) ≥ det(M) for our specific class of matrices?
M is a {0,1,2}-matrix with 1s on diagonal, 0 or 2 off-diagonal, where
exactly one of M_{ij}, M_{ji} is 2 and the other is 0 for i ≠ j.

This is a TOURNAMENT MATRIX — a very special structure.
For tournament matrices: perm ≥ det ⟺ H² ≥ Pf² ⟺ H ≥ |Pf|.
""")

# Check perm ≥ det for all tournament matrices
print("Checking perm(I+2A) ≥ det(I+2A) for tournament matrices:")
for n in [3, 4, 5, 6]:
    m = n*(n-1)//2
    violations = 0
    total = 0
    max_test = min(1 << m, 2**15)

    for bits in range(max_test):
        b = [(bits >> i) & 1 for i in range(m)]
        A_np = np.zeros((n,n), dtype=int)
        idx = 0
        for i in range(n):
            for j in range(i+1, n):
                if b[idx]: A_np[i][j] = 1
                else: A_np[j][i] = 1
                idx += 1

        M = np.eye(n, dtype=int) + 2 * A_np
        d = int(round(np.linalg.det(M.astype(float))))
        p = permanent(M.tolist())

        total += 1
        if p < d:
            violations += 1

    print(f"  n={n}: {violations} violations in {total} tournaments {'✓ ALL PASS' if violations == 0 else '✗ VIOLATIONS'}")

print("""
════════════════════════════════════════════════════════════════════════════════
SYNTHESIS
════════════════════════════════════════════════════════════════════════════════

TWO CONJECTURES:

CONJECTURE A: perm(I+2A) = H(T)² for all tournaments T.
  If true, this gives H as a permanent, connecting to Valiant's theory.

CONJECTURE B: perm(M) ≥ det(M) for all tournament matrices M = I+2A.
  Combined with Conj A and THM-174: H² = perm ≥ det = Pf² → H ≥ |Pf|.

Note: perm ≥ det for nonneg matrices is related to the van der Waerden
conjecture (proved by Egorychev and Falikman). For doubly stochastic
matrices, perm ≥ n!/n^n. For our {0,1,2} tournament matrices, there may
be a direct proof using the tournament structure.

THE 2-3 CONNECTION:
  Matrix entries are 0, 1, 2 — the TRINITY numbers.
  perm counts unsigned weighted permutations with weights 0, 1, 2.
  det counts signed weighted permutations with weights 0, 1, 2.
  The difference perm - det = Q involves CROSS TERMS between + and - signs.
""")
