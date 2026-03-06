#!/usr/bin/env python3
"""
COMPLETE INDUCTIVE PROOF that M[a,b] has only even r-powers.

THEOREM: For a c-tournament on vertex set W with t_{ij} = r + s_{ij},
s_{ij} = -s_{ji}, the transfer matrix
  M_W[a,b] = sum_S (-1)^|S| E_a(S+{a}) B_b(R+{b})
has only even powers of r.

PROOF STRUCTURE:
  Key Identity: odd_r(B_b(W)) = r * col_sum_W(b)
  where col_sum_W(b) = sum_{v != b} M_W[v,b].

  This is proved by strong induction on |W|, using:
  1. Column recurrence for M
  2. Definite r-parity of total Ham path weight T(W)
  3. The inductive hypothesis at size |W|-1

  The proof closes without circularity.

kind-pasteur-2026-03-06-S25
"""

from itertools import permutations
from sympy import symbols, expand, Poly, Symbol

# ============================================================
# Core functions
# ============================================================

def setup(n):
    r = symbols('r')
    sv = {}
    for i in range(n):
        for j in range(i+1, n):
            sv[(i,j)] = symbols(f's{i}{j}')
    def s(i, j):
        if i == j: return 0
        if i < j: return sv[(i,j)]
        return -sv[(j,i)]
    def t(i, j):
        if i == j: return 0
        return r + s(i, j)
    return r, sv, s, t

def transfer_M(t_fn, vertex_set, a, b):
    """M_W[a,b] = sum_S (-1)^|S| E_a(S+a) B_b(R+b)."""
    V = sorted(vertex_set)
    U = [v for v in V if v != a and v != b]
    result = 0
    for mask in range(1 << len(U)):
        S = [U[i] for i in range(len(U)) if mask & (1 << i)]
        R = [U[i] for i in range(len(U)) if not (mask & (1 << i))]
        sign = (-1)**len(S)
        S_set = set(S) | {a}
        R_set = set(R) | {b}
        ea = 0
        for p in permutations(sorted(S_set)):
            if p[-1] != a: continue
            prod = 1
            for i in range(len(p)-1): prod *= t_fn(p[i], p[i+1])
            ea += prod
        if len(S_set) == 1: ea = 1
        bb = 0
        for p in permutations(sorted(R_set)):
            if p[0] != b: continue
            prod = 1
            for i in range(len(p)-1): prod *= t_fn(p[i], p[i+1])
            bb += prod
        if len(R_set) == 1: bb = 1
        result += sign * ea * bb
    return expand(result)

def total_ham(vertex_set, t_fn):
    """T(W) = sum over ALL permutations of W of product of t-weights."""
    vs = sorted(vertex_set)
    if len(vs) <= 1: return 1
    total = 0
    for perm in permutations(vs):
        prod = 1
        for i in range(len(perm)-1):
            prod *= t_fn(perm[i], perm[i+1])
        total += prod
    return expand(total)

def B_from(vertex_set, start, t_fn):
    """B_v(W) = Ham paths from v through W."""
    vs = sorted(vertex_set)
    if len(vs) <= 1: return 1
    total = 0
    for perm in permutations(vs):
        if perm[0] != start: continue
        prod = 1
        for i in range(len(perm)-1):
            prod *= t_fn(perm[i], perm[i+1])
        total += prod
    return expand(total)

def E_to(vertex_set, end, t_fn):
    """E_v(W) = Ham paths ending at v through W."""
    vs = sorted(vertex_set)
    if len(vs) <= 1: return 1
    total = 0
    for perm in permutations(vs):
        if perm[-1] != end: continue
        prod = 1
        for i in range(len(perm)-1):
            prod *= t_fn(perm[i], perm[i+1])
        total += prod
    return expand(total)

def odd_part(expr, r):
    if expr == 0: return 0
    p = Poly(expand(expr), r)
    d = p.as_dict()
    return expand(sum(c * r**deg[0] for deg, c in d.items() if deg[0] % 2 == 1))

def even_part(expr, r):
    if expr == 0: return 0
    p = Poly(expand(expr), r)
    d = p.as_dict()
    return expand(sum(c * r**deg[0] for deg, c in d.items() if deg[0] % 2 == 0))


# ============================================================
# STEP 0: Verify prerequisites
# ============================================================
print("=" * 70)
print("STEP 0: Verify prerequisites")
print("=" * 70)

# (a) T(W) has definite r-parity (-1)^{|W|-1}
print("\n(a) T(W) has definite r-parity (-1)^{|W|-1}:")
for m in range(1, 7):
    r_sym, sv, s, t = setup(m)
    W = set(range(m))
    T = total_ham(W, t)
    if T == 0:
        print(f"  |W|={m}: T=0")
        continue
    p = Poly(T, r_sym)
    powers = [d[0] for d in p.as_dict().keys()]
    expected_parity = (m - 1) % 2
    actual_parities = set(pw % 2 for pw in powers)
    ok = actual_parities == {expected_parity}
    print(f"  |W|={m}: r-powers = {sorted(powers)}, parity {actual_parities}, expected {{{expected_parity}}}: {'OK' if ok else 'FAIL'}")

# (b) B_b(-r) = (-1)^{m-1} E_b(r) [proved algebraically, verify here]
print("\n(b) B_b(W; -r) = (-1)^{m-1} E_b(W; r):")
for m in range(2, 6):
    r_sym, sv, s, t = setup(m)
    W = set(range(m))
    b = 0
    Bb = B_from(W, b, t)
    Eb = E_to(W, b, t)
    from sympy import sympify
    Bb_neg = expand(sympify(Bb).subs(r_sym, -r_sym))
    check = expand(Bb_neg - (-1)**(m-1) * Eb)
    print(f"  |W|={m}: B_b(-r) = (-1)^{{m-1}} E_b(r)? {check == 0}")


# ============================================================
# STEP 1: Verify column recurrence
# ============================================================
print("\n" + "=" * 70)
print("STEP 1: Column recurrence")
print("M_W[a,b] = (-1)^{m-2} E_a(W') + sum_w t(b,w) M_{W'}[a,w]")
print("where W' = W\\{b}")
print("=" * 70)

for m in range(3, 7):
    r_sym, sv, s, t = setup(m)
    W = set(range(m))
    all_ok = True
    for a in sorted(W):
        for b in sorted(W):
            if a == b: continue
            W_prime = W - {b}
            # Direct
            M_direct = transfer_M(t, W, a, b)
            # Column recurrence
            Ea = E_to(W_prime, a, t)
            rec = expand((-1)**(m-2) * Ea)
            for w in sorted(W_prime - {a}):
                rec = expand(rec + t(b, w) * transfer_M(t, W_prime, a, w))
            if expand(M_direct - rec) != 0:
                all_ok = False
                print(f"  m={m}, a={a}, b={b}: FAIL")
    print(f"  m={m}: Column recurrence holds for ALL (a,b) pairs: {all_ok}")


# ============================================================
# STEP 2: Verify Key Identity base cases
# ============================================================
print("\n" + "=" * 70)
print("STEP 2: Key Identity base cases")
print("odd(B_b(W)) = r * col_sum_W(b)")
print("=" * 70)

# m=1
print("  m=1: B_b({b}) = 1, odd(1) = 0, col_sum = 0, r*0 = 0. OK")

# m=2
r_sym, sv, s, t = setup(2)
W = {0, 1}
for b in [0, 1]:
    v = 1 - b
    Bb = B_from(W, b, t)
    Mvb = transfer_M(t, W, v, b)
    odd_Bb = odd_part(Bb, r_sym)
    print(f"  m=2, b={b}: B_b = {Bb}, odd = {odd_Bb}, M[{v},{b}] = {Mvb}, r*M = {expand(r_sym*Mvb)}, match: {expand(odd_Bb - r_sym*Mvb) == 0}")


# ============================================================
# STEP 3: Verify inductive step components
# ============================================================
print("\n" + "=" * 70)
print("STEP 3: Inductive step verification")
print("=" * 70)

for m in range(3, 7):
    r_sym, sv, s, t = setup(m)
    W = set(range(m))
    print(f"\n--- m = {m} ---")

    for b in [0]:  # Check one b value per m (symmetry ensures all work)
        W_prime = W - {b}
        n_prime = m - 1  # |W'|

        # Compute col_sum_W(b) directly
        col_sum = 0
        for v in sorted(W - {b}):
            col_sum = expand(col_sum + transfer_M(t, W, v, b))

        # Compute T(W')
        T_Wp = total_ham(W_prime, t)
        even_T = even_part(T_Wp, r_sym)

        # Compute col_sum_{W'}(w) for each w
        cs_Wp = {}
        for w in sorted(W_prime):
            cs = 0
            for v in sorted(W_prime - {w}):
                cs = expand(cs + transfer_M(t, W_prime, v, w))
            cs_Wp[w] = cs

        # Compute Sigma_{W'} = sum_w cs_{W'}(w)
        Sigma = sum(cs_Wp.values())
        Sigma = expand(Sigma)

        # VERIFY: r * Sigma = odd(T)
        odd_T = odd_part(T_Wp, r_sym)
        check_sigma = expand(r_sym * Sigma - odd_T)
        print(f"  r*Sigma = odd(T)? {check_sigma == 0}")

        # VERIFY: (-1)^{m-2} T + odd(T) = even(T)
        lhs_combine = expand((-1)**(m-2) * T_Wp + odd_T)
        check_combine = expand(lhs_combine - even_T)
        print(f"  (-1)^{{m-2}}*T + odd(T) = even(T)? {check_combine == 0}")

        # VERIFY: col_sum = even(T) + sum s_{bw} cs_w
        s_weighted = 0
        for w in sorted(W_prime):
            s_weighted = expand(s_weighted + s(b, w) * cs_Wp[w])

        target = expand(even_T + s_weighted)
        check_colsum = expand(col_sum - target)
        print(f"  col_sum = even(T) + sum s*cs? {check_colsum == 0}")

        # VERIFY THE FULL KEY IDENTITY: odd(B_b) = r * col_sum
        Bb = B_from(W, b, t)
        odd_Bb = odd_part(Bb, r_sym)
        check_key = expand(odd_Bb - r_sym * col_sum)
        print(f"  KEY IDENTITY: odd(B_b) = r*col_sum? {check_key == 0}")

        # VERIFY the decomposition proof:
        # B_b = r*T + sum s_{bv} B_v
        sum_sBv = 0
        for v in sorted(W_prime):
            Bv = B_from(W_prime, v, t)
            sum_sBv = expand(sum_sBv + s(b, v) * Bv)

        check_decomp = expand(Bb - r_sym * T_Wp - sum_sBv)
        print(f"  B_b = r*T + sum s*B_v? {check_decomp == 0}")

        # VERIFY: odd(sum s*B_v) = r * sum s*cs  (using B_v = even(B_v) + r*cs_v)
        # sum s*B_v = sum s*even(B_v) + r*sum s*cs_v
        # odd of first term = 0 (even*const = even)
        # odd of second term = r*sum s*cs (since sum s*cs is already even)
        sum_s_even_Bv = 0
        sum_s_cs = 0
        for v in sorted(W_prime):
            Bv = B_from(W_prime, v, t)
            even_Bv = even_part(Bv, r_sym)
            sum_s_even_Bv = expand(sum_s_even_Bv + s(b, v) * even_Bv)
            sum_s_cs = expand(sum_s_cs + s(b, v) * cs_Wp[v])

        # Check even(B_v) has even r-powers (and hence s*even(B_v) does too)
        check_seven = odd_part(sum_s_even_Bv, r_sym)
        print(f"  odd(sum s*even(B_v)) = 0? {check_seven == 0}")

        # Check sum s*cs has even r-powers
        check_scs = odd_part(sum_s_cs, r_sym)
        print(f"  odd(sum s*cs) = 0? {check_scs == 0}")

        # Full verification: odd(sum s*B_v) = r * sum s*cs
        odd_sBv = odd_part(sum_sBv, r_sym)
        check_odd_sBv = expand(odd_sBv - r_sym * sum_s_cs)
        print(f"  odd(sum s*B_v) = r*sum(s*cs)? {check_odd_sBv == 0}")

        # FINAL: odd(B_b) = r*even(T) + r*sum(s*cs) = r*[even(T) + sum(s*cs)] = r*col_sum
        final = expand(r_sym * even_T + r_sym * sum_s_cs - r_sym * col_sum)
        print(f"  r*even(T) + r*sum(s*cs) = r*col_sum? {final == 0}")


# ============================================================
# STEP 4: Verify M has even r-powers (consequence)
# ============================================================
print("\n" + "=" * 70)
print("STEP 4: Consequence - M[a,b] has only even r-powers")
print("=" * 70)

for m in range(2, 7):
    r_sym, sv, s, t = setup(m)
    W = set(range(m))
    all_even = True
    for a in sorted(W):
        for b in sorted(W):
            if a == b: continue
            M_ab = transfer_M(t, W, a, b)
            if M_ab == 0: continue
            p = Poly(M_ab, r_sym)
            odd_powers = [d[0] for d, c in p.as_dict().items() if d[0] % 2 == 1 and c != 0]
            if odd_powers:
                all_even = False
                print(f"  m={m}, M[{a},{b}]: ODD r-powers {odd_powers}!")
    print(f"  m={m}: All M[a,b] have even r-powers: {all_even}")

# Also verify M[a,b] = M[b,a] (symmetry, consequence of even r-powers + M(-r)=M^T)
print("\n  Symmetry M[a,b] = M[b,a] (consequence):")
for m in range(2, 7):
    r_sym, sv, s, t = setup(m)
    W = set(range(m))
    all_sym = True
    for a in sorted(W):
        for b in sorted(W):
            if a >= b: continue
            M_ab = transfer_M(t, W, a, b)
            M_ba = transfer_M(t, W, b, a)
            if expand(M_ab - M_ba) != 0:
                all_sym = False
    print(f"  m={m}: M symmetric: {all_sym}")


# ============================================================
# WRITE-UP
# ============================================================
print("\n" + "=" * 70)
print("COMPLETE PROOF (WRITE-UP)")
print("=" * 70)

print("""
THEOREM: For a c-tournament with t_{ij} = r + s_{ij}, s_{ij} = -s_{ji},
the transfer matrix M_W[a,b] has only even powers of r.
Equivalently, M_W[a,b] = M_W[b,a] (symmetry).

PROOF by strong induction on m = |W|.

We prove the KEY IDENTITY:
  odd_r(B_b(W)) = r * col_sum_W(b)
where col_sum_W(b) = sum_{v != b} M_W[v,b].

This implies M has even r-powers via the row recurrence
(see Step 4 explanation below).

--- Base cases ---
m = 1: B_b({b}) = 1 (degree 0), odd = 0. col_sum = 0. 0 = r*0.  [OK]
m = 2: B_b({b,v}) = r + s_{bv}. odd = r. M[v,b] = 1. r*1 = r.   [OK]

--- Inductive step (m >= 3) ---
Assume Key Identity for all vertex sets of size < m.

INGREDIENTS:

(I) Decomposition of B_b(W):
    B_b(W) = sum_{v in W'} t(b,v) B_v(W')
           = r * T(W') + sum_v s_{bv} B_v(W')
    where W' = W\\{b}, T(W') = sum_v B_v(W').

(II) Column recurrence for M:
    M_W[a,b] = (-1)^{m-2} E_a(W') + sum_{w in W'\\{a}} t(b,w) M_{W'}[a,w]

    Summing over a in W':
    col_sum_W(b) = (-1)^{m-2} T(W') + r * Sigma + sum_w s_{bw} cs_{W'}(w)
    where Sigma = sum_{a!=w in W'} M_{W'}[a,w] (total off-diagonal sum).

(III) By induction at size |W'| = m-1:
    odd(B_w(W')) = r * cs_{W'}(w) for each w in W'.

    Summing over w: odd(T(W')) = r * Sigma.

    Since T(W') has definite r-parity (-1)^{m-2}:
      (-1)^{m-2} T(W') + odd(T(W')) = even_r(T(W'))

    Therefore:
      col_sum_W(b) = even_r(T(W')) + sum_w s_{bw} cs_{W'}(w)   ...(*)

(IV) From the inductive hypothesis, M_{W'} has even r-powers, so:
    - Each cs_{W'}(w) has even r-powers
    - Each even_r(B_v(W')) has even r-powers
    - B_v(W') = even_r(B_v) + r * cs_{W'}(v)

PROOF OF KEY IDENTITY:

odd(B_b(W)) = odd(r * T + sum_v s_{bv} B_v)
            = r * even_r(T) + odd(sum_v s_{bv} B_v)

For the second term:
  sum_v s_{bv} B_v = sum_v s_{bv} even_r(B_v) + r * sum_v s_{bv} cs_v

  - s_{bv} * even_r(B_v) has only even r-powers  =>  odd part = 0
  - s_{bv} * cs_v has only even r-powers          =>  sum is even in r
  - r * (even in r) contributes entirely to odd part

  Therefore: odd(sum s B_v) = r * sum_v s_{bv} cs_v

Combining:
  odd(B_b) = r * even_r(T) + r * sum_v s_{bv} cs_v
           = r * [even_r(T) + sum_v s_{bv} cs_v]
           = r * col_sum_W(b)     [by (*)]

QED.

--- Why this implies M has even r-powers ---

Row recurrence: M_W[a,b] = B_b(W\\{a}) - sum_u t(u,a) M_{W\\{a}}[u,b]

By induction, M_{W\\{a}} has even r-powers.
So sum_u s_{ua} M^{(a)}[u,b] has even r-powers.
And r * col_sum^{(a)} has odd r-powers.

odd(M_W[a,b]) = odd(B_b(W\\{a})) - r * col_sum_{W\\{a}}(b) = 0
by the Key Identity at size m-1 (applied to B_b on W\\{a}).

Therefore M_W[a,b] has even r-powers. QED.
""")

print("=" * 70)
print("ALL VERIFICATIONS COMPLETE")
print("=" * 70)
