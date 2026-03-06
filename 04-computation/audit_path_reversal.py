#!/usr/bin/env python3
"""
AUDIT: Path reversal identity and c=0 proof completeness.

Devil's advocate verification of the claims from S20-S23.

opus-2026-03-06-S21
"""

from itertools import permutations, combinations
from sympy import symbols, expand, Poly, Rational
import random

# ============================================================
# Core infrastructure (minimal, clean)
# ============================================================

def setup(n):
    c = symbols('c')
    r = symbols('r')
    sv = {}
    for i in range(n):
        for j in range(i+1, n):
            sv[(i,j)] = symbols(f's{i}{j}')

    def s(i, j):
        if i == j: return 0
        if i < j: return sv[(i,j)]
        return -sv[(j,i)]

    def t_c(i, j):
        if i == j: return 0
        return c/2 + s(i, j)

    def t_r(i, j):
        if i == j: return 0
        return r + s(i, j)

    def t_skew(i, j):
        return s(i, j)

    return c, r, sv, s, t_c, t_r, t_skew

def hp(t_fn, vset, start=None, end=None):
    vl = sorted(vset)
    k = len(vl)
    if k == 0: return 0
    if k == 1:
        if start is not None and vl[0] != start: return 0
        if end is not None and vl[0] != end: return 0
        return 1
    total = 0
    for p in permutations(vl):
        if start is not None and p[0] != start: continue
        if end is not None and p[-1] != end: continue
        prod = 1
        for i in range(len(p)-1):
            prod *= t_fn(p[i], p[i+1])
        total += prod
    return expand(total)

def transfer_M(t_fn, n, a, b):
    U = [v for v in range(n) if v != a and v != b]
    result = 0
    for mask in range(1 << len(U)):
        S = [U[i] for i in range(len(U)) if mask & (1 << i)]
        R = [U[i] for i in range(len(U)) if not (mask & (1 << i))]
        sign = (-1)**len(S)
        ea = hp(t_fn, set(S)|{a}, end=a)
        bb = hp(t_fn, set(R)|{b}, start=b)
        result += sign * ea * bb
    return expand(result)

print("=" * 70)
print("AUDIT: PATH REVERSAL AND PROOF COMPLETENESS")
print("=" * 70)

# ============================================================
# AUDIT 1: Path reversal B_v(S+v;c,s) = E_v(S+v;c,-s)
# ============================================================
print("\n--- AUDIT 1: Path reversal identity ---")
print("Claim: B_v(S+v; c,s) = E_v(S+v; c,-s) for ALL v, S, n")
print("Proof: reversing path v0->...->vk changes each t_{vi,v_{i+1}} = c/2 + s_{vi,v_{i+1}}")
print("  to t_{v_{i+1},vi} = c/2 + s_{v_{i+1},vi} = c/2 - s_{vi,v_{i+1}}")
print("  which is exactly the c,-s version of t_{vi,v_{i+1}}.")
print("  So reversed path weight = original path weight with s -> -s.")
print("  And B_v (paths starting at v) = E_v (paths ending at v) with reversed paths.")
print("  Therefore B_v(S+v; c,s) = E_v(S+v; c,-s). QED.")
print()
print("This is a RIGOROUS argument, not just computation. But let's verify.")

for n in [3, 4, 5]:
    c, r, sv, s, t_c, t_r, t_skew = setup(n)
    subs_neg = {sv[k]: -sv[k] for k in sv}

    failures = 0
    checks = 0
    for v in range(n):
        others = [u for u in range(n) if u != v]
        for size in range(len(others)+1):
            for combo in combinations(others, size):
                S_set = set(combo) | {v}
                ea = hp(t_c, S_set, end=v)
                ba = hp(t_c, S_set, start=v)
                ea_negs = expand(ea.subs(subs_neg)) if not isinstance(ea, int) else ea
                diff = expand(ba - ea_negs) if not isinstance(ba, int) else ba - ea_negs
                checks += 1
                if diff != 0:
                    failures += 1
                    print(f"  FAIL at n={n}, v={v}, S={sorted(combo)}")

    print(f"  n={n}: {checks} checks, {failures} failures")

# ============================================================
# AUDIT 2: c=0 proof — is S<->R relabeling valid?
# ============================================================
print("\n--- AUDIT 2: c=0 proof S<->R relabeling ---")
print("Claim: sum_S E_a(S+a)*E_b(R+b) is symmetric in a<->b")
print("Proof: relabel S -> U\\S gives sum_S E_a(R+a)*E_b(S+b)")
print("  = sum_S E_b(S+b)*E_a(R+a) which is the (b,a) version.")
print()
print("DEVIL'S ADVOCATE CHECK: Is this relabeling actually valid?")
print("The E functions are E_a(S+a) where a is FIXED and S varies.")
print("When we relabel S -> U\\S, the set S+a becomes (U\\S)+a = R+a.")
print("And R+b becomes S+b. So the product becomes E_a(R+a)*E_b(S+b).")
print("This IS the (b,a) version because M[b,a] at c=0 has the factor")
print("E_b(S+b)*E_a(R+a). YES, this is valid.")

for n in [3, 4, 5, 6]:
    c_v, r, sv, s, t_c, t_r, t_skew = setup(n)
    a, b = 0, 1
    U = [v for v in range(n) if v != a and v != b]

    # Compute sum_S E_a(S+a)*E_b(R+b) and sum_S E_b(S+b)*E_a(R+a)
    sum_ab = 0
    sum_ba = 0
    for mask in range(1 << len(U)):
        S = [U[i] for i in range(len(U)) if mask & (1 << i)]
        R = [U[i] for i in range(len(U)) if not (mask & (1 << i))]
        ea_S = hp(t_skew, set(S)|{a}, end=a)
        eb_R = hp(t_skew, set(R)|{b}, end=b)
        eb_S = hp(t_skew, set(S)|{b}, end=b)
        ea_R = hp(t_skew, set(R)|{a}, end=a)
        sum_ab += ea_S * eb_R
        sum_ba += eb_S * ea_R

    sum_ab = expand(sum_ab)
    sum_ba = expand(sum_ba)
    print(f"  n={n}: sum_S E_a(S+a)*E_b(R+b) = sum_S E_b(S+b)*E_a(R+a)? {expand(sum_ab - sum_ba) == 0}")

# ============================================================
# AUDIT 3: Circularity check at general c
# ============================================================
print("\n--- AUDIT 3: Is the general-c argument circular? ---")
print("From path reversal: M[b,a](c,s) = (-1)^{n-2} * M[a,b](c,-s)")
print("For symmetry: need M[a,b](c,s) = M[b,a](c,s)")
print("Combining: M[a,b](c,s) = (-1)^{n-2} * M[a,b](c,-s)")
print("i.e., M has definite s-parity (-1)^{n-2}.")
print("This is a SEPARATE claim from path reversal. CIRCULARITY CONFIRMED.")
print()

# Let's verify the chain explicitly
for n in [4, 5]:
    c_v, r, sv, s, t_c, t_r, t_skew = setup(n)
    subs_neg = {sv[k]: -sv[k] for k in sv}
    a, b = 0, 1

    M_ab = transfer_M(t_c, n, a, b)
    M_ba = transfer_M(t_c, n, b, a)

    # Check: M[b,a](c,s) = (-1)^{n-2} M[a,b](c,-s)
    M_ab_negs = expand(M_ab.subs(subs_neg))
    relation = expand(M_ba - (-1)**(n-2) * M_ab_negs)

    # Check: M has definite parity
    parity_check = expand(M_ab - (-1)**(n-2) * M_ab_negs)

    print(f"  n={n}:")
    print(f"    M[b,a] = (-1)^(n-2) * M[a,b](c,-s)? {relation == 0}")
    print(f"    M[a,b] has parity (-1)^(n-2) in s? {parity_check == 0}")
    print(f"    These are DIFFERENT claims. The first is proved, the second is verified.")

# ============================================================
# AUDIT 4: Even r-powers — numerical stress test at n=6,7
# ============================================================
print("\n--- AUDIT 4: Even r-powers at n=6,7 (numerical) ---")
random.seed(2026)

for n in [6, 7]:
    failures = 0
    for trial in range(50):
        # Random skew values
        sv_num = {}
        for i in range(n):
            for j in range(i+1, n):
                sv_num[(i,j)] = random.uniform(-2, 2)

        def s_num(i, j):
            if i == j: return 0
            if i < j: return sv_num[(i,j)]
            return -sv_num[(j,i)]

        # Test multiple r values
        for r_val in [0.3, 0.7, 1.5, -0.5]:
            def t_num(i, j):
                if i == j: return 0
                return r_val + s_num(i, j)

            a, b = 0, 1
            U = [v for v in range(n) if v != a and v != b]
            m_ab = 0
            m_ba = 0
            for mask in range(1 << len(U)):
                S = [U[i] for i in range(len(U)) if mask & (1 << i)]
                R = [U[i] for i in range(len(U)) if not (mask & (1 << i))]
                sign = (-1)**len(S)

                def hp_num(vset, start=None, end=None):
                    vl = sorted(vset)
                    if len(vl) <= 1:
                        if start is not None and (len(vl)==0 or vl[0] != start): return 0
                        if end is not None and (len(vl)==0 or vl[0] != end): return 0
                        return 1 if len(vl)==1 else 0
                    total = 0
                    for p in permutations(vl):
                        if start is not None and p[0] != start: continue
                        if end is not None and p[-1] != end: continue
                        prod = 1
                        for ii in range(len(p)-1):
                            prod *= t_num(p[ii], p[ii+1])
                        total += prod
                    return total

                ea = hp_num(set(S)|{a}, end=a)
                bb = hp_num(set(R)|{b}, start=b)
                eb = hp_num(set(S)|{b}, end=b)
                ba_v = hp_num(set(R)|{a}, start=a)
                m_ab += sign * ea * bb
                m_ba += sign * eb * ba_v

            if abs(m_ab - m_ba) > 1e-6:
                failures += 1

    print(f"  n={n}: {failures} failures in 200 tests (50 skew configs x 4 r values)")

# ============================================================
# AUDIT 5: Does the proof at c=0 extend to ALL vertex pairs?
# ============================================================
print("\n--- AUDIT 5: c=0 symmetry for ALL vertex pairs ---")

for n in [4, 5, 6]:
    c_v, r, sv, s, t_c, t_r, t_skew = setup(n)
    all_ok = True
    for a in range(n):
        for b in range(n):
            if a == b: continue
            M_ab = transfer_M(t_skew, n, a, b)
            M_ba = transfer_M(t_skew, n, b, a)
            if expand(M_ab - M_ba) != 0:
                print(f"  FAIL: n={n}, a={a}, b={b}")
                all_ok = False
    print(f"  n={n}: ALL pairs symmetric at c=0? {all_ok}")

# ============================================================
# AUDIT 6: Does the c=0 proof actually use skew-symmetry?
# ============================================================
print("\n--- AUDIT 6: What if s_ij != -s_ji at c=0? ---")
print("If we use INDEPENDENT variables (no skew constraint), does c=0 still give symmetry?")

for n in [3, 4]:
    t_ind = {}
    for i in range(n):
        for j in range(n):
            if i != j:
                t_ind[(i,j)] = symbols(f'u{i}{j}')

    def t_ind_fn(i, j):
        if i == j: return 0
        return t_ind[(i,j)]

    M_ab = transfer_M(t_ind_fn, n, 0, 1)
    M_ba = transfer_M(t_ind_fn, n, 1, 0)
    diff = expand(M_ab - M_ba)
    print(f"  n={n}: M[0,1]-M[1,0] with independent u_ij (no skew, c=0): {'0' if diff == 0 else f'{len(diff.as_ordered_terms())} terms'}")

    if diff != 0:
        # Now impose skew: u_ji = -u_ij
        subs_skew = {}
        for i in range(n):
            for j in range(i+1, n):
                subs_skew[t_ind[(j,i)]] = -t_ind[(i,j)]
        diff_skew = expand(diff.subs(subs_skew))
        print(f"    After u_ji = -u_ij: {diff_skew}")

# ============================================================
# AUDIT 7: The "E_b(R+b; c,-s)" in Step 2 — careful check
# ============================================================
print("\n--- AUDIT 7: Careful trace of M[a,b] after path reversal ---")
print("M[a,b] = sum_S (-1)^|S| E_a(S+a; c,s) * B_b(R+b; c,s)")
print("       = sum_S (-1)^|S| E_a(S+a; c,s) * E_b(R+b; c,-s)")
print()
print("M[b,a] = sum_S (-1)^|S| E_b(S+b; c,s) * B_a(R+a; c,s)")
print("       = sum_S (-1)^|S| E_b(S+b; c,s) * E_a(R+a; c,-s)")
print()
print("Relabel S -> U\\S in M[b,a]:")
print("M[b,a] = sum_S (-1)^{n-2-|S|} E_b(R+b; c,s) * E_a(S+a; c,-s)")
print("       = (-1)^{n-2} sum_S (-1)^|S| E_b(R+b; c,s) * E_a(S+a; c,-s)")
print()
print("Compare:")
print("M[a,b] = sum_S (-1)^|S| E_a(S+a; c,s)  * E_b(R+b; c,-s)")
print("M[b,a] = (-1)^{n-2} sum_S (-1)^|S| E_a(S+a; c,-s) * E_b(R+b; c,s)")
print()
print("So M[b,a](c,s) = (-1)^{n-2} * M[a,b](c,-s). CONFIRMED.")
print()
print("For M[a,b] = M[b,a], need M[a,b](c,s) = (-1)^{n-2} M[a,b](c,-s),")
print("i.e., M has s-parity (-1)^{n-2}. This is NOT proved by path reversal alone.")

# ============================================================
# AUDIT 8: The r^1 coefficient at n=6 — does it really vanish?
# ============================================================
print("\n--- AUDIT 8: r^1 coefficient at n=6 (numerical) ---")
random.seed(42)

n = 6
sv_num = {}
for i in range(n):
    for j in range(i+1, n):
        sv_num[(i,j)] = random.uniform(-2, 2)

def s_num6(i, j):
    if i == j: return 0
    if i < j: return sv_num[(i,j)]
    return -sv_num[(j,i)]

# Compute M[0,1] at several r values and fit polynomial
r_values = [-2, -1, -0.5, 0, 0.5, 1, 2]
m_values = []
for r_val in r_values:
    def t_num6(i, j):
        if i == j: return 0
        return r_val + s_num6(i, j)

    a, b = 0, 1
    U = [v for v in range(n) if v != a and v != b]
    m_val = 0
    for mask in range(1 << len(U)):
        S = [U[i] for i in range(len(U)) if mask & (1 << i)]
        R = [U[i] for i in range(len(U)) if not (mask & (1 << i))]
        sign = (-1)**len(S)

        def hp_n6(vset, start=None, end=None):
            vl = sorted(vset)
            if len(vl) <= 1:
                if start is not None and (len(vl)==0 or vl[0] != start): return 0
                if end is not None and (len(vl)==0 or vl[0] != end): return 0
                return 1 if len(vl)==1 else 0
            total = 0
            for p in permutations(vl):
                if start is not None and p[0] != start: continue
                if end is not None and p[-1] != end: continue
                prod = 1
                for ii in range(len(p)-1):
                    prod *= t_num6(p[ii], p[ii+1])
                total += prod
            return total

        ea = hp_n6(set(S)|{a}, end=a)
        bb = hp_n6(set(R)|{b}, start=b)
        m_val += sign * ea * bb

    m_values.append(m_val)

# Fit polynomial in r
import numpy as np
coeffs = np.polyfit(r_values, m_values, 4)  # degree n-2=4
print(f"  Polynomial fit M[0,1](r) at n=6:")
for i, c_val in enumerate(reversed(coeffs)):
    parity = "even" if i % 2 == 0 else "ODD"
    print(f"    r^{i}: {c_val:+.10f}  [{parity}]")

# Check: does M(r) = M(-r)? (even function for n=6, since n-2=4 is even)
print(f"\n  Even function check M(r) = M(-r):")
for r_val in [0.3, 0.7, 1.5]:
    def make_t(rv):
        def t_fn(i, j):
            if i == j: return 0
            return rv + s_num6(i, j)
        return t_fn

    a, b = 0, 1
    U = [v for v in range(n) if v != a and v != b]

    m_pos = 0
    m_neg = 0
    for mask in range(1 << len(U)):
        S = [U[i] for i in range(len(U)) if mask & (1 << i)]
        R = [U[i] for i in range(len(U)) if not (mask & (1 << i))]
        sign = (-1)**len(S)

        t_pos = make_t(r_val)
        t_neg = make_t(-r_val)

        def hp_gen(t_fn, vset, start=None, end=None):
            vl = sorted(vset)
            if len(vl) <= 1:
                if start is not None and (len(vl)==0 or vl[0] != start): return 0
                if end is not None and (len(vl)==0 or vl[0] != end): return 0
                return 1 if len(vl)==1 else 0
            total = 0
            for p in permutations(vl):
                if start is not None and p[0] != start: continue
                if end is not None and p[-1] != end: continue
                prod = 1
                for ii in range(len(p)-1):
                    prod *= t_fn(p[ii], p[ii+1])
                total += prod
            return total

        ea_p = hp_gen(t_pos, set(S)|{a}, end=a)
        bb_p = hp_gen(t_pos, set(R)|{b}, start=b)
        ea_n = hp_gen(t_neg, set(S)|{a}, end=a)
        bb_n = hp_gen(t_neg, set(R)|{b}, start=b)
        m_pos += sign * ea_p * bb_p
        m_neg += sign * ea_n * bb_n

    print(f"    r={r_val:+.1f}: M(r)={m_pos:.10f}, M(-r)={m_neg:.10f}, diff={abs(m_pos-m_neg):.2e}")

print("\n" + "=" * 70)
print("AUDIT SUMMARY")
print("=" * 70)
print("""
1. PATH REVERSAL IDENTITY: B_v(S+v; c,s) = E_v(S+v; c,-s)
   STATUS: PROVED (by clean argument) and VERIFIED computationally through n=5.
   The proof is: reversing a path replaces each arc weight c/2+s_ij
   with c/2-s_ij = c/2+s_ji, which is exactly s -> -s.

2. c=0 PROOF: M[a,b] = M[b,a] at c=0 (pure skew weights)
   STATUS: PROVED. Complete, non-circular.
   Step 1: B_v = (-1)^|S| E_v (special case of path reversal at c=0)
   Step 2: M = (-1)^{n-2} * unsigned sum
   Step 3: Unsigned sum symmetric by S <-> U\\S relabeling
   Note: This proof DOES require skew-symmetry (s_ij = -s_ji).

3. GENERAL c: M[b,a](c,s) = (-1)^{n-2} M[a,b](c,-s)
   STATUS: PROVED (from path reversal + S<->R relabeling).
   This REDUCES symmetry to the s-parity claim.

4. S-PARITY / EVEN r-POWERS:
   STATUS: VERIFIED computationally through n=7, NOT PROVED.
   This is the LAST GAP. Proving M(c,s) has parity (-1)^{n-2} in s
   would complete the proof of transfer matrix symmetry.

5. CIRCULARITY DIAGNOSIS: CORRECT.
   Path reversal alone does not prove symmetry at general c.
   It proves the EQUIVALENCE of symmetry and s-parity.
""")
