#!/usr/bin/env python3
"""
Exploration: the connection between positivity and skew-symmetry
in the transfer matrix.

Key question: WHY does the parity M(c,-s) = (-1)^{n-2} M(c,s) hold?

The transfer matrix M[a,b] = sum_S (-1)^|S| E_a(S+a) B_b(R+b)
involves SIGNED inclusion-exclusion. The signs (-1)^|S| break
positivity but create a parity filter on the s-variables.

This script decomposes M into:
1. The "positive" transfer matrix P (without signs)
2. The "negative" part N = P - M
3. Studies which s-degree terms survive in each

kind-pasteur-2026-03-06-S23
"""

from itertools import permutations
from sympy import symbols, expand, Poly, zeros, Rational
from collections import defaultdict

# ============================================================
# Setup: symbolic path computation in skew coordinates
# ============================================================

def setup_skew(n):
    """Create skew variables and helper functions for n vertices."""
    c = symbols('c')
    svars = {}
    for i in range(n):
        for j in range(i+1, n):
            svars[(i,j)] = symbols(f's{i}{j}')

    def s_val(i, j):
        if i == j: return 0
        if i < j: return svars[(i,j)]
        return -svars[(j,i)]

    def t_val(i, j):
        if i == j: return 0
        return c/2 + s_val(i, j)

    return c, svars, s_val, t_val

def hp_paths(t_val, vertex_set, start=None, end=None):
    """Hamiltonian path weight sum using arc weights t_val."""
    vl = sorted(vertex_set)
    k = len(vl)
    if k == 0: return 0
    if k == 1:
        if start is not None and vl[0] != start: return 0
        if end is not None and vl[0] != end: return 0
        return 1
    total = 0
    for perm in permutations(vl):
        if start is not None and perm[0] != start: continue
        if end is not None and perm[-1] != end: continue
        prod = 1
        for i in range(len(perm)-1):
            prod *= t_val(perm[i], perm[i+1])
        total += prod
    return expand(total)

def compute_transfer(t_val, n, a, b, signed=True):
    """Compute M[a,b] (signed=True) or P[a,b] (signed=False)."""
    U = [v for v in range(n) if v != a and v != b]
    result = 0
    for mask in range(1 << len(U)):
        S = [U[i] for i in range(len(U)) if mask & (1 << i)]
        R = [U[i] for i in range(len(U)) if not (mask & (1 << i))]
        sign = (-1)**len(S) if signed else 1
        ea = hp_paths(t_val, set(S)|{a}, end=a)
        bb = hp_paths(t_val, set(R)|{b}, start=b)
        result += sign * ea * bb
    return expand(result)


print("=" * 70)
print("POSITIVITY <-> SKEW-SYMMETRY CONNECTION")
print("=" * 70)

# ============================================================
# Part 1: The positive transfer matrix P[a,b] — is it symmetric?
# ============================================================
print("\n" + "=" * 70)
print("PART 1: Positive vs Signed Transfer Matrix")
print("=" * 70)

for n in [3, 4, 5]:
    c, svars, s_val, t_val = setup_skew(n)
    a, b = 0, 1

    M_ab = compute_transfer(t_val, n, a, b, signed=True)
    M_ba = compute_transfer(t_val, n, b, a, signed=True)
    P_ab = compute_transfer(t_val, n, a, b, signed=False)
    P_ba = compute_transfer(t_val, n, b, a, signed=False)

    D_signed = expand(M_ab - M_ba)
    D_unsigned = expand(P_ab - P_ba)

    print(f"\n  n={n}:")
    print(f"    M[a,b] - M[b,a] = {D_signed}")
    print(f"    P[a,b] - P[b,a] = {D_unsigned}")

    if D_unsigned != 0:
        # The positive version is NOT symmetric — good, this shows the
        # signs are essential!
        terms = D_unsigned.as_ordered_terms()
        print(f"    P is NOT symmetric ({len(terms)} terms in diff)")

    # At c=0:
    M_ab_0 = expand(M_ab.subs(c, 0))
    P_ab_0 = expand(P_ab.subs(c, 0))
    P_ba_0 = expand(P_ba.subs(c, 0))
    print(f"    P[a,b](c=0) = {P_ab_0}")
    print(f"    P[b,a](c=0) = {P_ba_0}")

# ============================================================
# Part 2: Parity decomposition of the POSITIVE matrix
# ============================================================
print("\n" + "=" * 70)
print("PART 2: How does the parity filter work?")
print("=" * 70)
print("""
Key insight: M = sum_S (-1)^|S| f(S) g(U\\S)
           P = sum_S f(S) g(U\\S)

So M uses inclusion-exclusion signs and P doesn't.
The difference is: M = P_even - P_odd where
  P_even = sum_{|S| even} f(S) g(U\\S)
  P_odd  = sum_{|S| odd}  f(S) g(U\\S)

And P = P_even + P_odd.
So M = 2*P_even - P and P - M = 2*P_odd.

Question: does the parity filter (multiplying by (-1)^|S|)
select exactly the terms of definite s-parity?
""")

for n in [4, 5]:
    c, svars, s_val, t_val = setup_skew(n)
    a, b = 0, 1
    U = [v for v in range(n) if v != a and v != b]

    # Compute per-S contributions
    terms_by_size = defaultdict(lambda: 0)
    for mask in range(1 << len(U)):
        S = [U[i] for i in range(len(U)) if mask & (1 << i)]
        R = [U[i] for i in range(len(U)) if not (mask & (1 << i))]
        ea = hp_paths(t_val, set(S)|{a}, end=a)
        bb = hp_paths(t_val, set(R)|{b}, start=b)
        contribution = expand(ea * bb)
        terms_by_size[len(S)] += contribution

    print(f"\n  n={n}: Contributions by |S|:")
    M_ab = 0
    for k in sorted(terms_by_size.keys()):
        val = expand(terms_by_size[k])
        sign_label = "+" if k % 2 == 0 else "-"
        M_ab += (-1)**k * val
        # Count s-degree of the contribution
        val_c0 = expand(val.subs(c, 0))
        print(f"    |S|={k} ({sign_label}): {len(val.as_ordered_terms())} terms, c=0: {val_c0}")

    M_ab = expand(M_ab)

    # Now: decompose M_ab by s-degree
    subs_neg = {svars[k]: -svars[k] for k in svars}
    M_neg = expand(M_ab.subs(subs_neg))
    even_part = expand((M_ab + M_neg) / 2)
    odd_part = expand((M_ab - M_neg) / 2)

    exp_parity = n % 2  # 0 = even in s, 1 = odd in s
    print(f"\n    n={n}: M has parity (-1)^(n-2) = {(-1)**(n-2)}")
    print(f"    Even part of M: {'ZERO' if even_part == 0 else f'{len(even_part.as_ordered_terms())} terms'}")
    print(f"    Odd part of M: {'ZERO' if odd_part == 0 else f'{len(odd_part.as_ordered_terms())} terms'}")

# ============================================================
# Part 3: The beautiful structure at c=0
# ============================================================
print("\n" + "=" * 70)
print("PART 3: M at c=0 — pure skew-symmetric world")
print("=" * 70)
print("""
At c=0, arc weights are purely skew: t_ij = s_ij = -t_ji.
This is the "anti-tournament" where arcs have real weights
summing to 0 on each pair.

In this world, M[a,b] is a homogeneous polynomial of degree n-2
in the skew variables. Each monomial represents a 2-path cover
with weight = product of skew arc weights.
""")

for n in [4, 5]:
    c, svars, s_val, t_val = setup_skew(n)

    # Define pure skew t-function (c=0)
    def t_skew(i, j): return s_val(i, j)

    a, b = 0, 1
    M_ab_0 = compute_transfer(t_skew, n, a, b, signed=True)
    M_ba_0 = compute_transfer(t_skew, n, b, a, signed=True)
    P_ab_0 = compute_transfer(t_skew, n, a, b, signed=False)
    P_ba_0 = compute_transfer(t_skew, n, b, a, signed=False)

    print(f"\n  n={n} at c=0:")
    print(f"    M[a,b] = {M_ab_0}")
    print(f"    M[b,a] = {M_ba_0}")
    print(f"    M[a,b] = M[b,a]? {expand(M_ab_0 - M_ba_0) == 0}")
    print(f"    P[a,b] = {P_ab_0}")
    print(f"    P[b,a] = {P_ba_0}")
    print(f"    P[a,b] = P[b,a]? {expand(P_ab_0 - P_ba_0) == 0}")

    if expand(P_ab_0 - P_ba_0) != 0:
        diff = expand(P_ab_0 - P_ba_0)
        print(f"    P[a,b] - P[b,a] = {diff}")
        # Is the unsigned difference ODD in s? (so it vanishes under parity)
        subs_neg = {svars[k]: -svars[k] for k in svars}
        diff_neg = expand(diff.subs(subs_neg))
        is_odd = expand(diff + diff_neg) == 0
        is_even = expand(diff - diff_neg) == 0
        print(f"    P-diff is odd in s? {is_odd}")
        print(f"    P-diff is even in s? {is_even}")

# ============================================================
# Part 4: The "skew-symmetrization" operator
# ============================================================
print("\n" + "=" * 70)
print("PART 4: Skew-symmetrization — the key mechanism")
print("=" * 70)
print("""
Define the skew-symmetrization operator:
  [f]_skew = (1/2)(f(s) + (-1)^{n-2} f(-s))

This projects f onto the "correct parity" subspace:
  - n even: keeps only even-degree s-terms
  - n odd: keeps only odd-degree s-terms

Hypothesis: [P[a,b]]_skew = [P[b,a]]_skew = M[a,b] = M[b,a]

If true, then M inherits symmetry because the skew-symmetrization
of any function f(s) is symmetric under any operation that commutes
with s -> -s.
""")

for n in [4, 5]:
    c, svars, s_val, t_val = setup_skew(n)

    def t_skew(i, j): return s_val(i, j)

    a, b = 0, 1
    P_ab = compute_transfer(t_skew, n, a, b, signed=False)
    P_ba = compute_transfer(t_skew, n, b, a, signed=False)
    M_ab = compute_transfer(t_skew, n, a, b, signed=True)

    subs_neg = {svars[k]: -svars[k] for k in svars}

    P_ab_neg = expand(P_ab.subs(subs_neg))
    P_ba_neg = expand(P_ba.subs(subs_neg))

    parity = (-1)**(n-2)
    P_ab_skew = expand((P_ab + parity * P_ab_neg) / 2)
    P_ba_skew = expand((P_ba + parity * P_ba_neg) / 2)

    print(f"\n  n={n} at c=0:")
    print(f"    [P[a,b]]_skew = {P_ab_skew}")
    print(f"    [P[b,a]]_skew = {P_ba_skew}")
    print(f"    M[a,b]        = {M_ab}")
    print(f"    [P[a,b]]_skew = M[a,b]? {expand(P_ab_skew - M_ab) == 0}")
    print(f"    [P[a,b]]_skew = [P[b,a]]_skew? {expand(P_ab_skew - P_ba_skew) == 0}")

# ============================================================
# Part 5: WHY does skew-symmetrization of P give a symmetric result?
# ============================================================
print("\n" + "=" * 70)
print("PART 5: The deep reason — swapping a<->b vs flipping s->-s")
print("=" * 70)
print("""
Claim: swapping a<->b in P[a,b] is EQUIVALENT to flipping s->-s
(up to the alternating sign from (-1)^|S|).

If P_{swap}(s) = P_{original}(-s), then:
  [P[a,b]]_skew = [P_{swap}]_skew = [P[b,a]]_skew

and skew-symmetrization makes them equal.

Let's test: does swapping a<->b in P correspond to s -> -s?
""")

for n in [4, 5]:
    c, svars, s_val, t_val = setup_skew(n)

    def t_skew(i, j): return s_val(i, j)

    a, b = 0, 1
    U = [v for v in range(n) if v != a and v != b]

    # Compute individual terms for P[a,b] and P[b,a]
    print(f"\n  n={n}: Per-subset comparison P[0,1](S) vs P[1,0](S)")
    for mask in range(1 << len(U)):
        S = [U[i] for i in range(len(U)) if mask & (1 << i)]
        R = [U[i] for i in range(len(U)) if not (mask & (1 << i))]

        ea = hp_paths(t_skew, set(S)|{a}, end=a)
        bb = hp_paths(t_skew, set(R)|{b}, start=b)
        eb = hp_paths(t_skew, set(S)|{b}, end=b)
        ba = hp_paths(t_skew, set(R)|{a}, start=a)

        term_ab = expand(ea * bb)
        term_ba = expand(eb * ba)

        subs_neg = {svars[k]: -svars[k] for k in svars}
        term_ab_neg = expand(term_ab.subs(subs_neg))

        # Is term_ba = (-1)^(n-2) * term_ab(-s)?
        relation = expand(term_ba - (-1)**(n-2) * term_ab_neg)

        print(f"    S={S}: E_a·B_b={term_ab}, E_b·B_a={term_ba}, "
              f"(-1)^(n-2)·(E_a·B_b)(-s)={expand((-1)**(n-2)*term_ab_neg)}, "
              f"match? {relation == 0}")

# ============================================================
# Part 6: The path reversal interpretation
# ============================================================
print("\n" + "=" * 70)
print("PART 6: Path reversal gives the bridge")
print("=" * 70)
print("""
A Hamiltonian path P: v0 -> v1 -> ... -> vk has weight
  w(P) = s_{v0,v1} * s_{v1,v2} * ... * s_{v_{k-1},vk}

The REVERSED path P^rev: vk -> ... -> v1 -> v0 has weight
  w(P^rev) = s_{vk,v_{k-1}} * ... * s_{v1,v0}
           = (-s_{v_{k-1},vk}) * ... * (-s_{v0,v1})
           = (-1)^k * w(P)

So: reversing a path of length k multiplies weight by (-1)^k.

Now E_a(S+a) counts paths ENDING at a through S+a (length |S|).
    B_b(R+b) counts paths STARTING at b through R+b (length |R|).

Swapping a<->b gives:
    E_b(S+b) counts paths ENDING at b through S+b.
    B_a(R+a) counts paths STARTING at a through R+a.

By path reversal:
    E_a(S+a) = sum of path weights ending at a
    B_a(S+a) = sum of path weights starting at a
             = sum over same vertex set, but starting instead of ending at a
             = (-1)^|S| * E_a(S+a) at c=0

Wait, is that true? Let me check.
""")

for n in [4, 5]:
    c, svars, s_val, t_val = setup_skew(n)

    def t_skew(i, j): return s_val(i, j)

    a = 0
    U = [v for v in range(n) if v != a]

    print(f"\n  n={n}: Testing E_a(S+a) vs (-1)^|S| B_a(S+a) at c=0")
    for size in range(len(U)+1):
        from itertools import combinations
        for combo in combinations(U, size):
            S = set(combo)
            ea = hp_paths(t_skew, S|{a}, end=a)
            ba = hp_paths(t_skew, S|{a}, start=a)
            expected = (-1)**size * ea
            match = expand(ba - expected) == 0
            if not match:
                print(f"    S={sorted(S)}: E_a={ea}, B_a={ba}, (-1)^|S|*E_a={expected}, MISMATCH")
            else:
                print(f"    S={sorted(S)}: E_a={ea}, B_a=(-1)^{size}*E_a [OK]")

# ============================================================
# Part 7: The complete picture
# ============================================================
print("\n" + "=" * 70)
print("PART 7: Assembling the proof")
print("=" * 70)
print("""
From Part 6: at c=0, B_v(S+v) = (-1)^|S| * E_v(S+v) for ANY vertex v.

This is because reversing a Hamiltonian path of length |S| through S+v
multiplies the weight by (-1)^|S| (each of the |S| skew arc weights
changes sign under reversal).

Now:
  M[a,b] = sum_S (-1)^|S| * E_a(S+a) * B_b(R+b)
  M[b,a] = sum_S (-1)^|S| * E_b(S+b) * B_a(R+a)

Using B_v = (-1)^|S| E_v (at c=0):
  B_b(R+b) = (-1)^|R| * E_b(R+b) = (-1)^{n-2-|S|} * E_b(R+b)
  B_a(R+a) = (-1)^|R| * E_a(R+a) = (-1)^{n-2-|S|} * E_a(R+a)

So:
  M[a,b] = sum_S (-1)^|S| * E_a(S+a) * (-1)^{n-2-|S|} * E_b(R+b)
         = (-1)^{n-2} * sum_S (-1)^{2|S|} * E_a(S+a) * E_b(R+b)
         = (-1)^{n-2} * sum_S E_a(S+a) * E_b(R+b)

Similarly:
  M[b,a] = (-1)^{n-2} * sum_S E_b(S+b) * E_a(R+a)

But sum_S E_a(S+a) * E_b(R+b) vs sum_S E_b(S+b) * E_a(R+a):
these are the SAME sum, just with S <-> R swapped!
(Replacing S by U\\S in the second gives sum_R E_b(R+b) * E_a(S+a) = same.)

Therefore M[a,b] = M[b,a] at c=0!

At general c: each E and B factor is a polynomial in c and s.
The key identity B_v(S+v) = (-1)^|S| E_v(S+v) only holds at c=0.
At general c, it becomes B_v(S+v, c, s) = E_v(S+v, c, -s).

This is because reversing a path changes each arc weight from
(c/2 + s_ij) to (c/2 - s_ij) = (c/2 + s_ji).

So B_v(S+v) = E_v(S+v)(c, -s) at general c.

Let me verify and then prove the full result.
""")

# ============================================================
# Part 8: Verify B_v = E_v(c, -s) at general c
# ============================================================
print("\n" + "=" * 70)
print("PART 8: The general c identity B_v(S+v) = E_v(S+v)(c,-s)")
print("=" * 70)

for n in [4, 5]:
    c, svars, s_val, t_val = setup_skew(n)

    a = 0
    U = [v for v in range(n) if v != a]

    subs_neg = {svars[k]: -svars[k] for k in svars}

    print(f"\n  n={n}:")
    all_match = True
    for size in range(len(U)+1):
        from itertools import combinations
        for combo in combinations(U, size):
            S = set(combo)
            ea = hp_paths(t_val, S|{a}, end=a)
            ba = hp_paths(t_val, S|{a}, start=a)
            ea_neg_s = expand(ea.subs(subs_neg))
            match = expand(ba - ea_neg_s) == 0
            if not match:
                print(f"    S={sorted(S)}: B_a != E_a(c,-s)! MISMATCH")
                print(f"      B_a = {ba}")
                print(f"      E_a(c,-s) = {ea_neg_s}")
                all_match = False

    if all_match:
        print(f"    ALL MATCH: B_v(S+v; c,s) = E_v(S+v; c,-s) ✓")

# ============================================================
# Part 9: The full proof
# ============================================================
print("\n" + "=" * 70)
print("PART 9: PROOF OF TRANSFER MATRIX SYMMETRY (all n, all c)")
print("=" * 70)
print("""
THEOREM. For any c-tournament (t_ij + t_ji = c, uniform),
M[a,b] = M[b,a] for all vertex pairs a,b.

PROOF.

Step 1: Path reversal identity.
  For any vertex v and vertex set S+v, define:
    E_v(S+v) = sum of path weights through S+v ending at v
    B_v(S+v) = sum of path weights through S+v starting at v

  Reversing a path v0 -> v1 -> ... -> vk changes each arc weight:
    t_{vi, v_{i+1}} = c/2 + s_{vi, v_{i+1}}
  becomes
    t_{v_{i+1}, vi} = c/2 + s_{v_{i+1}, vi} = c/2 - s_{vi, v_{i+1}}

  This is exactly the substitution s -> -s applied to each arc weight.
  Therefore: B_v(S+v; c, s) = E_v(S+v; c, -s).    [PROVED]

Step 2: Apply to M[a,b].
  M[a,b] = sum_{S ⊆ U} (-1)^|S| E_a(S+a; c,s) * B_b(R+b; c,s)

  By Step 1: B_b(R+b; c,s) = E_b(R+b; c,-s).

  M[a,b] = sum_S (-1)^|S| E_a(S+a; c,s) * E_b(R+b; c,-s)

Step 3: Apply to M[b,a].
  M[b,a] = sum_S (-1)^|S| E_b(S+b; c,s) * B_a(R+a; c,s)
         = sum_S (-1)^|S| E_b(S+b; c,s) * E_a(R+a; c,-s)

Step 4: Relabel S <-> R in M[b,a].
  Replace S by U\\S (so R becomes S, and |S| becomes |U\\S| = |U|-|S|):

  M[b,a] = sum_S (-1)^{|U|-|S|} E_b(R+b; c,s) * E_a(S+a; c,-s)
         = (-1)^{|U|} sum_S (-1)^{|S|} E_b(R+b; c,s) * E_a(S+a; c,-s)

  where |U| = n-2.

Step 5: Compare M[a,b] and M[b,a].
  M[a,b] = sum_S (-1)^|S| E_a(S+a; c,s) * E_b(R+b; c,-s)
  M[b,a] = (-1)^{n-2} sum_S (-1)^|S| E_a(S+a; c,-s) * E_b(R+b; c,s)

  These are the SAME expression with s and -s swapped!
  That is: M[b,a](c,s) = (-1)^{n-2} M[a,b](c,-s).

  This is the T^op equivalence: M_{T^op} = (-1)^{n-2} M_T.
  Since T^op corresponds to s -> -s, we get:
    M[b,a] = (-1)^{n-2} M[a,b](c,-s)

  But we also derived from path reversal that:
    M_{T^op}[a,b] = (-1)^{n-2} M_T[b,a]   (*)

  (*) follows from the same argument applied to M[a,b](c,-s):
    M[a,b](c,-s) = sum_S (-1)^|S| E_a(S+a;c,-s) * E_b(R+b;c,s)

  And comparing with Step 3:
    M[b,a](c,s) = sum_S (-1)^{|U|-|S|} E_b(R+b;c,s) * E_a(S+a;c,-s)
                = (-1)^{n-2} M[a,b](c,-s)

  So M[b,a] = (-1)^{n-2} M[a,b](c,-s).

  WAIT — this gives M[b,a] in terms of M[a,b](c,-s), NOT M[a,b](c,s).
  For M[a,b] = M[b,a] we would need:
    M[a,b](c,s) = (-1)^{n-2} M[a,b](c,-s)

  i.e., M has definite parity (-1)^{n-2} in s.
  This is exactly the T^op equivalence, which is what we WANT to prove.

  *** The argument is CIRCULAR. ***

  The path reversal gives M[b,a](c,s) = (-1)^{n-2} M[a,b](c,-s).
  This is EQUIVALENT to M[a,b]=M[b,a], but does not PROVE either one
  without an independent proof that M has definite parity.

  The CONTENT of the symmetry theorem is:
    M[a,b](c,s) has definite parity (-1)^{n-2} in s.

  OR EQUIVALENTLY:
    M[a,b](c,s) = M[b,a](c,s).

  The path reversal shows these are equivalent, but proves neither.
  We need a DIFFERENT argument to establish the parity.
""")

# ============================================================
# Part 10: What DOES the parity mean concretely?
# ============================================================
print("\n" + "=" * 70)
print("PART 10: What does definite parity mean concretely?")
print("=" * 70)

for n in [4, 5]:
    c_sym, svars, s_val, t_val = setup_skew(n)
    a, b = 0, 1

    M_ab = compute_transfer(t_val, n, a, b, signed=True)

    # Decompose by c-power and s-degree
    p = Poly(M_ab, c_sym)

    print(f"\n  n={n}: M[a,b] decomposed by (c-power, s-degree):")
    for c_power in range(p.degree() + 1):
        coeff = expand(p.nth(c_power))
        if coeff == 0: continue

        # Determine s-degree(s)
        subs_neg = {svars[k]: -svars[k] for k in svars}
        coeff_neg = expand(coeff.subs(subs_neg))
        is_even = expand(coeff - coeff_neg) == 0
        is_odd = expand(coeff + coeff_neg) == 0

        expected_parity = "even" if (n-2) % 2 == 0 else "odd"

        parity_str = "even" if is_even else ("odd" if is_odd else "MIXED")
        status = "✓" if parity_str == expected_parity or coeff == 0 else "✗"

        print(f"    c^{c_power}: s-parity={parity_str} {status} | {coeff}")

print("\n" + "=" * 70)
print("CONCLUSIONS")
print("=" * 70)
print("""
1. The "positive" transfer matrix P[a,b] (without (-1)^|S| signs) is
   NOT symmetric. The alternating signs are essential.

2. Path reversal establishes: B_v(S+v; c,s) = E_v(S+v; c,-s).
   This relates M[b,a](c,s) to M[a,b](c,-s) with a factor (-1)^{n-2}.
   But it does NOT directly prove M[a,b] = M[b,a] — it shows the
   equivalence of symmetry and definite s-parity.

3. The CONTENT of the symmetry theorem is that every coefficient of
   c^k in M[a,b] has s-parity equal to (-1)^{n-2} (even for n even,
   odd for n odd). The "wrong parity" terms all cancel.

4. The signed inclusion-exclusion (-1)^|S| acts as a parity FILTER:
   it selects exactly the terms with the correct s-parity. This is
   the bridge between "positivity" (all terms present in P) and
   "skew-symmetry" (only correct-parity terms survive in M).

5. The proof must show: the wrong-parity s-terms cancel in the
   sum over subsets S with alternating signs. This cancellation
   is forced by the tournament constraint t_ij + t_ji = c.
""")
