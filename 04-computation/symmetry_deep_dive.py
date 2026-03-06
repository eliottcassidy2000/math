"""
Deep dive into WHY M[a,b] = M[b,a].

Key insight from top_equivalence_proof.py:
  M_{T^op}[a,b] = (-1)^{n-2} * M_T[a,b]  (verified n=4,...,7)

But there's a TRIVIAL re-indexing identity:
  M_T[b,a] = (-1)^{n-2} * M_{T^op}[a,b]   (always true, by definition)

Combining: M_T[b,a] = (-1)^{n-2} * (-1)^{n-2} * M_T[a,b] = M_T[a,b].

So T^op equivalence <=> symmetry. They're the SAME claim.

This script explores WHY it holds by analyzing the combinatorial pairing structure.
"""

from itertools import permutations, combinations
from sympy import symbols, expand, Poly, degree
from collections import defaultdict

# ==========================================================
# Part 1: Verify the re-indexing identity (trivial direction)
# ==========================================================
print("=" * 70)
print("SYMMETRY DEEP DIVE")
print("=" * 70)

print("\n--- Part 1: Re-indexing identity verification ---")
print("""
THEOREM (trivial): M_T[b,a] = (-1)^{n-2} * M_{T^op}[a,b].

PROOF:
  M_T[b,a] = sum_{S' ⊆ U} (-1)^|S'| E_b^T(S'∪{b}) B_a^T(R'∪{a})

  Setting S' = R (complement), R' = S:
  = sum_S (-1)^{n-2-|S|} E_b^T(R∪{b}) B_a^T(S∪{a})
  = (-1)^{n-2} sum_S (-1)^|S| B_a^T(S∪{a}) E_b^T(R∪{b})
  = (-1)^{n-2} M_{T^op}[a,b]    [since E_a^{T^op} = B_a^T, B_b^{T^op} = E_b^T]  ∎

So proving M_{T^op}[a,b] = (-1)^{n-2} M_T[a,b] is EQUIVALENT to proving symmetry.
""")

# ==========================================================
# Part 2: The "EB vs BE" formulation
# ==========================================================
print("--- Part 2: EB vs BE formulation ---")
print("""
Symmetry means:
  sum_S (-1)^|S| E_a(S∪{a}) B_b(R∪{b}) = sum_S (-1)^|S| E_b(S∪{b}) B_a(R∪{a})

Define F(a,b) = sum_S (-1)^|S| E_a(S∪{a}) B_b(R∪{b}).
Define G(a,b) = sum_S (-1)^|S| B_a(S∪{a}) E_b(R∪{b}).

The T^op equivalence says F(a,b) = (-1)^{n-2} G(a,b).
The re-indexing says F(b,a) = (-1)^{n-2} G(a,b).

So F(a,b) = F(b,a) ⟺ F(a,b) = (-1)^{n-2} G(a,b).

The CONTENT of the symmetry is: the "E-first-B-second" pairing equals
(-1)^{n-2} times the "B-first-E-second" pairing.
""")

# ==========================================================
# Part 3: Symbolic n=4 — trace the pairing structure
# ==========================================================
print("--- Part 3: n=4 symbolic pairing analysis ---")

n = 4
verts = list(range(n))
# Create symbolic tournament variables
tvars = {}
for i in range(n):
    for j in range(i+1, n):
        tvars[(i,j)] = symbols(f't{i}{j}')

def t_sym(i, j):
    if i == j: return 0
    if i < j: return tvars[(i,j)]
    return 1 - tvars[(min(i,j), max(i,j))]

def ham_paths_through(T_func, vertex_set, start=None, end=None):
    """Count Hamiltonian paths in T[vertex_set] with optional start/end constraints."""
    vlist = sorted(vertex_set)
    k = len(vlist)
    if k == 0:
        return 0
    if k == 1:
        if start is not None and vlist[0] != start: return 0
        if end is not None and vlist[0] != end: return 0
        return 1

    total = 0
    for perm in permutations(vlist):
        if start is not None and perm[0] != start: continue
        if end is not None and perm[-1] != end: continue
        prod = 1
        for i in range(len(perm)-1):
            prod *= T_func(perm[i], perm[i+1])
        total += prod
    return expand(total)

a, b = 0, 1
U = [v for v in verts if v != a and v != b]

print(f"  n={n}, a={a}, b={b}, U={U}")
print(f"  |U| = {len(U)}, n-2 = {n-2}, (-1)^(n-2) = {(-1)**(n-2)}")

# Compute F(a,b) = M[a,b] term by term
print(f"\n  Individual terms of M[{a},{b}]:")
M_ab = 0
terms_F = {}
for mask in range(1 << len(U)):
    S = [U[i] for i in range(len(U)) if mask & (1 << i)]
    R = [U[i] for i in range(len(U)) if not (mask & (1 << i))]
    sign = (-1) ** len(S)

    E_a_val = ham_paths_through(t_sym, set(S) | {a}, end=a)
    B_b_val = ham_paths_through(t_sym, set(R) | {b}, start=b)

    term = expand(sign * E_a_val * B_b_val)
    terms_F[tuple(sorted(S))] = (sign, E_a_val, B_b_val, term)
    M_ab = expand(M_ab + term)

    print(f"    S={str(S):12s}  sign={sign:+d}  E_{a}={E_a_val}  B_{b}={B_b_val}  term={term}")

print(f"\n  M[{a},{b}] = {M_ab}")

# Now compute G(a,b) = sum (-1)^|S| B_a(S∪{a}) E_b(R∪{b})
print(f"\n  Individual terms of G({a},{b}) [B-first, E-second]:")
G_ab = 0
terms_G = {}
for mask in range(1 << len(U)):
    S = [U[i] for i in range(len(U)) if mask & (1 << i)]
    R = [U[i] for i in range(len(U)) if not (mask & (1 << i))]
    sign = (-1) ** len(S)

    B_a_val = ham_paths_through(t_sym, set(S) | {a}, start=a)
    E_b_val = ham_paths_through(t_sym, set(R) | {b}, end=b)

    term = expand(sign * B_a_val * E_b_val)
    terms_G[tuple(sorted(S))] = (sign, B_a_val, E_b_val, term)
    G_ab = expand(G_ab + term)

    print(f"    S={str(S):12s}  sign={sign:+d}  B_{a}={B_a_val}  E_{b}={E_b_val}  term={term}")

print(f"\n  G({a},{b}) = {G_ab}")
print(f"  (-1)^(n-2) * G({a},{b}) = {expand((-1)**(n-2) * G_ab)}")
print(f"  F = (-1)^(n-2)*G? {expand(M_ab - (-1)**(n-2) * G_ab) == 0}")

# ==========================================================
# Part 4: Term-by-term pairing
# ==========================================================
print(f"\n--- Part 4: Term-by-term pairing (S <-> complement) ---")
print("  For each S, compare F's term at S with G's term at complement(S):")

for mask in range(1 << len(U)):
    S = tuple(sorted(U[i] for i in range(len(U)) if mask & (1 << i)))
    R = tuple(sorted(U[i] for i in range(len(U)) if not (mask & (1 << i))))

    f_sign, f_E, f_B, f_term = terms_F[S]
    g_sign, g_B, g_E, g_term = terms_G[R]  # G at complement

    # F(S) has sign (-1)^|S|, G(R) has sign (-1)^|R| = (-1)^{n-2-|S|}
    # Ratio of signs: (-1)^|R| / (-1)^|S| = (-1)^{n-2-2|S|} = (-1)^{n-2} (since (-1)^{-2|S|}=1)

    print(f"  S={list(S)}, R={list(R)}:")
    print(f"    F term (S): {f_sign:+d} * E_{a}({set(S)|{a}}) * B_{b}({set(R)|{b}}) = {f_term}")
    print(f"    G term (R): {g_sign:+d} * B_{a}({set(R)|{a}}) * E_{b}({set(S)|{b}}) = {g_term}")

    # The E and B functions:
    # F(S): E_a(S∪{a}) * B_b(R∪{b})   — paths ending at a through S∪{a}, paths starting from b through R∪{b}
    # G(R): B_a(R∪{a}) * E_b(S∪{b})   — paths starting from a through R∪{a}, paths ending at b through S∪{b}

    # These are the SAME sets! E_a(S∪{a})*B_b(R∪{b}) vs B_a(R∪{a})*E_b(S∪{b})
    # The vertex sets are identical (S∪{a} and R∪{b}), but the roles flip:
    # "end at a in S∪{a}, start from b in R∪{b}" vs "start from a in R∪{a}, end at b in S∪{b}"

    # So the question reduces to: for each partition S∪{a} | R∪{b}:
    # Does E_a(S∪{a}) * B_b(R∪{b}) = B_a(R∪{a}) * E_b(S∪{b})?
    # (up to sign correction)

    ratio_match = expand(f_term - (-1)**(n-2) * g_term) == 0
    print(f"    F(S) = (-1)^(n-2) * G(R)? {ratio_match}")
    if not ratio_match:
        print(f"    DIFF: {expand(f_term - (-1)**(n-2) * g_term)}")

# ==========================================================
# Part 5: The key question — individual term matching?
# ==========================================================
print(f"\n--- Part 5: Does each (S, complement(S)) pair match individually? ---")
print("""
For symmetry, we need: sum_S (-1)^|S| [E_a(S∪a)B_b(R∪b) - E_b(S∪b)B_a(R∪a)] = 0

This could happen because:
(A) Each term vanishes: E_a(S∪a)B_b(R∪b) = E_b(S∪b)B_a(R∪a) for each S
(B) Cancellation across different S values
""")

# Check (A)
individual_match = True
for mask in range(1 << len(U)):
    S = [U[i] for i in range(len(U)) if mask & (1 << i)]
    R = [U[i] for i in range(len(U)) if not (mask & (1 << i))]

    E_a_val = ham_paths_through(t_sym, set(S) | {a}, end=a)
    B_b_val = ham_paths_through(t_sym, set(R) | {b}, start=b)
    E_b_val = ham_paths_through(t_sym, set(S) | {b}, end=b)
    B_a_val = ham_paths_through(t_sym, set(R) | {a}, start=a)

    diff = expand(E_a_val * B_b_val - E_b_val * B_a_val)
    print(f"  S={S}: E_a*B_b - E_b*B_a = {diff}")
    if diff != 0:
        individual_match = False

print(f"\n  Individual terms match? {individual_match}")
if not individual_match:
    print("  => Symmetry arises from CANCELLATION across terms, not individual matching!")

# ==========================================================
# Part 6: The cancellation structure
# ==========================================================
print(f"\n--- Part 6: Cancellation structure ---")
print("  Computing D(S) = E_a(S∪a)*B_b(R∪b) - E_b(S∪b)*B_a(R∪a) for each S:")

D_values = {}
for mask in range(1 << len(U)):
    S = tuple(sorted(U[i] for i in range(len(U)) if mask & (1 << i)))
    R = [U[i] for i in range(len(U)) if not (mask & (1 << i))]

    E_a_val = ham_paths_through(t_sym, set(S) | {a}, end=a)
    B_b_val = ham_paths_through(t_sym, set(R) | {b}, start=b)
    E_b_val = ham_paths_through(t_sym, set(S) | {b}, end=b)
    B_a_val = ham_paths_through(t_sym, set(R) | {a}, start=a)

    D = expand(E_a_val * B_b_val - E_b_val * B_a_val)
    D_values[S] = D
    print(f"  S={str(list(S)):12s} |S|={len(S)}: D(S) = {D}")

print(f"\n  Weighted sum: sum (-1)^|S| D(S) = {expand(sum((-1)**len(S) * D for S, D in D_values.items()))}")

# Check: does D(S) + D(complement(S)) = 0?
print("\n  Complement pairing D(S) + D(R):")
seen = set()
for mask in range(1 << len(U)):
    S = tuple(sorted(U[i] for i in range(len(U)) if mask & (1 << i)))
    R = tuple(sorted(U[i] for i in range(len(U)) if not (mask & (1 << i))))
    if S in seen: continue
    seen.add(S)
    seen.add(R)
    d_sum = expand(D_values[S] + D_values.get(R, 0))
    d_diff = expand(D_values[S] - D_values.get(R, 0))
    print(f"  D({list(S)}) + D({list(R)}) = {d_sum}")
    print(f"  D({list(S)}) - D({list(R)}) = {d_diff}")

# ==========================================================
# Part 7: n=5 numerical verification of cancellation pattern
# ==========================================================
print(f"\n--- Part 7: n=5 cancellation pattern (numerical) ---")
import random
random.seed(42)

n5 = 5
for trial in range(3):
    # Random tournament
    T = [[0]*n5 for _ in range(n5)]
    for i in range(n5):
        for j in range(i+1, n5):
            if random.random() < 0.5:
                T[i][j] = 1
            else:
                T[j][i] = 1

    a5, b5 = 0, 1
    U5 = [v for v in range(n5) if v != a5 and v != b5]

    def hp_num(T, vset, start=None, end=None):
        vl = sorted(vset)
        count = 0
        for p in permutations(vl):
            if start is not None and p[0] != start: continue
            if end is not None and p[-1] != end: continue
            ok = True
            for i in range(len(p)-1):
                if T[p[i]][p[i+1]] != 1:
                    ok = False
                    break
            if ok: count += 1
        return count

    print(f"\n  Trial {trial+1}: T = {[sum(row) for row in T]} (score seq)")

    M_ab5 = 0
    M_ba5 = 0
    for mask in range(1 << len(U5)):
        S = [U5[i] for i in range(len(U5)) if mask & (1 << i)]
        R = [U5[i] for i in range(len(U5)) if not (mask & (1 << i))]
        sign = (-1) ** len(S)

        ea = hp_num(T, set(S)|{a5}, end=a5)
        bb = hp_num(T, set(R)|{b5}, start=b5)
        eb = hp_num(T, set(S)|{b5}, end=b5)
        ba = hp_num(T, set(R)|{a5}, start=a5)

        D = ea*bb - eb*ba
        M_ab5 += sign * ea * bb
        M_ba5 += sign * eb * ba

        if D != 0:
            print(f"    S={S}: E_a={ea}, B_b={bb}, E_b={eb}, B_a={ba}, D={D}, sign*D={sign*D}")

    print(f"    M[{a5},{b5}]={M_ab5}, M[{b5},{a5}]={M_ba5}, symmetric={M_ab5==M_ba5}")

# ==========================================================
# Part 8: Path-level analysis — what pairs to create symmetry?
# ==========================================================
print(f"\n--- Part 8: Path-level pairing at n=4 ---")
print("""
M[a,b] counts SIGNED PRODUCTS of path pairs.
Each product is: (-1)^|S| * (path ending at a through S∪{a}) * (path starting from b through R∪{b})

The two paths together cover all of V (since S∪{a} and R∪{b} partition V, with S∪R = U, {a}∩{b}=∅).
So each term corresponds to a PAIR of paths that together visit every vertex exactly once.

This is equivalent to: choose a set S, then an ordering of S∪{a} ending at a,
and an ordering of R∪{b} starting at b, with all arcs going the right direction.

The pair (path_1_ending_at_a, path_2_starting_from_b) is essentially a
"2-path cover" of V where path_1 ends at a and path_2 starts at b.
""")

# Enumerate all 2-path covers at n=4
n = 4
a, b = 0, 1
U = [v for v in range(n) if v != a and v != b]

# For a random n=4 tournament
random.seed(123)
T4 = [[0]*n for _ in range(n)]
for i in range(n):
    for j in range(i+1, n):
        if random.random() < 0.5:
            T4[i][j] = 1
        else:
            T4[j][i] = 1

print(f"  Tournament: {[(i, [j for j in range(n) if T4[i][j]]) for i in range(n)]}")

covers_ab = []  # (S, path_ending_a, path_starting_b, sign)
covers_ba = []  # (S, path_ending_b, path_starting_a, sign)

for mask in range(1 << len(U)):
    S = [U[i] for i in range(len(U)) if mask & (1 << i)]
    R = [U[i] for i in range(len(U)) if not (mask & (1 << i))]
    sign = (-1) ** len(S)

    # M[a,b] terms: paths ending at a in S∪{a}, paths starting from b in R∪{b}
    for p1 in permutations(sorted(set(S)|{a})):
        if p1[-1] != a: continue
        ok1 = all(T4[p1[i]][p1[i+1]] for i in range(len(p1)-1))
        if not ok1: continue
        for p2 in permutations(sorted(set(R)|{b})):
            if p2[0] != b: continue
            ok2 = all(T4[p2[i]][p2[i+1]] for i in range(len(p2)-1))
            if not ok2: continue
            covers_ab.append((tuple(S), p1, p2, sign))

    # M[b,a] terms: paths ending at b in S∪{b}, paths starting from a in R∪{a}
    for p1 in permutations(sorted(set(S)|{b})):
        if p1[-1] != b: continue
        ok1 = all(T4[p1[i]][p1[i+1]] for i in range(len(p1)-1))
        if not ok1: continue
        for p2 in permutations(sorted(set(R)|{a})):
            if p2[0] != a: continue
            ok2 = all(T4[p2[i]][p2[i+1]] for i in range(len(p2)-1))
            if not ok2: continue
            covers_ba.append((tuple(S), p1, p2, sign))

print(f"\n  M[{a},{b}] has {len(covers_ab)} 2-path covers:")
for S, p1, p2, sign in covers_ab:
    print(f"    sign={sign:+d}  [{' -> '.join(map(str,p1))}] | [{' -> '.join(map(str,p2))}]  (S={list(S)})")

print(f"\n  M[{b},{a}] has {len(covers_ba)} 2-path covers:")
for S, p1, p2, sign in covers_ba:
    print(f"    sign={sign:+d}  [{' -> '.join(map(str,p1))}] | [{' -> '.join(map(str,p2))}]  (S={list(S)})")

val_ab = sum(sign for _, _, _, sign in covers_ab)
val_ba = sum(sign for _, _, _, sign in covers_ba)
print(f"\n  M[{a},{b}] = {val_ab}, M[{b},{a}] = {val_ba}")

# ==========================================================
# Part 9: The involution — is there a sign-reversing involution?
# ==========================================================
print(f"\n--- Part 9: Looking for an involution ---")
print("""
For M[a,b] = M[b,a], we need:
  sum covers_ab(sign) = sum covers_ba(sign)

Equivalently: sum (covers_ab - covers_ba)(sign) = 0.

Is there a BIJECTION phi: covers_ab -> covers_ba preserving signs?
Or a sign-reversing involution on the "excess" terms?
""")

# A natural bijection attempt: given a 2-path cover (p1 ending at a, p2 starting from b),
# can we rearrange to get (p1' ending at b, p2' starting from a)?

# The two paths together form a linear ordering of all vertices,
# split into two segments. The "break point" is between the two paths.
# Path 1: x_1 -> x_2 -> ... -> a  (vertices S∪{a})
# Path 2: b -> y_1 -> y_2 -> ... (vertices R∪{b})

# If we concatenate: x_1 -> ... -> a | b -> y_1 -> ... -> y_k
# This is almost a Hamiltonian path from x_1 to y_k, except there may
# be no arc from a to b.

# A natural map: REVERSE both paths:
# Path 1 reversed: a -> x_{|S|} -> ... -> x_1  (starts from a through S∪{a})
# Path 2 reversed: y_k -> ... -> y_1 -> b  (ends at b through R∪{b})
# This gives a cover for M[b,a]!

# But wait — reversing a path in a tournament means following REVERSE arcs.
# So the reversed path exists in T^op, not in T.
# This doesn't directly give a bijection.

# What about: swap the roles of a,b in the partition?
# Given S with cover (p1 through S∪{a} ending at a, p2 through R∪{b} starting from b),
# consider S' = ? with cover (p1' through S'∪{b} ending at b, p2' through R'∪{a} starting from a).

# The simplest attempt: keep the same paths but relabel.
# But a and b are specific vertices, not interchangeable.

# Let me look at what happens with CONCATENATION.
# If there IS an arc a->b, then p1|p2 = x1->...->a->b->y1->...->yk is a Hamiltonian path.
# This has sign (-1)^|S| in M[a,b].
# The SAME path, viewed as a cover for M[b,a] with S' = {y1,...,yk}:
# Path ending at b: x1->...->a->b, through {x1,...,a,b} wait, this doesn't work cleanly.

# Let me try a different approach: the "cut" bijection.
# A Hamiltonian path P = v1->v2->...->vn can be "cut" at any point to give a 2-path cover.
# Cut after position k: (v1->...->vk) and (v_{k+1}->...->vn).
# For M[a,b], we need vk = a and v_{k+1} = b, i.e., the cut is at the arc a->b.
# So the contribution of P to M[a,b] is (-1)^{k-1} if a->b is an arc of P.

# For M[b,a], we need vk = b and v_{k+1} = a, i.e., cut at arc b->a.

# Hmm, but M[a,b] also gets contributions from 2-path covers where a->b is NOT an arc!

# Actually, let me reconsider. M[a,b] = sum_S (-1)^|S| E_a(S∪{a}) B_b(R∪{b}).
# This involves TWO SEPARATE paths: one ending at a, one starting from b.
# They don't need to be connected by an arc a->b.

# The total contribution is from all (path1, path2) pairs where:
# - path1 goes through some subset S∪{a} and ends at a
# - path2 goes through the complement R∪{b} and starts from b
# - sign is (-1)^|S|

# Key observation: |S| = |path1| - 1 (since path1 has |S|+1 vertices)
# So sign = (-1)^{|path1|-1}

print("  Bijection analysis:")
print(f"  covers_ab sum = {val_ab}")
print(f"  covers_ba sum = {val_ba}")

# ==========================================================
# Part 10: The Hamiltonian path connection
# ==========================================================
print(f"\n--- Part 10: Hamiltonian path decomposition ---")
print("""
Every Hamiltonian path P = v1->v2->...->vn contributes to M[a,b] for specific (a,b).
If a = v_i and b = v_j with i < j, the path can be cut to contribute to M[a,b]:
  Path1 = v_1->...->v_i (ending at a), which is through {v_1,...,v_i}
  Path2 = v_j->...->v_n (starting from b), which is through {v_j,...,v_n}
  BUT: we also need {v_{i+1},...,v_{j-1}} to be in S or R consistently.

Actually, S = {v_1,...,v_{i-1}} ∪ {intermediate vertices assigned to a's side}
This is MORE COMPLEX than simple cutting.

Let me instead enumerate ALL Hamiltonian paths and see how they relate to M.
""")

ham_paths = []
for p in permutations(range(n)):
    ok = all(T4[p[i]][p[i+1]] for i in range(n-1))
    if ok:
        ham_paths.append(p)

print(f"  Hamiltonian paths in T: {len(ham_paths)}")
for p in ham_paths:
    print(f"    {' -> '.join(map(str, p))}")

# Compute M[a,b] directly
M_full = [[0]*n for _ in range(n)]
for aa in range(n):
    for bb in range(n):
        if aa == bb: continue
        UU = [v for v in range(n) if v != aa and v != bb]
        val = 0
        for mask in range(1 << len(UU)):
            S = [UU[i] for i in range(len(UU)) if mask & (1 << i)]
            R = [UU[i] for i in range(len(UU)) if not (mask & (1 << i))]
            sign = (-1) ** len(S)
            ea = hp_num(T4, set(S)|{aa}, end=aa)
            bbb = hp_num(T4, set(R)|{bb}, start=bb)
            val += sign * ea * bbb
        M_full[aa][bb] = val

print(f"\n  Transfer matrix M:")
for i in range(n):
    print(f"    {M_full[i]}")

print(f"\n  Symmetric? {all(M_full[i][j] == M_full[j][i] for i in range(n) for j in range(n) if i != j)}")
print(f"  Trace = {sum(M_full[i][i] for i in range(n))}")
print(f"  H(T) = {len(ham_paths)}")
print(f"  (-1)^n = {(-1)**n}")
print(f"  Expected trace: H if n odd, 0 if n even = {len(ham_paths) if n % 2 else 0}")

# ==========================================================
# Part 11: The matrix E^T Lambda B decomposition
# ==========================================================
print(f"\n--- Part 11: Cauchy-Binet structure ---")
print("""
M[a,b] = sum_S (-1)^|S| E_a(S∪{a}) B_b(R∪{b})

Think of this as a dot product:
  M[a,b] = <e_a, Lambda * b_b>

where e_a[S] = E_a(S∪{a}), b_b[S] = B_b(S∪{b}), Lambda[S,S'] = (-1)^|S| delta(S, U\S').

Actually: Lambda = (-1)^|S| matching S with complement R=U\S.

M[a,b] = sum_S (-1)^|S| e_a[S] * b_b[U\S]
        = sum_S e_a[S] * ((-1)^|S| b_b[U\S])

Define c_b[S] = (-1)^|U\S| b_b[S] = (-1)^{|U|-|S|} b_b[S].
Then (-1)^|S| b_b[U\S] = (-1)^|S| * (-1)^{|U|-|U\S|} c_b[U\S] * (-1)^{|U\S|-|U|+|S|}...

Let me just be concrete.
""")

# For n=4, |U|=2
print(f"  n=4, |U|=2, subsets of U = {U}:")
for mask in range(1 << len(U)):
    S = tuple(sorted(U[i] for i in range(len(U)) if mask & (1 << i)))
    R = tuple(sorted(U[i] for i in range(len(U)) if not (mask & (1 << i))))

    ea = hp_num(T4, set(S)|{a}, end=a)
    bb_val = hp_num(T4, set(R)|{b}, start=b)
    eb = hp_num(T4, set(S)|{b}, end=b)
    ba = hp_num(T4, set(R)|{a}, start=a)

    print(f"    S={list(S)}: E_0={ea}, B_1={bb_val}, E_1={eb}, B_0={ba}")

# The matrix formulation: index subsets by size
# E[size][a] = vector of E_a values for subsets of that size
# Similarly for B

print(f"\n  Organized by |S|:")
by_size = defaultdict(list)
for mask in range(1 << len(U)):
    S = tuple(sorted(U[i] for i in range(len(U)) if mask & (1 << i)))
    by_size[len(S)].append(S)

for sz in sorted(by_size.keys()):
    print(f"\n  |S| = {sz}, sign = {(-1)**sz}:")
    for S in by_size[sz]:
        R = tuple(sorted(set(U) - set(S)))
        ea = hp_num(T4, set(S)|{a}, end=a)
        bb_val = hp_num(T4, set(R)|{b}, start=b)
        eb = hp_num(T4, set(S)|{b}, end=b)
        ba = hp_num(T4, set(R)|{a}, start=a)
        print(f"    S={list(S)}: E_a*B_b={ea}*{bb_val}={ea*bb_val}, E_b*B_a={eb}*{ba}={eb*ba}, diff={ea*bb_val-eb*ba}")

print(f"\n" + "=" * 70)
print("CONCLUSIONS")
print("=" * 70)
print("""
1. M_T[b,a] = (-1)^{n-2} * M_{T^op}[a,b] is TRIVIAL (re-indexing).
2. M_{T^op}[a,b] = (-1)^{n-2} * M_T[a,b] is the NONTRIVIAL claim (= symmetry).
3. Individual terms do NOT match: E_a*B_b != E_b*B_a in general.
4. Symmetry arises from CANCELLATION in the alternating sum over subsets S.
5. The cancellation mechanism involves subtle interactions between:
   - Path counts in complementary vertex subsets
   - The alternating sign (-1)^|S|
   - The tournament structure T
""")
