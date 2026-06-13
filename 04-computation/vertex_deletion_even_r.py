#!/usr/bin/env python3
"""
Vertex deletion recurrence for M[a,b] — does it preserve even r-powers?

KEY APPROACH: Delete vertex v from U = V\{a,b}. Express:
M_n[a,b] = f(M_{n-1}[a,b], M_{n-1}[v,b], M_{n-1}[a,v], ...) + corrections

If ALL terms in the recurrence are even in r when the inductive hypothesis
says M_{n-1} entries are even in r, then M_n is even in r.

ALSO: The transfer recurrence from transfer_recurrence.py:
M_n[a,b] = B_b(V\{a}) - sum_v t(v,a) * M_{n-1}^{(a)}[v,b]

where M_{n-1}^{(a)} is the transfer matrix on V\{a}.

Instance: opus-2026-03-06-S24
"""

from itertools import permutations
from sympy import symbols, expand, Symbol, Poly, diff

def make_symbols(n):
    r = Symbol('r')
    s = {}
    for i in range(n):
        for j in range(n):
            if i < j:
                s[(i,j)] = Symbol(f's{i}{j}')
                s[(j,i)] = -s[(i,j)]
    return r, s

def edge_weight(i, j, r, s):
    if i == j: return 0
    return r + s[(i,j)]

def ham_paths_ending_at(vertex_set, target, r, s):
    vs = list(vertex_set)
    if len(vs) == 1: return 1
    total = 0
    for perm in permutations([v for v in vs if v != target]):
        path = list(perm) + [target]
        w = 1
        for k in range(len(path)-1): w *= edge_weight(path[k], path[k+1], r, s)
        total += w
    return expand(total)

def ham_paths_beginning_at(vertex_set, source, r, s):
    vs = list(vertex_set)
    if len(vs) == 1: return 1
    total = 0
    for perm in permutations([v for v in vs if v != source]):
        path = [source] + list(perm)
        w = 1
        for k in range(len(path)-1): w *= edge_weight(path[k], path[k+1], r, s)
        total += w
    return expand(total)

def compute_M(a, b, vertices, r, s):
    V = sorted(vertices)
    U = [v for v in V if v != a and v != b]
    total = 0
    for mask in range(2**len(U)):
        S = [U[i] for i in range(len(U)) if (mask >> i) & 1]
        R = [u for u in U if u not in S]
        sign = (-1)**len(S)
        Sa = set(S + [a])
        Rb = set(R + [b])
        Ea = ham_paths_ending_at(Sa, a, r, s)
        Bb = ham_paths_beginning_at(Rb, b, r, s)
        total += sign * Ea * Bb
    return expand(total)


def test_vertex_deletion_recurrence(n, a=0, b=1):
    """
    Express M_n[a,b] in terms of M_{n-1} entries.

    Delete vertex v from the tournament. The 2-path-covers of V
    can be partitioned by which path contains v:
    - v is in the E_a path (S-side): contributes to some M_{n-1} entries
    - v is in the B_b path (R-side): contributes to other M_{n-1} entries

    For each S ⊆ U:
    Case 1: v ∈ S (v is in E_a path through S∪{a})
      E_a(S∪{a}) = sum over paths ...→v→...→a through S∪{a}
      When we delete v, each path splits into:
      - A path through (S\{v})∪{a} ending at a, preceded by an edge to v
      - OR v is the starting vertex, contributing edge v→next

    Case 2: v ∉ S (v is in B_b path through R∪{b})
      Similar decomposition.

    Actually, let's just verify the ROW recurrence:
    M_n[a,b] = B_b(V\{a}) - sum_{v∈U} t(v,a) * M_{n-1}^{(a)}[v,b]
    """
    print(f"\n{'='*60}")
    print(f"VERTEX DELETION: n={n}, a={a}, b={b}")
    print(f"{'='*60}")

    r, s = make_symbols(n)
    V = list(range(n))
    U = [v for v in V if v != a and v != b]

    # Direct computation
    M_direct = compute_M(a, b, V, r, s)

    # Row recurrence: M_n[a,b] = B_b(V\{a}) - sum_{v∈U} t(v,a) * M_{n-1}^{(a)}[v,b]
    V_minus_a = [v for v in V if v != a]
    Bb_full = ham_paths_beginning_at(set(V_minus_a), b, r, s)

    recurrence_sum = 0
    for v in U:
        M_sub = compute_M(v, b, V_minus_a, r, s)
        recurrence_sum += edge_weight(v, a, r, s) * M_sub

    M_recurrence = expand(Bb_full - recurrence_sum)

    match = expand(M_direct - M_recurrence) == 0
    print(f"\n  Row recurrence verified: {match}")

    # Now analyze r-parity of each component
    print(f"\n  --- R-parity analysis of recurrence components ---")

    # B_b(V\{a}) has degree n-2 (n-1 vertices, n-2 edges)
    Bb_poly = Poly(Bb_full, r)
    Bb_odd = [d[0] for d, c in Bb_poly.as_dict().items() if d[0] % 2 == 1]
    Bb_even = [d[0] for d, c in Bb_poly.as_dict().items() if d[0] % 2 == 0]
    print(f"  B_b(V\\{{a}}): odd r-powers = {Bb_odd}, even r-powers = {Bb_even}")

    # Each t(v,a) * M_{n-1}^{(a)}[v,b] term
    for v in U:
        M_sub = compute_M(v, b, V_minus_a, r, s)
        term = expand(edge_weight(v, a, r, s) * M_sub)

        # Check if M_sub is even in r (inductive hypothesis)
        M_sub_neg = expand(M_sub.subs(r, -r))
        M_sub_even = expand(M_sub - M_sub_neg) == 0

        # Check r-parity of the full term t(v,a) * M_sub
        if term != 0:
            term_poly = Poly(term, r)
            term_odd = [d[0] for d, c in term_poly.as_dict().items() if d[0] % 2 == 1]
            term_even = [d[0] for d, c in term_poly.as_dict().items() if d[0] % 2 == 0]
        else:
            term_odd, term_even = [], []

        print(f"  t({v},{a}) * M^(a)[{v},{b}]: M_sub even? {M_sub_even}, "
              f"term odd r = {term_odd}, even r = {term_even}")

    # KEY OBSERVATION:
    # If M_{n-1}^{(a)}[v,b] is even in r (inductive hypothesis),
    # then t(v,a) * M_{n-1} has BOTH odd and even r-powers
    # (because t(v,a) = r + s_{va} shifts r-parity).
    #
    # So the recurrence sum has BOTH odd and even r-powers,
    # and B_b(V\{a}) also has both.
    # The cancellation happens between B_b and the sum.

    # What is the ODD part of B_b vs ODD part of sum?
    if Bb_full != 0:
        Bb_dict = Poly(Bb_full, r).as_dict()
        Bb_odd_part = sum(c * r**d[0] for d, c in Bb_dict.items() if d[0] % 2 == 1)
    else:
        Bb_odd_part = 0

    rec_full = expand(recurrence_sum)
    if rec_full != 0:
        rec_dict = Poly(rec_full, r).as_dict()
        rec_odd_part = sum(c * r**d[0] for d, c in rec_dict.items() if d[0] % 2 == 1)
    else:
        rec_odd_part = 0

    odd_cancel = expand(Bb_odd_part - rec_odd_part)
    print(f"\n  Odd-part(B_b) - Odd-part(sum) = {odd_cancel}")
    print(f"  Cancellation works: {odd_cancel == 0}")

    # DEEPER: Can we understand WHY the odd parts cancel?
    # Decompose: sum_{v∈U} t(v,a) * M^{(a)}[v,b]
    # = r * sum_v M^{(a)}[v,b] + sum_v s_{va} * M^{(a)}[v,b]
    #
    # The r * sum_v M^{(a)}[v,b] term shifts everything by r^1.
    # So odd parts of M^{(a)} become even parts (times r), and vice versa.
    #
    # The s_{va} * M^{(a)} term doesn't change r-parity.

    print(f"\n  --- Decomposition of recurrence sum ---")

    # r-part: r * sum_v M^{(a)}[v,b]
    col_sum = 0
    for v in U:
        M_sub = compute_M(v, b, V_minus_a, r, s)
        col_sum = expand(col_sum + M_sub)
    r_part = expand(r * col_sum)

    # s-part: sum_v s_{va} * M^{(a)}[v,b]
    s_part = 0
    for v in U:
        M_sub = compute_M(v, b, V_minus_a, r, s)
        s_part = expand(s_part + s[(v, a)] * M_sub)

    print(f"  sum_v M^(a)[v,{b}] = {col_sum}")
    print(f"  r * (sum) = {r_part}")
    print(f"  sum_v s_va * M^(a)[v,{b}] = {s_part}")
    print(f"  r_part + s_part = sum? {expand(r_part + s_part - recurrence_sum) == 0}")

    # So: M_n[a,b] = B_b - r * col_sum - s_part
    # The odd r-powers of M_n come from:
    # odd(B_b) - r * even(col_sum) - odd(s_part) = 0
    # (since r * even becomes odd, and s_part preserves parity)

    # If M_{n-1} is even in r, then col_sum is even in r.
    # So r * col_sum is ODD in r.
    # And s_part (sum of s_va * M_sub) is EVEN in r (s_va doesn't change parity).

    # Therefore: odd(M_n) = odd(B_b) - r * col_sum = 0
    # IFF: odd parts of B_b(V\{a}) = r * col_sum

    # This is a CONCRETE condition! Can we verify it?

    if col_sum != 0:
        r_col_dict = Poly(r * col_sum, r).as_dict()
        r_col_odd = sum(c * r**d[0] for d, c in r_col_dict.items() if d[0] % 2 == 1)
    else:
        r_col_odd = 0

    print(f"\n  KEY TEST: odd(B_b) = r * col_sum?")
    print(f"  odd(B_b) = {expand(Bb_odd_part)}")
    print(f"  r * col_sum (odd part only) = {expand(r_col_odd)}")
    print(f"  Match: {expand(Bb_odd_part - r_col_odd) == 0}")

    # If this works, then odd(M_n) = 0 follows from:
    # 1. M_{n-1} is even in r (inductive hypothesis)
    # 2. odd(B_b(V\{a})) = r * sum_v M^{(a)}[v,b]
    #
    # Condition 2 is an identity about boundary Hamiltonian path sums
    # and column sums of the transfer matrix!


def main():
    for n in [4, 5]:
        test_vertex_deletion_recurrence(n)

    print(f"\n\n{'='*60}")
    print("INDUCTIVE PROOF STRUCTURE")
    print(f"{'='*60}")
    print("""
CLAIM: If M_{n-1} is even in r, then M_n is even in r.

PROOF SKETCH (if the key identity holds):
  1. Row recurrence: M_n[a,b] = B_b(V\\{a}) - sum_v t(v,a) M^{(a)}[v,b]
  2. Decompose: sum_v t(v,a) M^{(a)}[v,b]
     = r * sum_v M^{(a)}[v,b] + sum_v s_{va} * M^{(a)}[v,b]
  3. By inductive hypothesis, each M^{(a)}[v,b] is even in r.
  4. So sum_v s_{va} * M^{(a)}[v,b] is even in r (s_{va} doesn't shift).
  5. And r * sum_v M^{(a)}[v,b] is ODD in r.
  6. The odd part of M_n = odd(B_b) - r * sum_v M^{(a)}[v,b].
  7. KEY IDENTITY: odd(B_b(V\\{a})) = r * sum_v M^{(a)}[v,b].
  8. If step 7 holds, then odd(M_n) = 0. QED.

So the proof reduces to proving the KEY IDENTITY in step 7.
""")


if __name__ == '__main__':
    main()
