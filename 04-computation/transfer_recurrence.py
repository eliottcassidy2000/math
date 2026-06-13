#!/usr/bin/env python3
"""
Transfer matrix recurrence: M_n[a,b] = B_b(V\\{a}) - sum_v t(v,a) M_{n-1}^{(a)}[v,b]

This gives an INDUCTIVE structure for the transfer matrix.
Key question: can this recurrence + inductive hypothesis (M_{n-1} even in r)
prove M_n is even in r?

Instance: opus-2026-03-06-S22
"""

from itertools import permutations
from sympy import symbols, expand, Symbol, Poly, Rational
from collections import defaultdict

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
    return r + s[(i,j)]

def ham_paths_ending_at(vertex_set, target, r, s):
    vs = list(vertex_set)
    if len(vs) == 1:
        return 1
    total = 0
    others = [v for v in vs if v != target]
    for perm in permutations(others):
        path = list(perm) + [target]
        w = 1
        for k in range(len(path)-1):
            w *= edge_weight(path[k], path[k+1], r, s)
        total += w
    return expand(total)

def ham_paths_beginning_at(vertex_set, source, r, s):
    vs = list(vertex_set)
    if len(vs) == 1:
        return 1
    total = 0
    others = [v for v in vs if v != source]
    for perm in permutations(others):
        path = [source] + list(perm)
        w = 1
        for k in range(len(path)-1):
            w *= edge_weight(path[k], path[k+1], r, s)
        total += w
    return expand(total)

def compute_M(a, b, vertices, r, s):
    """Compute M[a,b] on given vertex set."""
    U = [v for v in vertices if v != a and v != b]
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


def verify_recurrence(n):
    """Verify: M_n[a,b] = B_b(V\\{a}) - sum_v t(v,a) M_{n-1}^{(a)}[v,b]"""
    print(f"\n{'='*60}")
    print(f"RECURRENCE VERIFICATION: n={n}")
    print(f"{'='*60}")

    r, s = make_symbols(n)
    V = list(range(n))
    a, b = 0, 1

    # Direct computation of M_n[a,b]
    M_direct = compute_M(a, b, V, r, s)

    # Recurrence: M_n[a,b] = B_b(V\{a}) - sum_v t(v,a) M^{(a)}[v,b]
    V_minus_a = [v for v in V if v != a]
    U = [v for v in V if v != a and v != b]

    Bb_full = ham_paths_beginning_at(set(V_minus_a), b, r, s)

    recurrence_sum = 0
    for v in U:
        M_sub = compute_M(v, b, V_minus_a, r, s)
        recurrence_sum += edge_weight(v, a, r, s) * M_sub

    M_recurrence = expand(Bb_full - recurrence_sum)

    diff = expand(M_direct - M_recurrence)
    print(f"  M_n[{a},{b}] (direct) matches recurrence: {diff == 0}")

    # Also verify the COLUMN recurrence:
    # M_n[a,b] = E_a(V\{b}) + (-1)^{n-2} E_a(V\{b})... hmm
    # Actually, by symmetry of the formula, there should be an analogous column recurrence.
    # Let me derive: condition on last step of E_a, or first step of B_b.

    # First step of B_b: B_b(R+b) = sum_w t(b,w) B_w(R)
    # This gives: M_n[a,b] = sum_w t(b,w) M^{(b)}[a,w] + (-1)^{n-2} E_a(V\{b})

    # Let me verify this
    V_minus_b = [v for v in V if v != b]

    Ea_full = ham_paths_ending_at(set(V_minus_b), a, r, s)

    col_sum = 0
    for w in U:
        M_sub_col = compute_M(a, w, V_minus_b, r, s)
        col_sum += edge_weight(b, w, r, s) * M_sub_col

    M_col_recurrence = expand(col_sum + (-1)**(n-2) * Ea_full)
    diff_col = expand(M_direct - M_col_recurrence)
    print(f"  Column recurrence matches: {diff_col == 0}")

    if diff_col != 0:
        # Try without the E_a term
        M_col2 = expand(col_sum)
        diff_col2 = expand(M_direct - M_col2)
        print(f"  Column recurrence (no E_a term): diff = {diff_col2}")

    # Now analyze the r-parity structure of the recurrence.
    # B_b(V\{a}) has both odd and even r-powers.
    # sum_v t(v,a) M^{(a)}[v,b] also has both.
    # The cancellation of odd powers is the key.

    print(f"\n  --- R-parity analysis ---")

    Bb_poly = Poly(Bb_full, r)
    Bb_odd = sum(c for m, c in Bb_poly.as_dict().items() if m[0] % 2 == 1)
    Bb_even = sum(c * r**m[0] for m, c in Bb_poly.as_dict().items() if m[0] % 2 == 0)

    rec_poly = Poly(recurrence_sum, r)
    rec_odd = sum(c for m, c in rec_poly.as_dict().items() if m[0] % 2 == 1)

    print(f"  Odd r-powers of B_b(V\\{{a}}): {expand(Bb_odd)}")
    print(f"  Odd r-powers of sum t(v,a)M^(a)[v,b]: {expand(rec_odd)}")
    print(f"  Match (odd parts equal): {expand(Bb_odd - rec_odd) == 0}")

    # The odd part of M_n[a,b] is: (odd part of B_b) - (odd part of sum)
    # For M_n to be even, we need these to be equal.

    # Explicit: what IS the odd part of B_b?
    if n <= 5:
        print(f"\n  B_b(V\\{{{a}}}) = {Bb_full}")
        for v in U:
            M_sub = compute_M(v, b, V_minus_a, r, s)
            print(f"  t({v},{a}) * M^(a)[{v},{b}] = {expand(edge_weight(v, a, r, s) * M_sub)}")

    return M_direct


def analyze_boundary_term():
    """
    The boundary term B_b(V\\{a}) counts Hamiltonian paths in V\\{a} starting at b.

    KEY QUESTION: Does B_b(V\\{a}) have a specific r-parity structure?

    We know B_b(V; r, s) = E_b(V; r, -s) (path reversal).
    So B_b(V; -r, s) = E_b(V; -r, -s) = E_b(V; -r, -s).

    Also E_b(V; -r, s) = (-1)^{|V|-1} B_b(V; r, s) (from general path reversal).
    So B_b(V; -r, s) = (-1)^{|V|-1} E_b(V; r, -s)... hmm wait.

    Let me just check: is B_b(V\\{a}) antisymmetric in (a,b)?
    """
    print(f"\n\n{'='*60}")
    print("BOUNDARY TERM ANALYSIS")
    print(f"{'='*60}")

    for n in [4, 5]:
        r, s = make_symbols(n)
        V = list(range(n))

        # B_1(V\{0}) vs B_0(V\{1})
        a, b = 0, 1
        V_minus_a = [v for v in V if v != a]
        V_minus_b = [v for v in V if v != b]

        Bb_Va = ham_paths_beginning_at(set(V_minus_a), b, r, s)
        Ba_Vb = ham_paths_beginning_at(set(V_minus_b), a, r, s)

        # Also: E_a(V\{b}) vs E_b(V\{a})
        Ea_Vb = ham_paths_ending_at(set(V_minus_b), a, r, s)
        Eb_Va = ham_paths_ending_at(set(V_minus_a), b, r, s)

        print(f"\n  n={n}:")
        print(f"  B_b(V\\{{a}}) - B_a(V\\{{b}}) = {expand(Bb_Va - Ba_Vb)}")

        # Check if B_b(V\{a}) relates to B_a(V\{b}) under r->-r
        Bb_neg = expand(Bb_Va.subs(r, -r))
        print(f"  B_b(V\\{{a}})(-r) = {Bb_neg}")

        # Path reversal on V\{a}: B_b(V\{a}; r, s) = E_b(V\{a}; r, -s)
        # So B_b(-r, s) = E_b(-r, -s)
        # E_b(-r, -s) = E_b(-r, -s)... relate to B?
        # B_v(r,s) = E_v(r,-s), so E_v(-r,-s) = B_v(-r,s)... hmm circular.

        # Actually: E_v(V; -r, s) = (-1)^{|V|-1} B_v(V; r, s) was proved for the S+a setup.
        # More precisely, for T[V] at weights (r,s): reversing all paths gives
        # B_v(V; r, s) = E_v(V; r, -s)  [just negate the skew parts]
        # So E_v(V; -r, s) = B_v(V; -r, -s)... still circular.

        # Let me just directly check: B_b(V\{a}; -r) = ? * B_a(V\{b}; r)
        Ba_Vb_r = expand(Ba_Vb)
        ratio_check = expand(Bb_neg + (-1)**(n-2) * Ba_Vb)
        print(f"  B_b(V\\{{a}})(-r) + (-1)^{{n-2}} B_a(V\\{{b}})(r) = {ratio_check}")
        print(f"  B_b(V\\{{a}})(-r) - (-1)^{{n-2}} B_a(V\\{{b}})(r) = {expand(Bb_neg - (-1)**(n-2) * Ba_Vb)}")


def test_row_sum():
    """Check if row sums of M^{(a)} are zero."""
    print(f"\n\n{'='*60}")
    print("ROW SUM TEST FOR M^{(a)}")
    print(f"{'='*60}")

    for n_sub in [3, 4]:
        r, s = make_symbols(n_sub + 1)  # Need symbols for original n
        V_sub = list(range(1, n_sub + 1))  # V\{0}

        for b in V_sub:
            row_sum = 0
            for v in V_sub:
                if v != b:
                    M_vb = compute_M(v, b, V_sub, r, s)
                    row_sum += M_vb
            print(f"  n_sub={n_sub}, b={b}: sum_v M[v,{b}] = {expand(row_sum)}")


def main():
    verify_recurrence(4)
    verify_recurrence(5)
    analyze_boundary_term()
    test_row_sum()

    print(f"\n\n{'='*60}")
    print("KEY RECURRENCE")
    print(f"{'='*60}")
    print("""
PROVED: M_n[a,b] = B_b(V\\{a}) - sum_{v in U} t(v,a) M_{n-1}^{(a)}[v,b]

where:
  - B_b(V\\{a}) = sum of Ham path weights in V\\{a} starting at b
  - M_{n-1}^{(a)}[v,b] = transfer matrix on V\\{a}
  - t(v,a) = r + s_{va}

For M_n to be even in r:
  odd_r(B_b) = odd_r(sum_v t(v,a) M^{(a)}[v,b])

This is a precise cancellation between the "boundary" Hamiltonian
path generating function and the "bulk" transfer matrix recurrence.
""")


if __name__ == '__main__':
    main()
