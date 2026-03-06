#!/usr/bin/env python3
"""
CLEAN PROOF: M[b,a](r) = M[a,b](-r) — Transfer matrix r-swap identity

Instance: opus-2026-03-06-S11

THEOREM: For the c-tournament transfer matrix M[a,b] defined by
  M[a,b](r) = sum_{S⊂U} (-1)^|S| E_a(S∪{a}; r) B_b((U\S)∪{b}; r)
we have M[b,a](r) = M[a,b](-r) for all a,b, where U = V \ {a,b}.

PROOF:
Step 1: Path reversal identity.
  B_v(V; r) = (-1)^{|V|-1} E_v(V; -r)

  Proof: A path v_1→...→v_k→v ending at v has weight
  product(r + s_{v_i, v_{i+1}}). The reversed path v→v_k→...→v_1
  starting at v has weight product(r + s_{v_{i+1}, v_i}) = product(r - s_{v_i, v_{i+1}}).
  Since product(r - s_j) = (-1)^k product(-r + s_j) = (-1)^k [original at r→-r],
  and k = |V|-1 edges, we get B_v(V;r) = (-1)^{|V|-1} E_v(V;-r). □

Step 2: Rewrite M using path reversal.
  B_b(R∪{b}; r) = (-1)^{|R|} E_b(R∪{b}; -r)  [since |R∪{b}|-1 = |R| = |U|-|S|]

  M[a,b](r) = sum_S (-1)^|S| E_a(S∪{a}; r) (-1)^{|U|-|S|} E_b(R∪{b}; -r)
             = (-1)^{|U|} sum_S (-1)^{2|S|} E_a(S∪{a}; r) E_b(R∪{b}; -r)
             = (-1)^{n-2} sum_S E_a(S∪{a}; r) E_b((U\S)∪{b}; -r)    ... (★)

Step 3: Compute M[a,b](-r) from (★).
  M[a,b](-r) = (-1)^{n-2} sum_S E_a(S∪{a}; -r) E_b((U\S)∪{b}; r)    ... (★★)

  [using E_b(V; -(-r)) = E_b(V; r) when applying the substitution to B_b]

Step 4: Compute M[b,a](r) and compare.
  M[b,a](r) = (-1)^{n-2} sum_S E_b(S∪{b}; r) E_a((U\S)∪{a}; -r)

  Substitute T = U\S (so S = U\T):
  = (-1)^{n-2} sum_T E_b((U\T)∪{b}; r) E_a(T∪{a}; -r)

  Comparing with (★★): M[a,b](-r) = (-1)^{n-2} sum_T E_a(T∪{a}; -r) E_b((U\T)∪{b}; r)

  These are IDENTICAL (multiplication commutes). □

COROLLARY: The following are equivalent:
  (a) M[a,b] = M[b,a]
  (b) M[a,b](r) = M[a,b](-r) (only even powers of r)
  (c) sum_S (-1)^|S| E_a(S∪{a}) B_b((U\S)∪{b}) has no odd r-powers

NOTE: This corrects the earlier claim M[b,a] = (-1)^{n-2} M[a,b](-r)
which had an incorrect sign factor due to omitting the (-1)^|R| from
the B→E conversion in Step 2.

FURTHER RESULTS:
- Leading coefficient [r^{n-2}] of M[a,b] = (n-2)! if n even, 0 if n odd.
  This follows from sum_{k=0}^m (-1)^k C(m,k) k! (m-k)! = m! if m even, 0 if m odd.
- At n=5: [r^2] = 2 * sum_{i∈{a,b}, u∈U} s_{iu} (only endpoint-to-internal edges).
- At n=6: [r^4] = 4! = 24 = (n-2)!.
"""

from itertools import permutations
from sympy import symbols, expand, Poly, Symbol

def setup(n):
    r = Symbol('r')
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

def transfer_M(t_fn, n, a, b):
    U = [v for v in range(n) if v != a and v != b]
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
            for i in range(len(p)-1):
                prod *= t_fn(p[i], p[i+1])
            ea += prod
        if len(S_set) == 1: ea = 1
        bb = 0
        for p in permutations(sorted(R_set)):
            if p[0] != b: continue
            prod = 1
            for i in range(len(p)-1):
                prod *= t_fn(p[i], p[i+1])
            bb += prod
        if len(R_set) == 1: bb = 1
        result += sign * ea * bb
    return expand(result)

def verify(n, pairs=None):
    r, sv, s, t = setup(n)
    if pairs is None:
        pairs = [(0, 1)]

    print(f"\nn={n}:")
    for a, b in pairs:
        Mab = transfer_M(t, n, a, b)
        Mba = transfer_M(t, n, b, a)
        Mab_negr = expand(Mab.subs(r, -r))

        sym = expand(Mab - Mba) == 0
        swap = expand(Mba - Mab_negr) == 0
        even = expand(Mab - Mab_negr) == 0

        p = Poly(Mab, r)
        deg = p.degree()

        print(f"  M[{a},{b}]: deg={deg}")
        print(f"    M[{a},{b}] = M[{b},{a}]? {sym}")
        print(f"    M[{b},{a}] = M[{a},{b}](-r)? {swap}")
        print(f"    M[{a},{b}](r) = M[{a},{b}](-r)? {even}")

        for k in range(deg + 1):
            c = expand(p.nth(k))
            if c == 0:
                parity = "EVEN" if k % 2 == 0 else "ODD"
                print(f"    [r^{k}] = 0 ({parity}{'→ vanishes!' if k%2==1 else ''})")

if __name__ == '__main__':
    for n in [3, 4, 5, 6]:
        verify(n, [(0, 1)] if n <= 5 else [(0, 1)])
