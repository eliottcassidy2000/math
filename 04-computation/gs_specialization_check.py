#!/usr/bin/env python3
"""
Verify the GS-OCF specialization: p_1 -> 1, p_{2k+1} -> 2 for k>=1, p_{2k} -> 0.

From Corollary 19/20 of Grinberg-Stanley (arXiv 2412.10572):
  U_T = sum_{sigma in S_V(T,T_bar)} (-1)^phi(sigma) p_{type(sigma)}

And ham(T_bar) = ps_1(U_T)(1) = sum sigma... but with specific specialization.

Actually, Corollary 20 says:
  ham(D_bar) = sum_{sigma in S(D), all cycles odd} 2^{psi(sigma)}

where psi = number of nontrivial cycles.

So the specialization is: keep only terms with ALL ODD parts, then replace
p_{lambda} -> 2^{number of parts >= 3} * 1^{number of parts = 1}
= 2^{l(lambda) - m_1(lambda)} where m_1 = number of 1's.

Wait, psi counts NONTRIVIAL cycles, and each contributes 2. So
the weight is 2^{(number of nontrivial cycles)} = 2^{l(lambda) - m_1(lambda)}.

But p_{type(sigma)} with type lambda = (lambda_1, ..., lambda_k)
specializes under p_j -> 1 for all j to: 1.
Under p_1 -> 1, p_{2k+1} -> 2, p_{2k} -> 0:
  p_lambda -> prod_{i} [p_{lambda_i}] = prod_{lambda_i=1} 1 * prod_{lambda_i odd >= 3} 2 * prod_{lambda_i even} 0
  = 0 if any lambda_i is even
  = 2^{#{lambda_i >= 3}} = 2^{l(lambda) - m_1(lambda)} if all lambda_i odd.

So H(T) = sum_{sigma, all parts odd} (-1)^phi(sigma) * 2^{l(lambda)-m_1(lambda)}

But phi = sum_{gamma X-cycle} (len(gamma) - 1). For tournaments,
X_bar = X^op, and (-1)^phi...

Let me just compute and check.

opus-2026-03-07-S36
"""
from itertools import permutations
from collections import defaultdict

def redei_berge_coeffs(n, edges):
    edge_set = set(edges)
    opp_edges = set((j, i) for (i, j) in edge_set)  # T^op = reverse all arcs

    coeffs = defaultdict(int)

    for sigma in permutations(range(n)):
        visited = [False] * n
        cycles = []
        for start in range(n):
            if visited[start]:
                continue
            cycle = []
            curr = start
            while not visited[curr]:
                visited[curr] = True
                cycle.append(curr)
                curr = sigma[curr]
            cycles.append(tuple(cycle))

        valid = True
        phi = 0
        for cyc in cycles:
            if len(cyc) == 1:
                continue
            is_X = all((cyc[i], cyc[(i+1) % len(cyc)]) in edge_set for i in range(len(cyc)))
            is_Xbar = all((cyc[i], cyc[(i+1) % len(cyc)]) in opp_edges for i in range(len(cyc)))
            if not is_X and not is_Xbar:
                valid = False
                break
            if is_X:
                phi += len(cyc) - 1

        if valid:
            cycle_type = tuple(sorted([len(c) for c in cycles], reverse=True))
            coeffs[cycle_type] += (-1)**phi

    return dict(coeffs)

def specialization(coeffs, spec):
    """Specialize: spec maps part size k to value."""
    total = 0
    for lam, coeff in coeffs.items():
        val = coeff
        for part in lam:
            val *= spec.get(part, 0)
        total += val
    return total

def main():
    # Cyclic tournament C_5
    n = 5
    cyc_edges = [(i, (i+1) % n) for i in range(n)] + [(i, (i+2) % n) for i in range(n)]
    coeffs = redei_berge_coeffs(n, cyc_edges)

    print(f"=== C_5: H = 15 ===")
    print(f"Coefficients: {coeffs}")

    # Try different specializations
    # spec1: p_k -> 1 for all k
    s1 = {k: 1 for k in range(1, n+1)}
    print(f"  p_k -> 1: {specialization(coeffs, s1)}")

    # spec2: p_1 -> 1, p_odd -> 2, p_even -> 0
    s2 = {1: 1, 2: 0, 3: 2, 4: 0, 5: 2}
    print(f"  p_1->1, p_odd->2, p_even->0: {specialization(coeffs, s2)}")

    # spec3: p_k -> k for all k (principal spec at m=1 is different)
    # Actually ps_1(p_k)(m) = m^k? No, ps_1(p_k)(m) = sum_{i=1}^m 1^k = m.
    # So ps_1(p_lambda)(m) = m^{l(lambda)}.
    # At m=1: ps_1(p_lambda)(1) = 1 for all lambda. So sum of coeffs = ps_1(U)(1).
    # But ham(D_bar) uses a DIFFERENT formula.

    # The issue: ps_1(U_T)(m) counts something else.
    # For tournament T, ham(T_bar) where T_bar = T^op:
    # ham(T^op) = ham(T) by path reversal.
    # But U_T evaluated at ps_1 gives u_T(m), the Redei-Berge polynomial.
    # u_T(1) = # of listings with no T-descents? No...
    # u_T(m) = # of (listing, coloring) pairs where coloring is weakly increasing
    # with strict increase at T-descents.
    # u_T(1) = # of listings with NO T-descent = # of acyclic orderings compatible with T.

    # For tournament, T-descent of listing (v1,...,vn) at position i means
    # (v_i, v_{i+1}) is an edge of T, i.e., v_i beats v_{i+1}.
    # A listing with no T-descent means v_i does NOT beat v_{i+1} for all i,
    # i.e., v_{i+1} beats v_i for all i. This is a Hamiltonian path in T^op!

    # So u_T(1) = ham(T^op) = ham(T).

    # But we computed sum of coeffs = 3 for C_5, while ham = 15. Why?

    # Because ps_1(p_lambda)(1) is NOT 1 for all lambda!
    # ps_1(p_k)(1) = 1, so ps_1(p_lambda)(1) = 1^{l(lambda)} = 1. Yes it is 1.
    # But the sum is 3, and H = 15. Contradiction!

    # Wait... the sign might be wrong. Let me recheck.
    # Theorem 2.2 says U_X = sum (-1)^phi p_type
    # But Equation (4) says omega(U_X) = U_X_bar.
    # And ps_1(omega(f))(m) = ps_1(f)(-m) * (-1)^n? No...
    # Actually omega(p_lambda) = (-1)^{|lambda|-l(lambda)} p_lambda.
    # So U_X_bar = omega(U_X) = sum (-1)^phi (-1)^{n-l(type)} p_type.

    # For tournament, X_bar = X^op.
    # U_{X^op} = omega(U_X).
    # ham(X) = ps_1(U_{X^op})(1) = ps_1(omega(U_X))(1).
    # But ps_1(omega(f))(m) = ???

    # omega(F_I) = F_{I^c} where I^c = [n-1] \ I.
    # ps_1(F_I)(m) = C(m + n - 1 - |I|, n) ... let me look this up.
    # Actually ps_1(F_I)(m) = C(m, n) for I = empty, C(m-1, n-1) for |I|=n-1.
    # In general ps_1(F_I)(m) = C(m - |I| + n - 1 - |I|, n) ... no.

    # Let me just compute it differently.

    # Method 2: directly count permutations with all odd nontrivial cycles
    # that are T-cycles, weighted by 2^{# nontrivial cycles}.
    edge_set = set(cyc_edges)

    H_ocf = 0
    for sigma in permutations(range(n)):
        visited = [False] * n
        cycles = []
        for start in range(n):
            if visited[start]:
                continue
            cycle = []
            curr = start
            while not visited[curr]:
                visited[curr] = True
                cycle.append(curr)
                curr = sigma[curr]
            cycles.append(tuple(cycle))

        # Check: all nontrivial cycles must be directed cycles of T
        # AND have odd length
        ok = True
        nontrivial = 0
        for cyc in cycles:
            if len(cyc) == 1:
                continue
            if len(cyc) % 2 == 0:
                ok = False
                break
            is_T = all((cyc[i], cyc[(i+1) % len(cyc)]) in edge_set for i in range(len(cyc)))
            if not is_T:
                ok = False
                break
            nontrivial += 1

        if ok:
            H_ocf += 2**nontrivial

    print(f"\n  Direct OCF computation (all odd T-cycles, 2^psi): {H_ocf}")
    print(f"  H(C_5) = 15")
    print(f"  Match: {H_ocf == 15}")

    # Paley T_7
    print(f"\n=== Paley T_7: H = 189 ===")
    n = 7
    QR = {1, 2, 4}
    paley_edges = [(i, j) for i in range(n) for j in range(n) if i != j and (j - i) % n in QR]
    edge_set = set(paley_edges)

    H_ocf = 0
    for sigma in permutations(range(n)):
        visited = [False] * n
        cycles = []
        for start in range(n):
            if visited[start]:
                continue
            cycle = []
            curr = start
            while not visited[curr]:
                visited[curr] = True
                cycle.append(curr)
                curr = sigma[curr]
            cycles.append(tuple(cycle))

        ok = True
        nontrivial = 0
        for cyc in cycles:
            if len(cyc) == 1:
                continue
            if len(cyc) % 2 == 0:
                ok = False
                break
            is_T = all((cyc[i], cyc[(i+1) % len(cyc)]) in edge_set for i in range(len(cyc)))
            if not is_T:
                ok = False
                break
            nontrivial += 1

        if ok:
            H_ocf += 2**nontrivial

    print(f"  Direct OCF: {H_ocf}")
    print(f"  H(T_7) = 189")
    print(f"  Match: {H_ocf == 189}")

if __name__ == "__main__":
    main()
