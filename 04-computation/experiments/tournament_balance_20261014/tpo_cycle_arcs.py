# Transitive-plus-one tournaments T = TT_m + z (z -> u_i iff w_i = 1), strong iff w_1 = 1 and w_m = 0.
# The unique Hamiltonian cycle is z u_1 u_2 ... u_m z.  Exact probabilities from the path <-> valid-set bijection of
# THM-4655 1(b): a path is (A increasing) z (B increasing); its class is an initial segment A = {1..k}, k in K (weight 1),
# or a crossing (beta, alpha) with w_beta = 1, w_alpha = 0, beta < alpha (weight 2^(alpha-beta-1)): indices < beta and
# alpha lie in A, beta and indices > alpha in B, indices strictly between are independent fair coins.
# Checks, for every strong word with m <= M:
#  (i)   the THM-4655 1(e) dichotomy: w_2 = 1 => (z, u_1) balanced, P(z before u_1) in [1/2, 2/3);
#                                     w_2 = 0 => (u_1, u_2) balanced, P(u_2 before u_1) in [1/3, 1/2]
#        (and the mirror statement at the other end), so the Hamiltonian cycle has >= 2 balanced arcs when m >= 3;
#  (ii)  every u_i has a balanced partner among z, u_(i-1), u_(i+1)  (FINITE-EXACT observation);
#  (iii) the minimum number of balanced arcs on the Hamiltonian cycle.
from fractions import Fraction as Fr
from itertools import product
import sys
def classes(w):
    m = len(w); out = []
    K = sorted(set([0, m] + [k for k in range(1, m) if w[k - 1] == 0 and w[k] == 1]))
    for k in K: out.append((1, [None] + [Fr(1) if i <= k else Fr(0) for i in range(1, m + 1)]))
    for b in range(1, m + 1):
        if w[b - 1] != 1: continue
        for a in range(b + 1, m + 1):
            if w[a - 1] != 0: continue
            out.append((2 ** (a - b - 1), [None] + [Fr(1) if (i < b or i == a) else (Fr(1, 2) if b < i < a else Fr(0))
                                                    for i in range(1, m + 1)]))
    return out
def probs(w):
    m = len(w); cl = classes(w); H = sum(c[0] for c in cl)
    pA = [None] + [sum(wt * q[i] for wt, q in cl) / H for i in range(1, m + 1)]                 # P(u_i before z)
    rev = {i: sum(wt * q[i + 1] * (1 - q[i]) for wt, q in cl) / H for i in range(1, m)}          # P(u_(i+1) before u_i)
    return H, pA, rev
bal = lambda p: Fr(1, 3) <= p <= Fr(2, 3)
M = int(sys.argv[1]) if len(sys.argv) > 1 else 14
for m in range(3, M + 1):
    dich_fail = 0; nb_fail = 0; minbal = 99
    for mid in product((0, 1), repeat=m - 2):
        w = (1,) + mid + (0,)
        H, pA, rev = probs(w)
        z_u1 = 1 - pA[1]                      # P(z before u_1)
        u2_u1 = rev[1]                        # P(u_2 before u_1)
        if w[1] == 1: ok = Fr(1, 2) <= z_u1 < Fr(2, 3)
        else:         ok = Fr(1, 3) <= u2_u1 <= Fr(1, 2)
        if not ok: dich_fail += 1
        arcs = [z_u1] + [1 - rev[i] for i in range(1, m)] + [pA[m]]   # p(tail before head) around z u_1 ... u_m z
        minbal = min(minbal, sum(bal(p) for p in arcs))
        for i in range(1, m + 1):
            c = [bal(pA[i])] + ([bal(rev[i - 1])] if i > 1 else []) + ([bal(rev[i])] if i < m else [])
            if not any(c): nb_fail += 1
    print(f"m={m} (n={m + 1}): {2 ** (m - 2)} words; dichotomy (i) failures {dich_fail}; "
          f"u_i without balanced partner in {{z, u_(i-1), u_(i+1)}}: {nb_fail}; min #balanced Ham-cycle arcs {minbal}")
