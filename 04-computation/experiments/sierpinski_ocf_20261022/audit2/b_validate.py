# Audit2 (B): validate the brute-force cluster machinery (b_partitions.py) against log Phi_T(z) computed directly from
# the linear-forest numbers N_L (a_forest_identity.forest_sums), exactly, on small tournaments; then report
#  - the sum over all excess levels equals [z^L] log Phi_T (Moebius inversion + exponential formula),
#  - tuples with an odd run cancel (odd=True and odd=False give the same totals),
#  - for regular T the degree-1 pruning changes nothing,
#  - non-regular T have nonzero tree terms (excess -1), regular T have none.
import sys, random
from fractions import Fraction
import numpy as np
sys.path.insert(0, '.')
from a_forest_identity import forest_sums, tour_from_bits
from b_partitions import cluster_sums
from common import rot_A, paley_A, skew

def log_series(psi, K):
    lam = {}
    for L in range(1, K + 1):
        s = psi.get(L, Fraction(0))
        for k in range(1, L):
            s -= Fraction(k, L) * lam[k] * psi.get(L - k, Fraction(0))
        lam[L] = s
    return lam

def direct_log(A, K):
    n = A.shape[0]
    S = skew(A).tolist()
    N, _, _ = forest_sums(S, n)
    psi = {L: Fraction(N[L], n ** L) for L in range(1, n)}
    return log_series(psi, K)

def run(A, label, K, odd_upto, regular):
    n = A.shape[0]
    S = skew(A)
    lam = direct_log(A, K)
    print(f"{label} (n={n}, regular={regular}):")
    cs_odd = cluster_sums(S, odd_upto, regular=False, emax=None, odd=True, verbose=False)
    cs_even = cluster_sums(S, K, regular=False, emax=None, odd=False, verbose=False)
    cs_reg = cluster_sums(S, K, regular=True, emax=None, odd=False, verbose=False) if regular else None
    allok = True
    for L in range(1, K + 1):
        tot_e = sum(cs_even[L].values(), Fraction(0))
        ok = (tot_e == lam[L])
        line = f"  L={L}: [z^L] log Phi = {float(lam[L]):+.12e}; cluster total (even runs) {'==' if ok else '!='}"
        if L <= odd_upto:
            tot_o = sum(cs_odd[L].values(), Fraction(0))
            ok2 = (tot_o == lam[L]) and all(cs_odd[L].get(e, 0) == cs_even[L].get(e, 0) for e in set(cs_odd[L]) | set(cs_even[L]))
            line += f"; all runs {'==' if ok2 else '!='} (per excess)"
            ok = ok and ok2
        if regular:
            ok3 = all(cs_reg[L].get(e, 0) == cs_even[L].get(e, 0) for e in set(cs_reg[L]) | set(cs_even[L]))
            line += f"; degree-1 pruning {'harmless' if ok3 else 'CHANGES RESULT'}"
            ok = ok and ok3
        tree = cs_even[L].get(-1, Fraction(0))
        line += f"; tree (e=-1) term = {float(tree):+.6e}"
        print(line)
        allok = allok and ok
    return allok

if __name__ == '__main__':
    random.seed(7)
    ok = True
    for n, reps in ((5, 3), (6, 3), (7, 2)):
        m = n * (n - 1) // 2
        for r in range(reps):
            A = np.array(tour_from_bits(n, [random.randint(0, 1) for _ in range(m)]), dtype=np.int64)
            ok &= run(A, f"random tournament #{r}", min(n - 1, 6), 5 if n <= 6 else 4, False)
    ok &= run(rot_A(7), "R_7", 6, 4, True)
    ok &= run(paley_A(7), "P_7", 6, 4, True)
    print("VALIDATION PASSED" if ok else "VALIDATION FAILED")
