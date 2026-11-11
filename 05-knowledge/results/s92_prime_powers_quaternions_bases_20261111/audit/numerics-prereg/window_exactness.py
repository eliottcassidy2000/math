# Auditor numerics-prereg: the window outputs of the session do not print the exactness flag ("every triple of type d
# lies in a class with d | c - c_sigma"), although THM-5134 says "exact, all runs".  This calls the session's own
# run_census / table (from the copy in session_copies/) on the registered window and on one in-sample window and
# prints the flags, the diagonal counts, and the largest class-level remainder of (B) on the window (not printed
# by the session either).
import sys
sys.path.insert(0, "session_copies")
import odd_modulus_census as OC  # noqa: E402

moduli = [9, 27, 81, 25, 125, 49, 45, 75]
for KHI, KLO in ((28, 27), (25, 24)):
    A, B = 1 << KLO, 1 << KHI
    q, divs, cnt, diag = OC.run_census(A, B, KHI, moduli)
    tab = OC.table(moduli, divs, cnt, diag)
    print("window [2^%d, 2^%d): triples %.0f / %.0f; diagonal triples %.0f / %.0f" % (
        KLO, KHI, tab[(0, 9)]["total"], tab[(1, 9)]["total"], tab[(0, 9)]["ndiag"], tab[(1, 9)]["ndiag"]))
    print("   exactness flag for every (sheet, modulus): %s" % all(tab[(sh, m)]["exact"] for sh in (0, 1) for m in moduli))
    print("   class-level remainder of (B) on the window, m max_c |#{N = c} - sum (d/m) #type d| / #T, sheet 0 / sheet 1:")
    for m in moduli:
        print("      m = %3d: %.5f / %.5f" % (m, tab[(0, m)]["rem"], tab[(1, m)]["rem"]))
