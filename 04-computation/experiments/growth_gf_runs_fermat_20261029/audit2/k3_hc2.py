"""HYP-9445 (1), k = 3: every strong class with >= 2 Hamiltonian cycles, n = 5..7, has torsion in coker M_3.
Certificate: nonsingular -> |det| > 1 (FLINT det); singular -> FLINT Smith form, then rank drop mod the smallest
prime of the smallest invariant factor > 1 (independent elimination)."""
import flint
from alib import gentourng, adj_from_ascii, path_table, ham_cycles, m3_matrix, snf_diag, rank_mod_p
for n in (5, 6, 7):
    tot = ok = sing = 0
    bad = []
    for s in gentourng(n):
        A = adj_from_ascii(n, s)
        H = path_table(A)
        if ham_cycles(A, H) < 2:
            continue
        tot += 1
        M, _, _ = m3_matrix(A, H)
        fm = flint.fmpz_mat([[int(x) for x in row] for row in M])
        r = int(fm.rank())
        if r == M.shape[0]:
            good = abs(int(fm.det())) > 1
        else:
            sing += 1
            d = sorted(x for x in snf_diag(M) if x > 1)
            good = False
            if d:
                q = int(min(f for f, e in flint.fmpz(d[0]).factor()))
                good = rank_mod_p(M, q) < r
        ok += good
        if not good:
            bad.append(s)
    print('n=%d k=3: strong HC>=2 classes %d, torsion certified %d (singular %d), counterexamples %s'
          % (n, tot, ok, sing, bad), flush=True)
