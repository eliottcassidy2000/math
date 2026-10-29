"""Two-run sectors versus Hamiltonian cycles: Smith forms of the plus (unordered covers) and minus (signed) sectors."""
import sys, collections
from fast_runs import m2_fast, ham_cycles, tuples2
from runs import tournaments, tour_from_bits, gp_matrix, gp_batch, parse_vec
for n in (5, 6, 7):
    tup = tuples2(n)
    m = len(tup) // 2
    recs, cmds = [], []
    for bits in tournaments(n):
        A = tour_from_bits(n, bits)
        T, hp = m2_fast(A)
        hc = ham_cycles(A, hp)
        if hc == 0:
            continue
        M = [[int(T[s[0], s[1], t[0], t[1]]) for t in tup] for s in tup]
        P = [[M[2*a][2*b] + M[2*a][2*b+1] for b in range(m)] for a in range(m)]
        S = [[M[2*a][2*b] - M[2*a][2*b+1] for b in range(m)] for a in range(m)]
        cmds += ['matsnf(%s)' % gp_matrix(P), 'matsnf(%s)' % gp_matrix(S)]
        recs.append(hc)
    out = gp_batch(cmds)
    c = collections.Counter()
    for i, hc in enumerate(recs):
        dP = [abs(x) for x in parse_vec(out[2*i])]
        dS = [abs(x) for x in parse_vec(out[2*i+1])]
        tP = any(x > 1 for x in dP); tS = any(x > 1 for x in dS)
        c[('HC=1' if hc == 1 else 'HC>=2', 'P torsion' if tP else 'P tf', 'S torsion' if tS else 'S tf')] += 1
    print(n, sorted(c.items()))
