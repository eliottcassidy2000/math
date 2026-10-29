"""Two-run observer and its sectors for the rotational tournaments R_n (i -> i+1..i+(n-1)/2), n = 3..13 odd,
and the Paley tournaments P_7, P_11.  Smith forms (PARI); Hamiltonian cycle counts by DP."""
from fast_runs import m2_fast, ham_cycles, tuples2
from runs import circulant, gp_matrix, gp_batch, parse_vec, snf_summary, fmt_snf
cases = [('R%d' % n, n, set(range(1, (n - 1) // 2 + 1))) for n in (3, 5, 7, 9, 11, 13)]
cases += [('P7', 7, {1, 2, 4}), ('P11', 11, {1, 3, 4, 5, 9})]
cmds, meta = [], []
for name, n, S in cases:
    A = circulant(n, S)
    T, hp = m2_fast(A)
    hc = ham_cycles(A, hp)
    tup = tuples2(n)
    m = len(tup) // 2
    M = [[int(T[s[0], s[1], t[0], t[1]]) for t in tup] for s in tup]
    P = [[M[2*a][2*b] + M[2*a][2*b+1] for b in range(m)] for a in range(m)]
    Sg = [[M[2*a][2*b] - M[2*a][2*b+1] for b in range(m)] for a in range(m)]
    cmds += ['matsnf(%s)' % gp_matrix(X) for X in (M, P, Sg)]
    meta.append((name, n, hc))
out = gp_batch(cmds)
for i, (name, n, hc) in enumerate(meta):
    f = [fmt_snf(snf_summary(parse_vec(out[3*i + j]))) for j in range(3)]
    print('%-4s HC=%-8d M2: %s | P2: %s | S2: %s' % (name, hc, f[0][:70], f[1][:60], f[2][:60]))
