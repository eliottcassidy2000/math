from runs import *
A = circulant(5, {1, 2})
tuples, M = run_matrix(A, 3)
nb = len(tuples) // 6
N, P, S, L0, L1 = unglued_k3(M, nb)
cmds = ['matsnf(%s)' % gp_matrix(X) for X in (M, N, P, S, L0, L1)] + ['matrank(%s)' % gp_matrix(M)]
out = gp_batch(cmds)
names = ['M3', 'unglued N', 'plus P', 'minus S', 'std col0', 'std col1']
for nm, line in zip(names, out[:6]):
    print(nm, fmt_snf(snf_summary(parse_vec(line))))
print('rank', out[6], 'total mass', sum(map(sum, M)))
