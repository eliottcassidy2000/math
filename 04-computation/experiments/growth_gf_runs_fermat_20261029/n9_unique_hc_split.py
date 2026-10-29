"""Order 9: split the 144 unique-Hamiltonian-cycle classes into unimodular / singular torsion-free / torsion,
re-deriving M_2 Smith forms from the class list recorded by stream_m2.py (out/stream_m2_n9_part*.out)."""
import ast, glob
from fast_runs import m2_fast, ham_cycles, tuples2
from runs import tour_from_bits, gp_matrix, gp_batch, parse_vec, snf_summary, fmt_snf
uniq = []
for f in sorted(glob.glob('out/stream_m2_n9_part*.out')):
    for line in open(f):
        if line.startswith('unique-HC list (bits, torsion): '):
            uniq += ast.literal_eval(line.split(': ', 1)[1])
assert len(uniq) == 144, len(uniq)
tup = tuples2(9)
cmds = []
for bits, _ in uniq:
    A = tour_from_bits(9, bits)
    T, hp = m2_fast(A)
    assert ham_cycles(A, hp) == 1
    M = [[int(T[s[0], s[1], t[0], t[1]]) for t in tup] for s in tup]
    cmds.append('matsnf(%s)' % gp_matrix(M))
out = []
for i in range(0, len(cmds), 24):
    out += gp_batch(cmds[i:i + 24])
split = {'unimodular': 0, 'singular torsion-free': 0, 'torsion': 0}
singular_tf = []
for (bits, torsion), line in zip(uniq, out):
    d = [abs(x) for x in parse_vec(line)]
    tor = any(x > 1 for x in d)
    sing = 0 in d
    assert tor == torsion
    if tor:
        split['torsion'] += 1
        split['torsion and singular' if sing else 'torsion and NONSINGULAR'] = split.get('torsion and singular' if sing else 'torsion and NONSINGULAR', 0) + 1
    elif sing:
        split['singular torsion-free'] += 1
        singular_tf.append((bits, fmt_snf(snf_summary(d))))
    else:
        split['unimodular'] += 1
print('order 9 unique-HC classes (144):', split)
print('singular torsion-free:', singular_tf)
