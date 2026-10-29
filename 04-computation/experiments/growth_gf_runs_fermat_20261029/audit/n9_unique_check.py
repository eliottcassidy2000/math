"""Independent check of the producer's order-9 unique-Hamiltonian-cycle split (HYP-9445 (2), post-launch commit
03c05b4762): reads only the 144 class strings recorded in ../out/stream_m2_n9_part*.out, recomputes HC by the
prepend DP, M_2 by the subset DP and Smith forms by FLINT."""
import ast, glob
from collections import Counter
from runlib import adj_from_bits, hp_table, ham_cycles, run_matrix, snf_diag, snf_summary, fmt
uniq = []
for f in sorted(glob.glob('../out/stream_m2_n9_part*.out')):
    for line in open(f):
        if line.startswith('unique-HC list (bits, torsion): '):
            uniq += ast.literal_eval(line.split(': ', 1)[1])
bits_list = [b for b, _ in uniq]
assert len(bits_list) == 144 and len(set(bits_list)) == 144
c = Counter()
agree = 0
for b, tor_prod in uniq:
    A = adj_from_bits(9, b)
    hp = hp_table(A)
    assert ham_cycles(A, hp) == 1
    _, M = run_matrix(A, 2, hp)
    d = snf_diag(M)
    tor = any(x > 1 for x in d)
    sing = 0 in d
    agree += tor == tor_prod
    c[('torsion' if tor else 'tf', 'singular' if sing else 'nonsingular')] += 1
print('144 distinct recorded classes, each with exactly one Hamiltonian cycle (F(12) = 144, Garey)')
print('split:', sorted(c.items()), '; torsion flags agree with the producer on', agree, 'of 144')
