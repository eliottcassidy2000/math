"""HYP-9445 (2): what separates the torsion-carrying unique-Hamiltonian-cycle classes (k = 2)?
For each unique-HC class of order 7 and 8: max number of Hamiltonian cycles over proper strong subtournaments
of each order, scores, and torsion status of M_2 (from out/fast_m2_n*.out)."""
import re, itertools, collections
from runs import tour_from_bits
from torsion_vs_hc import ham_cycles_sub
for n in (7, 8):
    rows = []
    for l in open('out/fast_m2_n%d.out' % n):
        m = re.match(r'^([01]+) HC=(\d+) M2=(.*?) unglued', l)
        if m and int(m.group(2)) == 1:
            snf = m.group(3).split()
            tf = all(t.startswith('1^') or t.startswith('0^') for t in snf)
            rows.append((m.group(1), tf, m.group(3)))
    table = collections.Counter()
    for bits, tf, snf in rows:
        A = tour_from_bits(n, bits)
        prof = []
        for k in range(4, n):
            best = 0
            for S in itertools.combinations(range(n), k):
                best = max(best, ham_cycles_sub(A, S))
            prof.append(best)
        sc = ''.join(map(str, sorted(sum(r) for r in A)))
        table[(tuple(prof), tf)] += 1
    print('order %d unique-HC classes: (max #HC over subtournaments of orders 4..%d, torsion-free) -> count' % (n, n - 1))
    for k, v in sorted(table.items()):
        print('   ', k, v)
