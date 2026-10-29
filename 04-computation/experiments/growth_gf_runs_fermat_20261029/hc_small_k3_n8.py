"""HYP-9445 (1) at k = 3, order 8, in the delicate region: every class with exactly 1, 2 or 3 Hamiltonian cycles
(55 + 60 + 124 classes, read from out/fast_m2_n8.out). Three-run observer (336 x 336), Smith form by PARI."""
import re
from runs import tour_from_bits, run_matrix, gp_matrix, gp_batch, parse_vec
sel = []
for l in open('out/fast_m2_n8.out'):
    m = re.match(r'^([01]+) HC=(\d+) ', l)
    if m and int(m.group(2)) in (1, 2, 3):
        sel.append((m.group(1), int(m.group(2))))
cmds = []
for bits, hc in sel:
    A = tour_from_bits(8, bits)
    _, M = run_matrix(A, 3)
    cmds.append('matsnf(%s)' % gp_matrix(M))
out = []
for i in range(0, len(cmds), 20):
    out += gp_batch(cmds[i:i + 20])
stats = {}
for (bits, hc), line in zip(sel, out):
    d = [abs(x) for x in parse_vec(line)]
    torsion = any(x > 1 for x in d)
    stats[(hc, torsion)] = stats.get((hc, torsion), 0) + 1
print('order 8, k = 3: (HC, has torsion) -> classes:', sorted(stats.items()))
