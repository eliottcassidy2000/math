"""Signed two-run sector S_2 versus Hamiltonian cycles at order 8 (streamed): is coker S_2 torsion-free iff HC = 1?"""
import sys, subprocess, collections
from fast_runs import m2_fast, ham_cycles, tuples2
from runs import tour_from_bits, gp_matrix, gp_batch, parse_vec
n = int(sys.argv[1]) if len(sys.argv) > 1 else 8
tup = tuples2(n)
m = len(tup) // 2
proc = subprocess.Popen(['gentourng', '-q', str(n)], stdout=subprocess.PIPE, text=True)
c = collections.Counter()
viol = []
chunk, meta = [], []
def flush():
    res = gp_batch(chunk) if chunk else []
    for (bits, hc), line in zip(meta, res):
        d = [abs(x) for x in parse_vec(line)]
        t = any(x > 1 for x in d)
        sing = 0 in d
        c[('HC=1' if hc == 1 else 'HC>=2', 'S torsion' if t else 'S tf', 'S singular' if sing else 'S nonsingular')] += 1
        if (hc == 1) == t:
            viol.append((bits, hc, line[:120]))
    chunk.clear(); meta.clear()
for line in proc.stdout:
    bits = line.strip()
    if not bits:
        continue
    A = tour_from_bits(n, bits)
    T, hp = m2_fast(A)
    hc = ham_cycles(A, hp)
    if hc == 0:
        continue
    S = [[int(T[tup[2*a][0], tup[2*a][1], tup[2*b][0], tup[2*b][1]] - T[tup[2*a][0], tup[2*a][1], tup[2*b+1][0], tup[2*b+1][1]]) for b in range(m)] for a in range(m)]
    chunk.append('matsnf(%s)' % gp_matrix(S)); meta.append((bits, hc))
    if len(chunk) >= 300:
        flush()
flush()
print('n=%d strong classes:' % n, sorted(c.items()))
print('violations of [S_2 torsion-free <=> unique Hamiltonian cycle]:', viol[:20] if viol else 'none', len(viol))
