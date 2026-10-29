"""Exhaustive two-run torsion test at a given order, streamed in chunks (memory stays small).

For every strong class (gentourng), counts (number of Hamiltonian cycles, torsion in coker M_2) and lists every class
with at least two Hamiltonian cycles and a torsion-free M_2 (a counterexample to HYP-9445 (1)), plus the unique-cycle
classes.  usage: python3 stream_m2.py n out_path
"""
import sys
import subprocess
from fast_runs import m2_fast, ham_cycles, tuples2
from runs import tour_from_bits, gp_matrix, gp_batch, parse_vec


def main():
    n = int(sys.argv[1])
    outp = sys.argv[2]
    tup = tuples2(n)
    proc = subprocess.Popen(['gentourng', '-q', str(n)], stdout=subprocess.PIPE, text=True)
    stats = {}
    counter = []
    uniq = []
    chunk, meta = [], []
    total = 0

    def flush():
        if not chunk:
            return
        res = gp_batch(chunk)
        for (bits, hc), line in zip(meta, res):
            d = [abs(x) for x in parse_vec(line)]
            torsion = any(x > 1 for x in d)
            key = ('HC=1' if hc == 1 else 'HC>=2', 'torsion' if torsion else 'torsion-free')
            stats[key] = stats.get(key, 0) + 1
            if hc >= 2 and not torsion:
                counter.append((bits, hc, line[:200]))
            if hc == 1:
                uniq.append((bits, torsion))
        chunk.clear()
        meta.clear()

    for line in proc.stdout:
        bits = line.strip()
        if not bits:
            continue
        total += 1
        A = tour_from_bits(n, bits)
        T, hp = m2_fast(A)
        hc = ham_cycles(A, hp)
        if hc == 0:
            stats[('HC=0', 'skipped')] = stats.get(('HC=0', 'skipped'), 0) + 1
            continue
        M = [[int(T[s[0], s[1], t[0], t[1]]) for t in tup] for s in tup]
        chunk.append('matsnf(%s)' % gp_matrix(M))
        meta.append((bits, hc))
        if len(chunk) >= 200:
            flush()
            if total % 10000 < 200:
                with open(outp, 'w') as f:
                    f.write('progress: %d classes read; %s; counterexamples so far %d\n' % (total, sorted(stats.items()), len(counter)))
    flush()
    with open(outp, 'w') as f:
        f.write('n=%d classes=%d\n' % (n, total))
        f.write('SUMMARY: %s\n' % sorted(stats.items()))
        f.write('unique-HC classes: %d, torsion-free: %d\n' % (len(uniq), sum(1 for _, t in uniq if not t)))
        f.write('counterexamples to HYP-9445 (1) at k = 2: %s\n' % (counter if counter else 'none'))


if __name__ == '__main__':
    main()
