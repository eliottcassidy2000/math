"""k = 3 observer on every unique-Hamiltonian-cycle class, n = 3..7 (FLINT Smith forms), and the P_2 remark."""
import csv
from alib import gentourng, adj_from_ascii, path_table, ham_cycles, m3_matrix, snf_diag, snf_summary, fmt_snf
for n in range(3, 8):
    tf = tors = 0
    out = []
    for s in gentourng(n):
        A = adj_from_ascii(n, s)
        H = path_table(A)
        if ham_cycles(A, H) != 1:
            continue
        M, _, _ = m3_matrix(A, H)
        d = snf_diag(M)
        t = any(x > 1 for x in d)
        tors += t
        tf += (not t)
        if t:
            out.append((s, fmt_snf(snf_summary(d))))
    print('n=%d k=3 unique-HC: torsion-free %d, torsion %d %s' % (n, tf, tors, out), flush=True)
for n in range(3, 9):
    rows = list(csv.DictReader(open('out/m2_classes_n%d.tsv' % n), delimiter='\t'))
    def tor(sn):
        return any(int(tok.split('^')[0]) > 1 for tok in sn.split())
    u = [r for r in rows if r['HC'] == '1']
    print('n=%d P2 torsion among unique-HC classes: %d/%d' % (n, sum(tor(r['snf_P2']) for r in u), len(u)))
