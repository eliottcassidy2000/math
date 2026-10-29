"""Audit of HYP-9445 (1)-(2) data (two-run observer M_2 and its signed sector S_2), every class n = 3..8.

Independent of the producer: M_2 from alib's prepend-vertex path table (alib.py), Hamiltonian cycles by the same
table but checked against brute force in test_alib.py, strong connectivity by BFS, Smith forms and determinants by
FLINT (python-flint fmpz_mat.snf / det), not PARI.
    S_2[{a1<a2},{b1<b2}] = M_2[(a1,a2),(b1,b2)] - M_2[(a1,a2),(b2,b1)],  P_2 = the same with +.
"Torsion" = some nonzero invariant factor > 1.
Writes out/m2_classes_n<n>.tsv (bits, strong, HC, Smith M2, Smith S2, Smith P2, det M2) and a summary.
"""
import sys
import os
from collections import Counter
import numpy as np
from alib import (gentourng, adj_from_ascii, path_table, ham_cycles, is_strong, m2_matrix, snf_diag, snf_summary,
                  fmt_snf, det_exact)

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, 'out')
PROD = os.path.join(HERE, '..', 'out')


def torsion(d):
    return any(x > 1 for x in d)


def only_odd(d):
    big = [x for x in d if x > 1]
    return bool(big) and all(x % 2 for x in big)


def sectors(M, n):
    m = M.shape[0] // 2
    P = np.array([[M[2 * a, 2 * b] + M[2 * a, 2 * b + 1] for b in range(m)] for a in range(m)], dtype=np.int64)
    S = np.array([[M[2 * a, 2 * b] - M[2 * a, 2 * b + 1] for b in range(m)] for a in range(m)], dtype=np.int64)
    return P, S


def fib(k):
    a, b = 0, 1
    for _ in range(k):
        a, b = b, a + b
    return a


summary = []
for n in range(3, int(sys.argv[1]) + 1 if len(sys.argv) > 1 else 9):
    cls = gentourng(n)
    rows = []
    for s in cls:
        A = adj_from_ascii(n, s)
        H = path_table(A)
        hc = ham_cycles(A, H)
        st = is_strong(A)
        assert st == (hc > 0)
        M, pr = m2_matrix(A, H)
        # pair ordering check: (b0,b1),(b1,b0) per base; sectors use rows 2a (sorted) and columns 2b, 2b+1
        assert all(pr[2 * k][0] < pr[2 * k][1] and pr[2 * k + 1] == pr[2 * k][::-1] for k in range(len(pr) // 2))
        # run symmetry of M_2 (swap both tuples)
        sw = [k ^ 1 for k in range(len(pr))]
        assert (M[np.ix_(sw, sw)] == M).all()
        P, S = sectors(M, n)
        dM, dS, dP = snf_diag(M), snf_diag(S), snf_diag(P)
        det = det_exact(M)
        rows.append(dict(bits=s, strong=st, hc=hc, dM=dM, dS=dS, dP=dP, det=det))
    with open(os.path.join(OUT, 'm2_classes_n%d.tsv' % n), 'w') as f:
        f.write('bits\tstrong\tHC\tsnf_M2\tsnf_S2\tsnf_P2\tdet_M2\n')
        for r in rows:
            f.write('%s\t%d\t%d\t%s\t%s\t%s\t%d\n' % (r['bits'], r['strong'], r['hc'], fmt_snf(snf_summary(r['dM'])),
                                                     fmt_snf(snf_summary(r['dS'])), fmt_snf(snf_summary(r['dP'])),
                                                     r['det']))
    strong2 = [r for r in rows if r['strong'] and r['hc'] >= 2]
    uniq = [r for r in rows if r['hc'] == 1]
    lines = []
    lines.append('n=%d: classes=%d strong=%d unique-HC=%d (F(2n-6)=%d) HC>=2=%d'
                 % (n, len(rows), sum(r['strong'] for r in rows), len(uniq), fib(2 * n - 6) if n >= 3 else 0,
                    len(strong2)))
    lines.append('   HC>=2: S2 torsion %d/%d ; M2 torsion %d/%d ; M2 torsion with only odd factors >1: %d ; '
                 'P2 torsion %d/%d'
                 % (sum(torsion(r['dS']) for r in strong2), len(strong2), sum(torsion(r['dM']) for r in strong2),
                    len(strong2), sum(only_odd(r['dM']) for r in strong2), sum(torsion(r['dP']) for r in strong2),
                    len(strong2)))
    uS_tf = [r for r in uniq if not torsion(r['dS'])]
    uS_t = [r for r in uniq if torsion(r['dS'])]
    lines.append('   unique-HC: S2 torsion-free %d/%d ; S2 torsion classes: %s'
                 % (len(uS_tf), len(uniq), [(r['bits'], fmt_snf(snf_summary(r['dS']))) for r in uS_t]))
    lines.append('   unique-HC: S2 singular %d (torsion-free and singular %d)'
                 % (sum(0 in r['dS'] for r in uniq), sum((0 in r['dS']) and not torsion(r['dS']) for r in uniq)))
    uM_tf = [r for r in uniq if not torsion(r['dM'])]
    nonsing = [r for r in uniq if r['det'] != 0]
    lines.append('   unique-HC: M2 torsion-free %d ; M2 nonsingular %d ; dets %s'
                 % (len(uM_tf), len(nonsing), dict(Counter(r['det'] for r in uniq))))
    lines.append('   unique-HC: [M2 torsion <=> M2 singular] violations: %d ; singular torsion-free: %s'
                 % (sum(torsion(r['dM']) != (r['det'] == 0) for r in uniq),
                    [r['bits'] for r in uniq if r['det'] == 0 and not torsion(r['dM'])]))
    lines.append('   unique-HC with M2 torsion: %s'
                 % [(r['bits'], fmt_snf(snf_summary(r['dM']))) for r in uniq if torsion(r['dM'])])
    # strong with HC >= 2 whose S2 or M2 is torsion-free (counterexamples)
    cex = [(r['bits'], r['hc'], fmt_snf(snf_summary(r['dS'])), fmt_snf(snf_summary(r['dM'])))
           for r in strong2 if not torsion(r['dS']) or not torsion(r['dM'])]
    lines.append('   COUNTEREXAMPLES (strong, HC>=2, S2 or M2 torsion-free): %s' % cex)
    # unimodular M2 only with a unique HC (strong classes)
    lines.append('   strong classes with unimodular M2 (|det|=1): %d, all unique-HC: %s'
                 % (sum(1 for r in rows if r['strong'] and abs(r['det']) == 1),
                    all(r['hc'] == 1 for r in rows if r['strong'] and abs(r['det']) == 1)))
    # comparison with the producer's PARI outputs, class by class
    prodf = os.path.join(PROD, 'fast_m2_n%d.out' % n)
    if os.path.exists(prodf):
        prod = {}
        for ln in open(prodf):
            parts = ln.split()
            if len(parts) >= 3 and parts[1].startswith('HC='):
                hcp = int(parts[1][3:])
                m2p = []
                for tok in parts[2:]:
                    if tok.startswith('unglued='):
                        break
                    m2p.append(tok.replace('M2=', ''))
                prod[parts[0]] = (hcp, ' '.join(m2p))
        agree = sum(1 for r in rows if r['bits'] in prod and prod[r['bits']] == (r['hc'], fmt_snf(snf_summary(r['dM']))))
        lines.append('   producer fast_m2_n%d.out: %d entries; HC and Smith(M2) agree for %d/%d classes'
                     % (n, len(prod), agree, len(rows)))
    for ln in lines:
        print(ln, flush=True)
    summary += lines
with open(os.path.join(OUT, 'm2_torsion_summary.out'), 'w') as f:
    f.write('\n'.join(summary) + '\n')
