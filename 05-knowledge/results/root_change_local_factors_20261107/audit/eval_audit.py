"""eval_audit.py - auditor's evaluation of the registered tests R2-R5 from the auditor's own census (census_2p24.txt)
and the auditor's own residue automaton (resdp_audit.c raw output), applying PREREGISTRATION.md exactly as written:
  R2: rho_p = mean(N = 0 mod p)/mean(N != 0 mod p) inside the REGISTERED band (table parsed from PREREGISTRATION.md);
      all seven must hold. Also shown against the auditor's own DP band (n = 18..23, x0.97 / x1.03).
  R3: rho_5 > rho_17 > rho_11 > max(rho_7, rho_13, rho_23, rho_31).
  R4: for each p of R2, max/min of the p-1 nonzero class means <= 1.10.
  R5: at least 25 of 27 classes r mod 27 have census mean(r)/overall mean inside the DP per-class band over
      n = 18..23 widened by 5% (multiplicatively, the convention the registered R2 table demonstrably uses);
      also reported: the unwidened band.
Also compares every class mean with the packet's lf_2p24.txt.
usage: python3 eval_audit.py <PREREGISTRATION.md> <census_2p24.txt> <lf_2p24.txt> <dp bands raw>
"""
import re
import sys

prereg, census, lf, dpraw = sys.argv[1:5]
txt = open(prereg).read()
row_p = re.search(r'\| p \| ([\d |]+)\|\s*\n\s*\| --- [^\n]*\n\s*\| band \| ([^\n]*) \|', txt)
ps = [int(x) for x in row_p.group(1).split('|') if x.strip()]
bands = [tuple(float(y) for y in x.strip().strip('[]').split(',')) for x in row_p.group(2).split('|')]
BAND = dict(zip(ps, bands))
print('registered R2 bands parsed from PREREGISTRATION.md:', BAND)


def means_of(fn, fmt):
    m = {}; overall = None
    for l in open(fn):
        g = re.match(r'mod (\d+) class means: (.*)$', l.strip())
        if g:
            m[int(g.group(1))] = [float(x) for x in g.group(2).split()]
        g = re.search(r'overall mean ([\d.]+)', l)
        if g:
            overall = float(g.group(1))
    return m, overall


mine, ov = means_of(census, 'audit')
theirs, ov2 = means_of(lf, 'packet')
dmax = max(abs(a - b) for q in mine for a, b in zip(mine[q], theirs[q]))
print('class means mod 27,5,7,11,13,17,23,31: auditor vs packet lf_2p24.txt, max |diff| = %.2e (packet prints 4 decimals); overall %.6f vs %.4f'
      % (dmax, ov, ov2))

raw = {}
for l in open(dpraw):
    g = re.match(r'q=(\d+) n=(\d+) model L: ([^|]+)\|', l)
    if g:
        raw[(int(g.group(1)), int(g.group(2)))] = [float(x) for x in g.group(3).split()]

rho = {}
ok2 = True
for p in ps:
    c = mine[p]; rho[p] = c[0] / (sum(c[1:]) / (p - 1))
    lo, hi = BAND[p]
    r = [raw[(p, n)][0] / (sum(raw[(p, n)][1:]) / (p - 1)) for n in range(18, 24)]
    alo, ahi = 0.97 * min(r), 1.03 * max(r)
    ins = lo <= rho[p] <= hi; ok2 &= ins
    print('R2 p=%2d rho=%.5f registered band [%.4f, %.4f] %-7s | auditor DP band [%.5f, %.5f] %s | margin to registered band %+.5f'
          % (p, rho[p], lo, hi, 'inside' if ins else 'OUTSIDE', alo, ahi, 'inside' if alo <= rho[p] <= ahi else 'OUTSIDE',
             min(rho[p] - lo, hi - rho[p])))
print('R2', 'CONFIRMED' if ok2 else 'FAILED')
r3 = rho[5] > rho[17] > rho[11] > max(rho[7], rho[13], rho[23], rho[31])
print('R3 %s: %.5f > %.5f > %.5f > %.5f' % ('CONFIRMED' if r3 else 'FAILED', rho[5], rho[17], rho[11], max(rho[7], rho[13], rho[23], rho[31])))
ok4 = True
for p in ps:
    c = mine[p][1:]; ok4 &= max(c) / min(c) <= 1.10
    print('R4 p=%2d nonzero max/min %.5f' % (p, max(c) / min(c)))
print('R4', 'CONFIRMED' if ok4 else 'FAILED')
lo27 = [min(raw[(27, n)][r] for n in range(18, 24)) for r in range(27)]
hi27 = [max(raw[(27, n)][r] for n in range(18, 24)) for r in range(27)]
ins5 = 0; ins0 = 0
for r in range(27):
    L = mine[27][r] / ov
    a = 0.95 * lo27[r] <= L <= 1.05 * hi27[r]; b = lo27[r] <= L <= hi27[r]
    ins5 += a; ins0 += b
    print('R5 class %2d census %.5f DP band [%.5f, %.5f] widened [%.5f, %.5f] %s%s' % (r, L, lo27[r], hi27[r], 0.95 * lo27[r], 1.05 * hi27[r],
          'inside' if a else 'OUTSIDE', '' if b else ' (outside the unwidened band)'))
print('R5 %d of 27 inside the 5%%-widened band -> %s; %d of 27 inside the unwidened DP band' % (ins5, 'CONFIRMED' if ins5 >= 25 else 'FAILED', ins0))
