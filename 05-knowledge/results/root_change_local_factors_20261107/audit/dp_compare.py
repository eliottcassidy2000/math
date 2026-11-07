"""dp_compare.py - compare the auditor's residue-automaton output (resdp_audit.c, raw files) with the packet's
dp_bands.txt, dp27_classes.txt, dp27_long.txt, dp_primes_n1000.txt and q3_long.txt, at the packet's printed precision,
and print the auditor's own R2 / R5 bands (n = 18..23, widened by 3% and 5%).
usage: python3 dp_compare.py <packet dir> <bands.raw> <q27long.raw> <n1000.raw> <q3long.raw>
"""
import re
import sys

pk, f_bands, f_long, f_1000, f_q3 = sys.argv[1:6]


def parse_raw(fn):
    out = {}
    for l in open(fn):
        m = re.match(r'q=(\d+) n=(\d+) model L: ([^|]+)\| level1 L: ([^|]+)\| level2 L: (.+)$', l.strip())
        if m:
            q, n = int(m.group(1)), int(m.group(2))
            out[(q, n)] = {'model': [float(x) for x in m.group(3).split()],
                           'l1': [float(x) for x in m.group(4).split()], 'l2': [float(x) for x in m.group(5).split()]}
    return out


def vgroup(q, p):
    kq = 0; t = q
    while t % p == 0:
        t //= p; kq += 1
    g = {}
    for r in range(q):
        v = kq if r == 0 else next(k for k in range(kq + 1) if r % p ** (k + 1))
        g.setdefault(min(v, kq), []).append(r)
    return [g[k] for k in sorted(g)]


def vagg(L, q, p):
    return [sum(L[r] for r in grp) / len(grp) for grp in vgroup(q, p)]


raw = parse_raw(f_bands); raw.update(parse_raw(f_long)); raw.update(parse_raw(f_1000)); raw.update(parse_raw(f_q3))
worst = 0.0; nchk = 0


def cmp(a, b, tol):
    global worst, nchk
    d = abs(a - b); nchk += 1
    if d > worst:
        worst = d
    return d <= tol


bad = []
# dp_bands.txt
for l in open(pk + '/dp_bands.txt'):
    m = re.match(r'q=27 n=(\d+) v3=0\.\.3: (.*)$', l.strip())
    if m:
        n = int(m.group(1)); want = [float(x) for x in m.group(2).split()]
        got = vagg(raw[(27, n)]['model'], 27, 3)
        if not all(cmp(a, b, 5.01e-5) for a, b in zip(got, want)):
            bad.append(('dp_bands q27', n, got, want))
        continue
    m = re.match(r'q=(\d+) n=(\d+) L\(0\)=([\d.]+) L\(nz\)avg=([\d.]+) ratio=([\d.]+)', l.strip())
    if m:
        q, n = int(m.group(1)), int(m.group(2)); L = raw[(q, n)]['model']
        nz = sum(L[1:]) / (q - 1)
        got = [L[0], nz, L[0] / nz]; want = [float(m.group(k)) for k in (3, 4, 5)]
        if not all(cmp(a, b, 5.01e-5) for a, b in zip(got, want)):
            bad.append(('dp_bands', q, n, got, want))
# dp27_classes.txt
for l in open(pk + '/dp27_classes.txt'):
    m = re.match(r'q=27 n=(\d+) .*by class: (.*)$', l.strip())
    if m:
        n = int(m.group(1)); want = [float(x) for x in m.group(2).split()]
        got = raw[(27, n)]['model']
        if not all(cmp(a, b, 5.01e-5) for a, b in zip(got, want)):
            bad.append(('dp27_classes', n))
# dp27_long.txt (v3 aggregates)
for l in open(pk + '/dp27_long.txt'):
    m = re.match(r'q=27 n=(\d+)\s+L by v_3\(N\)=0\.\.3: ([\d. ]+)', l.strip())
    if m:
        n = int(m.group(1)); want = [float(x) for x in m.group(2).split()[:4]]
        got = vagg(raw[(27, n)]['model'], 27, 3) if n > 0 else None
        if n == 0:
            continue
        if not all(cmp(a, b, 5.01e-5) for a, b in zip(got, want)):
            bad.append(('dp27_long', n, got, want))
# dp_primes_n1000.txt
for l in open(pk + '/dp_primes_n1000.txt'):
    m = re.match(r'q=(\d+) n=1000\s+L by v_\d+\(N\)=0\.\.1: ([\d.]+) ([\d.]+)\s+\| by class: (.*)$', l.strip())
    if m:
        q = int(m.group(1)); L = raw[(q, 1000)]['model']
        want_cls = [float(x) for x in m.group(4).split() if re.match(r'^\d+\.\d+$', x)]
        ok = cmp(sum(L[1:]) / (q - 1), float(m.group(2)), 5.01e-5) and cmp(L[0], float(m.group(3)), 5.01e-5)
        ok &= all(cmp(a, b, 5.01e-5) for a, b in zip(L, want_cls))
        if not ok:
            bad.append(('n1000', q))
# q3_long.txt
for l in open(pk + '/q3_long.txt'):
    m = re.match(r'n=(\d+) ([\d.]+) ([\d.]+) ([\d.]+)', l.strip())
    if m:
        n = int(m.group(1)); want = [float(m.group(k)) for k in (2, 3, 4)]
        if not all(cmp(a, b, 5.01e-5) for a, b in zip(raw[(3, n)]['model'], want)):
            bad.append(('q3_long', n, raw[(3, n)]['model'], want))
print('compared %d printed numbers with the auditor DP; largest |difference| %.2e; mismatches beyond print rounding: %d'
      % (nchk, worst, len(bad)))
for b in bad:
    print('  MISMATCH', b)

# R2 bands (ratio rho = L(0)/mean nonzero L over n = 18..23, widened 3%) and R5 per-class bands (widened 5%)
print('auditor DP bands (n = 18..23) and registered R2 bands:')
for p in (5, 7, 11, 13, 17, 23, 31):
    rh = []
    for n in range(18, 24):
        L = raw[(p, n)]['model']; rh.append(L[0] / (sum(L[1:]) / (p - 1)))
    print('  p=%2d DP rho band [%.4f, %.4f] -> widened [%.4f, %.4f]' % (p, min(rh), max(rh), 0.97 * min(rh), 1.03 * max(rh)))
lo = [min(raw[(27, n)]['model'][r] for n in range(18, 24)) for r in range(27)]
hi = [max(raw[(27, n)]['model'][r] for n in range(18, 24)) for r in range(27)]
print('auditor R5 per-class DP band mod 27 (n = 18..23): ' + ' '.join('%d:[%.4f,%.4f]' % (r, lo[r], hi[r]) for r in range(27)))
# large-n extrapolation for q = 3
a, b, c = raw[(3, 4000)]['model'][0], raw[(3, 8000)]['model'][0], raw[(3, 16000)]['model'][0]
print('q=3 L_n(0) at n=4000/8000/16000: %.10f %.10f %.10f; successive differences ratio %.4f; Richardson 2L(16000)-L(8000) = %.10f (34/27 = %.10f)'
      % (a, b, c, (a - b) / (b - c), 2 * c - b, 34 / 27))
print('q=3 parity check: |L_4001(0) - L_4000(0)| = %.2e, |L_16001(0) - L_16000(0)| = %.2e'
      % (abs(raw[(3, 4001)]['model'][0] - a), abs(raw[(3, 16001)]['model'][0] - c)))

# explanation of the dp_bands.txt column differences: recompute L(nz)avg and ratio from 4-decimal per-class values
ok = 0; tot = 0
for l in open(pk + '/dp_bands.txt'):
    m = re.match(r'q=(\d+) n=(\d+) L\(0\)=([\d.]+) L\(nz\)avg=([\d.]+) ratio=([\d.]+)', l.strip())
    if m:
        q, n = int(m.group(1)), int(m.group(2)); L = [round(x, 4) for x in raw[(q, n)]['model']]
        nz = sum(L[1:]) / (q - 1); r = L[0] / nz
        tot += 1; ok += abs(L[0] - float(m.group(3))) < 1e-9 and abs(nz - float(m.group(4))) <= 5.01e-5 and abs(r - float(m.group(5))) <= 5.01e-5
print('dp_bands.txt prime rows reproduced exactly when L(nz)avg and ratio are computed from the 4-decimal per-class values: %d of %d' % (ok, tot))
