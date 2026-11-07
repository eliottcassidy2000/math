"""Evaluate S15 registered tests R2-R5 (PREREGISTRATION.md) from the census class means (lf_2p24.txt) and the
residue-automaton prediction (dp_bands.txt for primes, dp27_classes.txt for classes mod 27)."""
import re, sys
lf = open(sys.argv[1]).read().splitlines()
bands_txt = open(sys.argv[2]).read().splitlines()
dp27 = open(sys.argv[3]).read().splitlines()
overall = float(re.search(r'overall mean ([\d.]+)', lf[0]).group(1))
means = {}
for l in lf[1:]:
    m = re.match(r'mod (\d+) class means: (.*)$', l)
    if m:
        means[int(m.group(1))] = [float(x) for x in m.group(2).split()]
BAND = {5: (1.4375, 1.5674), 7: (1.1093, 1.2409), 11: (1.2604, 1.3569), 13: (1.1001, 1.1707),
        17: (1.3379, 1.4857), 23: (1.0650, 1.1486), 31: (1.0514, 1.1227)}
rho = {}
ok2 = True
for p, (lo, hi) in BAND.items():
    c = means[p]; nz = sum(c[1:]) / (p - 1); rho[p] = c[0] / nz
    inside = lo <= rho[p] <= hi; ok2 &= inside
    print('R2 p=%d rho=%.4f band [%.4f, %.4f] %s' % (p, rho[p], lo, hi, 'inside' if inside else 'OUTSIDE'))
print('R2', 'CONFIRMED' if ok2 else 'FAILED')
r3 = rho[5] > rho[17] > rho[11] > max(rho[7], rho[13], rho[23], rho[31])
print('R3 order rho_5 > rho_17 > rho_11 > max(rho_7, rho_13, rho_23, rho_31): %s (%.4f, %.4f, %.4f; others max %.4f)'
      % ('CONFIRMED' if r3 else 'FAILED', rho[5], rho[17], rho[11], max(rho[7], rho[13], rho[23], rho[31])))
ok4 = True
for p in BAND:
    c = means[p][1:]; r = max(c) / min(c); ok4 &= r <= 1.10
    print('R4 p=%d nonzero-class max/min %.4f' % (p, r))
print('R4', 'CONFIRMED' if ok4 else 'FAILED')
# R5: per-class band from the DP, n = 18..23
lo27 = [9e9] * 27; hi27 = [-9e9] * 27
for l in dp27:
    c = [float(x) for x in l.split('by class:')[1].split()]
    for r in range(27):
        lo27[r] = min(lo27[r], c[r]); hi27[r] = max(hi27[r], c[r])
inside = 0
for r in range(27):
    L = means[27][r] / overall
    ok = lo27[r] * 0.95 <= L <= hi27[r] * 1.05
    inside += ok
    print('R5 class %2d census %.4f DP band [%.4f, %.4f] %s' % (r, L, lo27[r], hi27[r], '' if ok else 'OUTSIDE'))
print('R5 %d of 27 classes inside -> %s (needs >= 25)' % (inside, 'CONFIRMED' if inside >= 25 else 'FAILED'))
