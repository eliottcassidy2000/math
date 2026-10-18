#!/usr/bin/env python3
"""audit F1: score the pre-registered items Q1-Q5, T1-T4, L1-L7 directly from the session's receipts (../*.out),
using the thresholds in the committed pre-registration JSON files, and print PASS/FAIL per item."""
import re, os, json
D = os.path.dirname(os.path.abspath(__file__)); R = os.path.join(D, '..')
RES = os.path.join(D, '..', '..', '..', '..', '05-knowledge', 'results')

def fap_octaves(fn):
    """octave -> dict from fap.c 'octave 2^k:' lines (and the two continuation lines when present)."""
    out = {}; L = open(os.path.join(R, fn)).read().splitlines()
    for i, l in enumerate(L):
        m = re.match(r'octave 2\^\s*(\d+):', l)
        if not m: continue
        k = int(m.group(1)); d = {}
        g = re.search(r'E_V max ([\d.]+) mean ([\d.]+)', l); d['EVmax'], d['EVmean'] = float(g.group(1)), float(g.group(2))
        g = re.search(r'E_clk max ([\d.]+) mean ([\d.]+)', l); d['G0mean'] = float(g.group(2))
        g = re.search(r'cross-port \(E_V-E_port\)/N mean ([+-][\d.]+) range \[([+-][\d.]+), ([+-][\d.]+)\]', l)
        if g: d['Xmean'], d['Xlo'], d['Xhi'] = map(float, g.groups())
        else:
            g = re.search(r'cross-port \(E_V-E_port\)/N mean ([+-][\d.]+) max\|\.\| ([\d.]+)', l); d['Xmean'] = float(g.group(1)); d['Xabs'] = float(g.group(2))
        g = re.search(r'G1/G0 mean ([+-][\d.]+) range \[([+-][\d.]+), ([+-][\d.]+)\]', l)
        if g: d['r1'], d['r1lo'], d['r1hi'] = map(float, g.groups())
        if 'Xlo' in d: d['Xabs'] = max(abs(d['Xlo']), abs(d['Xhi']))
        out[k] = d
    return out

def lean_rows(fn):
    out = {}
    for l in open(os.path.join(R, fn)):
        m = re.match(r'\s*(\d+) \| (.*)$', l)
        if m: out[int(m.group(1))] = [float(x) for x in re.findall(r'[+-]?\d+\.\d+', m.group(2))]
    return out

res = []
def item(name, ok, obs): res.append((name, 'PASS' if ok else 'FAIL', obs))

# Q1-Q5: octaves 23-25 (receipts committed after the Q pre-registration)
v1k = fap_octaves('fap_V1000_k25_full.out'); v100 = fap_octaves('fap_V100_k25_full.out')
mx = [v1k[k]['Xabs'] for k in (23, 24, 25)]
item('Q1 V_1000 max|X|/N <= 0.10 each octave', all(x <= 0.10 for x in mx), mx)
item('Q2 mean of the three maxima < 0.08', sum(mx) / 3 < 0.08, round(sum(mx) / 3, 4))
ev = [v1k[k]['EVmax'] for k in (23, 24, 25)]
item('Q3 V_1000 max E_V/N <= 2.5', all(x <= 2.5 for x in ev), ev)
m100 = [v100[k]['Xabs'] for k in (23, 24, 25)]
item('Q4 V_100 max|X|/N <= 0.35 each octave', all(x <= 0.35 for x in m100), m100)
xm = [v1k[k]['Xmean'] for k in (23, 24, 25)]
item('Q5 V_1000 octave means of X/N in [-0.06, 0.06]', all(-0.06 <= x <= 0.06 for x in xm), xm)
# T1-T4: octave 26
t = fap_octaves('fap_V1000_k26_full2.out')[26]
item('T1 mean G1/G0 in [-0.290, -0.265]', -0.290 <= t['r1'] <= -0.265, t['r1'])
item('T2 every N G1/G0 in [-0.300, -0.250]', -0.300 <= t['r1lo'] and t['r1hi'] <= -0.250, (t['r1lo'], t['r1hi']))
item('T3 mean G0/N in [2.85, 3.20]', 2.85 <= t['G0mean'] <= 3.20, t['G0mean'])
item('T4 max E_V/N <= 2.5 and mean X/N in [-0.15, 0.10]', t['EVmax'] <= 2.5 and -0.15 <= t['Xmean'] <= 0.10, (t['EVmax'], t['Xmean']))
# L1-L7: octave 27 (fap_lean)
r = lean_rows('fap_lean_V1000_k27.out')[27]
# columns: 0 EVmean 1 EVmax 2 G0 3 r1 4 r1lo 5 r1hi 6 z0m 7 z0x 8 gam 9 gamx 10 hlm 11 hlx 12 E+ 13 dip 14 e0lo 15 e0hi 16 sblo 17 sbhi
item('L1 mean G1/G0 in [-0.295,-0.265], every N in [-0.300,-0.250]', -0.295 <= r[3] <= -0.265 and -0.300 <= r[4] and r[5] <= -0.250, (r[3], r[4], r[5]))
item('L2 (G0-E0)/N < 0 every N, mean in [-1.95,-1.70]', r[7] < 0 and -1.95 <= r[6] <= -1.70, (r[6], r[7]))
item('L3 Gam1/N < 0 every N, mean in [-0.25,-0.05]', r[9] < 0 and -0.25 <= r[8] <= -0.05, (r[8], r[9]))
item('L4 HL/N <= 0.07 every N', r[11] <= 0.07, r[11])
item('L5 E_V/N <= 2.0 every N, mean in [0.9,1.5]', r[1] <= 2.0 and 0.9 <= r[0] <= 1.5, (r[0], r[1]))
item('L6 mean G0/N in [2.92,3.10]', 2.92 <= r[2] <= 3.10, r[2])
item('L7 mean E_V/E+ in [0.18,0.35]', 0.18 <= r[13] <= 0.35, r[13])
item('L-check E0/N within 1e-4 of 34/7, SIB_V/N within 1e-3 of 5/7', abs(r[14] - 34 / 7) < 1e-4 and abs(r[15] - 34 / 7) < 1e-4 and abs(r[16] - 5 / 7) < 1e-3 and abs(r[17] - 5 / 7) < 1e-3, (r[14], r[15], r[16], r[17]))
for name, v, obs in res: print('%-62s %s   observed %s' % (name, v, obs))
# the thresholds above are copied from the committed JSON files; print them for the record
for fn in ('first_arrival_port_orthogonality_prereg_20261018.json', 'first_arrival_one_tick_prereg_20261018.json', 'first_arrival_pair_ledger_prereg_20261018.json'):
    j = json.load(open(os.path.join(RES, fn)))
    print('\n%s:' % fn); [print('  %s: %s' % kv) for kv in j['predictions'].items()]
