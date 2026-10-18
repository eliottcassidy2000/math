#!/usr/bin/env python3
"""audit F1: compare the independent every-N scan (aud_scan_*_K21.out, octaves 2^8..2^21) with the session's receipts.

Receipts compared (all in ../):
  ledger_table.out            (three sections: V_1000 k26, V_100 k25, H_1000 k25)
  fap_lean_V1000_k27.out      (V_1000, incl. E0/N and SIB/N ranges)
  fap_V1000_k26_full2.out, fap_V100_k25_full.out, fap_H1000_k25.out (cross-port term X = E_V - E_port: mean and range)
Every printed number of an octave row must agree digit for digit.
"""
import re, sys, os
D = os.path.dirname(os.path.abspath(__file__))
R = os.path.join(D, '..')

def mine(fn):
    rows = {}
    for l in open(os.path.join(D, fn)):
        m = re.match(r'\s*(\d+) \| (.*)$', l)
        if not m: continue
        f = re.findall(r'[+-]?\d+\.\d+', m.group(2))
        rows[int(m.group(1))] = f
    return rows

def ledger_sections():
    secs, cur = {}, None
    for l in open(os.path.join(R, 'ledger_table.out')):
        if l.startswith('=='): cur = l.split()[1]; secs[cur] = {}; continue
        m = re.match(r'\s*(\d+) \| (.*)$', l)
        if m and cur: secs[cur][int(m.group(1))] = re.findall(r'[+-]?\d+\.\d+', m.group(2))
    return secs

def lean():
    rows = {}
    for l in open(os.path.join(R, 'fap_lean_V1000_k27.out')):
        m = re.match(r'\s*(\d+) \| (.*)$', l)
        if m: rows[int(m.group(1))] = re.findall(r'[+-]?\d+\.\d+', m.group(2))
    return rows

def cross(fn):
    rows = {}
    for l in open(os.path.join(R, fn)):
        m = re.match(r'octave 2\^\s*(\d+):.*cross-port \(E_V-E_port\)/N mean ([+-][\d.]+) (?:range \[([+-][\d.]+), ([+-][\d.]+)\]|max\|\.\| ([\d.]+))', l)
        if m: rows[int(m.group(1))] = m.groups()[1:]
    return rows

bad = 0; ncmp = 0
secs = ledger_sections()
# my columns: 0 EVmean 1 EVmax 2 G0 3 r1 4 r1lo 5 r1hi 6 z0m 7 z0x 8 gam 9 gamx 10 hlm 11 hlx 12 E+ 13 dip 14 e0lo 15 e0hi 16 svlo 17 svhi 18 Xm 19 Xlo 20 Xhi 21 Xmax 22 gamout
pairs = [('aud_scan_V_closed_Y1000_K21.out', 'fap_V1000_k26_ledger.out', 'fap_V1000_k26_full2.out'),
         ('aud_scan_V_closed_Y100_K21.out', 'fap_V100_k25_ledger.out', 'fap_V100_k25_full.out'),
         ('aud_scan_V_height_Y1000_K21.out', 'fap_H1000_k25_ledger.out', 'fap_H1000_k25.out')]
for mf, sec, cf in pairs:
    me = mine(mf); th = secs[sec]; cr = cross(cf)
    for k in sorted(me):
        if k not in th: continue
        a = me[k][:14]; b = th[k]                    # ledger_table prints 14 numbers per row
        b = [x if x.startswith(('+', '-')) or i in (0, 1, 2, 12, 13) else x for i, x in enumerate(b)]
        for i, (x, y) in enumerate(zip(a, b)):
            ncmp += 1
            if float(x) != float(y): bad += 1; print('MISMATCH', mf, 'octave', k, 'col', i, 'mine', x, 'receipt', y)
        if k in cr:
            c = cr[k]
            ncmp += 1
            if float(me[k][18]) != float(c[0]): bad += 1; print('MISMATCH X mean', mf, k, me[k][18], c[0])
            if c[1] is not None:
                ncmp += 2
                if float(me[k][19]) != float(c[1]) or float(me[k][20]) != float(c[2]): bad += 1; print('MISMATCH X range', mf, k, me[k][19:21], c[1:3])
    print('%s vs %s / %s: octaves %s compared' % (mf, sec, cf, sorted(set(me) & set(th))))
# fap_lean (V_1000) incl. E0 and SIB_V ranges
me = mine('aud_scan_V_closed_Y1000_K21.out'); ln = lean()
for k in sorted(me):
    if k not in ln: continue
    a = me[k][:18]; b = ln[k]
    for i, (x, y) in enumerate(zip(a, b)):
        ncmp += 1
        if float(x) != float(y): bad += 1; print('MISMATCH lean', k, i, x, y)
print('fap_lean_V1000_k27.out: octaves %s compared (incl. E0/N and SIB_V/N ranges)' % sorted(set(me) & set(ln)))
print('numbers compared: %d ; mismatches: %d' % (ncmp, bad))
sys.exit(1 if bad else 0)
