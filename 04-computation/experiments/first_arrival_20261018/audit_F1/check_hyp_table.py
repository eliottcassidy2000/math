#!/usr/bin/env python3
"""audit F1: check every number of the HYP-9325 V_1000 evidence table against the receipts
(../ledger_table.out section fap_V1000_k26_ledger.out for octaves 11..26, ../fap_lean_V1000_k27.out for 27), and against
the independent scan (aud_scan_V_closed_Y1000_K21.out) where available (octaves 11..21).
A table entry passes when it equals the receipt value rounded to the printed number of decimals (|diff| <= half a unit)."""
import re, os, sys
D = os.path.dirname(os.path.abspath(__file__)); R = os.path.join(D, '..')
HYP = os.path.join(D, '..', '..', '..', '..', '05-knowledge', 'hypotheses', 'HYP-9325-one-tick-law-of-first-arrivals-into-the-627-seed-frontier.md')
num = r'[+-]?\d+\.\d+'

def rec_rows(path, start_marker=None):
    rows, on = {}, start_marker is None
    for l in open(path):
        if start_marker is not None and l.startswith('=='): on = (start_marker in l); continue
        m = re.match(r'\s*(\d+) \| (.*)$', l)
        if m and on: rows[int(m.group(1))] = [float(x) for x in re.findall(num, m.group(2))]
    return rows

led = rec_rows(os.path.join(R, 'ledger_table.out'), 'fap_V1000_k26_ledger.out')
lean = rec_rows(os.path.join(R, 'fap_lean_V1000_k27.out'))
mine = rec_rows(os.path.join(D, 'aud_scan_V_closed_Y1000_K21.out'))
# receipt columns: 0 EVmean 1 EVmax 2 G0 3 r1 4 r1lo 5 r1hi 6 z0m 7 z0x 8 gam 9 gamx 10 hlm 11 hlx 12 E+ 13 dip
hyprows = {}
for l in open(HYP, encoding='utf-8'):
    m = re.match(r'\|\s*(\d+) \|(.*)\|\s*$', l)
    if not m: continue
    k = int(m.group(1)); cells = m.group(2)
    vals = re.findall(num, cells.replace('−', '-'))
    hyprows[k] = vals
order = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13]   # HYP columns map to these receipt columns
bad = 0; n = 0
for k in sorted(hyprows):
    src = lean if k == 27 else led
    if k not in src: print('octave', k, 'missing in receipts'); bad += 1; continue
    rec = src[k]
    for hv, ci in zip(hyprows[k], order):
        dec = len(hv.split('.')[1]); h = float(hv); r = rec[ci]
        n += 1
        if abs(h - r) > 0.5 * 10 ** (-dec) + 1e-12:
            bad += 1; print('MISMATCH octave %d col %d: HYP %s receipt %s' % (k, ci, hv, r))
        if k in mine and abs(mine[k][ci] - r) > 1e-12:
            bad += 1; print('independent scan differs: octave %d col %d mine %s receipt %s' % (k, ci, mine[k][ci], r))
print('HYP-9325 table: %d octave rows, %d numbers checked against receipts (and octaves 11-21 against the independent scan); mismatches %d' % (len(hyprows), n, bad))
sys.exit(1 if bad else 0)
