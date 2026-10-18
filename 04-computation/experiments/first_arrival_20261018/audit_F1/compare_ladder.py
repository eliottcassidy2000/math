#!/usr/bin/env python3
"""audit F1: compare the independent ladder E^(k)/N (aud_point_V1000.out) with fap_levels receipts (../fap_levels_V1000*.out)."""
import re, os, sys
D = os.path.dirname(os.path.abspath(__file__)); R = os.path.join(D, '..')
mine = {}
cur = None
for l in open(os.path.join(D, 'aud_point_V1000.out')):
    m = re.match(r'==== N = (\d+) \(M', l)
    if m: cur = int(m.group(1))
    if l.startswith('ladder E^(k)/N:') and cur:
        mine[cur] = {int(a): b for a, b in re.findall(r'\[(\d+)\]([\d.]+)', l)}
theirs = {}
for fn in ('fap_levels_V1000.out', 'fap_levels_V1000_scales.out'):
    cur = None
    for l in open(os.path.join(R, fn)):
        m = re.match(r'N=(\d+) ', l)
        if m: cur = int(m.group(1)); theirs.setdefault(cur, {}); continue
        m = re.match(r'\s+(\d+): ([\d.]+)\s', l)
        if m and cur: theirs[cur][int(m.group(1))] = m.group(2)
n = bad = 0
for N in sorted(set(mine) & set(theirs)):
    for k in sorted(set(mine[N]) & set(theirs[N])):
        n += 1
        if mine[N][k] != theirs[N][k]: bad += 1; print('MISMATCH N=%d k=%d mine %s receipt %s' % (N, k, mine[N][k], theirs[N][k]))
    print('N=%d: %d ladder levels compared' % (N, len(set(mine[N]) & set(theirs[N]))))
print('ladder values compared: %d ; mismatches: %d' % (n, bad))
sys.exit(1 if bad else 0)
