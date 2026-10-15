#!/usr/bin/env python3
"""Audit C: compare b(p) from the independent eta/theta construction (audit_newform_bp.json) with
(i) PARI mfcoefs from this audit's own gp script (audit_gp_bp.txt), (ii) the audited session's
g_trace_a2_half.txt (Re(a_p^2) for primes 5..30000), (iii) the a_p list printed in the audited g_coeffs.out."""
import json
D = '/Users/e/Documents/GitHub/math-wt-wilson-20261014/04-computation/experiments/'
mine = {int(k): tuple(v) for k, v in json.load(open(D + 'audit_C_20261014/audit_newform_bp.json')).items()}
gp = {}
for line in open(D + 'audit_C_20261014/audit_gp_bp.txt'):
    p, a, b = map(int, line.split())
    gp[p] = (a, b)
bad_gp = [p for p in gp if mine[p] != gp[p]]
print(f"(i) PARI (own script) vs eta/theta construction: {len(gp)} primes 3..30000, mismatches {bad_gp[:10]}")
txt = open(D + 'paley_k6p6_20261014/g_trace_a2_half.txt').read().strip().strip('[]')
re_a2 = [int(s) for s in txt.split(',')]
ps = sorted(p for p in mine if 5 <= p <= 30000)
assert len(ps) == len(re_a2)
bad_s = [p for p, r in zip(ps, re_a2) if mine[p][0] ** 2 - mine[p][1] ** 2 != r]
print(f"(ii) session g_trace_a2_half.txt (Re a_p^2) vs eta/theta construction: {len(ps)} primes, mismatches {bad_s[:10]}")
bad_c = []
n = 0
for line in open(D + 'paley_k6p6_20261014/g_coeffs.out'):
    parts = line.split()
    if len(parts) == 3 and parts[0].isdigit():
        p, a, b = map(int, parts)
        n += 1
        if mine[p] != (a, b):
            bad_c.append(p)
print(f"(iii) session g_coeffs.out a_p (p < 100): {n} primes, mismatches {bad_c}")
