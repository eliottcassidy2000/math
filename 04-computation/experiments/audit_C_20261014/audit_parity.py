#!/usr/bin/env python3
"""Audit C: unconditional parity laws for M_q (proof in REPORT.md, finding m6):
   q = 1 mod 4  =>  M_q = 2 (mod 8)    [so v_2(M_q) = 1 and M_q != 0, with no use of HYP-9276]
   q = 3 mod 4  =>  M_q = 0 (mod 16)
Proof idea: E_a has full 2-torsion so t(a) = -(q+1-#E_a) = q+1 (mod 4) for a not in {0,1}; the a = 1 terms give -(1+eps).
Checked on every prime value in audit_identity_M.json and on the prime powers of audit_prime_powers.out."""
import json, re
D = '/Users/e/Documents/GitHub/math-wt-wilson-20261014/04-computation/experiments/audit_C_20261014/'
M = {int(k): v for k, v in json.load(open(D + 'audit_identity_M.json'))['fft'].items()}
for line in open(D + 'audit_prime_powers.out'):
    m = re.search(r"q = \d+\^\d+ = +(\d+) +eps = [+-]1 +M_q = (-?\d+)", line)
    if m:
        M[int(m.group(1))] = int(m.group(2))
b1 = [q for q, v in M.items() if q % 4 == 1 and v % 8 != 2]
b3 = [q for q, v in M.items() if q % 4 == 3 and v % 16 != 0]
print(f"{sum(q % 4 == 1 for q in M)} fields q = 1 mod 4: M_q = 2 mod 8 fails at {b1[:10]}")
print(f"{sum(q % 4 == 3 for q in M)} fields q = 3 mod 4: M_q = 0 mod 16 fails at {b3[:10]}")
