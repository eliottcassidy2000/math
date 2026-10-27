"""Auditor's driver: compares type tables and recomputes every ledger number of THM-4816 / THM-4817.

usage: python3 audit_main.py [use13=mine|producer]
"""
import os, sys
from fractions import Fraction as Fr
from math import factorial, comb
from collections import Counter
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from audit_ledger import *

use13 = sys.argv[1] if len(sys.argv) > 1 else ('mine' if os.path.exists(os.path.join(HERE, 'my_nr13.txt')) and os.path.getsize(os.path.join(HERE, 'my_nr13.txt')) > 0 else 'producer')
NB = 16
A, SC = burnside_table(NB)
print("== Burnside (auditor) ==")
for n in range(1, NB + 1):
    print(f"n={n:2d} A={A[n]:,} SC={SC[n]:,}")

klein = {k: parse_types(os.path.join(KLEIN, f'prime_types_n{k}.txt')) for k in range(3, 11)}
mine = {}
for k in range(3, 14):
    fn = os.path.join(HERE, f'my_nr{k}.txt')
    if os.path.exists(fn) and os.path.getsize(fn) > 0 or k == 4:
        mine[k] = parse_types(fn) if os.path.exists(fn) else []
prod_nr = {k: parse_types(os.path.join(PROD, f'nr{k}.txt')) for k in range(3, 14)}
prod_t = {k: parse_types(os.path.join(PROD, f't{k}.txt')) for k in range(3, 14)}

print("\n== non-rigid rows: auditor census vs producer nr{n} vs klein census (non-rigid rows) ==")
for k in range(3, 14):
    if k not in mine:
        print(f"n={k}: auditor census not available")
        continue
    dm = table_dict(mine[k])
    dp = table_dict(prod_nr[k])
    line = f"n={k}: auditor {sum(dm.values())} classes in {len(dm)} rows; producer {sum(dp.values())} in {len(dp)} rows; identical={dm == dp}"
    if k <= 10:
        dk = table_dict([r for r in klein[k] if r[2] > 1])
        line += f"; klein non-rigid {sum(dk.values())} in {len(dk)} rows, identical={dk == dm}"
    print(line)

print("\n== producer full tables t{n} vs klein census files (n=3..10) ==")
for k in range(3, 11):
    print(f"n={k}: identical={table_dict(prod_t[k]) == table_dict(klein[k])}")

# ---------------- ledger ----------------
N = 14
L = Ledger(A, SC, N)
tower = {}
src = {}
for n in range(2, N + 1):
    L.close_order(n)
    if n < 3:
        continue
    if n <= 10:
        nr = [r for r in mine[n]] if n in mine else []
        rig, rn, rs = L.rigid_rows(n, nr)
        mytab = table_dict(nr + rig)
        same_klein = (mytab == table_dict(klein[n]))
        print(f"n={n}: auditor tower table (own census + own ledger rigid rows) == klein census: {same_klein}")
        L.types[n] = klein[n]          # use the census itself below order 11
        src[n] = 'klein census'
    elif n <= 13:
        if n == 13 and use13 == 'producer':
            nr = prod_nr[13]
            src[n] = 'producer nr13 + auditor ledger'
        else:
            nr = mine[n]
            src[n] = 'auditor census + auditor ledger'
        rig, rn, rs = L.rigid_rows(n, nr)
        L.types[n] = nr + rig
        print(f"n={n}: rigid non-SC {rn:,}, rigid SC {rs:,}; table == producer t{n}: {table_dict(L.types[n]) == table_dict(prod_t[n])}")
    tower[n] = L.types.get(n)

cs, c = intransitivity(L.types, L.P, L.SCP, N)
slab, npl, Plab, lab = labelled(N)

print("\n== ledger (auditor) ==")
qty = [
    ("A(n)", lambda n: A[n]), ("non-strong", lambda n: L.NS[n]), ("strong non-prime", lambda n: L.SNP[n]),
    ("P(n)", lambda n: L.P[n]), ("S(n)", lambda n: L.S[n]), ("SC(n)", lambda n: SC[n]),
    ("SC non-strong", lambda n: L.SCns[n]), ("SC strong non-prime", lambda n: L.SCsnp[n]),
    ("SCP(n)", lambda n: L.SCP[n]), ("c_s(n)", lambda n: cs[n]), ("c(n)", lambda n: c[n]),
    ("s_lab(n)", lambda n: slab[n]), ("labelled strong non-prime", lambda n: npl[n]), ("P_lab(n)", lambda n: Plab[n]),
    ("l(n)", lambda n: lab[n]),
]
for name, f in qty:
    print(f"{name:28s} " + " | ".join(f"n={n}: {fmt(f(n))}" for n in range(3, N + 1)))

print("\n== labelled-count identity  sum count*n!/|Aut| over the table == P_lab(n) ==")
for n in range(3, 14):
    tab = L.types[n]
    s = sum(Fr(r[1] * factorial(n), r[2]) for r in tab)
    print(f"n={n}: sum={fmt(s)}  P_lab={fmt(Plab[n])}  equal={s == Plab[n]}  ({src[n]})")

print("\n== Burnside defects (n! P(n) - P_lab(n))/n! ==")
for n in range(3, N + 1):
    d = (factorial(n) * L.P[n] - Plab[n]) / factorial(n)
    line = f"n={n}: {fmt(d)} = {float(d):.6f}"
    if n in L.types and L.types[n] is not None:
        cen = sum(Fr(r[1]) * (1 - Fr(1, r[2])) for r in L.types[n] if r[2] > 1)
        line += f"   census sum (1-1/|Aut|) = {fmt(cen)}  equal={cen == d}"
    print(line)

print("\n== self-checks ==")
for n in range(3, N + 1):
    print(f"n={n}: P-SCP even: {(L.P[n] - L.SCP[n]) % 2 == 0};  P/A={float(L.P[n] / A[n]):.4f}  P_lab/2^C={float(Plab[n] / 2 ** comb(n, 2)):.4f}  c/A={float(c[n] / A[n]):.4f}")

# machine-readable dump
with open(os.path.join(HERE, 'audit_ledger_values.txt'), 'w') as fo:
    for name, f in qty:
        for n in range(3, N + 1):
            v = f(n)
            fo.write(f"{name}\t{n}\t{v}\n")
print("\nsources:", src)
