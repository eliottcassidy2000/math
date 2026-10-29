#!/usr/bin/env python3
"""Results-note inversion-pair claim: compile a11_cycles_all.c, run it for every odd b <= 401 with starts <= 2*10^6,
reduce every cycle point n/b to p/q, keep p, q <= 401, count distinct periodic rationals and mutual pairs p/q, q/p."""
import os, subprocess, sys, tempfile
from math import gcd
here = os.path.dirname(os.path.abspath(__file__))
tmp = tempfile.mkdtemp()
exe = os.path.join(tmp, "a11")
subprocess.run(["cc", "-O2", "-o", exe, os.path.join(here, "a11_cycles_all.c")], check=True)
S = set()
for b in range(1, 402, 2):
    out = subprocess.run([exe, str(b), "2000000"], capture_output=True, text=True).stdout
    for line in out.splitlines():
        if "cycle min" not in line:
            continue
        for n in map(int, line.split("points:")[1].split()):
            g = gcd(n, b); p, q = n // g, b // g
            if p <= 401 and q <= 401:
                S.add((p, q))
mut = {frozenset(((p, q), (q, p))) for (p, q) in S if (q, p) in S and p != q}
print("distinct periodic reduced rationals p/q (p, q <= 401):", len(S), "; mutual pairs (p != q):", len(mut))
