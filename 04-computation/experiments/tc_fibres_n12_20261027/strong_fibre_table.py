"""strong_fibre_table.py -- mac-mini-2026-10-27-potentials
Assemble the complete strong presentation fibres F_s(v) = {strong classes with tc = v} for odd v <= 23 (THM-4812 (g)):
orders <= 10 from the strong census (tc_strong_fibres.cpp outputs tcstrong_n.txt for n = 3..10), orders >= 11 empty by
fibre_scan_v5_21.out (all 567 jobs, no hit) together with THM-4812 (d) for v in {7, 21} and the order bounds
(type_lemma style counting) for every v <= 21.  Usage: python3 strong_fibre_table.py DIR_WITH_tcstrong_FILES"""
import sys
d = sys.argv[1]
tab = {}
for n in range(3, 11):
    for line in open(f"{d}/tcstrong_{n}.txt"):
        if line.startswith('#'):
            continue
        v, c = map(int, line.split())
        tab.setdefault(v, {})[n] = c
print("strong presentation fibres F_s(v), v odd <= 23 (K1 counted separately for v = 1):")
for v in range(1, 24, 2):
    row = tab.get(v, {})
    total = sum(row.values()) + (1 if v == 1 else 0)
    print(f"  v = {v:2d}: |F_s(v)| = {total:3d}   by order: {dict(sorted(row.items()))}" + ("  (+ K1)" if v == 1 else ""))
