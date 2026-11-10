#!/usr/bin/env python3
"""AUDIT check (d), spot checks on the auditor's own deletion-pair tables (d_window.c):
 - every pair listed as "missing n" in the packet's out/window_le10.out: at which base orders (2..10) does it occur?
 - the universal floor n = h or n >= (3h+1)/2 (Session214) on every realised h < 1024;
 - the pairs quoted in the README text ((13,27), least ratio (257,387)).
usage: d_spotchecks.py DIR"""
import sys, os, glob, re
D = sys.argv[1] if len(sys.argv) > 1 else "/tmp/claude-audit-ensemble"
HM, NM = 1024, 16384; ROW = NM // 16
def load(fn):
    b = open(fn, "rb").read(); return [int.from_bytes(b[r * ROW:(r + 1) * ROW], "little") for r in range(HM // 2)]
per = {}; A = [0] * (HM // 2)
for K in range(2, 11):
    T = [0] * (HM // 2)
    for fn in sorted(glob.glob(f"{D}/win_{K}.A.bin") + glob.glob(f"{D}/win_{K}_*.A.bin")):
        for r, x in enumerate(load(fn)): T[r] |= x
    per[K] = T
    for r in range(HM // 2): A[r] |= T[r]
def has(T, h, n): return T[h >> 1] >> (n >> 1) & 1
here = os.path.dirname(os.path.abspath(__file__))
claimed = {}
for l in open(os.path.join(here, "..", "out", "window_le10.out")):
    mm = re.match(r"\s*h=\s*(\d+) missing n: \[(.*)\]", l)
    if mm: claimed[int(mm.group(1))] = [int(x) for x in mm.group(2).split(",")]
pairs = [(h, n) for h in sorted(claimed) for n in claimed[h]]
occ = [(h, n, [K for K in range(2, 11) if has(per[K], h, n)]) for h, n in pairs]
print(f"packet's missing pairs: {len(pairs)} pairs over {len(claimed)} values of h; pairs that DO occur at some base order <= 10 in the auditor's tables: {[(h, n, k) for h, n, k in occ if k]}")
print("the entries named in the audit brief: " + ", ".join(f"({h},{n}) absent at orders 2..8: {not any(has(per[K], h, n) for K in range(2, 9))}, absent at 9, 10: {not (has(per[9], h, n) or has(per[10], h, n))}" for h, n in [(3, 11), (5, 21), (9, 35)]))
print("first 14 further entries, all absent at every base order <= 10:", all(not k for h, n, k in occ[:17]), [(h, n) for h, n, k in occ[3:17]])
viol = 0; onfloor = []
for h in range(1, HM, 2):
    if not has(A, h, h): continue
    n = h + 2
    while n < NM and not has(A, h, n): n += 2
    if n < NM:
        if 2 * n < 3 * h + 1: viol += 1
        if 2 * n in (3 * h + 1, 3 * h + 3): onfloor.append((h, n))
print(f"realised h < 1024: least proper extension below (3h+1)/2: {viol} cases; least proper extension on the floor (2n = 3h+1 or 3h+3): {onfloor}")
print("pair (13,27) occurs with a base of order <= 10:", bool(has(A, 13, 27)), "; extension values of the realisers of 13 up to 45:", [n for n in range(13, 46, 2) if has(A, 13, n)])
print("pair (257,387) [least ratio 1.5058 for h >= 237]: base orders at which it occurs:", [K for K in range(2, 11) if has(per[K], 257, 387)])
print("pair (39,183) [missing with bases of order <= 8, present in the packet's order <= 10 census]: base orders:", [K for K in range(2, 11) if has(per[K], 39, 183)])
