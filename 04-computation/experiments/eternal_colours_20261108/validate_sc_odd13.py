#!/usr/bin/env python3
# mac-mini-2026-11-08-firstrealizers: validation record of sc_odd13.c. Compares the unpruned validation runs with the
# complete inventories: -DM=2 (orders 6, 7) and -DM=3 (orders 8, 9) with the census; -DM=4 (orders 10, 11, all H) with
# sc_even10.out and audit K1's complete order 11 (audit_K1/out/o11_full.txt), in full, with SHA-256 digests of the
# order-11 value/strong lists. Then the orders 12, 13 run below 1125 with audit K1's order 12 below 625.
# usage: validate_sc_odd13.py WORK   (WORK holds cls_n.txt and v13_m2.out, v13_m3.out, v13_m4.out; run in the packet dir)
import sys, os, hashlib
WORK = sys.argv[1]
def census(n):
    sc = {}
    for line in open(os.path.join(WORK, f"cls_{n}.txt")):
        hd = line.split("|")[0].split(); H, s, c = int(hd[1]), int(hd[3]), int(hd[4])
        if c: sc[H] = sc.get(H, 0) | s
    return sc
def run(fn):
    E, O = {}, {}
    for line in open(os.path.join(WORK, fn)):
        t, h, s = line.split(); (E if t == "E" else O)[int(h)] = int(s)
    return E, O
def digest(d): return hashlib.sha256("".join(f"{h} {s}\n" for h, s in sorted(d.items())).encode()).hexdigest()[:16]
for m, (ne, no) in ((2, (6, 7)), (3, (8, 9))):
    E, O = run(f"v13_m{m}.out")
    for n, got in ((ne, E), (no, O)):
        ref = census(n)
        print(f"-DM={m}, order {n}: {len(got)} values; census {len(ref)} self-converse values; identical with strong flags: {got == ref}")
E, O = run("v13_m4.out")
ref10 = {}
for line in open("sc_even10.out"):
    if line.startswith("S "): _, h, s = line.split(); ref10[int(h)] = int(s)
print(f"-DM=4, order 10: {len(E)} values; sc_even10.out {len(ref10)}; identical with strong flags: {E == ref10}")
ref11 = {}
for line in open("audit_K1/out/o11_full.txt"):
    h, c, s = map(int, line.split()); ref11[h] = int(s > 0)
print(f"-DM=4, order 11: {len(O)} values; audit K1 complete order 11 {len(ref11)}; identical with strong flags: {O == ref11}; "
      f"sha256[:16] {digest(O)} vs {digest(ref11)}")
E12 = {}
for line in open("sc_odd13_1124.out"):
    t, h, s = line.split()
    if t == "E" and int(h) < 625: E12[int(h)] = int(s)
a12 = {}
for line in open("audit_K1/out/o12_lt625.txt"):
    h, c, s = map(int, line.split()); a12[h] = int(s > 0)
print(f"-DM=5, order 12 below 625: {len(E12)} values; audit K1 {len(a12)}; identical with strong flags: {E12 == a12}")
