#!/usr/bin/env python3
# audit K1: order-10 self-converse values, three independent sources:
#   (a) auditor's aud_sc.c (reduced anti-invariant model, sigma = (0 1)(2 3)...(8 9), pull DP, Landau strongness);
#   (b) the packet's sc_even10.out (all 2^25 labellings, sigma = (0 5)...(4 9));
#   (c) THM-4922's audit class census at order 10 (aud_o10_byH.txt: per h, all classes and SC classes; canonical
#       augmentation, independent code), which carries no strong flag.
import sys, os
W, PKT, AUD = sys.argv[1], sys.argv[2], sys.argv[3]
a = {}
for line in open(os.path.join(W, "anti_10_red.txt")):
    h, c, s = map(int, line.split()); a[h] = int(s > 0)
b = {}
for line in open(os.path.join(PKT, "sc_even10.out")):
    t, h, s = line.split(); b[int(h)] = int(s)
c = set(); allv = set()
for line in open(AUD):
    if line.startswith("#"): continue
    f = line.split(); h, allc, scc = int(f[0]), int(f[1]), int(f[2]); allv.add(h)
    if scc: c.add(h)
print(f"(a) auditor: {len(a)} SC values, {sum(a.values())} with a strong SC realisation")
print(f"(b) packet : {len(b)} SC values, {sum(b.values())} strong")
print(f"(c) class census: {len(c)} SC values; all values at order 10: {len(allv)}")
print("(a) == (b) including strong flags:", a == b)
print("set(a) == (c):", set(a) == c)
if a != b:
    print(" only a:", sorted(set(a) - set(b))[:30], " only b:", sorted(set(b) - set(a))[:30],
          " flag diffs:", [h for h in a if h in b and a[h] != b[h]][:30])
if set(a) != c:
    print(" only a:", sorted(set(a) - c)[:30], " only c:", sorted(c - set(a))[:30])
print("all order-10 SC values = 1 mod 4 (THM-791 at n = 10):", all(h % 4 == 1 for h in a))
print("least strong SC value at order 10:", min(h for h in a if a[h]), "(Busch f(10) = 125)")
print("non-strong SC values below 625:", sorted(h for h in a if not a[h] and h < 625))
