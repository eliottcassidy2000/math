"""e1_summary.py - E1 / F1 of the S23 pre-registration and the tables of the pair spectrum.

input: the output of pairspec (lines "m code H Pi2") for every tournament class of order <= 9
       (class lists by classgen.c:  classgen 9 cls ;  pairspec m < cls_m.txt  for m = 1..9).
Pi_2(X) = number of ordered pairs of disjoint possibly empty paths covering X."""
import sys, collections


def odd(x):
    v = 0
    while x % 2 == 0:
        x //= 2
        v += 1
    return x, v


rows = [tuple(map(int, l.split())) for l in open(sys.argv[1])]
by = collections.defaultdict(dict)
for (m, code, h, p) in rows:
    o, e = odd(p)
    cur = by[h].get(o)
    if cur is None or (m, e) < cur[:2]:
        by[h][o] = (m, e, code)
print("classes:", len(rows), "| distinct H:", len(by), "| distinct pairs (H, odd part of Pi_2):", sum(len(v) for v in by.values()))
print("Pi_2 = (-1)^n H (mod 3) [THM-4756, step 2]: violations", sum(1 for (m, code, h, p) in rows if (p - (-1) ** m * h) % 3))
print()
print("F1. levels s with a class of order <= 9, H = 2^s + 1 and odd part of Pi_2 dividing 9:")
hit = 0
for s in range(1, 12):
    h = 2 ** s + 1
    d = by.get(h, {})
    good = sorted((o, d[o]) for o in d if 9 % o == 0)
    if 3 <= s <= 11 and good:
        hit += 1
    print("  s = %2d, H = %4d: %3d odd parts; least odd part %d; dividing 9: %s"
          % (s, h, len(d), min(d) if d else 0, ["odd part %d, order %d, Pi_2 = %d * 2^%d, code %d" % (o, v[0], o, v[1], v[2]) for o, v in good]))
print("  levels 3..11 with such a class: %d of 9  (prediction F1: at least 7)" % hit)
print()
print("classes with Pi_2 a power of two (odd part 1), least witness per H:")
for h in sorted(by):
    if 1 in by[h]:
        m, e, code = by[h][1]
        print("  H = %d: order %d, Pi_2 = 2^%d, code %d" % (h, m, e, code))
print()
print("odd parts of Pi_2 for H <= 33:", {h: sorted(by[h]) for h in sorted(by) if h <= 33})
print()
for mm in range(1, 10):
    r = [(p / h, h, p) for (m, code, h, p) in rows if m == mm]
    lo = min(r)
    print("order %d: least Pi_2 / H = %.4f (H = %d, Pi_2 = %d), greatest = %d (transitive)" % (mm, lo[0], lo[1], lo[2], 2 ** mm))
