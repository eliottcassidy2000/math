"""aud_numbers.py - audit task 3: recompute the quoted run statistics from the data files of the packet
(runs/*.txt, runs/*.err, certificate files).  Reads only; prints the numbers next to the quoted ones."""
import os, re, collections
HERE = os.path.dirname(os.path.abspath(__file__))
R = os.path.join(HERE, "..", "runs")
EXC33 = [245, 853, 909, 1049, 1407, 1467, 1515, 1577, 1629, 1829, 1887, 2023, 2037, 2085, 2343, 2393, 2703, 3001, 3033,
         3123, 3313, 3411, 3427, 3545, 3713, 3733, 3753, 3773, 3907, 3909, 3985, 3991, 4023]


def load(fn):
    """lines 'S sign c v2 depth' -> dict sign -> {c: (v2, depth)} (first occurrence kept, all occurrences listed)."""
    d = {1: collections.defaultdict(list), -1: collections.defaultdict(list)}
    other = []
    for l in open(os.path.join(R, fn)):
        p = l.split()
        if p and p[0] == "S":
            d[int(p[1])][int(p[2])].append((int(p[3]), int(p[4])))
        else:
            other.append(l.rstrip())
    return d, other


for fn in ("e1_d2.txt", "e1_d4.txt", "e1x_d3.txt", "e2_d2.txt", "e2_d3.txt", "e2o8_d2.txt", "e3_d6.txt"):
    d, other = load(fn)
    odd_ok = all(c % 2 == 1 and c <= 4097 for s in d for c in d[s])
    print("%-12s  +1: %4d odd c ; -1: %4d odd c ; all c odd and <= 4097: %s ; other lines: %d"
          % (fn, len(d[1]), len(d[-1]), odd_ok, len(other)))
    for o in other[:3]:
        print("      ", o[:160])
    # per depth (cumulative) counts
    for s in (1, -1):
        bydepth = collections.Counter()
        for c, lst in d[s].items():
            bydepth[min(x[1] for x in lst)] += 1
        cum = 0
        row = []
        for dep in sorted(bydepth):
            cum += bydepth[dep]
            row.append("d<=%d:%d" % (dep, cum))
        print("       sign %+d cumulative by least depth: %s" % (s, " ".join(row)))

d4, _ = load("e1_d4.txt")
d2, _ = load("e1_d2.txt")
print("E1: N_2 =", sum(1 for c, l in d4[1].items() if min(x[1] for x in l) <= 2), " N_4 =", len(d4[1]),
      " (quoted 805, 1612); offset -1 at depth <= 3:", len(d4[-1]), "(quoted 1,330)")
unc = [c for c in range(3, 130, 2) if c not in d4[1]]
print("E1: odd c in [3,129] without a +1 state of depth <= 4:", unc, " (quoted 5, 11, 23, 29, 73, 85, 91, 95, 127)")
print("E1: least uncovered c >= 3:", min(c for c in range(3, 4098, 2) if c not in d4[1]), "(quoted 5)")
pw2 = [(c, l) for c, l in d4[-1].items() if c == 1]
print("E1: -1 states with c = 1 (slope a power of two):", pw2, "(quoted: 2^23, 2^24 or 2^25 at depth 3)")
print("E1: c = 19 with +1:", d4[1].get(19), "(quoted: slope 19 * 2^6 at depth 2)")
dx, _ = load("e1x_d3.txt")
union = set(d4[1]) | set(dx[1])
print("E1 + extra transitive pair starts: union of +1 values:", len(union), "(quoted in HYP-9662: 1,649); 29 in extra run:", dx[1].get(29))
e2, _ = load("e2_d2.txt")
print("E2: +1:", len(e2[1]), " -1:", len(e2[-1]), "(quoted 1,872 and 1,782); exceptional values with a seed:",
      sum(1 for c in EXC33 if c in e2[1]), "of 33 ; missing:", [c for c in EXC33 if c not in e2[1]])
e2o8, _ = load("e2o8_d2.txt")
print("E2 with order-8 starts: +1:", len(e2o8[1]), " -1:", len(e2o8[-1]), "; 3773:", e2o8[1].get(3773), "; missing of the 33:",
      [c for c in EXC33 if c not in e2o8[1]])
e2d3, _ = load("e2_d3.txt")
print("E2 three levels: +1:", len(e2d3[1]), " -1:", len(e2d3[-1]), "; missing of the 33:", [c for c in EXC33 if c not in e2d3[1]])
e3, _ = load("e3_d6.txt")
n2 = sum(1 for c, l in e3[1].items() if min(x[1] for x in l) <= 2)
n4 = sum(1 for c, l in e3[1].items() if min(x[1] for x in l) <= 4)
print("E3: N_2' =", n2, " N_4' =", n4, " N_6' =", len(e3[1]), "(quoted 614, 1357, 1592); -1 at depth <= 5:", len(e3[-1]), "(quoted 1,554)")
for c in (5, 11, 23, 73, 85, 91, 95, 127, 29, 125):
    print("     E3 c = %d: %s" % (c, e3[1].get(c)))
print("E3: odd c <= 129 not occurring with +1:", [c for c in range(3, 130, 2) if c not in e3[1]])
print("E3: depth six adds", len(e3[1]) - n4, "values to the", n4, "of depth four (quoted 235 and 1,357)")

# states files
for fn in ("e2_states_order_le7.txt", "e2_states_order8.txt"):
    rows = [tuple(map(int, l.split())) for l in open(os.path.join(R, fn)) if l.strip()]
    print(fn, ": %d states (r g M E); signs: %s" % (len(rows), collections.Counter(r[0] for r in rows)))
    hits = [r for r in rows if r[0] == 1 and r[1] == 16 and r[2] == 64]
    print("     states with (r, h, M) = (1, 15, 64):", hits)
    pw = sorted(set(r[1] - r[0] for r in rows if r[0] == 1 and r[2] & (r[2] - 1) == 0 and r[2] > 0))
    print("     H(B) of seeds (r = 1) with a power of two as slope:", pw)
    pwm = sorted(set((r[1] - r[0], r[2]) for r in rows if r[0] == -1 and r[2] & (r[2] - 1) == 0 and r[2] > 0))
    print("     (H(B), M) of anti-seeds with a power of two as slope:", pwm)

# src9grid / src43grid
for fn, d in (("src9grid.txt", 9), ("src43grid.txt", 43)):
    rows = [tuple(map(int, l.split())) for l in open(os.path.join(R, fn)) if l.strip() and not l.startswith("#")]
    lev = collections.defaultdict(list)
    bad = 0
    for (m, code, Om, r, h, M) in rows:
        s = (r - 1).bit_length() - 1
        if r - 1 != 2 ** s:
            bad += 1
            continue
        x = M
        v = 0
        while x % 2 == 0:
            x //= 2
            v += 1
        if x != d:
            bad += 1
            continue
        lev[s].append((m, code, Om, r, h, M, v, bin(Om).count("1")))
    print(fn, ": %d lines (m code Omask offset H(B) slope); malformed for the %d-grid: %d" % (len(rows), d, bad))
    for s in sorted(lev):
        L = lev[s]
        ok_v = [x for x in L if x[6] <= s]
        single = [x for x in ok_v if x[7] == 1]
        best = min(ok_v, key=lambda x: (x[0], x[7])) if ok_v else None
        bs = min(single, key=lambda x: x[0]) if single else None
        print("     level %2d: %3d lines, %3d with v <= s, %3d of them one-vertex; least order with v <= s: %s ; least one-vertex: %s"
              % (s, len(L), len(ok_v), len(single), best[:7] if best else None,
                 (bs[:7], "excess U-H = %d" % (bs[5] - bs[4])) if bs else None))
# twin file
print("twin_le8.txt lines:", sum(1 for l in open(os.path.join(R, "twin_le8.txt")) if l.strip()))
