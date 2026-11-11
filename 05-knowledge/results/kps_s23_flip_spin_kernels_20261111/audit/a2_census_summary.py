"""a2_census_summary.py (AUDIT) - census table of section 2.4 from the audit's own walk enumeration (a2_census.c),
the list of discriminants of symmetric kernels of order <= 6, carriers found at orders <= 7, and a line-by-line
comparison with the author's runs/ker_m.txt (comparison of OUTPUTS only)."""
import sys, os, collections
HERE = os.path.dirname(os.path.abspath(__file__))
out = []


def P(s):
    print(s)
    out.append(s)


def ispow2(x):
    return x > 0 and x & (x - 1) == 0


def odd(x):
    w = 0
    while x % 2 == 0:
        x //= 2
        w += 1
    return x, w


rows = {}
for m in range(3, 8):
    fn = os.path.join(HERE, "a2_census_%d.out" % m)
    if not os.path.exists(fn):
        continue
    rows[m] = [tuple(int(t) for t in l.split()) for l in open(fn) if l.strip()]
P("order | all | <= 2 visits | join-like | elliptic | parabolic (not join-like) | hyperbolic | det = 0 in all (parabolic by the DEFINITION) | max visits seen (bound m - 1)")
disc_neg, disc_pos, disc_zero = set(), set(), 0
hyper5 = []
plus, minus = {}, {}
for m in sorted(rows):
    c = collections.Counter()
    maxv = 0
    for (mm, code, x, y, maxi, maxj, n11, n12, n21, n22, det, two) in rows[m]:
        c["all"] += 1
        maxv = max(maxv, maxi, maxj)
        if not two:
            continue
        c["two"] += 1
        if n12 == 0 and n21 == 0 and n22 == 0:
            c["join"] += 1
        elif det > 0:
            c["ell"] += 1
        elif det == 0:
            c["par"] += 1
        else:
            c["hyp"] += 1
            if m == 5:
                hyper5.append((code, x, y, n11, n12, n21, n22))
        if det == 0:
            c["det0"] += 1
        if m <= 6 and n12 == n21 and n22 > 0:
            d = 4 * (n12 * n12 - n11 * n22)
            if d < 0:
                disc_neg.add(d)
            elif d > 0:
                disc_pos.add(d)
            else:
                disc_zero += 1
        al = n11 - 2 * n12 - 2 * n21 + 4 * n22
        be = n21 - 2 * n22
        ga = n12 - 2 * n22
        de = n22
        if de > 0 and ga == 0 and al == 1 and ispow2(be):
            plus.setdefault(odd(de)[0], []).append((m, code, x, y, al, be, ga, de))
        if de > 0 and ga == 0 and al == -1 and be < 0 and ispow2(-be):
            minus.setdefault(odd(de)[0], []).append((m, code, x, y, al, be, ga, de))
        if de > 0 and be == 0 and al == 1 and ispow2(ga):
            plus.setdefault(odd(de)[0], []).append((m, code, x, y, al, be, ga, de))
        if de > 0 and be == 0 and al == -1 and ga < 0 and ispow2(-ga):
            minus.setdefault(odd(de)[0], []).append((m, code, x, y, al, be, ga, de))
    P("%d | %d | %d | %d | %d | %d | %d | %d | %d" % (m, c["all"], c["two"], c["join"], c["ell"], c["par"], c["hyp"], c["det0"], maxv))
    if c["two"]:
        P("    share of hyperbolic among at-most-two-visit kernels: %.2f per cent" % (100.0 * c["hyp"] / c["two"]))
P("hyperbolic kernels of order 5: %s" % hyper5)
P("discriminants 4 (N12^2 - N11 N22) of symmetric two-visit kernels (N12 = N21, N22 > 0), orders <= 6:")
P("   negative: %s" % sorted(disc_neg, reverse=True))
P("   zero: %d kernels" % disc_zero)
P("   positive: %s" % sorted(disc_pos))
P("carriers (1, 2^u, 0, d 2^w) at orders <= 7 (either port order): %s" % {d: v[:2] for d, v in sorted(plus.items())})
P("anti-carriers (-1, -2^u, 0, d 2^w) at orders <= 7: %s" % {d: v[:2] for d, v in sorted(minus.items())})
# compare with the author's runs/ker_m.txt
for m in sorted(rows):
    fn = os.path.join(HERE, "..", "runs", "ker_%d.txt" % m)
    if not os.path.exists(fn):
        P("author file ker_%d.txt not present" % m)
        continue
    theirs = {}
    for l in open(fn):
        t = l.split("|")[0].split()
        if len(t) >= 12:
            v = [int(z) for z in t[:12]]
            theirs[(v[1], v[2], v[3])] = tuple(v[4:12])
    mine = {(r[1], r[2], r[3]): tuple(r[4:12]) for r in rows[m]}
    diff = sum(1 for k in mine if theirs.get(k) != mine[k])
    P("order %d: my rows %d, author's rows %d, rows that differ (maxi, maxj, N11, N12, N21, N22, det, two): %d" % (m, len(mine), len(theirs), diff))
open(os.path.join(HERE, "a2_census_summary.out"), "w").write("\n".join(out) + "\n")
