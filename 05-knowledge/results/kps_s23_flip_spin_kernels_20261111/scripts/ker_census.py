"""ker_census.py - signature census and grid carriers from the output of kernels_sat.c (stdin or files).

line: m code x y  N11 N12 N21 N22  det  two      (ports x < y; the pair (y, x) has the transposed matrix)
In the coordinates (h, pi) = (H, Pi_2) of the modules A (at x) and B (at y) the kernel form is
    al h h' + be pi h' + ga h pi' + de pi pi',
    de = N22, be = N21 - 2 N22, ga = N12 - 2 N22, al = N11 - 2 N12 - 2 N21 + 4 N22,   al de - be ga = det.
Transitive modules: H(K[TT_a, TT_b]) = al + be 2^a + ga 2^b + de 2^(a+b).

plus carrier : (al, be, ga, de) = (1, 2^u, 0, d 2^w) up to the order of the ports:  1 + 2^(a+u) + d 2^(a+b+w)
minus carrier: (-1, -2^u, 0, d 2^w):                                              d 2^(a+b+w) - 2^(a+u) - 1"""
import sys, collections


def odd(x):
    v = 0
    while x % 2 == 0:
        x //= 2
        v += 1
    return x, v


def ispow2(x):
    return x > 0 and x & (x - 1) == 0


sig = collections.Counter()
plus = {}
minus = {}
seeds = collections.Counter()
total = 0
files = sys.argv[1:]
streams = [open(f) for f in files] if files else [sys.stdin]
for st in streams:
    for l in st:
        t = l.split()
        if len(t) < 10:
            continue
        m, code, x, y, n11, n12, n21, n22, det, two = map(int, t[:10])
        total += 1
        if not two:
            sig[(m, "more than two visits")] += 2
            continue
        if n12 == 0 and n21 == 0 and n22 == 0:
            sig[(m, "join-like")] += 2
        else:
            sig[(m, "elliptic" if det > 0 else ("parabolic" if det == 0 else "hyperbolic"))] += 2
        de = n22
        if de <= 0:
            continue
        for (p12, p21, px, py) in ((n12, n21, x, y), (n21, n12, y, x)):
            be = p21 - 2 * de
            ga = p12 - 2 * de
            al = n11 - 2 * n12 - 2 * n21 + 4 * de
            if ga == 0 and al == 1 and ispow2(be):
                d, w = odd(de)
                u = be.bit_length() - 1
                key = d
                val = (u, w, m, code, px, py)
                if key not in plus or val[:3] < plus[key][:3]:
                    plus[key] = val
            if ga == 0 and al == -1 and be < 0 and ispow2(-be):
                d, w = odd(de)
                u = (-be).bit_length() - 1
                val = (u, w, m, code, px, py)
                if d not in minus or val[:3] < minus[d][:3]:
                    minus[d] = val
            if be == 0 and ga == 0 and abs(al) == 1:
                seeds[(m, al, odd(de)[0])] += 1
print("kernel-port pairs read (x < y):", total)
for m in sorted({k[0] for k in sig}):
    print("order %d (ordered port pairs):" % m, {k[1]: v for k, v in sorted(sig.items()) if k[0] == m})
print("plus carriers  1 + 2^(a+u) + d 2^(a+b+w)  [module TT_a at the first port listed]:")
for d in sorted(plus):
    u, w, m, code, px, py = plus[d]
    print("   d = %d: u = %d, w = %d, order %d, code %d, ports (%d, %d): N+(d, r, s) for s >= %d, r >= %d" % (d, u, w, m, code, px, py, u + 1, max(1, w - u + 1)))
print("minus carriers  d 2^(a+b+w) - 2^(a+u) - 1:")
for d in sorted(minus):
    u, w, m, code, px, py = minus[d]
    print("   d = %d: u = %d, w = %d, order %d, code %d, ports (%d, %d): N-(d, r, s) for s >= %d, r >= %d" % (d, u, w, m, code, px, py, u + 1, max(1, w - u + 1)))
dg = collections.defaultdict(set)
for (m, al, d), c in seeds.items():
    dg[(m, al)].add(d)
for k in sorted(dg):
    v = sorted(dg[k])
    print("diagonal kernels of order %d with offset %+d: %d odd slopes; first %s" % (k[0], k[1], len(v), v[:25]))
