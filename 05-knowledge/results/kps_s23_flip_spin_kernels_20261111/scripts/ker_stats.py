"""ker_stats.py - summary of the visit matrices written by kernels.c (E2).
Line: m code x y  maxi maxj  N11 N12 N21 N22  det  two [ | full matrix ]"""
import sys, os, collections
here = sys.argv[1] if len(sys.argv) > 1 else os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "runs")
for m in (3, 4, 5, 6, 7):
    fn = os.path.join(here, "ker_%d.txt" % m)
    if not os.path.exists(fn):
        continue
    rows = []
    for l in open(fn):
        head = l.split("|")[0].split()
        mm, code, x, y, maxi, maxj, n11, n12, n21, n22, det, two = map(int, head)
        rows.append((code, x, y, maxi, maxj, n11, n12, n21, n22, det, two))
    two = [r for r in rows if r[10] == 1]
    sig = collections.Counter()
    for r in two:
        det = r[9]
        rank1 = (r[6] == 0 and r[7] == 0 and r[8] == 0)
        sig["join-like (N12=N21=N22=0)" if rank1 else ("elliptic det>0" if det > 0 else ("parabolic det=0" if det == 0 else "hyperbolic det<0"))] += 1
    print("order %d: %d (class, ordered port pair); at most two visits: %d; more: %d" % (m, len(rows), len(two), len(rows) - len(two)))
    print("   signature among at-most-two-visit kernels:", dict(sig))
    dets = collections.Counter(r[9] for r in two if not (r[6] == 0 and r[7] == 0 and r[8] == 0))
    print("   determinants (value: count):", sorted(dets.items()))
    hyp = sorted((r for r in two if r[9] < 0), key=lambda r: (r[5], r[0]))[:6]
    for r in hyp:
        print("   hyperbolic example: code %d ports (%d,%d) G = [[%d,%d],[%d,%d]] det %d" % (r[0], r[1], r[2], r[5], r[6], r[7], r[8], r[9]))
    par = sorted((r for r in two if r[9] == 0 and (r[6] or r[7] or r[8])), key=lambda r: (r[5], r[0]))[:6]
    for r in par:
        print("   parabolic (not join-like) example: code %d ports (%d,%d) G = [[%d,%d],[%d,%d]]" % (r[0], r[1], r[2], r[5], r[6], r[7], r[8]))
    sym = [r for r in two if r[6] == r[7] and r[8] > 0]
    disc = collections.Counter(4 * (r[6] * r[6] - r[5] * r[8]) for r in sym)
    print("   symmetric kernels (N12 = N21, N22 > 0): %d; discriminants 4(N12^2 - N11 N22): %s" % (len(sym), sorted(disc.items(), reverse=True)))
    d4 = [r for r in sym if r[5] * r[8] - r[6] * r[6] == 1][:5]
    for r in d4:
        print("   discriminant -4 example: code %d ports (%d,%d) G = [[%d,%d],[%d,%d]]" % (r[0], r[1], r[2], r[5], r[6], r[7], r[8]))
    more = [r for r in rows if r[10] == 0]
    mx = collections.Counter((r[3], r[4]) for r in more)
    print("   more than two visits: (max i, max j) counts:", sorted(mx.items()))
