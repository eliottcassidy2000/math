# audit-5129, step 9: smaller checks.
#  (1) rooted classes of order 8 with the pairs (257, 288), (129, 288); orders at which (129,144), (513,576),
#      (257,288) first occur;
#  (2) least order of the least-exponent slope of each odd part, offsets <= 73 (is it in the census of order <= 9?);
#  (3) the four head triples at non-sink roots, with the base orders and one explicit base each;
#  (4) the family (2^s + 1, 9 * 2^(s-3)) for s = 3..9 in the censuses of order <= 9.
import numpy as np
import acore as A
import supply as S
from s2_classes import canon

Z = np.load("classes_le8.npz")
p9 = np.load("pairs9.npy")
print("== (1) rooted classes ==")
for (D, U) in ((257, 288), (129, 288), (129, 144), (513, 576), (65, 72), (65, 288), (33, 36), (17, 18), (9, 9), (5, 9)):
    for n in range(1, 9):
        H = Z["H%d" % n]; Uv = Z["U%d" % n]; codes = Z["codes%d" % n]
        hits = np.argwhere((H[:, None] == D) & (Uv == U))
        if len(hits):
            rc = set()
            adj = np.zeros(n, np.int64)
            for (r, y) in hits:
                A.adj_from_code(codes[r], n, adj)
                rc.add(int(canon(adj, n, int(y))))
            print("pair (%d, %d): order %d: %d (class, vertex) hits, %d rooted isomorphism classes" %
                  (D, U, n, len(hits), len(rc)))
    key = (D << 32) | U
    print("pair (%d, %d): present at order 9: %s" % (D, U, bool(np.any(p9 == key))))

print()
print("== (2) least orders ==")
strong = []
for n in range(3, 9):
    H = Z["H%d" % n]; st = Z["strong%d" % n]; P = Z["Pi2%d" % n]; U = Z["U%d" % n]
    for r in range(H.size):
        if st[r]:
            strong.append((n, int(H[r]), int(P[r]), [int(x) for x in U[r]]))
Root, PI, Aord, Suf, S0, minord = S.build(strong, 73)
late = []
hist = {}
for D in sorted(S0):
    DL = S.delta(S0[D])
    for d, j in sorted(DL.items()):
        mo = S.minorder(S0[D], d << j)
        hist[mo] = hist.get(mo, 0) + 1
        if mo >= 10:
            late.append((D, d, j, mo))
print("least order of the least-exponent slope d 2^j, histogram over the 1066 (offset, odd part) pairs:",
      dict(sorted(hist.items())))
print("least-exponent slopes that need order >= 10:", late)
# every pair of the census of order 9 with D <= 73 is in the product-law set with predicted least order <= 9, and
# conversely
sel = p9[(p9 >> 32) <= 73]
bad = 0
for k in sel:
    D = int(k >> 32); U = int(k & 0xFFFFFFFF)
    mo = S.minorder(S0[D], U) if D in S0 else None
    if mo is None or mo > 9:
        bad += 1
pred = set()
for D in S0:
    for s0, o in S0[D].items():
        j = 0
        while o + j <= 9:
            pred.add((D << 32) | (s0 << j))
            j += 1
print("census of order 9, D <= 73: %d pairs; outside the product-law prediction: %d; predicted pairs of least order"
      " <= 9: %d; predicted but absent: %d" % (sel.size, bad, len(pred), len(pred - set(int(x) for x in sel))))

print()
print("== (3) heads of 11, 13, 47, 67 at a root that is not a sink ==")
want = {(5, 1, 7): 11, (5, 1, 9): 13, (17, 1, 31): 47, (33, 7, 41): 67}
seen = {}
for n in range(1, 9):
    H = Z["H%d" % n]; Uv = Z["U%d" % n]; Ev = Z["E%d" % n]; Ue = Z["Ue%d" % n]; codes = Z["codes%d" % n]
    adj = np.zeros(n, np.int64)
    for r in range(H.size):
        for y in range(n):
            t = (int(H[r]), int(Ev[r, y]), int(Uv[r, y]))
            if t in want and t not in seen:
                A.adj_from_code(codes[r], n, adj)
                seen[t] = (n, [int(a) for a in adj], y, int(Ue[r, y]), A.popcount(int(adj[y])))
for t in want:
    n, masks, y, ue, od = seen[t]
    print("N = %d = %d - %d + %d: first at order %d, out-masks %s, root %d (out-degree %d), U_even = %d, "
          "D + 2 U_even = %d" % (want[t], t[0], t[1], t[2], n, masks, y, od, ue, t[0] + 2 * ue))
# all head values at non-sink roots with D <= N/3 would be tails; list all triples giving 11, 13, 47, 67 as heads
for target in (11, 13, 47, 67):
    reps = set()
    for n in range(1, 9):
        H = Z["H%d" % n]; Uv = Z["U%d" % n]; Ev = Z["E%d" % n]
        m = (H[:, None] - Ev + Uv) == target
        for r, y in np.argwhere(m):
            if Ev[r, y] < H[r]:
                reps.add((int(H[r]), int(Ev[r, y]), int(Uv[r, y])))
    print("head triples (D, E, U) with D - E + U = %d and E < D, orders <= 8: %s" % (target, sorted(reps)))

print()
print("== (4) the family (2^s + 1, 9 * 2^(s-3)) ==")
for s in range(3, 12):
    D = 2 ** s + 1; U = 9 * 2 ** (s - 3)
    key = (D << 32) | U
    first = None
    for n in range(1, 9):
        H = Z["H%d" % n]; Uv = Z["U%d" % n]
        if np.any((H[:, None] == D) & (Uv == U)):
            first = n
            break
    if first is None and np.any(p9 == key):
        first = 9
    print("s = %d: (%d, %d): %s" % (s, D, U, "least order %d" % first if first else "absent at orders <= 9"))
