# audit-5129, step 3: the lists of (B), the omissions of (C), the windows of (D), from the strong classes of
# order <= 8 (classes_le8.npz of step 2).  Complete for D <= 73 once f(9) = 75 is known (step 4).
import numpy as np
import supply as S

Z = np.load("classes_le8.npz")
strong = []
nstrong_with_H = {}
for n in range(3, 9):
    H = Z["H%d" % n]; st = Z["strong%d" % n]; P = Z["Pi2%d" % n]; U = Z["U%d" % n]
    for r in range(H.size):
        if st[r]:
            strong.append((n, int(H[r]), int(P[r]), [int(x) for x in U[r]]))
            nstrong_with_H[int(H[r])] = nstrong_with_H.get(int(H[r]), 0) + 1
DMAX = 73
Root, PI, Aord, Suf, S0, minord = S.build(strong, DMAX)
DL = {D: S.delta(S0[D]) for D in S0}

print("== attained odd numbers <= 73 ==")
att = sorted(Aord)
print("attained:", att)
print("not attained:", [d for d in range(1, DMAX + 1, 2) if d not in Aord])
print("least orders of a tournament with H = d:", {d: Aord[d] for d in att})

print()
print("== (B) odd parts of the slopes with least exponents, offsets 3, 5, 9, 17, 33, 65 ==")
for D in (3, 5, 9, 17, 33, 65):
    dl = DL[D]
    print("offset %d: %d odd parts: %s" % (D, len(dl), " ".join("%d:%d" % (d, dl[d]) for d in sorted(dl))))
print()
print("== number of odd parts |Delta(D)| for every attained D <= 73 ==")
print(" ".join("%d:%d" % (D, len(DL[D])) for D in sorted(DL)))

print()
print("== dyadic slopes (odd part 1) at offsets <= 73 ==")
for D in sorted(DL):
    if 1 in DL[D]:
        print("offset %d: least dyadic slope 2^%d = %d" % (D, DL[D][1], 2 ** DL[D][1]))

print()
print("== exponent of the odd part 9 at the Fermat-type offsets (9-grid sources) ==")
for s, D in enumerate((3, 5, 9, 17, 33, 65), start=1):
    print("level s=%d offset %d: 9 listed with least exponent %s; kps Prop.4 lower bound max(0, s-3) = %d"
          % (s, D, DL[D].get(9), max(0, s - 3)))
print("(65, 72) is a slope:", S.is_slope(DL[65], 72), "; (65, 288):", S.is_slope(DL[65], 288),
      "; (65,144):", S.is_slope(DL[65], 144))

print()
print("== (C) r(N) for odd N in [3, 225) ==")
r = {}
for N in range(3, 225, 2):
    c = 0
    for D in range(1, N // 3 + 1, 2):
        if D in DL and S.is_slope(DL[D], (N - D) // 2):
            c += 1
    r[N] = c
zeros = [N for N in r if r[N] == 0]
print("r(N) = 0 exactly for N in", zeros)
blk = [r[N] for N in range(129, 225, 2)]
print("on [129, 225): min r = %d (at N = %s), mean = %.4f, median = %s" %
      (min(blk), [N for N in range(129, 225, 2) if r[N] == min(blk)], float(np.mean(blk)), float(np.median(blk))))
print("r(N), N = 3..223:", " ".join("%d:%d" % (N, r[N]) for N in sorted(r)))
# max D used: (C) needs D <= N/3 < 75
print("largest offset needed: floor(223/3) =", 223 // 3)

print()
print("== (D) window counts: slopes in [h, 3h/2] / integers in [h, 3h/2], attained h <= 73 ==")
share = {}
line = []
for h in att:
    lo, hi = h, (3 * h) // 2
    ints = hi - lo + 1
    k = sum(1 for u in range(lo, hi + 1) if S.is_slope(DL[h], u))
    share[h] = k / ints
    line.append("%d:%d/%d" % (h, k, ints))
print(" ".join(line))
print("whole window present for h in", [h for h in att if share[h] == 1.0])
hs = [h for h in att if 37 <= h <= 73]
sh = [share[h] for h in hs]
print("37 <= h <= 73: mean share %.4f, min %.4f (h = %d), max %.4f (h = %d)" %
      (float(np.mean(sh)), min(sh), hs[int(np.argmin(sh))], max(sh), hs[int(np.argmax(sh))]))
hs2 = [h for h in att if 23 <= h <= 73]
x = np.array([share[h] for h in hs2]); y = np.array([nstrong_with_H.get(h, 0) for h in hs2], dtype=float)
print("23 <= h <= 73: correlation(share, number of strong classes with H = h) = %.4f (n = %d)" %
      (float(np.corrcoef(x, y)[0, 1]), len(hs2)))
print("number of strong classes (order <= 8) with H = h, h <= 73:",
      " ".join("%d:%d" % (h, nstrong_with_H.get(h, 0)) for h in att))
print("offsets h <= 73 with exactly two strong classes:", [h for h in att if nstrong_with_H.get(h, 0) == 2])
print("offsets with at most two strong classes:", [h for h in att if nstrong_with_H.get(h, 0) <= 2])
order_by_share = sorted(hs2, key=lambda h: share[h])
print("poorest shares on 23 <= h <= 73:", [(h, round(share[h], 3), nstrong_with_H.get(h, 0)) for h in order_by_share[:8]])

print()
print("== census of orders <= 8 against the product law (both directions) ==")
census = {m: set() for m in range(1, 9)}
for n in range(1, 9):
    H = Z["H%d" % n]; U = Z["U%d" % n]
    for rr in range(H.size):
        if H[rr] <= DMAX:
            for yv in range(n):
                census[n].add((int(H[rr]), int(U[rr, yv])))
bad1 = 0
bad2 = 0
tot = 0
for m in range(1, 9):
    for (D, U) in census[m]:
        tot += 1
        mo = S.minorder(S0[D], U) if D in S0 else None
        if mo is None or mo > m:
            bad1 += 1
    # predicted pairs of least order exactly m must be in the census of order m
    for D in S0:
        for s0, o in S0[D].items():
            if o <= m:
                Uv = s0 << (m - o)
                if S.minorder(S0[D], Uv) <= m and (D, Uv) not in census[m]:
                    bad2 += 1
allp = set().union(*census.values())
print("distinct census pairs with D <= 73, orders <= 8: %d; census pairs outside the product-law set or below the"
      " predicted least order: %d; predicted pairs missing from the census: %d" % (len(allp), bad1, bad2))

print()
print("== least order needed for each odd part (any exponent), offsets <= 73 ==")
late = []
for D in sorted(S0):
    for d in sorted(DL[D]):
        best = min(o for s0, o in S0[D].items() if S.odd_part(s0)[0] == d)
        if best >= 10:
            late.append((D, d, best))
print("odd parts whose least realisation has order >= 10:", late)
late9 = []
for D in sorted(S0):
    for d in sorted(DL[D]):
        best = min(o for s0, o in S0[D].items() if S.odd_part(s0)[0] == d)
        if best == 9:
            late9.append((D, d))
print("number of (offset, odd part) first realised at order 9:", len(late9))
print("total number of (offset, odd part) pairs, D <= 73:", sum(len(DL[D]) for D in DL))

# full lists, for the comparison with the committed output
with open("s3_slope_lists.txt", "w") as f:
    for D in sorted(DL):
        f.write("%d | %s\n" % (D, " ".join("%d:%d" % (d, DL[D][d]) for d in sorted(DL[D]))))
print("wrote s3_slope_lists.txt")
