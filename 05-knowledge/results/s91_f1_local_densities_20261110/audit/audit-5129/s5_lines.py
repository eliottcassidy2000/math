# audit-5129, step 5: line counts.  r_m(N) = number of lattice points (D, U) on the line D + 2U = N that are
# pairs of a rooted tournament of order <= m.  Two readings of "(with appended sinks)":
#   (a) pairs (D, U 2^j), (D, U) of order <= m, j >= 0 arbitrary  (sinks appended beyond order m);
#   (b) pairs of order <= m only.
# Exact r(N) (all orders) from the product law for N < 225 (offsets <= 73, step 3) as the reference.
import sys
import numpy as np
import supply as S

NMAX = 1 << 15
Z = np.load("classes_le8.npz")
pairs = {}
for n in range(1, 9):
    H = Z["H%d" % n]; U = Z["U%d" % n]
    k = ((H[:, None].astype(np.int64) << 32) | U.astype(np.int64)).ravel()
    pairs[n] = np.unique(k)
pairs[9] = np.load("pairs9.npy")
try:
    p7 = np.load("pairs7.npy")
    print("order 7: pairs from the class list %d, from all 2^21 labelled tournaments %d, equal: %s" %
          (pairs[7].size, p7.size, np.array_equal(pairs[7], p7)))
except FileNotFoundError:
    print("pairs7.npy not present yet")
for n in range(1, 9):
    print("order %d: %d distinct pairs (D,U); contained in the pairs of order %d: %s" %
          (n, pairs[n].size, n + 1, bool(np.all(np.isin(pairs[n], pairs[n + 1])))))
print("order 9: %d distinct pairs" % pairs[9].size)
extra = sys.argv[1:]
if extra:
    pairs[10] = np.load(extra[0])
    print("order 10: %d distinct pairs (from %s); contains the pairs of order 9: %s" %
          (pairs[10].size, extra[0], bool(np.all(np.isin(pairs[9], pairs[10])))))


def line_counts(keys, sinks):
    D = (keys >> 32).astype(np.int64)
    U = (keys & 0xFFFFFFFF).astype(np.int64)
    mark = np.zeros((NMAX // 2, NMAX // 6 + 2), dtype=bool)   # [(N-1)/2, (D-1)/2]
    j = 0
    while True:
        N = D + 2 * (U << j)
        ok = N < NMAX
        if not ok.any():
            break
        mark[(N[ok] - 1) // 2, (D[ok] - 1) // 2] = True
        if not sinks:
            break
        j += 1
    return mark.sum(axis=1)      # index (N-1)/2


def lattice(N):
    return (N // 3 + 1) // 2


Ls = np.array([lattice(N) for N in range(1, NMAX, 2)])
res = {}
for m in sorted(pairs):
    if m < 7:
        continue
    res[(m, "a")] = line_counts(pairs[m], True)
    res[(m, "b")] = line_counts(pairs[m], False)

# exact r(N), N < 225
strong = []
for n in range(3, 9):
    H = Z["H%d" % n]; st = Z["strong%d" % n]; P = Z["Pi2%d" % n]; U = Z["U%d" % n]
    for r in range(H.size):
        if st[r]:
            strong.append((n, int(H[r]), int(P[r]), [int(x) for x in U[r]]))
Root, PI, Aord, Suf, S0, minord = S.build(strong, 73)
DL = {D: S.delta(S0[D]) for D in S0}
rex = {}
for N in range(3, 225, 2):
    rex[N] = sum(1 for D in range(1, N // 3 + 1, 2) if D in DL and S.is_slope(DL[D], (N - D) // 2))

print()
print("== r_m(N) against the exact r(N), odd N in [3, 225) ==")
for key in sorted(res):
    rm = res[key]
    above = [N for N in rex if rm[(N - 1) // 2] > rex[N]]
    eq = sum(1 for N in rex if rm[(N - 1) // 2] == rex[N])
    zeros = [N for N in range(3, 129, 2) if rm[(N - 1) // 2] == 0]
    print("r_%d (%s): exceeds exact r at %d values; equals exact r at %d of %d; zeros below 129: %s" %
          (key[0], key[1], len(above), eq, len(rex), zeros))

print()
print("== blocks: min r_m (at N), median, share = sum r / sum lattice points ==")
for key in sorted(res):
    rm = res[key]
    print("r_%d (%s):" % key)
    for k in range(6, 15):
        idx = np.arange((1 << k) // 2, (1 << (k + 1)) // 2)      # N = 2 idx + 1 in [2^k, 2^(k+1))
        blk = rm[idx]
        Nmin = 2 * int(idx[int(np.argmin(blk))]) + 1
        print("   [2^%d, 2^%d): min %d (N = %d), median %s, share %.4f, mean of r/L %.4f, #N with r = 0: %d" %
              (k, k + 1, int(blk.min()), Nmin, float(np.median(blk)), float(blk.sum() / Ls[idx].sum()),
               float(np.mean(blk / Ls[idx])), int(np.sum(blk == 0))))

print()
print("== L1: odd N in [129, 2^14) with r_m(N) < N/64 ==")
for key in sorted(res):
    rm = res[key]
    Ns = np.arange(129, 1 << 14, 2)
    v = Ns[rm[(Ns - 1) // 2] * 64 < Ns]
    print("r_%d (%s): %d violations%s" % (key[0], key[1], v.size,
                                         "" if v.size == 0 else " (least %d, largest %d)" % (v.min(), v.max())))
print()
print("== r_m(N) >= 1 on odd N in [69, 32767] ==")
for key in sorted(res):
    rm = res[key]
    Ns = np.arange(69, NMAX, 2)
    z = Ns[rm[(Ns - 1) // 2] == 0]
    print("r_%d (%s): %d values with r_m = 0%s" % (key[0], key[1], z.size, "" if z.size == 0 else " (least %d)" % z.min()))

print()
print("== values at the N named in the table of (F) ==")
for N in (131, 259, 521, 1027, 2101, 4109, 8203, 16447):
    print("N = %5d: N/64 = %7.2f ; " % (N, N / 64) + " ; ".join("r_%d(%s) = %d" % (k[0], k[1], res[k][(N - 1) // 2])
                                                                for k in sorted(res)))
np.savez_compressed("line_counts.npz", **{"r%d%s" % k: v for k, v in res.items()})

if 10 in pairs:
    rm = res[(10, "a")]
    print()
    print("== deciles of 3D/N on [4096, 8192), order 10, reading (a) ==")
    keys = pairs[10]
    D = (keys >> 32).astype(np.int64); U = (keys & 0xFFFFFFFF).astype(np.int64)
    mark = np.zeros((NMAX // 2, NMAX // 6 + 2), dtype=bool)
    j = 0
    while True:
        N = D + 2 * (U << j)
        ok = N < NMAX
        if not ok.any():
            break
        mark[(N[ok] - 1) // 2, (D[ok] - 1) // 2] = True
        j += 1
    num = np.zeros((2, 10)); den = np.zeros((2, 10))
    for N in range(4097, 8192, 2):
        Ds = np.arange(1, N // 3 + 1, 2)
        Us = (N - Ds) // 2
        dec = np.minimum((3 * Ds * 10) // N, 9)
        par = Us & 1          # 1 = odd U
        src = mark[(N - 1) // 2, (Ds - 1) // 2]
        for p in (0, 1):
            sel = par == p
            den[p] += np.bincount(dec[sel], minlength=10)
            num[p] += np.bincount(dec[sel], weights=src[sel], minlength=10)
    print("odd U :", " ".join("%.2f" % x for x in num[1] / den[1]))
    print("even U:", " ".join("%.2f" % x for x in num[0] / den[0]))
