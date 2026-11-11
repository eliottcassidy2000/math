# audit-5129, step 8: analysis of the independent order-10 census of step 7 (reading (a): sinks appended).
import numpy as np
import supply as S

NMAX = 1 << 15
codes = np.load("o10_codes_0.npy"); strong = np.load("o10_strong_0.npy"); H = np.load("o10_H_0.npy")
u, first = np.unique(codes, return_index=True)
print("kept candidates: %d ; isomorphism classes of order 10: %d (known 9733056) ; strong: %d (known 9355949)" %
      (codes.size, u.size, int(strong[first].sum())))
Hs = H[first]; ss = strong[first] == 1
print("f(10) = least H of a strong tournament of order 10: %d (classes: %d) ; max H at order 10: %d" %
      (int(Hs[ss].min()), int(np.sum(Hs[ss] == Hs[ss].min())), int(Hs.max())))
print("strong classes of order 10 with H <= 223: %d" % int(np.sum(Hs[ss] <= 223)))
mark = np.unpackbits(np.load("o10_mark_0.npy"), axis=1)[:, :NMAX // 6 + 2].astype(bool)
r10 = mark.sum(axis=1)
LC = np.load("line_counts.npz")
r9 = LC["r9a"]
print("r_10 >= r_9 at every odd N < 2^15: %s ; lines where r_10 > r_9: %d (least N = %d)" %
      (bool(np.all(r10 >= r9)), int(np.sum(r10 > r9)), 2 * int(np.argmax(r10 > r9)) + 1))
# the order-9 marks must be contained in the order-10 marks
p9 = np.load("pairs9.npy")
D = (p9 >> 32).astype(np.int64); U = (p9 & 0xFFFFFFFF).astype(np.int64)
miss = 0
j = 0
while True:
    N = D + 2 * (U << j)
    ok = N < NMAX
    if not ok.any():
        break
    miss += int(np.sum(~mark[(N[ok] - 1) // 2, (D[ok] - 1) // 2]))
    j += 1
print("order-9 sources absent from the order-10 marks:", miss)

Ns_all = np.arange(1, NMAX, 2)
L = (Ns_all // 3 + 1) // 2
print()
print("== blocks ==")
print("block | min r_10 (at N) | median | sum r / sum L | mean r/L | mean r/(N/6)  || same three shares for r_9")
for k in range(6, 15):
    idx = np.arange((1 << k) // 2, (1 << (k + 1)) // 2)
    blk = r10[idx]; b9 = r9[idx]
    Nb = Ns_all[idx]
    amin = int(np.argmin(blk))
    print("[2^%d, 2^%d) | %d (%d) | %s | %.4f | %.4f | %.4f || %.4f | %.4f | %.4f" %
          (k, k + 1, int(blk.min()), int(Nb[amin]), float(np.median(blk)),
           blk.sum() / L[idx].sum(), np.mean(blk / L[idx]), np.mean(blk / (Nb / 6)),
           b9.sum() / L[idx].sum(), np.mean(b9 / L[idx]), np.mean(b9 / (Nb / 6))))
    allmin = Nb[blk == blk.min()]
    if allmin.size > 1:
        print("      (minimum attained at N = %s)" % allmin.tolist())

print()
Ns = np.arange(129, 1 << 14, 2)
v = Ns[r10[(Ns - 1) // 2] * 64 < Ns]
print("L1: odd N in [129, 2^14) with r_10(N) < N/64: %d ; min of 64 r_10(N)/N on the range: %.4f at N = %d" %
      (v.size, float(np.min(64 * r10[(Ns - 1) // 2] / Ns)), int(Ns[np.argmin(64 * r10[(Ns - 1) // 2] / Ns)])))
Ns2 = np.arange(1 << 14, NMAX, 2)
v2 = Ns2[r10[(Ns2 - 1) // 2] * 64 < Ns2]
print("    (beyond the claim) odd N in [2^14, 2^15) with r_10(N) < N/64: %d" % v2.size)
Ns3 = np.arange(69, NMAX, 2)
print("r_10(N) >= 1 on all odd N in [69, 32767]: %s ; odd N < 69 with r_10 = 0: %s" %
      (bool(np.all(r10[(Ns3 - 1) // 2] >= 1)), [N for N in range(3, 69, 2) if r10[(N - 1) // 2] == 0]))
for N in (131, 259, 521, 1027, 2101, 4109, 8203, 16447):
    print("   r_10(%d) = %d, r_9 = %d" % (N, r10[(N - 1) // 2], r9[(N - 1) // 2]))

print()
print("== deciles of 3D/N on [4096, 8192) ==")
for (lo, hi) in ((4096, 8192), (2048, 4096)):
    num = np.zeros((2, 10)); den = np.zeros((2, 10))
    for N in range(lo + 1, hi, 2):
        Ds = np.arange(1, N // 3 + 1, 2)
        Us = (N - Ds) // 2
        dec = np.minimum((3 * Ds * 10) // N, 9)
        par = Us & 1
        src = mark[(N - 1) // 2, (Ds - 1) // 2]
        for p in (0, 1):
            sel = par == p
            den[p] += np.bincount(dec[sel], minlength=10)
            num[p] += np.bincount(dec[sel], weights=src[sel].astype(float), minlength=10)
    print("[%d, %d) odd U :" % (lo, hi), " ".join("%.2f" % x for x in num[1] / den[1]))
    print("[%d, %d) even U:" % (lo, hi), " ".join("%.2f" % x for x in num[0] / den[0]))

# exact r(N) for N < 225 against r_10
Z = np.load("classes_le8.npz")
st = []
for n in range(3, 9):
    Hn = Z["H%d" % n]; sg = Z["strong%d" % n]; P = Z["Pi2%d" % n]; Un = Z["U%d" % n]
    for r in range(Hn.size):
        if sg[r]:
            st.append((n, int(Hn[r]), int(P[r]), [int(x) for x in Un[r]]))
Root, PI, Aord, Suf, S0, minord = S.build(st, 73)
DL = {Dd: S.delta(S0[Dd]) for Dd in S0}
neq = 0
for N in range(3, 225, 2):
    ex = sum(1 for Dd in range(1, N // 3 + 1, 2) if Dd in DL and S.is_slope(DL[Dd], (N - Dd) // 2))
    if ex != r10[(N - 1) // 2]:
        neq += 1
print()
print("odd N in [3, 225) where r_10(N) differs from the exact r(N) of the product law:", neq)
np.save("r10a.npy", r10)
