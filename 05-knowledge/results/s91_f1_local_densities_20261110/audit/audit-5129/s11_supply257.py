# audit-5129, step 11: complete slope lists for every offset D <= 257 (all strong tournaments with H <= 257, from
# step 10), the table of (B'), and the EXACT r(N) for odd N <= 773 against the order-10 line counts.
import pickle
import numpy as np
import supply as S

rows = pickle.load(open("strong_le257.pkl", "rb"))
DMAX = 257
Root, PI, Aord, Suf, S0, minord = S.build(rows, DMAX)
DL = {D: S.delta(S0[D]) for D in S0}
na = [d for d in range(1, DMAX + 1, 2) if d not in Aord]
print("odd numbers <= 257 that are not path counts:", na)

print()
print("== (B') ==")
for s, D in ((7, 129), (8, 257)):
    dl = DL[D]
    lo, hi = D, (3 * D) // 2
    k = sum(1 for u in range(lo, hi + 1) if S.is_slope(dl, u))
    print("offset %d: %d odd parts; with j <= s: %d; least exponent of d = 9: %s; dyadic slope: %s; "
          "slopes in [h, 3h/2]: %d of %d (%.3f)" %
          (D, len(dl), sum(1 for d in dl if dl[d] <= s), dl.get(9), dl.get(1), k, hi - lo + 1, k / (hi - lo + 1)))
    with open("s11_slopes_%d.txt" % D, "w") as f:
        f.write(" ".join("%d:%d" % (d, dl[d]) for d in sorted(dl)) + "\n")
print("numbers of odd parts at 3, 5, 9, 17, 33, 65, 129, 257:", [len(DL[D]) for D in (3, 5, 9, 17, 33, 65, 129, 257)])
print("least exponents of d = 9 at levels 2..8:", [DL[2 ** s + 1].get(9) for s in range(2, 9)])
print("(129, 144) slope: %s ; (129, 72): %s ; (257, 288): %s ; (257, 144): %s" %
      (S.is_slope(DL[129], 144), S.is_slope(DL[129], 72), S.is_slope(DL[257], 288), S.is_slope(DL[257], 144)))
print("least order of (129,144): %s ; (129,288): %s ; (257,288): %s" %
      (S.minorder(S0[129], 144), S.minorder(S0[129], 288), S.minorder(S0[257], 288)))

print()
print("== HYP-9661 (A) shares at h = 17, 33, 65, 129, 257 ==")
for h in (17, 33, 65, 129, 257):
    def share(lo, hi):          # integers lo <= u <= hi
        return sum(1 for u in range(lo, hi + 1) if S.is_slope(DL[h], u)) / (hi - lo + 1)
    print("h = %d: [h, 9h/8]: %.2f ; [h, 2h): %.2f ; [2h, 4h): %.2f" %
          (h, share(h, (9 * h) // 8), share(h, 2 * h - 1), share(2 * h, 4 * h - 1)))

print()
print("== exact r(N), odd N in [3, 773] (offsets D <= N/3 <= 257) ==")
rex = {}
for N in range(3, 774, 2):
    rex[N] = sum(1 for D in range(1, N // 3 + 1, 2) if D in DL and S.is_slope(DL[D], (N - D) // 2))
print("r(N) = 0 exactly at:", [N for N in rex if rex[N] == 0])
r10 = np.load("r10a.npy")
LC = np.load("line_counts.npz")
r9 = LC["r9a"]
d10 = [N for N in rex if rex[N] != r10[(N - 1) // 2]]
d9 = [N for N in rex if rex[N] != r9[(N - 1) // 2]]
print("r_10 > exact r anywhere: %s ; r_10 < exact r at %d values (least %s) ; r_9 < exact r at %d values (least %s)" %
      (any(r10[(N - 1) // 2] > rex[N] for N in rex), len(d10), d10[:1], len(d9), d9[:1]))
for lo, hi in ((129, 225), (128, 256), (256, 512), (512, 774)):
    Ns = [N for N in rex if lo <= N < hi]
    v = [rex[N] for N in Ns]
    w = [int(r10[(N - 1) // 2]) for N in Ns]
    print("[%d, %d): exact r: min %d (N = %d), mean %.2f, median %s, mean share r/(N/6) %.4f ; r_10: min %d, mean %.2f ; "
          "min of 64 r/N: %.3f" % (lo, hi, min(v), Ns[int(np.argmin(v))], float(np.mean(v)), float(np.median(v)),
                                  float(np.mean([6 * rex[N] / N for N in Ns])), min(w), float(np.mean(w)),
                                  min(64 * rex[N] / N for N in Ns)))
# total size of the supply
print("total number of (offset, odd part) pairs, D <= 257:", sum(len(DL[D]) for D in DL))
# dyadic slopes
dy = [D for D in sorted(DL) if 1 in DL[D]]
print("offsets <= 257 with a dyadic slope: %d of %d attained ; <= 73: %d" % (len(dy), len(DL), sum(1 for D in dy if D <= 73)))
