# Auditor recompute-5128: machine comparison of the numbers quoted in THM-5128 / README (sections 1, 2, 5) with my
# recomputed values (census_28.npz from a2_census.py; rates from a1_rates.build_T; A_p partial sums from a5_Ap.block_sums;
# predictions from a6_window.predict).  A line is flagged DIFF when the quoted number does not round-trip at its own
# precision (half a unit of the last quoted digit, plus 1e-12).
import math
import numpy as np
from a0_defs import psi_table32
from a1_rates import rates
from a5_Ap import block_sums
from a6_window import predict

Z = np.load("census_28.npz")
all_p = [int(p) for p in Z["all_p"]]
small_p = [int(p) for p in Z["small_p"]]
tot, hist, dist = Z["tot"], Z["hist"], Z["dist"]
Dsm, Dlg, Dboth, Ddiag = Z["Dsm"], Z["Dlg"], Z["Dboth"], Z["Ddiag"]
ndiff = 0
nok = 0


def chk(label, claimed, mine, digits):
    global ndiff, nok
    tol = 0.5 * 10 ** (-digits) + 1e-12
    ok = abs(claimed - mine) <= tol
    if ok:
        nok += 1
    else:
        ndiff += 1
    print(f"  {'OK  ' if ok else 'DIFF'} {label}: quoted {claimed:.{digits}f}, recomputed {mine:.{digits + 2}f}")


print("== (B) rates max_k rho(T_k)/2 ==")
for m, c in {3: 0.8660, 5: 0.7418, 7: 0.9320, 11: 0.8587, 13: 0.8272, 17: 0.8289, 19: 0.8053, 23: 0.7931, 31: 0.9002,
             43: 0.8342, 73: 0.8516, 127: 0.8754, 29: 0.7894}.items():
    r = rates(m)
    chk(f"m = {m}", c, max(r[k][0] for k in range(1, m)), 4)
def rate_of(m):
    r = rates(m)
    return max(r[k][0] for k in range(1, m))


others = {m: rate_of(m) for m in [29, 37, 41, 47, 53, 59, 61, 67, 71, 79, 83, 89]}
lo = min(others.values())
hi = max(others.values())
print(f"  other primes up to 89: min {lo:.6f} (m = {min(others, key=others.get)}), max {hi:.6f} (m = {max(others, key=others.get)}); quoted range [0.776, 0.806]: "
      f"{'OK' if lo >= 0.776 and hi <= 0.8065 else 'DIFF'}")
if not (lo >= 0.776 and hi <= 0.8065):
    ndiff += 1
for p, c in {3: 0.073, 5: 0.116, 7: 0.042}.items():
    rate = rate_of(p)
    chk(f"theta_{p}/2 (with rho_p = the spectral rate)", c, 0.5 * math.log(1 / rate) / math.log(2 / rate ** 2), 3)


def cum(sheet, b, arr, j=None):
    return int(arr[sheet, :b + 1].sum()) if j is None else int(arr[sheet, :b + 1, j].sum())


def rem(sheet, b, p):
    """(p max_c |R_sym|, p max_c |R_small|, F_large p^2)"""
    js = small_p.index(p)
    j = all_p.index(p)
    T = cum(sheet, b, tot)
    h = hist[sheet, :b + 1, js, :p].sum(axis=0) / T
    cs = (p - 1) if sheet == 0 else 0
    Fs = cum(sheet, b, Dsm, j) / T
    f = (cum(sheet, b, Dsm, j) / 2 + cum(sheet, b, Dlg, j) / 2 + cum(sheet, b, Ddiag, j)) / T
    out = []
    for D in (2 * f - 1 / p ** 2, Fs):
        pred = np.full(p, (1 - D) / p)
        pred[cs] += D
        out.append(p * np.abs(h - pred).max())
    return out[0], out[1], cum(sheet, b, Dlg, j) / T * p * p


print("== (C) VERIFIED parenthetical: p * max_c |remainder| of the form with 2f - 1/p^2 ('symmetric') ==")
for b, table in [(26, {5: (0.0003, 0.0002), 17: (0.002, 0.002), 13: (0.004, 0.004), 3: (0.009, 0.007), 31: (0.011, 0.008), 7: (0.059, 0.044)}),
                 (28, {5: (0.00013, 0.00013), 17: (0.0013, 0.0013), 13: (0.002, 0.002), 3: (0.014, 0.008), 31: (0.009, 0.006), 7: (0.053, 0.049)})]:
    for p, (c0, c1) in table.items():
        for sheet, c in ((0, c0), (1, c1)):
            digits = len(str(c).split(".")[1])
            chk(f"X = 2^{b}, p = {p}, sheet {sheet}", c, rem(sheet, b, p)[0], digits)
chk("X = 2^26, p = 7: form with F_small ('one-sided'), sheet 0", 0.117, rem(0, 26, 7)[1], 3)
chk("X = 2^26, p = 7: F_large p^2, sheet 0", 1.47, rem(0, 26, 7)[2], 2)
print("   the form of (C2) as stated (F_small), p * max_c |remainder|, sheets 0 / 1:")
for b in (26, 28):
    print(f"     X = 2^{b}: " + "; ".join(f"p = {p}: {rem(0, b, p)[1]:.4f} / {rem(1, b, p)[1]:.4f}" for p in (5, 17, 13, 3, 31, 7)))

print("== (D4) census values (F_small - 1/p^2) #T / (2X) at X = 2^28 ==")
for p, (c0, c1) in {5: (0.0300, 0.0434), 17: (0.0138, 0.0197), 11: (0.0118, 0.0155), 13: (-0.0035, -0.0044), 31: (-0.0017, -0.0021), 7: (-0.005, -0.010)}.items():
    j = all_p.index(p)
    for sheet, c in ((0, c0), (1, c1)):
        T = cum(sheet, 28, tot)
        v = (cum(sheet, 28, Dsm, j) - T / p ** 2) / (2 * 2.0 ** 28)
        digits = 3 if p == 7 else 4
        chk(f"p = {p}, sheet {sheet}", c, v, digits)

print("== (D4) A_p, partial sums over the states below 2^26 ==")
psi = psi_table32(1 << 26)
primes = np.array([3, 5, 7, 11, 13, 17, 31, 43], dtype=np.int64)
SW, DV, CN = block_sums(psi, 26, primes)
quoted = {3: (0.0557, 0.0789, 4), 5: (0.02943, 0.04250, 5), 7: (0.006, 0.007, 3), 11: (0.01096, 0.01449, 5), 13: (-0.00377, -0.00488, 5),
          17: (0.01355, 0.01931, 5), 31: (-0.00194, -0.00261, 5), 43: (0.00235, 0.00323, 5)}
for j, p in enumerate(primes):
    p = int(p)
    for sheet in (0, 1):
        A = float(np.sum(DV[sheet, j] - SW[sheet] / p ** 2))
        chk(f"A_{p}^{sheet} (partial sum at 2^26)", quoted[p][sheet], A, quoted[p][2])

print("== (D1) counts: two-point slopes at the top of the range against kappa ==")
c = np.cumsum(tot, axis=1)
for sheet, k in ((0, 0.1615389173), (1, 0.2125139092)):
    est = (c[sheet, 28] / 2.0 ** 28 - c[sheet, 27] / 2.0 ** 27) / math.log(2)
    chk(f"kappa_{sheet} from T(2^28), T(2^27)", k, est, 5)
    est = (c[sheet, 24] / 2.0 ** 24 - c[sheet, 23] / 2.0 ** 23) / math.log(2)
    chk(f"kappa_{sheet} from T(2^24), T(2^23)", k, est, 4)

print("== (G) window [2^24, 2^27): E measured ==")
wb = [25, 26, 27]
quotedE = {5: (0.3130, 0.3280, 4), 11: (0.1690, 0.1707, 4), 17: (0.2107, 0.2296, 4), 19: (0.0510, 0.0495, 4), 29: (0.0051, 0.0020, 4),
           43: (0.0885, 0.0930, 4), 257: (0.1859, 0.2023, 4), 2731: (0.063, 0.075, 3), 65537: (0.144, 0.168, 3), 3: (0.539, 0.544, 3),
           7: (0.233, 0.230, 3), 13: (0.0441, 0.0446, 4), 31: (0.0041, -0.0038, 4)}
for p, (c0, c1, digits) in quotedE.items():
    j = all_p.index(p)
    for sheet, cq in ((0, c0), (1, c1)):
        T = int(tot[sheet, wb].sum())
        d = int(dist[sheet, wb, j].sum())
        chk(f"E measured, p = {p}, sheet {sheet}", cq, d / ((T - d) / (p - 1)) - 1, digits)
T0, T1 = int(tot[0, wb].sum()), int(tot[1, wb].sum())
nodd = (2 ** 27 - 2 ** 24) // 2
chk("mass per odd N, measured", 14.1405, (T0 + T1) / nodd, 4)

print("== (G) window [2^24, 2^27): E predicted by the rule W ==")
psiw = psi_table32((1 << 27) // 2 + 4)
pp = np.array([3, 5, 7, 11, 13, 17, 19, 29, 31, 43, 257, 2731, 65537], dtype=np.int64)
Wtot, Wdiv, Wtrk, nterms = predict(psiw, pp, float(2 ** 24), float(2 ** 27), 28, np.array([0, 1, 2, 3, 8, 16, 64, 128, 1024], dtype=np.int64))
chk("mass per odd N, predicted", 14.1417, float(Wtot.sum()) / nodd, 4)
quotedP = {5: (0.3120, 0.3284, 4), 11: (0.1705, 0.1749, 4), 17: (0.2081, 0.2288, 4), 19: (0.0511, 0.0497, 4), 29: (0.0032, 0.0007, 4),
           43: (0.0877, 0.0941, 4), 257: (0.1815, 0.2044, 4), 2731: (0.061, 0.069, 3), 65537: (0.183, 0.206, 3), 3: (0.515, 0.532, 3),
           7: (0.174, 0.175, 3), 13: (0.0459, 0.0448, 4), 31: (-0.0061, -0.0089, 4)}
for j, p in enumerate(pp):
    p = int(p)
    for sheet in (0, 1):
        f = Wdiv[sheet, j] / Wtot[sheet]
        D = 2 * f - 1 / p ** 2
        chk(f"E predicted, p = {p}, sheet {sheet}", quotedP[p][sheet], p * D / (1 - D), quotedP[p][2])
print(f"   per sheet: predicted {Wtot[0] / nodd:.4f}, {Wtot[1] / nodd:.4f}; measured {T0 / nodd:.4f}, {T1 / nodd:.4f}; predicted/measured - 1 = {100 * (Wtot[0] / T0 - 1):+.2f} %, {100 * (Wtot[1] / T1 - 1):+.2f} %")

print("== README section 1 table: generic limit p (1/p^2) / (1 - 1/p^2) ==")
for p, cq, digits in [(3, 0.375, 3), (5, 0.208, 3), (7, 0.146, 3), (11, 0.092, 3), (13, 0.077, 3), (17, 0.059, 3), (19, 0.053, 3), (29, 0.034, 3),
                      (31, 0.032, 3), (43, 0.023, 3), (257, 0.004, 3), (65537, 0.00002, 5)]:
    chk(f"generic limit, p = {p}", cq, (1 / p) / (1 - 1 / p ** 2), digits)

print(f"\n{nok} quoted numbers reproduced at their own precision; {ndiff} flagged DIFF")
