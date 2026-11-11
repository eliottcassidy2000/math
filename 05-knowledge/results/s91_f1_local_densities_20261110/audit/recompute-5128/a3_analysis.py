# Auditor recompute-5128, tasks 2, 3, 6 (analysis of the exact census in census_B.npz produced by a2_census.py).
import sys
import numpy as np
import mpmath as mp

B = int(sys.argv[1]) if len(sys.argv) > 1 else 28
Z = np.load(f"census_{B}.npz")
small_p = [int(p) for p in Z["small_p"]]
all_p = [int(p) for p in Z["all_p"]]
tracked = [int(h) for h in Z["tracked"]]
tot, tots, diag, hist, dist = Z["tot"], Z["tots"], Z["diag"], Z["hist"], Z["dist"]
Dsm, Dlg, Dboth, Ddiag, trk = Z["Dsm"], Z["Dlg"], Z["Dboth"], Z["Ddiag"], Z["trk"]
J = {p: j for j, p in enumerate(all_p)}

mp.mp.dps = 30
l = [-mp.polylog(2, -mp.mpf(2) ** (-j)) for j in range(0, 40)]
kappa0 = l[0] - 2 * l[1] + l[2]
kappa1 = l[1] - l[2]
print("== constants (mpmath, 30 digits) ==")
print(f"kappa_0 = l0 - 2 l1 + l2 = {mp.nstr(kappa0, 12)}   (claimed 0.1615389173)")
print(f"kappa_1 = l1 - l2        = {mp.nstr(kappa1, 12)}   (claimed 0.2125139092)")
print(f"closed forms: pi^2/12 + 2 Li2(-1/2) - Li2(-1/4) = {mp.nstr(mp.pi ** 2 / 12 + 2 * mp.polylog(2, -0.5) - mp.polylog(2, -0.25), 12)};  Li2(-1/4) - Li2(-1/2) = {mp.nstr(mp.polylog(2, -0.25) - mp.polylog(2, -0.5), 12)}")
gbar = 2 * (kappa0 + kappa1)
print(f"g_bar = 2 (kappa_0 + kappa_1) = {mp.nstr(gbar, 12)}; pi^2/6 + 2 Li2(-1/2) = {mp.nstr(mp.pi ** 2 / 6 + 2 * mp.polylog(2, -0.5), 12)}   (claimed 0.7481056530)")
print(f"w_0 = {mp.nstr(kappa0 / (kappa0 + kappa1), 10)}, w_1 = {mp.nstr(kappa1 / (kappa0 + kappa1), 10)}   (claimed 0.43186124, 0.56813876)")
print(f"word pairs: sum (s+1)(l_s - 2 l_(s+1) + l_(s+2)), s < 38 = {mp.nstr(sum((s + 1) * (l[s] - 2 * l[s + 1] + l[s + 2]) for s in range(38)), 12)}; pi^2/12 = {mp.nstr(mp.pi ** 2 / 12, 12)}")
print(f"word-pair shares: kappa_0/(pi^2/12) = {mp.nstr(kappa0 / (mp.pi ** 2 / 12), 10)}, 1 - that = {mp.nstr(1 - kappa0 / (mp.pi ** 2 / 12), 10)}   (claimed 0.19640777, 0.80359223)")
print(f"8/pi^2 = {mp.nstr(8 / mp.pi ** 2, 8)}")
k0, k1 = float(kappa0), float(kappa1)

print()
print("== task 2: ordered triples with N <= X = 2^b (N is odd, so N <= 2^b <=> N < 2^b) ==")
c = np.cumsum(tot, axis=1)  # c[sheet, b] = #{N < 2^b}
print(" b      T_0(X)        T_1(X)     T_0/(X ln X)  T_1/(X ln X)   (T_0 - k0 X ln X)/X  (T_1 - k1 X ln X)/X   share sheet 0")
for b in range(16, B + 1):
    X = 2.0 ** b
    print(f"{b:2d} {c[0, b]:13d} {c[1, b]:13d}   {c[0, b] / (X * np.log(X)):.6f}      {c[1, b] / (X * np.log(X)):.6f}       {(c[0, b] - k0 * X * np.log(X)) / X:+.5f}           {(c[1, b] - k1 * X * np.log(X)) / X:+.5f}        {c[0, b] / (c[0, b] + c[1, b]):.5f}")
print("two-point estimates of kappa from consecutive X: (T(2X)/(2X) - T(X)/X)/ln 2")
for b in range(16, B):
    e0 = (c[0, b + 1] / 2.0 ** (b + 1) - c[0, b] / 2.0 ** b) / np.log(2)
    e1 = (c[1, b + 1] / 2.0 ** (b + 1) - c[1, b] / 2.0 ** b) / np.log(2)
    print(f"  b = {b:2d} -> {b + 1:2d}: kappa_0 est = {e0:.6f} ({e0 - k0:+.6f})   kappa_1 est = {e1:.6f} ({e1 - k1:+.6f})")
for (lo, hi) in [(16, 24), (20, B), (16, B)]:
    bs = np.arange(lo, hi + 1)
    Xs = 2.0 ** bs
    A = np.vstack([np.log(Xs), np.ones(len(bs))]).T
    f0 = np.linalg.lstsq(A, c[0, lo:hi + 1] / Xs, rcond=None)[0]
    f1 = np.linalg.lstsq(A, c[1, lo:hi + 1] / Xs, rcond=None)[0]
    print(f"least squares T/X = kappa ln X + c on X = 2^{lo}..2^{hi}: sheet 0 kappa = {f0[0]:.6f} (c = {f0[1]:+.4f}), sheet 1 kappa = {f1[0]:.6f} (c = {f1[1]:+.4f});  exact {k0:.6f}, {k1:.6f}")
print("per padding s (X = 2^b, b = B): count, count/(X ln X), two-point slope, claimed constant l_s - 2 l_(s+1) + l_(s+2)")
cs = np.cumsum(tots, axis=1)
for s in range(0, 8):
    X = 2.0 ** B
    const = float(l[s] - 2 * l[s + 1] + l[s + 2])
    slope = (cs[s, B] / 2.0 ** B - cs[s, B - 1] / 2.0 ** (B - 1)) / np.log(2)
    print(f"  s = {s}: {cs[s, B]:12d}  {cs[s, B] / (X * np.log(X)):.6f}  two-point {slope:.6f}   claimed {const:.6f}")
print(f"diagonal triples (h1 = h2) with N < 2^{B}: sheet 0 {diag[0].sum()}, sheet 1 {diag[1].sum()}")


def stats(sel_bins, sheet, p):
    """class counts and divisible-state counts over a set of bins"""
    j = J[p]
    T = int(tot[sheet, sel_bins].sum())
    cs_ = (p - 1) if sheet == 0 else 0
    if p in small_p:
        h = hist[sheet][sel_bins][:, small_p.index(p), :p].sum(axis=0)
    else:
        h = None
    d = int(dist[sheet, sel_bins, j].sum())
    dsm = int(Dsm[sheet, sel_bins, j].sum())
    dlg = int(Dlg[sheet, sel_bins, j].sum())
    dbo = int(Dboth[sheet, sel_bins, j].sum())
    ddg = int(Ddiag[sheet, sel_bins, j].sum())
    return dict(T=T, cs=cs_, h=h, d=d, dsm=dsm, dlg=dlg, dbo=dbo, ddg=ddg)


def report_classes(label, sel_bins, primes):
    print(f"-- {label} --")
    for p in primes:
        for sheet in (0, 1):
            st = stats(sel_bins, sheet, p)
            T = st["T"]
            h = st["h"]
            shares = h / T
            fav = int(np.argmax(h))
            cs_ = st["cs"]
            others = (T - h[cs_]) / (p - 1)
            E = h[cs_] / others - 1
            Fs = st["dsm"] / T
            f = (st["dsm"] / 2 + st["dlg"] / 2 + st["ddg"]) / T
            pred_sym = np.full(p, 1 / p - Fs / p)
            pred_sym[cs_] += Fs
            Done = 2 * f - 1 / p ** 2
            pred_one = np.full(p, 1 / p - Done / p)
            pred_one[cs_] += Done
            rem_sym = (shares - pred_sym) * p
            rem_one = (shares - pred_one) * p
            E_sym = p * Fs / (1 - Fs)
            E_one = p * Done / (1 - Done)
            E_gen = p * (1 / p ** 2) / (1 - 1 / p ** 2)
            nonfav = np.delete(shares, cs_) * p
            print(f" p = {p:2d} sheet {sheet}: T = {T}; favoured class = {fav} (c_s = {cs_}); p*share(c_s) = {shares[cs_] * p:.5f}; other classes p*share in [{nonfav.min():.5f}, {nonfav.max():.5f}]")
            print(f"        excess E (c_s over mean of others) = {E:+.5f}; first-order: symmetric (F_small = {Fs:.6f}) {E_sym:+.5f}, one-sided (2f - 1/p^2 = {Done:.6f}) {E_one:+.5f}; generic limit {E_gen:+.5f}")
            print(f"        p*remainder at c_s: symmetric {rem_sym[cs_]:+.5f}, one-sided {rem_one[cs_]:+.5f}; max_c |p*remainder|: symmetric {np.abs(rem_sym).max():.5f}, one-sided {np.abs(rem_one).max():.5f}")
            print(f"        F_small*p^2 = {Fs * p * p:.4f}, F_large*p^2 = {st['dlg'] / T * p * p:.4f}, f*p^2 = {f * p * p:.4f}")


print()
print("== task 3: class distribution of N mod p ==")
report_classes("N <= 2^22 (cumulative)", list(range(0, 23)), [3, 5, 7, 13])
report_classes("window [2^21, 2^23)", [22, 23], [3, 5, 7, 13])

print()
print("== (C) check: p * remainder of the first-order law, cumulative X = 2^b (symmetric form with F_small; one-sided form with 2f - 1/p^2) ==")
print("   columns: p, b, sheet0: sym@c_s  sym max_c  one@c_s  one max_c | sheet1: the same")
for p in [3, 5, 7, 11, 13, 17, 31]:
    for b in [18, 20, 22, 24, 26, 28]:
        if b > B:
            continue
        row = f" p = {p:2d} X = 2^{b}:"
        for sheet in (0, 1):
            st = stats(list(range(0, b + 1)), sheet, p)
            T = st["T"]
            h = st["h"]
            shares = h / T
            cs_ = st["cs"]
            Fs = st["dsm"] / T
            f = (st["dsm"] / 2 + st["dlg"] / 2 + st["ddg"]) / T
            Done = 2 * f - 1 / p ** 2
            ps = np.full(p, 1 / p - Fs / p)
            ps[cs_] += Fs
            po = np.full(p, 1 / p - Done / p)
            po[cs_] += Done
            rs = (shares - ps) * p
            ro = (shares - po) * p
            row += f"  s{sheet}: {rs[cs_]:+.5f} {np.abs(rs).max():.5f} {ro[cs_]:+.5f} {np.abs(ro).max():.5f}  (F_large p^2 = {st['dlg'] / T * p * p:.3f})"
        print(row)

print()
print("== (D4) check: (F_small - 1/p^2) #T / (2X) at cumulative X = 2^b, sheets 0 / 1 ==")
for p in [3, 5, 7, 11, 13, 17, 31, 43]:
    row = f" p = {p:2d}:"
    for b in [22, 24, 26, 28]:
        if b > B:
            continue
        vals = []
        for sheet in (0, 1):
            st = stats(list(range(0, b + 1)), sheet, p)
            T = st["T"]
            vals.append((st["dsm"] - T / p ** 2) / (2 * 2.0 ** b))
        row += f"   2^{b}: {vals[0]:+.5f}, {vals[1]:+.5f}"
    print(row)

if B >= 27:
    print()
    print("== task 6: window [2^24, 2^27) ==")
    wb = [25, 26, 27]
    T0 = int(tot[0, wb].sum())
    T1 = int(tot[1, wb].sum())
    nodd = (2 ** 27 - 2 ** 24) // 2
    print(f"triples: sheet 0 {T0}, sheet 1 {T1}, total {T0 + T1}; odd N in window {nodd}; mass per odd N = {(T0 + T1) / nodd:.5f}  (claimed 14.1405; predicted 14.1417)")
    print(f"sub-window [2^24, 2^25): sheet 0 {int(tot[0, 25])}, sheet 1 {int(tot[1, 25])}, mass per odd N = {(int(tot[0, 25]) + int(tot[1, 25])) / 2 ** 23:.5f}")
    print(" p      E sheet0   E sheet1    | E_sym(F_small) s0  s1 | E_one(2f-1/p^2) s0  s1 | E_nd measured s0 s1 | stderr(E) s0 s1")
    for p in all_p:
        row = f" {p:6d}"
        Es, Esym, Eone, End, se = [], [], [], [], []
        for sheet in (0, 1):
            st = stats(wb, sheet, p)
            T = st["T"]
            d = st["d"]
            E = d / ((T - d) / (p - 1)) - 1
            Fs = st["dsm"] / T
            f = (st["dsm"] / 2 + st["dlg"] / 2 + st["ddg"]) / T
            Done = 2 * f - 1 / p ** 2
            anydiv = st["dsm"] + st["dlg"] - st["dbo"] + st["ddg"]
            Tnd = T - anydiv
            dnd = d - anydiv
            Es.append(E)
            Esym.append(p * Fs / (1 - Fs))
            Eone.append(p * Done / (1 - Done))
            End.append(dnd / ((Tnd - dnd) / (p - 1)) - 1)
            se.append(np.sqrt(p / T))  # Poisson: relative sd of the distinguished count ~ sqrt(p/T)
        row += f"  {Es[0]:+.4f}    {Es[1]:+.4f}    |  {Esym[0]:+.4f} {Esym[1]:+.4f}  |  {Eone[0]:+.4f} {Eone[1]:+.4f}  |  {End[0]:+.4f} {End[1]:+.4f}  |  {se[0]:.4f} {se[1]:.4f}"
        print(row)
    print(" sub-window [2^24, 2^25): E measured sheets 0 / 1")
    for p in [3, 5, 7, 11, 13, 17]:
        vals = []
        for sheet in (0, 1):
            st = stats([25], sheet, p)
            vals.append(st["d"] / ((st["T"] - st["d"]) / (p - 1)) - 1)
        print(f"   p = {p:2d}: {vals[0]:+.4f}  {vals[1]:+.4f}")
    print(" tracked first states: ordered triples with first state h in the window (all s), and by sheet")
    for t, h in enumerate(tracked):
        w0 = int(trk[t, 0, wb].sum())
        w1 = int(trk[t, 1:, :][:, wb].sum())
        print(f"   h1 = {h:5d}: sheet 0 {w0:10d}  sheet 1 {w1:10d}  total {w0 + w1:10d}")
