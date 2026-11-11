# Auditor recompute-5128: the excess at large primes bin by bin (from census_28.npz), including the bin [2^27, 2^28)
# that lies outside the session's pre-registered window: E (distinguished class over the mean of the others) and
# E_nd (after removal of the triples with a p-divisible state), per sheet.
import math
import numpy as np

Z = np.load("census_28.npz")
all_p = [int(p) for p in Z["all_p"]]
tot = Z["tot"]


def EE(p, bins, sheet):
    j = all_p.index(p)
    T = int(tot[sheet, bins].sum())
    d = int(Z["dist"][sheet, bins, j].sum())
    dsm = int(Z["Dsm"][sheet, bins, j].sum())
    dlg = int(Z["Dlg"][sheet, bins, j].sum())
    dbo = int(Z["Dboth"][sheet, bins, j].sum())
    ddg = int(Z["Ddiag"][sheet, bins, j].sum())
    anyd = dsm + dlg - dbo + ddg
    others = (T - d) / (p - 1)
    return d / others - 1, (d - anyd) / others - 1, T, anyd, d


# index of dispersion measured in a7_misses.out: about 2.78 (sheet 0), 3.08 (sheet 1)
FANO = (2.78, 3.08)
for p in [257, 2731, 65537, 65521, 65539, 8191]:
    print(f"== p = {p} ==")
    print("   bins            sheet   T            E        E_nd      se (sqrt(Fano p / T))   z(E_nd)   triples with a p-divisible state / distinguished class count")
    for (label, bins) in [("[2^24,2^25)", [25]), ("[2^25,2^26)", [26]), ("[2^26,2^27)", [27]), ("[2^24,2^27)", [25, 26, 27]), ("[2^27,2^28)", [28]), ("[2^20,2^24)", [21, 22, 23, 24])]:
        for sheet in (0, 1):
            E, End, T, anyd, d = EE(p, bins, sheet)
            se = math.sqrt(FANO[sheet] * p / T)
            print(f"   {label:14s}  {sheet}   {T:11d}   {E:+.4f}   {End:+.4f}   {se:.4f}                  {End / se:+5.2f}     {anyd} / {d}")
