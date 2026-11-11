# Auditor numerics-prereg: (1) additivity checks of my censuses (cumulative = cumulative + windows);
# (2) the exact remainder of THM-5134 (B) at every octave 2^22 .. 2^28 (my census), to see what the registered
#     P4 comparison (2^28 against 2^26) looks like among the neighbouring comparisons;
# (3) the check of the table "Numerical support" of THM-5134 against my values;
# (4) first states of the types 9, 27, 81, 25, 125, 49.
import json
import math
import re
from aud_common import divisors, remainder, groups, psi_table
from analyze import load

MODS12 = [9, 27, 81, 25, 125, 49, 45, 75, 63, 225, 15, 21]


def add(a, b, sign=1):
    # a, b: (T, mods) ; returns a + sign * b
    Ta, ma = a
    Tb, mb = b
    T = [Ta[0] + sign * Tb[0], Ta[1] + sign * Tb[1]]
    out = {}
    for m in ma:
        cls = [[x + sign * y for x, y in zip(ma[m][0][sh], mb[m][0][sh])] for sh in range(2)]
        typ = [{d: ma[m][1][sh][d] + sign * mb[m][1][sh][d] for d in ma[m][1][sh]} for sh in range(2)]
        diag = [{d: ma[m][2][sh][d] + sign * mb[m][2][sh][d] for d in ma[m][2][sh]} for sh in range(2)]
        out[m] = (cls, typ, diag)
    return T, out


def get(fn):
    J, T, mods = load(fn)
    return T, mods


cum = {k: get("my_cum_%d.json" % k) for k in (22, 24, 26, 28)}
win = {(24, 23): get("my_win_24_23.json"), (25, 24): get("my_win_25_24.json"), (26, 25): get("my_win_26_25.json"),
       (27, 26): get("my_win_27_26.json"), (28, 27): get("my_win_28_27.json")}

print("(1) additivity of my censuses")
c26 = add(add(cum[24], win[(25, 24)]), win[(26, 25)])
print("    cum 2^24 + [2^24,2^25) + [2^25,2^26) == cum 2^26:", c26[0] == cum[26][0] and c26[1] == cum[26][1])
c28 = add(add(cum[26], win[(27, 26)]), win[(28, 27)])
print("    cum 2^26 + [2^26,2^27) + [2^27,2^28) == cum 2^28:", c28[0] == cum[28][0] and c28[1] == cum[28][1])
cum[23] = add(cum[24], win[(24, 23)], -1)
cum[25] = add(cum[24], win[(25, 24)])
cum[27] = add(cum[26], win[(27, 26)])
for k in (23, 25, 27):
    print("    cumulative 2^%d by additivity: triples %s" % (k, cum[k][0]))

print()
print("(2) exact remainder of THM-5134 (B), m max_c |#{N=c} - sum_{d | gcd(m, c - c_s)} (d/m) D^(d)| / #T, sheet 0 / sheet 1, every octave")
print("      m  " + "  ".join("2^%d              " % k for k in range(22, 29)))
rem = {}
for m in MODS12:
    cells = []
    for k in range(22, 29):
        T, mods = cum[k]
        r0 = remainder(m, *mods[m], 0, False)
        r1 = remainder(m, *mods[m], 1, False)
        rem[(m, k)] = (r0, r1)
        cells.append("%.5f / %.5f" % (r0, r1))
    print("   %4d  " % m + "  ".join(cells))

print()
print("    number of the 12 pairs (m in 9, 27, 81, 25, 125, 49; sheet) whose remainder is smaller at the later size:")
P4 = (9, 27, 81, 25, 125, 49)
for a, b in ((22, 24), (24, 26), (26, 28), (23, 25), (25, 27), (22, 23), (23, 24), (24, 25), (25, 26), (26, 27), (27, 28), (24, 28), (22, 28)):
    n = sum(1 for m in P4 for sh in range(2) if rem[(m, b)][sh] < rem[(m, a)][sh])
    which = ", ".join("%d/s%d" % (m, sh) for m in P4 for sh in range(2) if not rem[(m, b)][sh] < rem[(m, a)][sh])
    print("      2^%d -> 2^%d: %2d of 12 smaller%s   not smaller: %s" % (a, b, n, "   <== the registered P4 comparison" if (a, b) == (26, 28) else "", which))

print()
print("(3) table 'Numerical support' of THM-5134 against my census (3 decimals; noise 4 decimals, sqrt(m/#T))")
thm = open("C:/tmp/wt-s92/01-canon/theorems/THM-5134-f1-local-densities-at-every-odd-modulus-prime-power-factor-and-the-four-square-singular-series.md", encoding="utf-8").read()
bad = 0
n = 0
for ln in thm.splitlines():
    mo = re.match(r"^\| (\d+) \| `([0-9.]+) / ([0-9.]+)` \| `([0-9.]+) / ([0-9.]+)` \| `([0-9.]+) / ([0-9.]+)` \| `([0-9.]+) / ([0-9.]+)` \| `([0-9.]+)` \|", ln)
    if not mo:
        continue
    m = int(mo.group(1))
    vals = [mo.group(i) for i in range(2, 11)]
    for j, k in enumerate((22, 24, 26, 28)):
        for sh in range(2):
            mine = rem[(m, k)][sh]
            pr = vals[2 * j + sh]
            n += 1
            okk = abs(mine - float(pr)) <= 0.00050001
            if not okk:
                bad += 1
                print("    MISMATCH m=%d 2^%d sheet %d: theorem %s, census %.5f" % (m, k, sh, pr, mine))
    T = cum[28][0]
    nz0 = math.sqrt(m / T[0])
    nz1 = math.sqrt(m / T[1])
    okn = abs(nz1 - float(vals[8])) <= 0.000050001 or abs(nz0 - float(vals[8])) <= 0.000050001
    n += 1
    if not okn:
        bad += 1
    print("    m = %3d: 8 remainders checked; noise in theorem %s, sqrt(m/#T) sheet 0 %.5f, sheet 1 %.5f %s" % (m, vals[8], nz0, nz1, "" if okn else "MISMATCH"))
print("    table entries checked: %d, mismatches: %d" % (n, bad))
G0 = groups(9, *cum[28][1][9], 0, False)
G1 = groups(9, *cum[28][1][9], 1, False)
print("    class factors at 9, X = 2^28, valuation 0 / 1 / >=2: sheet 0 %.4f / %.4f / %.4f, sheet 1 %.4f / %.4f / %.4f (theorem: 0.850 / 1.225 / 1.451 and 0.847 / 1.229 / 1.462)" % (
    G0[1][0], G0[3][0], G0[9][0], G1[1][0], G1[3][0], G1[9][0]))
G0 = groups(25, *cum[28][1][25], 0, False)
G1 = groups(25, *cum[28][1][25], 1, False)
print("    class factors at 25, X = 2^28, valuation 1 / >=2: sheet 0 %.4f / %.4f, sheet 1 %.4f / %.4f (Remark 2 of the theorem: 'nearly the same factor (1.25)')" % (
    G0[5][0], G0[25][0], G1[5][0], G1[25][0]))

print()
print("(4) first states h of each type (gcd(g(h), m) = d), h < 2^22")
psi = psi_table(1 << 22)
first = {}
for h in range(1 << 22):
    g = math.gcd(2 * int(psi[h]) + 1, 2 * h + 1)
    if g == 1:
        continue
    for d in (3, 9, 27, 81, 5, 25, 125, 7, 49):
        if g % d == 0 and d not in first:
            first[d] = (h, 2 * int(psi[h]) + 1, 2 * h + 1)
    if len(first) == 9:
        break
for d in (3, 9, 27, 81, 5, 25, 125, 7, 49):
    print("    d = %3d: first h with d | g(h): %s" % (d, first.get(d)))
