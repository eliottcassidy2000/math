#!/usr/bin/env python3
"""audit_tables.py -- independent audit: sections 4, 5 and 6.3 of the note from my own class census (audit_census.c).

Input: out/census_n.txt (n = 3..9) and out/census_10_part*.txt, lines
    sc sym strong Hmod8 p3 p5 p7 p9  classes  labelled_weight
Optional: out/cube_census_10_part*.out (exact labelled H mod 8 counts at order 11 from audit_cube.c).
"""
import glob, os, sys
from fractions import Fraction as F
from math import comb, factorial

here = os.path.dirname(os.path.abspath(__file__))
out = os.path.join(here, "out")

def dfact(k):
    r = 1
    while k > 1:
        r *= k; k -= 2
    return r

def load(n):
    files = [os.path.join(out, f"census_{n}.txt")] if n < 10 else sorted(glob.glob(os.path.join(out, "census_10_part*.txt")))
    rows = {}
    for fn in files:
        if not os.path.exists(fn):
            return None
        for line in open(fn):
            p = line.split()
            if len(p) != 10:
                continue
            key = tuple(int(x) for x in p[:8])
            c, w = int(p[8]), int(p[9])
            a = rows.setdefault(key, [0, 0]); a[0] += c; a[1] += w
    return rows

data = {}
for n in range(3, 11):
    r = load(n)
    if r:
        data[n] = r

def chi(h8):
    return 1 if h8 % 4 == 1 else -1

print("Checks of totals: classes and labelled weight per order")
A000568 = {3: 2, 4: 4, 5: 12, 6: 56, 7: 456, 8: 6880, 9: 191536, 10: 9733056}
for n, rows in data.items():
    cl = sum(v[0] for v in rows.values()); w = sum(v[1] for v in rows.values())
    print(f"  n={n}: classes {cl} (A000568: {A000568[n]}) {'OK' if cl == A000568[n] else 'MISMATCH'} ; weight {w} = 2^{comb(n,2)}: {'OK' if w == 2**comb(n,2) else 'MISMATCH'}")

print("\nSection 4 table: SC_n, bound, SC_n/2^C(n,2), D_n, share")
Dn = {}
for n, rows in data.items():
    tot = 2 ** comb(n, 2)
    SC = sum(v[1] for k, v in rows.items() if k[0] == 1)
    SCs = sum(v[1] * chi(k[3]) for k, v in rows.items() if k[0] == 1)
    D = sum(v[1] * chi(k[3]) for k, v in rows.items())
    Dn[n] = F(D, tot)
    bound = dfact(n - 1) * 2 ** (n * n // 4) if n % 2 == 0 else dfact(n) * 2 ** ((n * n - 1) // 4)
    low = 2 ** (n * n // 4) if n % 2 == 0 else 2 ** ((n * n - 1) // 4)
    share = F(SCs, tot) / Dn[n]
    print(f"  n={n}: SC_n = {SC} ; bound {bound} ; lower bound {low} ; within bounds {low <= SC <= bound} ; SC_n/2^C = {F(SC, tot)} ; signed sum over sc = {SCs} ; "
          f"D_n = {Dn[n]} ; D_n 4^n = {float(Dn[n] * 4**n):.1f} ; sc share of D_n = {float(share)*100:.1f}%")

print("\nSelf-converse classes by H mod 8 (1,3,5,7), and by H mod 4")
for n, rows in data.items():
    c8 = {r: sum(v[0] for k, v in rows.items() if k[0] == 1 and k[3] == r) for r in (1, 3, 5, 7)}
    print(f"  n={n}: sc classes {sum(c8.values())} ; mod 8: {c8} ; mod 4: =1: {c8[1]+c8[5]}, =3: {c8[3]+c8[7]}")

print("\nSection 5: three strata in units of 4^(-n)  (sc / sym non-sc / rigid non-sc)")
for n, rows in data.items():
    tot = 2 ** comb(n, 2); u = F(4 ** n, tot)
    sc = sum(v[1] * chi(k[3]) for k, v in rows.items() if k[0] == 1)
    sym = sum(v[1] * chi(k[3]) for k, v in rows.items() if k[0] == 0 and k[1] == 1)
    rig = sum(v[1] * chi(k[3]) for k, v in rows.items() if k[0] == 0 and k[1] == 0)
    print(f"  n={n}: D 4^n = {float(Dn[n]*4**n):.1f} ; sc {float(sc*u):.1f} ; sym {float(sym*u):.1f} ; rigid {float(rig*u):.1f}   (exact signed sums {sc}, {sym}, {rig})")

print("\nClass race: classes with H = 1 minus classes with H = 3 (mod 4); decomposition for even n")
for n, rows in data.items():
    c1 = sum(v[0] for k, v in rows.items() if chi(k[3]) == 1); c3 = sum(v[0] for k, v in rows.items() if chi(k[3]) == -1)
    scc = sum(v[0] for k, v in rows.items() if k[0] == 1)
    ns1 = sum(v[0] for k, v in rows.items() if k[0] == 0 and chi(k[3]) == 1); ns3 = sum(v[0] for k, v in rows.items() if k[0] == 0 and chi(k[3]) == -1)
    print(f"  n={n}: =1: {c1} =3: {c3} lead {c1-c3} ; sc classes {scc} ; non-sc classes =1: {ns1} =3: {ns3} (lead {ns1-ns3}; converse pairs lead {(ns1-ns3)//2})")

print("\nStrong race: signed counts of strong labelled tournaments, two routes")
t = {0: 1, 1: 1, 2: 2}
for n in data:
    t[n] = int(Dn[n] * 2 ** comb(n, 2))
# D_11 from my own exact vertex-cube census of the bases of order 10, if present
cube10 = sorted(glob.glob(os.path.join(out, "cube_census_10_part*.out")))
h8_11 = None
if cube10:
    h8_11 = [0, 0, 0, 0]; num = 0; ok = True
    for fn in cube10:
        txt = open(fn).read().split("\n")
        try:
            l2 = txt[1].split(); h = [int(l2[l2.index(str(r) + ":") + 1]) for r in (1, 3, 5, 7)]
            num += int(txt[2].split()[-1])
            h8_11 = [a + b for a, b in zip(h8_11, h)]
        except Exception as e:
            ok = False
    if ok and sum(h8_11) == 2 ** 55:
        t[11] = num
        print(f"  [own exact D_11 = {F(num, 2**55)} = {num}/2^55 ; quoted 260546929/2^39: {'OK' if F(num, 2**55) == F(260546929, 2**39) else 'MISMATCH'}]")
    else:
        h8_11 = None
if 11 not in t:
    t[11] = 260546929 * 2 ** 16
    print("  [D_11 taken from the note: 260546929/2^39]")
s = {}
for n in range(1, 12):
    s[n] = t[n] - sum(comb(n, k) * s[k] * t[n - k] for k in range(1, n))
# unsigned strong counts by the same inversion
tu = {n: 2 ** comb(n, 2) for n in range(0, 12)}
su = {}
for n in range(1, 12):
    su[n] = tu[n] - sum(comb(n, k) * su[k] * tu[n - k] for k in range(1, n))
for n in range(3, 12):
    line = f"  n={n}: by inversion: signed {s[n]}, strong labelled {su[n]}, strong-only lead {s[n]/su[n]:.6f}"
    if n in data:
        rows = data[n]
        sflag = sum(v[1] * chi(k[3]) for k, v in rows.items() if k[2] == 1)
        nflag = sum(v[1] for k, v in rows.items() if k[2] == 1)
        line += f" ; census flags: signed {sflag}, strong {nflag} : {'agree' if (sflag, nflag) == (s[n], su[n]) else 'DISAGREE'}"
        cs1 = sum(v[0] for k, v in rows.items() if k[2] == 1 and chi(k[3]) == 1); cs3 = sum(v[0] for k, v in rows.items() if k[2] == 1 and chi(k[3]) == -1)
        line += f" ; strong classes =1: {cs1}, =3: {cs3}"
    line += f" ; decomposable share of the signed sum {1 - s[n]/t[n]:.4f}"
    print(line)

print("\nRelative labelled masses of H = 1, 3, 5, 7 (mod 8)")
for n, rows in data.items():
    tot = 2 ** comb(n, 2)
    mm = [F(4 * sum(v[1] for k, v in rows.items() if k[3] == r), tot) for r in (1, 3, 5, 7)]
    print(f"  n={n}: " + "  ".join(f"{float(x):.5f}" for x in mm) + "    exact counts " + " ".join(str(sum(v[1] for k, v in rows.items() if k[3] == r)) for r in (1, 3, 5, 7)))
if h8_11:
    print("  n=11: " + "  ".join(f"{4*x/2**55:.5f}" for x in h8_11) + "    exact counts " + " ".join(map(str, h8_11)) + "   (own vertex-cube census of the bases of order 10)")
    print("        fractions: " + "  ".join(f"{x/2**55:.7f}" for x in h8_11))

print("\nSection 6.3: E_n[K] for every set K of odd lengths (exact)")
names = [(1, 0, 0, 0), (0, 1, 0, 0), (1, 1, 0, 0), (0, 0, 1, 0), (1, 0, 1, 0), (0, 1, 1, 0), (1, 1, 1, 0), (0, 0, 0, 1), (1, 0, 0, 1), (0, 1, 0, 1), (1, 1, 0, 1),
         (0, 0, 1, 1), (1, 0, 1, 1), (0, 1, 1, 1), (1, 1, 1, 1)]
def nm(K):
    return "{" + ",".join(str(l) for l, b in zip((3, 5, 7, 9), K) if b) + "}"
EK = {}
for n, rows in data.items():
    tot = 2 ** comb(n, 2)
    line = []
    for K in names:
        if any(b and l > n for l, b in zip((3, 5, 7, 9), K)):
            continue
        sg = sum(v[1] * (-1) ** sum(b * p for b, p in zip(K, k[4:8])) for k, v in rows.items())
        EK[(n, K)] = F(sg, tot)
        line.append(f"{nm(K)}: {F(sg, tot)}")
    print(f"  n={n}: " + " ; ".join(line))
    allK = tuple(1 if l <= n else 0 for l in (3, 5, 7, 9))
    print(f"        all odd lengths = D_n: {EK[(n, allK)] == Dn[n]}")

def chk(label, cond):
    print(f"  {'OK      ' if cond else 'MISMATCH'} {label}")

print("\nComparison with the numbers quoted in the note")
q = {4: (48, F(3, 4), F(1, 2)), 6: (7040, F(55, 256), F(1, 4)), 8: (6737920, F(1645, 65536), F(1361, 32768)), 10: (31596429312, F(1928493, 2147483648), F(650471, 268435456))}
for n, (SCq, frq, Dq) in q.items():
    if n in data:
        SC = sum(v[1] for k, v in data[n].items() if k[0] == 1)
        chk(f"n={n}: SC_n = {SCq}, fraction {frq}, D_n = {Dq}", SC == SCq and F(SC, 2 ** comb(n, 2)) == frq and Dn[n] == Dq)
for n, SCq, sg in ((3, 8, 4), (5, 704, 256), (7, 380416, 37824), (9, 958578688, 20782592)):
    if n in data:
        SC = sum(v[1] for k, v in data[n].items() if k[0] == 1); SCs = sum(v[1] * chi(k[3]) for k, v in data[n].items() if k[0] == 1)
        chk(f"n={n}: SC_n = {SCq}, signed sum {sg}", SC == SCq and SCs == sg)
Dq = {3: F(1, 2), 4: F(1, 2), 5: F(13, 32), 6: F(1, 4), 7: F(853, 8192), 8: F(1361, 32768), 9: F(63103, 2 ** 23), 10: F(650471, 2 ** 28)}
for n in data:
    chk(f"D_{n} = {Dq[n]}", Dn[n] == Dq[n])
lead = {3: 0, 4: 0, 5: 0, 6: 8, 7: 36, 8: 264, 9: 1268, 10: 22600}
for n, rows in data.items():
    c1 = sum(v[0] for k, v in rows.items() if chi(k[3]) == 1); c3 = sum(v[0] for k, v in rows.items() if chi(k[3]) == -1)
    chk(f"class race at n={n}: lead {lead[n]}", c1 - c3 == lead[n])
sq = {3: -2, 4: 24, 5: 176, 6: 4080, 7: 123952, 8: 8135680, 9: 332488064, 10: 75937800960, 11: 15259291222528}
for n in range(3, 12):
    chk(f"signed strong count at n={n}: {sq[n]}", s[n] == sq[n])
tab = {
    (5, (1, 0, 0, 0)): F(1, 4), (5, (0, 1, 0, 0)): F(7, 32), (5, (1, 1, 0, 0)): F(13, 32),
    (6, (1, 0, 0, 0)): F(1, 4), (6, (0, 1, 0, 0)): F(1, 16), (6, (1, 1, 0, 0)): F(1, 4),
    (7, (1, 0, 0, 0)): F(1, 8), (7, (0, 1, 0, 0)): F(11, 256), (7, (1, 1, 0, 0)): F(67, 256), (7, (0, 0, 1, 0)): F(811, 8192), (7, (1, 0, 1, 0)): F(1553, 8192),
    (7, (0, 1, 1, 0)): F(461, 8192), (7, (1, 1, 1, 0)): F(853, 8192),
    (8, (1, 0, 0, 0)): F(1, 8), (8, (0, 1, 0, 0)): F(11, 256), (8, (1, 1, 0, 0)): F(67, 256), (8, (0, 0, 1, 0)): F(1277, 32768), (8, (1, 0, 1, 0)): F(1361, 32768),
    (8, (0, 1, 1, 0)): F(437, 32768), (8, (1, 1, 1, 0)): F(1361, 32768),
    (9, (1, 0, 0, 0)): F(1, 16), (9, (0, 1, 0, 0)): F(23, 2048), (9, (1, 1, 0, 0)): F(263, 2048), (9, (0, 0, 1, 0)): F(16343, 2 ** 21), (9, (1, 0, 1, 0)): F(-1453, 2 ** 21),
    (9, (0, 1, 1, 0)): F(9019, 2 ** 18), (9, (1, 1, 1, 0)): F(35881, 2 ** 20), (9, (1, 1, 1, 1)): F(63103, 2 ** 23),
    (10, (1, 0, 0, 0)): F(1, 16), (10, (0, 1, 0, 0)): F(1, 256), (10, (1, 1, 0, 0)): F(1, 16), (10, (0, 0, 1, 0)): F(35881, 2 ** 24), (10, (1, 0, 1, 0)): F(9019, 2 ** 22),
    (10, (0, 1, 1, 0)): F(9019, 2 ** 18), (10, (1, 1, 1, 0)): F(35881, 2 ** 20), (10, (1, 1, 1, 1)): F(650471, 2 ** 28),
    (10, (0, 0, 0, 1)): F(1050641, 2 ** 28), (10, (1, 0, 1, 1)): F(533081, 2 ** 28),
}
for (n, K), v in sorted(tab.items()):
    if (n, K) in EK:
        chk(f"E_{n}[{nm(K)}] = {v}", EK[(n, K)] == v)
if 10 in data:
    chk("E_10[{7}] = E_10[{3,5}] E_10[{3,5,7}]", EK[(10, (0, 0, 1, 0))] == EK[(10, (1, 1, 0, 0))] * EK[(10, (1, 1, 1, 0))])
    chk("E_10[{3,7}] = E_10[{3,5}] E_10[{5,7}]", EK[(10, (1, 0, 1, 0))] == EK[(10, (1, 1, 0, 0))] * EK[(10, (0, 1, 1, 0))])
    # "and with no other sign involving c_7 or c_9": list every K, K' with E[K] E[K'] = E[K xor K'], K = {3,5}
    K35 = (1, 1, 0, 0); hits = []
    for K in names:
        Kx = tuple(a ^ b for a, b in zip(K, K35))
        if K == K35 or Kx == (0, 0, 0, 0):
            continue
        if EK[(10, K)] * EK[(10, K35)] == EK[(10, Kx)]:
            hits.append(nm(K) + " -> " + nm(Kx))
    print("        sets K with E_10[K xor {3,5}] = E_10[{3,5}] E_10[K]:", hits)
    chk("E[{5,7}], E[{3,5,7}] equal at n = 9 and n = 10", EK[(9, (0, 1, 1, 0))] == EK[(10, (0, 1, 1, 0))] and EK[(9, (1, 1, 1, 0))] == EK[(10, (1, 1, 1, 0))])
if 6 in data:
    print("        n = 6 analogue (c_7 = 0): E_6[{}]=1 vs E_6[{3,5}]^2 =", EK[(6, (1, 1, 0, 0))] ** 2, "; E_6[{3}] =", EK[(6, (1, 0, 0, 0))], "vs E_6[{3,5}] E_6[{5}] =", EK[(6, (1, 1, 0, 0))] * EK[(6, (0, 1, 0, 0))])
