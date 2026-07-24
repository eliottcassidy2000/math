#!/usr/bin/env python3
"""
snippet_eq27_lrc_families_opus_S4.py   opus-2026-07-23-S4

Test the strongest in-repo structural lead: the snippet compares two log-quantities
A=1285/896, B=8847357/2974400, and the repo compares the two tight LRC(14) families
AP={1..13}, GW={1..11,13,24} via power-sum invariants (J=S2*S6/S4^2, CONSTANTS-INDEX).
The coefficient 2457=3*F2(13)=3*819 where F2(13)=819 is the GW "acceleration".

Do A,B,896,1285,2974400,8847357,389,5872957, or 2457/6592 arise from the power sums
S_k of these LRC families (and neighbors {1..12}, {1..11,24}, {1..12,26}, {1..11,13,36})?
"""
from fractions import Fraction as F
from sympy import factorint

families = {
    "AP13={1..13}": list(range(1, 14)),
    "GW={1..11,13,24}": list(range(1, 12)) + [13, 24],
    "n13floor={1..12}": list(range(1, 13)),
    "second2/25={1..11,24}": list(range(1, 12)) + [24],
    "K2={1..12,26}": list(range(1, 13)) + [26],
    "3/41={1..11,13,36}": list(range(1, 12)) + [13, 36],
    "covermin={1..12,182}": list(range(1, 13)) + [182],
}
targets = {896: "Sn-1(A)", 1285: "Sn(A)", 389: "x(A)", 2181: "den tA",
           2974400: "Sn-1(B)", 8847357: "Sn(B)", 5872957: "x(B)", 11821757: "den tB",
           2457: "coef num", 6592: "coef den", 819: "F2(13)"}

def Sk(v, k): return sum(x**k for x in v)

print("=== power sums S_k and simple invariants per family ===")
hits = []
for name, v in families.items():
    S = {k: Sk(v, k) for k in range(0, 9)}
    print(f"\n{name}:  S1={S[1]} S2={S[2]} S3={S[3]} S4={S[4]} S5={S[5]} S6={S[6]}")
    # J separator and other ratios
    invs = {
        "S2": S[2], "S3": S[3], "S4": S[4], "S6": S[6],
        "S2*S6": S[2]*S[6], "S4^2": S[4]**2,
        "S1*S3": S[1]*S[3], "S2^2": S[2]**2, "S3^2": S[3]**2,
        "S2*S4": S[2]*S[4], "prod": 1,
    }
    p = 1
    for x in v: p *= x
    invs["prod(v)"] = p
    # scan invariants and pairwise ratios/diffs against targets
    for iname, iv in invs.items():
        if iv in targets:
            hits.append(f"{name}: {iname}={iv} == {targets[iv]}")
    # pairwise differences of power sums
    for a in range(1, 9):
        for b in range(0, a):
            d = S[a] - S[b]
            if d in targets:
                hits.append(f"{name}: S{a}-S{b}={d} == {targets[d]}")
            s = S[a] + S[b]
            if s in targets:
                hits.append(f"{name}: S{a}+S{b}={s} == {targets[s]}")

print("\n=== direct hits (power-sum invariant == a target integer) ===")
for h in hits:
    print("  ", h)
if not hits:
    print("   NONE — targets are not simple power-sum invariants of these families")

# J separator invariant (from CONSTANTS-INDEX: J(AP)=98827/83167, J(GW)=3083942821/1978381441)
print("\n=== J = S2*S6/S4^2 separator (verify against CONSTANTS-INDEX) ===")
for name, v in families.items():
    S = {k: Sk(v, k) for k in range(0, 9)}
    J = F(S[2]*S[6], S[4]**2)
    print(f"  {name}: J=S2S6/S4^2 = {J} = {float(J):.6f}")

# Does log(B)/log(A) or the ratio structure match any family-pair log-separation?
import mpmath as mp
mp.mp.dps = 40
logA = mp.log(mp.mpf(1285)/896); logB = mp.log(mp.mpf(8847357)/2974400)
print("\n=== do A,B match ratios J(fam1)/J(fam2) or S_k(fam1)/S_k(fam2)? ===")
fam_items = list(families.items())
for i in range(len(fam_items)):
    for j in range(len(fam_items)):
        if i == j: continue
        vi, vj = fam_items[i][1], fam_items[j][1]
        for k in range(1, 7):
            r = F(Sk(vi, k), Sk(vj, k))
            if r == F(1285, 896) or r == F(8847357, 2974400):
                print(f"  MATCH: S{k}({fam_items[i][0]})/S{k}({fam_items[j][0]}) = {r}")

# targets' factorizations reminder
print("\n=== target factorizations (for eyeballing family structure) ===")
for t, lbl in targets.items():
    print(f"  {lbl:10s} {t} = {dict(factorint(t))}")
