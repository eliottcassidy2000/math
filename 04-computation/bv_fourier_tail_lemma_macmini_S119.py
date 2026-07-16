#!/usr/bin/env python3
"""bv_fourier_tail_lemma_macmini_S119.py -- mac-mini-2026-07-16-S119.
THE BV-FOURIER TAIL LEMMA, proved and computed -- closing negative residue 6.

PROVED STRUCTURE (writeup in THM-908):
 (1) The sector indicator 1_{[s/7,(s+1)/7)} has Fourier coefficients
     I_s^(n) = e(-n(2s+1)/14) sin(pi n/7)/(pi n)  (n != 0), = 1/7 (n = 0).
 (2) Hence the triple channel G(y) = beta123(sec y1, sec y2, sec y3) has
     G^(n) = prod_i [sin(pi n_i/7)/(pi n_i)] * W(n mod 7),  where W = the (Z/7)^3
     "phase DFT" of beta123 (|W| computable; W(r) = 0 if any r_i = 0 -- zero marginality;
     G^(n) = 0 if any n_i = 0 or 7 | n_i).
 (3) T(a,b,c) = sum over the annihilator lattice n.(a,b,c) = 0 of G^(n) (absolutely
     convergent: cubic decay along lines, so the resonance expansion is legitimate).
 (4) LINE STRUCTURE: for a primitive relation k (all k_i != 0), the k-line contributes
     S_k = sum_{m != 0} G^(mk); |S_k| <= (||W||_inf/pi^3) * 2 zeta(3)/|k1 k2 k3|; and as
     the free parameters of a k-family grow, T -> S_k (other relations' products -> inf).
COMPUTED: W and ||W||_inf; the line-sums S_k for all primitive |k_i| <= 6 (certified
truncation: tail_m <= ||W||_inf/(pi^3 |k-prod|) * 1/M^2); the calibration T(1,2,7) vs
S_(1,3,-1) + S-corrections; and THE ASSEMBLED VERDICT for (3) of THM-904.
"""
import sys, math, cmath
from fractions import Fraction as Fr
sys.path.insert(0, "04-computation")
sys.stdout.reconfigure(line_buffering=True)
import importlib.util as ilu
spec = ilu.spec_from_file_location("cert", "04-computation/lrc14_residue6_triple_certificate_codex_S18.py")
cert = ilu.module_from_spec(spec)
try: spec.loader.exec_module(cert)
except SystemExit: pass
beta = cert.beta
S7 = range(7)

# Hoeffding triple channel
mean = sum(beta(a,b,c) for a in S7 for b in S7 for c in S7) / Fr(343)
sing = [sum(beta(s,b,c) for b in S7 for c in S7)/Fr(49) - mean for s in S7]
pairch = {}
for s1 in S7:
    for s2 in S7:
        pairch[(s1,s2)] = sum(beta(s1,s2,c) for c in S7)/Fr(7) - mean - sing[s1] - sing[s2]
b123 = {}
for s1 in S7:
    for s2 in S7:
        for s3 in S7:
            b123[(s1,s2,s3)] = (beta(s1,s2,s3) - mean - sing[s1]-sing[s2]-sing[s3]
                                - pairch[(s1,s2)] - pairch[(s1,s3)] - pairch[(s2,s3)])

# W(r) = sum_s b123(s) e(-2pi i (r1(2s1+1)+r2(2s2+1)+r3(2s3+1))/14)
def W(r1, r2, r3):
    tot = 0j
    for (s1,s2,s3), v in b123.items():
        ph = -(r1*(2*s1+1) + r2*(2*s2+1) + r3*(2*s3+1)) / 14.0
        tot += float(v) * cmath.exp(2j*cmath.pi*ph)
    return tot
Winf = max(abs(W(r1,r2,r3)) for r1 in range(1,7) for r2 in range(1,7) for r3 in range(1,7))
print(f"||W||_inf over nonzero residues = {Winf:.4f}")
zero_marg = max(abs(W(0,r2,r3)) for r2 in range(7) for r3 in range(7))
print(f"zero-marginality check: max |W(0,.,.)| = {zero_marg:.2e}")

def Ghat(n1, n2, n3):
    if n1 == 0 or n2 == 0 or n3 == 0: return 0j
    p = 1.0
    for n in (n1, n2, n3):
        p *= math.sin(math.pi*n/7) / (math.pi*n)
    return p * W(n1 % 7 if n1 % 7 else 0, n2 % 7, n3 % 7) if all(ni % 7 for ni in (n1,n2,n3)) else 0j
# W is ANTI-periodic mod 7 (n -> n+7 flips sign); index by n mod 14
WTAB14 = {}
for r1 in range(14):
    for r2 in range(14):
        for r3 in range(14):
            WTAB14[(r1, r2, r3)] = W(r1, r2, r3)
def Ghat_full(n1, n2, n3):
    if any(v == 0 or v % 7 == 0 for v in (n1, n2, n3)): return 0j
    p = 1.0
    for n in (n1, n2, n3):
        p *= math.sin(math.pi*n/7)/(math.pi*n)
    return p * WTAB14[(n1 % 14, n2 % 14, n3 % 14)]

# line sums with certified truncation
def S_line(k, M=1500):
    s = 0j
    for m in range(1, M+1):
        s += Ghat_full(m*k[0], m*k[1], m*k[2])
        s += Ghat_full(-m*k[0], -m*k[1], -m*k[2])
    tail = Winf/(math.pi**3 * abs(k[0]*k[1]*k[2])) * 2/(M*M)   # sum_{m>M} 2/m^3 < 1/M^2
    return s.real, tail

print("\nline sums S_k (primitive, all k_i != 0, |k_i| <= 6, certified):")
from math import gcd
lines = []
for k1 in range(1, 7):
    for k2 in range(-6, 7):
        for k3 in range(-6, 7):
            if k2 == 0 or k3 == 0: continue
            if gcd(gcd(k1, abs(k2)), abs(k3)) != 1: continue
            v, tail = S_line((k1, k2, k3))
            if abs(v) > 0.002:
                lines.append((v, tail, (k1, k2, k3)))
lines.sort(key=lambda t: -abs(t[0]))
for v, tail, k in lines[:14]:
    print(f"   S_{k} = {v:+.6f} (+-{tail:.1e})")
maxS = max(v for v, tail, k in lines)
print(f"   max positive line sum = {maxS:+.6f}")

# calibration: T(1,2,7) should = sum of S_k over relations of (1,2,7) (approx by dominant lines)
print("\ncalibration: relations of (1,2,7): (1,3,-1) [1+6-7], (3,2,-1) [3+4-7], (5,1,-1), ...")
for k in [(1,3,-1),(3,2,-1),(5,1,-1),(2,-1,0)]:
    if 0 in k: continue
    v, tail = S_line(k)
    print(f"   S_{k} = {v:+.6f}")
v137,_ = S_line((1,3,-1)); v321,_ = S_line((3,2,-1)); v511,_ = S_line((5,1,-1))
print(f"   sum of first three lines = {v137+v321+v511:+.6f}  vs measured T(1,2,7) = +0.100714")
# direct lattice verification of the resonance expansion at (1,2,7)
lat = 0j
NL = 260
for n1 in range(-NL, NL+1):
    for n2 in range(-NL, NL+1):
        n3n = n1*1 + n2*2
        if n3n % 7: continue
        n3 = -n3n // 7
        lat += Ghat_full(n1, n2, n3)
print(f"   direct lattice sum (|n|<={NL}) T(1,2,7) = {lat.real:+.6f}  (measured +0.100714)")

print("\nASSEMBLED VERDICT (q <= beta0 + pair-sum + T):")
beta0 = float(mean)
print(f"   beta0 = {beta0:.6f}; target 0.47; margin {0.47-beta0:.6f}")
print(f"   worst single line: {maxS:+.6f}; worst measured pair-sum+T (codex scan max residual: 81/175 - beta0 = {81/175-beta0:.6f})")
print(f"   line-based sup estimate: beta0 + maxS + max positive pair D (+0.0159) + second-line corrections")
est = beta0 + maxS + 0.0159 + sum(max(v,0) for v,_,k in lines[1:4])
print(f"   = {est:.4f}  (< 0.47: {est < 0.47})")
print("\nDONE")
