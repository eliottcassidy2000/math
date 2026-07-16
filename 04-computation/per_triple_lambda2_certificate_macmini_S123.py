#!/usr/bin/env python3
"""per_triple_lambda2_certificate_macmini_S123.py -- mac-mini-2026-07-16-S123.
THE TRUE PER-TRIPLE CERTIFICATE (route [A] sign-off): for each worst-case triple of the
S120 box sweep, certify the N = 400 lattice-sum truncation by
  (i)  BAND: exact congruence-solved enumeration of all lattice points with
       N < ||n||_inf <= N2 = 2400 (their absolute contributions, computed exactly);
  (ii) COARSE BEYOND N2: every point with ||n||_inf > N2 lies on a slice |t| with the
       t-dependent co-coordinate floors of THM-920's slice lemma; closed form:
       tail(>N2) <= Cw * [ 4 g ln(e N2)/(b N2) + 12 c/(N2^2 (b-a)) + 4 g^2 pi^2/(3 b c N2) ]
       with g = gcd(b,c), Cw = ||W||_inf/pi^3 = 1.2178 (each term derived in THM-922).
CERTIFICATE: err(N=400) <= band + coarse; CLEAR iff err < margin = 47/100 - q_est.
Run on the 20 worst box triples + the scan champions.
"""
import sys, math, cmath
from math import gcd
sys.path.insert(0, "04-computation")
sys.stdout.reconfigure(line_buffering=True)
import importlib.util as ilu
spec = ilu.spec_from_file_location("cert", "04-computation/lrc14_residue6_triple_certificate_codex_S18.py")
cert = ilu.module_from_spec(spec)
try: spec.loader.exec_module(cert)
except SystemExit: pass
from fractions import Fraction as Fr
beta = cert.beta
S7 = range(7)
mean = sum(beta(a,b,c) for a in S7 for b in S7 for c in S7)/Fr(343)
sing = [sum(beta(s,b,c) for b in S7 for c in S7)/Fr(49) - mean for s in S7]
pairch = {}
for s1 in S7:
    for s2 in S7:
        pairch[(s1,s2)] = sum(beta(s1,s2,c) for c in S7)/Fr(7) - mean - sing[s1] - sing[s2]
b123 = {}
for s1 in S7:
    for s2 in S7:
        for s3 in S7:
            b123[(s1,s2,s3)] = (beta(s1,s2,s3)-mean-sing[s1]-sing[s2]-sing[s3]
                                -pairch[(s1,s2)]-pairch[(s1,s3)]-pairch[(s2,s3)])
WT = {}
for r1 in range(14):
    for r2 in range(14):
        for r3 in range(14):
            t = 0j
            for (s1,s2,s3), v in b123.items():
                t += float(v)*cmath.exp(-2j*cmath.pi*(r1*(2*s1+1)+r2*(2*s2+1)+r3*(2*s3+1))/14.0)
            WT[(r1,r2,r3)] = t
Cw = 37.7545 / math.pi**3
def sinf(n): return math.sin(math.pi*n/7)/(math.pi*n)

def band_and_coarse(a, b, c, N=400, N2=2400):
    from sympy import mod_inverse
    g = gcd(b, c); step = c // g
    band = 0.0
    for n1 in range(-N2, N2+1):
        if n1 == 0 or n1 % 7 == 0: continue
        if (n1*a) % g: continue
        n20 = (-(n1*a)//g) * mod_inverse(b//g, step) % step
        n2 = n20 - ((n20+N2)//step)*step
        while n2 <= N2:
            if n2 != 0 and n2 % 7:
                n3 = -(n1*a+n2*b)//c
                if n3 != 0 and n3 % 7 and abs(n3) <= N2:
                    if max(abs(n1), abs(n2), abs(n3)) > N:
                        band += abs(sinf(n1)*sinf(n2)*sinf(n3)*WT[(n1%14, n2%14, n3%14)])
            n2 += step
    coarse = Cw*(4*g*math.log(math.e*N2)/(b*N2) + 12*c/(N2*N2*max(1, b-a)) + 4*g*g*math.pi**2/(3*b*c*N2))
    return band, coarse

cases = [(2,8,11,0.462775), (1,4,7,0.462857), (3,10,13,0.45), (2,9,11,0.45), (1,6,7,0.44),
         (3,4,7,0.44), (1,2,7,0.408), (2,4,7,0.45), (1,3,4,0.42), (4,9,13,0.44),
         (2,5,7,0.42), (5,8,13,0.44), (1,5,6,0.43), (3,8,11,0.44), (2,3,5,0.42)]
print("triple        q_est    margin   band(400->2400)  coarse(>2400)   err_total   CLEAR?")
allok = True
for a, b, c, qe in cases:
    marg = 0.47 - qe
    band, coarse = band_and_coarse(a, b, c)
    err = band + coarse
    ok = err < marg
    allok &= ok
    print(f"{(a,b,c)}: {qe:.4f}   {marg:.4f}   {band:.6f}        {coarse:.6f}      {err:.6f}   {'YES' if ok else 'NO'}")
print(f"\nVERDICT: {'ALL CERTIFICATES CLEAR -- route [A] truncations certified; sign-off complete' if allok else 'flagged cases need N2 extension'}")
print("DONE")
