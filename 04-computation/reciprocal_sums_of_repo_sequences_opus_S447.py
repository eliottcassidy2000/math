#!/usr/bin/env python3
"""
SUPERSEDED AS A THEOREM SOURCE by MISTAKE-209 and THM-2000/2005.
This file is retained as a historical exploratory prefix calculator.  It
mixes indexed multiplicity with literal support, labels uncertified census
prefixes as constants, uses C(n,3) for the false c3 maximum, and treats the
conjectural full H-value spectrum as proved.  Do not cite those rows.

RECIPROCAL SUMS of the repo's integer sequences (opus-S447).
Owner: the reciprocal of an integer sequence is a subset of the harmonic numbers; study sum 1/a_n
for as many of our sequences as possible; think figurate reciprocals (triangular -> 2), Abel-Dini,
Bertrand series.  The surviving result is the pure-binomial figurate ladder;
the broader poly/#P and arithmetic interpretations were retracted.
"""
from fractions import Fraction as F
from math import comb, log, e, pi

print("SUPERSEDED HISTORICAL PREFIX CALCULATOR -- USE THM-2000/2005")
print("All census rows below are indexed finite prefixes unless stated otherwise.")
print("="*72)
print("PART 1 -- FIGURATE / BINOMIAL invariant-SIZE sequences: EXACT rational sums")
print("  general law:  sum_{n>=k} 1/C(n,k) = k/(k-1)  (telescoping / hockey-stick)")
print("="*72)
def recip_binom_sum(k, N=4000):
    return sum(F(1, comb(n,k)) for n in range(k, N))
for k in range(2,7):
    s=recip_binom_sum(k, 3000)
    print(f"  sum 1/C(n,k), k={k}: partial(3000) = {float(s):.10f}   closed form k/(k-1) = {F(k,k-1)} = {k/(k-1):.10f}")
print("\n  TOURNAMENT identities (char_S coefficients / object sizes, THM-1920):")
print(f"    arc count = C(n,2) = triangular T_(n-1)        -> sum 1/arc  = 2          (Downey-Ong-Sellers triangular)")
print(f"    # tiles   = C(n-1,2) = triangular              -> sum 1/tiles= 2")
print(f"    triple slots C(n,3) = tetrahedral               -> sum 1/slots= 3/2")
print(f"    var-max   = 2*C(n,3) (transitive, THM-1930)    -> sum 1/varmax= 3/4")
print(f"    realizable H=1+2^(n-2) neighbor support (n>=3): sum = {float(sum(F(1,1+2**(n-2)) for n in range(3,60))):.8f}")

print("\n" + "="*72)
print("PART 2 -- COUNTING sequences: indexed finite prefixes (not constants)")
print("="*72)
SEQS = {
 "A000568 #tournaments (n>=1)":   [1,1,2,4,12,56,456,6880,191536,9733056,903753248,154108311168,48542114686912],
 "A038375 max Ham paths (n>=1)":  [1,1,3,5,15,45,189,661,3357],
 "A051337 strong tournaments (n>=3)":[1,1,6,35,353,6008],
 "A002854 even graphs E_n (n>=3)":[2,3,7,16,54,243,2038],
 "seed census modular-prime(n>=1)":[1,1,1,3,15],           # n=4=0 skipped (no reciprocal)
 "A000255 succession W (n>=1)":   [1,1,3,11,53,309,2119,16687,148329],
 "A000571 score seqs (n>=3)":     [2,4,9,22,59,167,490],
 "A002620 quarter-square (n>=2)": [1,2,4,6,9,12,16,20,25,30,36,42,49],  # slow (n^2/4)
}
for name,seq in SEQS.items():
    s=sum(1.0/a for a in seq if a>0)
    print(f"  {name:<34}: sum(1/a) over {len(seq)} terms = {s:.8f}")

# inverse-symbolic hints for the fast ones (compare to e, e-1, rationals via CF)
def cf(x, n=8):
    a=[]
    for _ in range(n):
        i=int(x//1); a.append(i); x-=i
        if x<1e-9: break
        x=1/x
    return a
print("\n  exploratory continued fractions of uncertified indexed prefixes:")
for name in ("A000568 #tournaments (n>=1)","A038375 max Ham paths (n>=1)","A002854 even graphs E_n (n>=3)"):
    s=sum(1.0/a for a in SEQS[name] if a>0)
    print(f"    {name:<34}: {s:.9f}  CF={cf(s)}   (e-1={e-1:.6f}, e={e:.6f}, 3/2={1.5})")

print("\n" + "="*72)
print("PART 3 -- CAYLEY-DICKSON levels n=2^k+1 (n=3,5,9,17,33,...): sum 1/(2^k+1)")
print("="*72)
cd=sum(1.0/(2**k+1) for k in range(1,60))
print(f"  sum_{{k>=1}} 1/(2^k+1) = {cd:.10f}   (cf. Erdos-Borwein sum 1/(2^k-1) = {sum(1.0/(2**k-1) for k in range(1,60)):.10f})")

print("\n" + "="*72)
print("PART 4 -- ABEL-DINI: exponent one is the divergent critical member")
print("  divergent S=sum a_n: sum a_n/S_n DIVERGES but sum a_n/S_n^(1+eps) CONVERGES (any eps>0)")
print("="*72)
# Demo on all odd numbers.  This is only a conditional model for the unknown
# global H-value support, not a theorem about that support.
odds=[2*n-1 for n in range(1,200000)]
terms=[1.0/a for a in odds]
S=0.0; ad0=0.0; ad1=0.0
partial=0.0
for t in terms:
    partial+=t
    ad0 += t/partial                      # sum t_n / S_n  -> diverges (~log log)
    ad1 += t/partial**1.05                 # sum t_n / S_n^1.05 -> converges
print(f"  CONDITIONAL all-odd model (not the proved H-spectrum); partial(2e5)={partial:.4f}")
print(f"    Abel-Dini boundary sum t/S     (diverges) partial = {ad0:.4f}")
print(f"    Abel-Dini sum t/S^1.05         (converges) partial = {ad1:.4f}")

print("\n" + "="*72)
print("PART 5 -- BERTRAND scale: where sequences sit (sum 1/(n (ln n)^a): div a<=1, conv a>1)")
print("="*72)
for a in (0.0,1.0,1.5,2.0):
    s=sum(1.0/(n*log(n)**a) for n in range(2,100000)) if a>0 else sum(1.0/n for n in range(2,100000))
    print(f"  sum 1/(n (ln n)^{a}) partial(1e5) = {s:.4f}  ({'DIVERGES' if a<=1 else 'converges'})")
print("\n READING: pure binomial rows have the proved rational ladder k/(k-1).")
print(" Census decimals above are indexed prefixes; growth proves convergence only with a tail bound,")
print(" and it does not prove an arithmetic type.  Global H-value support divergence remains OPEN.")
print(" Use THM-2000/2005 for support semantics, Abel-Dini, and certified rows.")
