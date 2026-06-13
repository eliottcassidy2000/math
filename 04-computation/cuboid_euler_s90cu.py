#!/usr/bin/env python3
"""
cuboid_euler_s90cu.py — Euler product for perfect cuboid existence.
Perfect cuboid: a,b,c, sqrt(a^2+b^2), sqrt(b^2+c^2), sqrt(a^2+c^2), sqrt(a^2+b^2+c^2) all integer.
Local factor at p = fraction of (a,b,c) mod p^k with all four sums being QR.
opus-2026-03-16-S90cu
"""
import math, sys, os
from fractions import Fraction
from collections import Counter
from math import isqrt

def local_density(m):
    """Count (a,b,c) mod m with a^2+b^2, b^2+c^2, a^2+c^2, a^2+b^2+c^2 all QR mod m."""
    qr = {(x*x) % m for x in range(m)}
    sq = [(a*a) % m for a in range(m)]
    count = 0
    for a in range(m):
        a2 = sq[a]
        for b in range(m):
            b2 = sq[b]
            if (a2+b2) % m not in qr: continue
            for c in range(m):
                c2 = sq[c]
                if (b2+c2) % m not in qr: continue
                if (a2+c2) % m not in qr: continue
                if (a2+b2+c2) % m not in qr: continue
                count += 1
    return count, m**3

def sieve_primes(N):
    s = [True]*(N+1); s[0]=s[1]=False
    for i in range(2, int(N**0.5)+1):
        if s[i]:
            for j in range(i*i, N+1, i): s[j]=False
    return [i for i in range(2, N+1) if s[i]]

def main():
    primes = sieve_primes(200)
    P = print
    P("="*72); P("EULER PRODUCT ANALYSIS FOR PERFECT CUBOID EXISTENCE"); P("="*72)

    # Part 1: Local densities mod p
    P("\n-- PART 1: Local density sigma_p = #{good (a,b,c) mod p} / p^3 --\n")
    P(f"{'p':>5}  {'count':>10}  {'p^3':>10}  {'sigma_p':>16}  {'decimal':>14}")
    P("-"*62)
    densities = {}
    for p in primes:
        cnt, tot = local_density(p)
        densities[p] = Fraction(cnt, tot)
        P(f"{p:>5}  {cnt:>10}  {tot:>10}  {str(densities[p]):>16}  {float(densities[p]):>14.10f}")

    # Part 2: Partial Euler product
    P("\n-- PART 2: Partial Euler product S(N) = prod_{p<=N} sigma_p --\n")
    partial = Fraction(1)
    for p in primes:
        partial *= densities[p]
        P(f"  p<={p:>3}: S = {float(partial):.12e}  (num has {len(str(partial.numerator))} digits)")
    P(f"\nFinal (p<=200): {float(partial):.15e}")
    P(f"  log2 = {math.log2(float(partial)):.4f}")

    # Part 3: Correlation ratio
    P("\n-- PART 3: Correlation R_p = sigma_p / ((p+1)/(2p))^4 (first 20 odd primes) --\n")
    P(f"{'p':>5}  {'sigma':>12}  {'indep':>12}  {'R_p':>10}")
    for p in primes[1:21]:
        indep = (Fraction(p+1, 2*p))**4
        P(f"{p:>5}  {float(densities[p]):>12.8f}  {float(indep):>12.8f}  {float(densities[p]/indep):>10.6f}")

    # Part 4: Higher-order local densities mod p^2
    P("\n-- PART 4: Local density mod p^2 for small primes --\n")
    for p in [2, 3, 5, 7]:
        c1, t1 = local_density(p)
        c2, t2 = local_density(p**2)
        s1, s2 = Fraction(c1,t1), Fraction(c2,t2)
        r = s2/s1 if s1 > 0 else Fraction(0)
        P(f"p={p}: mod {p}: {c1}/{t1}={float(s1):.10f}  "
          f"mod {p**2}: {c2}/{t2}={float(s2):.10f}  ratio={r} ({float(r):.8f})")

    # Part 5: Convergence analysis
    P("\n-- PART 5: Convergence analysis --\n")
    log_sum = sum(math.log(float(densities[p])) for p in primes if float(densities[p]) > 0)
    zero_primes = [p for p in primes if densities[p] == 0]
    P(f"Sum log(sigma_p) for p<=200: {log_sum:.6f}")
    P(f"Partial product: {math.exp(log_sum):.15e}")
    if zero_primes:
        P(f"\n*** LOCAL OBSTRUCTION at primes: {zero_primes} ***")
    else:
        P("No local obstruction for p <= 200.")
        tail = [p for p in primes if p > 50]
        c_vals = [-(math.log(float(densities[p])))*p for p in tail]
        c_avg = sum(c_vals)/len(c_vals)
        P(f"Tail exponent c (sigma_p ~ 1-c/p): {c_avg:.4f}")
        P(f"  => product ~ exp(-c*sum 1/p) -> 0 logarithmically (c>0)")

    # Part 6: p=7 analysis
    P("\n-- PART 6: Special analysis at p=7 --\n")
    cnt7 = densities[7]
    P(f"mod 7:  {cnt7} = {float(cnt7):.10f}")
    c49, t49 = local_density(49)
    s49 = Fraction(c49, t49)
    P(f"mod 49: {c49}/{t49} = {float(s49):.10f}  ratio={s49/cnt7}")
    qr7 = {(x*x)%7 for x in range(7)}
    P(f"QR mod 7 = {sorted(qr7)}")
    zc = Counter()
    sq7 = [(a*a)%7 for a in range(7)]
    for a in range(7):
        for b in range(7):
            if (sq7[a]+sq7[b])%7 not in qr7: continue
            for c in range(7):
                if (sq7[b]+sq7[c])%7 not in qr7: continue
                if (sq7[a]+sq7[c])%7 not in qr7: continue
                if (sq7[a]+sq7[b]+sq7[c])%7 not in qr7: continue
                zc[sum(1 for x in [a,b,c] if x==0)] += 1
    P(f"Survivors by #coords=0 mod 7: {dict(sorted(zc.items()))}  total={sum(zc.values())}")

    # Part 7: Primorial sieve
    P("\n-- PART 7: Primorial sieve mod 2310 = 2*3*5*7*11 --\n")
    pd = Fraction(1)
    for p in [2,3,5,7,11]:
        pd *= densities[p]
        P(f"  After p={p}: {float(pd):.12f}")
    P(f"Combined: {pd} = {float(pd):.15f}")
    P(f"  => 1 in {float(1/pd):.0f} triples survive mod-2310 sieve")

    # Part 8: Near-miss search
    P("\n-- PART 8: Near-miss search (a<=b<=c<=200, >=3 of 4 diagonals perfect square) --\n")
    nm = []
    for a in range(1, 201):
        a2 = a*a
        for b in range(a, 201):
            b2 = b*b
            for c in range(b, 201):
                c2 = c*c
                sqs = [isqrt(v)**2==v for v in [a2+b2, b2+c2, a2+c2, a2+b2+c2]]
                if sum(sqs) >= 3:
                    labs = ['ab','bc','ac','abc']
                    P(f"  ({a},{b},{c}): {sum(sqs)}/4 good={[l for l,s in zip(labs,sqs) if s]}")
                    nm.append(1)
    P(f"Total: {len(nm)} triples with >=3 perfect-square diagonals")

    # Summary
    P("\n"+"="*72); P("SUMMARY"); P("="*72)
    if zero_primes:
        P(f"LOCAL OBSTRUCTION at p = {zero_primes} => cuboids impossible")
    else:
        P(f"No local obstruction for p<=200. Partial product = {float(partial):.6e}")
        P("Product -> 0 as N->inf (heuristic: cuboids extremely rare, not locally blocked)")
        P("The perfect cuboid problem remains OPEN.")

if __name__ == "__main__":
    outpath = os.path.abspath(os.path.join(os.path.dirname(__file__),
                              "..", "05-knowledge", "results", "cuboid_euler_s90cu.out"))
    class Tee:
        def __init__(self, *s): self.s = s
        def write(self, d):
            for x in self.s: x.write(d)
        def flush(self):
            for x in self.s: x.flush()
    f = open(outpath, "w")
    sys.stdout = Tee(sys.__stdout__, f)
    main()
    f.close(); sys.stdout = sys.__stdout__
    print(f"\nOutput saved to {outpath}")
