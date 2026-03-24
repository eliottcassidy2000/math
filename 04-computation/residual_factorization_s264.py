#!/usr/bin/env python3
"""
residual_factorization_s264.py -- opus-2026-03-23-S264

UNDERSTAND THE RESIDUALS THROUGH FACTORIZATION

From S263: the exact decomposition
  T_n = twin_SL + 2*E + residual
  residual = complex_SL + MW = cut-cycle interaction

Residual sequence: 0, 2, 16, 76, 444, 2926, 31442

THE FACTORIZATION THEOREM (S263):
  Fix_tournament(σ) = Fix_even(σ) × 2^{k-1}
  twin_SL(σ) = 2C(f,2) × Fix_even(σ)

So: Fix_tournament = Fix_even × 2^{k-1} factors into:
  CYCLE PART: Fix_even = 2^{arc_orbits - k + 1}
  CUT PART: 2^{k-1}

The twin_SL captures C(f,2)/2^{k-1} of the cut part (twin constraint).
The residual = what remains after twin + simple edges are accounted for.

Let's factorize the residual sequence and understand its structure.

Author: opus-2026-03-23-S264
"""

import sys
from math import comb, factorial, gcd
from fractions import Fraction
from collections import Counter
from functools import reduce
sys.stdout.reconfigure(line_buffering=True)

def banner(title):
    print("\n" + "="*72)
    print("  " + title)
    print("="*72 + "\n")

# Known values
T = {3:4, 4:16, 5:88, 6:704, 7:8912, 8:188288, 9:6847200}
E = {3:1, 4:5, 5:30, 6:290, 7:4086, 8:91161, 9:3380751}
twin_SL = {3:2, 4:4, 5:12, 6:48, 7:296, 8:3040, 9:54256}

# Residual = T - twin_SL - 2E
residual = {n: T[n] - twin_SL[n] - 2*E[n] for n in T if n in E}

banner("THE RESIDUAL SEQUENCE")

print("  n | T_n      | twin_SL | 2E       | residual | res/T")
print("  --|----------|---------|----------|----------|------")
for n in range(3, 10):
    if n not in residual: continue
    r = residual[n]
    print("  %d | %8d | %7d | %8d | %8d | %.6f" %
          (n, T[n], twin_SL[n], 2*E[n], r, r/T[n]))

res_seq = [residual[n] for n in range(3, 10)]
print("\n  Residual sequence: %s" % res_seq)

banner("FACTORIZATIONS OF RESIDUALS")

for n in range(3, 10):
    r = residual[n]
    if r == 0:
        print("  n=%d: residual = 0" % n)
        continue
    temp = abs(r); factors = []
    for p in range(2, min(temp+1, 100000)):
        if p*p > temp and temp > 1: factors.append(temp); break
        while temp % p == 0: factors.append(p); temp //= p
    print("  n=%d: residual = %d = %s" %
          (n, r, " × ".join(str(f) for f in factors)))

banner("RESIDUAL RATIOS AND PATTERNS")

print("  Consecutive ratios residual(n+1)/residual(n):")
for i in range(1, len(res_seq)):
    if res_seq[i-1] != 0:
        print("    res(%d)/res(%d) = %d/%d = %.4f" %
              (i+3, i+2, res_seq[i], res_seq[i-1], res_seq[i]/res_seq[i-1]))

print("\n  Residual / (n-2)!:")
for n in range(3, 10):
    r = residual[n]
    print("    n=%d: %d / %d = %s = %.6f" %
          (n, r, factorial(n-2), Fraction(r, factorial(n-2)), r/factorial(n-2)))

print("\n  Residual / V_n:")
V = {3:2, 4:4, 5:12, 6:56, 7:456, 8:6880, 9:191536}
for n in range(3, 10):
    r = residual[n]
    print("    n=%d: %d / %d = %.4f" % (n, r, V[n], r/V[n]))

print("\n  Residual / n!:")
for n in range(3, 10):
    r = residual[n]
    print("    n=%d: %d / %d = %s" % (n, r, factorial(n), Fraction(r, factorial(n))))

banner("RESIDUAL AS CUT-CYCLE INTERACTION")

print("""
  From the factorization theorem (S263):

  Fix_tournament(σ) = Fix_even(σ) × 2^{k-1}

  where:
  - Fix_even(σ) = 2^{arc_orbits(σ) - k + 1} = the CYCLE part
  - 2^{k-1} = the CUT part (k = number of cycles of σ)

  TOURNAMENT = CYCLE × CUT

  The twin_SL captures the part where the CUT part has a specific
  structure (twin pair among fixed points):
    twin(σ) = 2C(f,2) × Fix_even(σ)

  The TOTAL = V_n uses Fix_even × 2^{k-1} (all of cut space).

  The difference (total - twin - 2E) = RESIDUAL:
  residual = #{orbits where the cut-cycle interaction creates
              complex self-loops or multi-edges that aren't twins}

  The residual MEASURES the interaction between the cycle space
  (Fix_even) and the cut space (2^{k-1}).

  When cycle and cut are INDEPENDENT: residual = 0.
  When they INTERACT: residual > 0.

  At n=3: residual = 0 (perfect independence at the triangle level!)
  At n≥4: residual > 0 (interaction emerges with the 4th vertex)
""")

# Decompose residual by cycle type
banner("RESIDUAL BY CYCLE TYPE")

def pair_orbits_from_type(ct):
    within = sum((c - 1) // 2 for c in ct)
    across = sum(gcd(ct[i], ct[j]) for i in range(len(ct)) for j in range(i+1, len(ct)))
    return within + across

def multinomial(n, ct):
    ct_counter = Counter(ct)
    denom = 1
    for c, a in ct_counter.items():
        denom *= (c ** a) * factorial(a)
    return factorial(n) // denom

def gen_odd_partitions(n):
    result = []
    def _gen(remaining, max_part, current):
        if remaining == 0:
            result.append(tuple(sorted(current)))
            return
        start = min(remaining, max_part)
        if start % 2 == 0: start -= 1
        for p in range(start, 0, -2):
            _gen(remaining - p, p, current + [p])
    _gen(n, n, [])
    return list(set(result))

for n in range(3, 10):
    parts = gen_odd_partitions(n)
    nf = factorial(n)

    print("  n=%d:" % n)

    total_V = 0; total_T = 0; total_twin = 0

    for ct in sorted(parts, key=lambda x: -len(x)):
        mult = multinomial(n, ct)
        po = pair_orbits_from_type(ct)
        k = len(ct)
        f = ct.count(1)
        nfc = k - f

        fix_even = 2 ** (po - k + 1)
        fix_tourn = fix_even * (2 ** (k - 1))
        twin_contrib = 2 * comb(f, 2) * fix_even if f >= 2 else 0

        V_contrib = mult * fix_tourn // nf
        T_contrib = mult * comb(f, 2) * fix_tourn // nf if f >= 2 else 0
        twin_c = mult * twin_contrib // nf

        total_V += mult * fix_tourn
        total_T += mult * comb(f, 2) * fix_tourn if f >= 2 else 0
        total_twin += mult * twin_contrib

        if n <= 7:
            # The "cut surplus" for this type: how much more than twin?
            # Fix_tourn = Fix_even × 2^{k-1}
            # twin = 2C(f,2) × Fix_even
            # ratio = Fix_tourn / twin = 2^{k-1} / (2C(f,2)) = 2^{k-2} / C(f,2)
            if f >= 2 and twin_contrib > 0:
                ratio = Fraction(fix_tourn, twin_contrib)
                print("    %s: k=%d, f=%d, Fix_even=2^%d, ratio=Fix_tourn/twin=%s" %
                      (ct, k, f, po-k+1, ratio))

    V_check = total_V // nf
    T_check = total_T // nf
    twin_check = total_twin // nf

    V_known = {3:2, 4:4, 5:12, 6:56, 7:456, 8:6880, 9:191536}
    print("    V=%d (known %d), twin=%d (known %d)" %
          (V_check, V_known[n], twin_check, twin_SL[n]))
    print()

banner("THE RATIO Fix_tourn/twin FOR EACH CYCLE TYPE")

print("""
  For cycle type with k cycles and f fixed points:

  Fix_tourn(σ) / twin(σ) = 2^{k-1} / (2C(f,2))
                          = 2^{k-2} / C(f,2)
                          = 2^{k-2} / (f(f-1)/2)
                          = 2^{k-1} / (f(f-1))

  For the identity (k=n, f=n):
    ratio = 2^{n-1} / (n(n-1)) = 2^{n-1} / 2m
    = 2^{n-2} / m

  This ratio GROWS exponentially with n.
  At n=5: 2^3/10 = 0.8 (twin captures 1/0.8 = 125% — wait, that's >1?)

  Actually: twin(σ) can be LARGER than Fix_tourn(σ) when f is large
  relative to k. Let me re-examine.

  twin(σ) = 2C(f,2) × Fix_even(σ)
  Fix_tourn(σ) = Fix_even(σ) × 2^{k-1}

  ratio = twin / Fix_tourn = 2C(f,2) / 2^{k-1} = f(f-1) / 2^{k-1}

  For identity (k=n, f=n):
    twin/Fix = n(n-1) / 2^{n-1}

  n=3: 6/4 = 1.5 (twin > Fix!)
  n=4: 12/8 = 1.5
  n=5: 20/16 = 1.25
  n=6: 30/32 = 0.9375
  n=7: 42/64 = 0.65625
  n=8: 56/128 = 0.4375

  The identity term has twin > Fix for n ≤ 5!
  This means the twin sum OVERESTIMATES at small n.
  The correction comes from non-identity cycle types where twin < Fix.
""")

print("  twin(identity) / Fix_tourn(identity) = n(n-1)/2^{n-1}:")
for n in range(3, 12):
    ratio = Fraction(n*(n-1), 2**(n-1))
    print("    n=%d: %s = %.6f" % (n, ratio, float(ratio)))

banner("THE RESIDUAL = BURNSIDE CORRECTION TO TWIN APPROXIMATION")

print("""
  The formula: T = twin_SL + 2E + residual

  Rearranging: 2E = T - twin_SL - residual

  The TWIN approximation: E ≈ (T - twin_SL) / 2
  The EXACT: E = (T - twin_SL - residual) / 2

  So: E_twin - E_exact = residual / 2

  The residual / 2 sequence: 0, 1, 8, 38, 222, 1463, 15721

  FACTORIZATIONS:
""")

for n in range(3, 10):
    r2 = residual[n] // 2
    if r2 == 0:
        print("  n=%d: residual/2 = 0" % n)
        continue
    temp = abs(r2); factors = []
    for p in range(2, min(temp+1, 100000)):
        if p*p > temp and temp > 1: factors.append(temp); break
        while temp % p == 0: factors.append(p); temp //= p
    print("  n=%d: residual/2 = %d = %s" %
          (n, r2, " × ".join(str(f) for f in factors)))

# The residual/2 sequence: 0, 1, 8, 38, 222, 1463, 15721
r2_seq = [residual[n]//2 for n in range(3, 10)]
print("\n  residual/2 = %s" % r2_seq)

print("\n  Ratios res/2(n+1) / res/2(n):")
for i in range(1, len(r2_seq)):
    if r2_seq[i-1] > 0:
        print("    %d/%d = %.4f" % (r2_seq[i], r2_seq[i-1], r2_seq[i]/r2_seq[i-1]))

print("\n  residual/2 / (n-2)!:")
for n in range(4, 10):
    r2 = residual[n] // 2
    print("    n=%d: %d / %d = %s = %.6f" %
          (n, r2, factorial(n-2), Fraction(r2, factorial(n-2)), r2/factorial(n-2)))

# Check: is residual/2 related to subfactorials?
# D_n = n! × sum_{k=0}^n (-1)^k/k! (derangements)
# D_1=0, D_2=1, D_3=2, D_4=9, D_5=44, D_6=265, D_7=1854
# res/2: 0, 1, 8, 38, 222, 1463, 15721
# Not matching derangements.

print("\n  residual/2 vs 2×D_{n-2} (D = derangements):")
for n in range(4, 10):
    r2 = residual[n] // 2
    n2 = n - 2
    # Derangement D_k = k! × sum (-1)^j/j! for j=0..k
    D_n2 = factorial(n2) * sum((-1)**j * Fraction(1, factorial(j)) for j in range(n2+1))
    D_n2 = int(D_n2)
    print("    n=%d: res/2=%d, 2×D_%d=%d, ratio=%.4f" %
          (n, r2, n2, 2*D_n2, r2/(2*D_n2) if D_n2 > 0 else 0))

banner("THE RESIDUAL THROUGH THE GF(2) LENS")

print("""
  From the factorization theorem:
  Fix_tournament = Fix_even × 2^{k-1}

  The 2^{k-1} is the CUT SPACE dimension.
  The CUT SPACE = GF(2)^{k-1} acts on the score sequence.

  twin_SL uses C(f,2) of the cut dimensions (twin pairs from fixed pts).
  The RESIDUAL comes from the REMAINING cut dimensions:
  2^{k-1} - 2C(f,2) of the cut space goes into non-twin self-loops + multi-edges.

  For the identity (k=n, f=n):
    Cut space = GF(2)^{n-1}
    Twin constraint uses C(n,2) = n(n-1)/2 dimensions... wait.
    Actually twin uses C(f,2) = C(n,2) fixed-point pairs.
    But the cut space is only (n-1)-dimensional!
    C(n,2) > n-1 for n ≥ 4, so the twin COUNT exceeds the cut DIMENSION.

  This is the source of the residual: the twin count OVERESTIMATES
  because it treats fixed-point pairs as independent, but they're
  constrained by the (n-1)-dimensional cut space.

  The residual = RANK DEFICIENCY of the twin constraint system.
  At n=3: C(3,2)=3, dim(cut)=2. 3-2=1 redundancy. But residual=0!
  At n=4: C(4,2)=6, dim(cut)=3. 6-3=3 redundancy. Residual=2.
  At n=5: C(5,2)=10, dim(cut)=4. 10-4=6 redundancy. Residual=16.

  The redundancy = C(n,2) - (n-1) = C(n-1,2) = dimension of the CYCLE SPACE.

  So: the residual is controlled by the CYCLE SPACE DIMENSION C(n-1,2).
  Let me check: residual/2 vs 2^{C(n-1,2)} ...

  residual/2: 0, 1, 8, 38, 222, 1463, 15721
  C(n-1,2):   1, 3, 6, 10, 15,   21,    28
  2^{C(n-1,2)}: 2, 8, 64, 1024, 32768, 2097152, ...

  residual/2 / 2^{C(n-1,2)}: -, 1/8, 1/8, 38/1024, 222/32768, ...
                             = -, 0.125, 0.125, 0.037, 0.0068, ...

  The fraction DECREASES — the residual is a SMALL fraction of the
  cycle space. It measures the INTERACTION (non-independence) between
  cut and cycle, which vanishes as the spaces grow apart.
""")

print("  residual/2 vs cycle space dimension C(n-1,2):")
for n in range(3, 10):
    r2 = residual[n] // 2
    cs = comb(n-1, 2)
    print("    n=%d: res/2=%d, C(n-1,2)=%d, 2^{C(n-1,2)}=%d, ratio=%.8f" %
          (n, r2, cs, 2**cs, r2 / 2**cs if 2**cs > 0 else 0))

print("\nDone. opus-2026-03-23-S264")
