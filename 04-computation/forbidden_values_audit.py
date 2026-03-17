"""
forbidden_values_audit.py — opus-S74

AUDIT: Which claims about forbidden H values in the repo are TRUE vs FALSE?

FINDING: Multiple writeups and hypotheses make FALSE claims about the
connection between k-nacci traces / Mersenne numbers and forbidden H values.

TRUTH:
  - ONLY H=7 and H=21 are permanently forbidden (proved: THM-029, THM-079)
  - 31 = 2^5-1 IS achievable (at n=6, bits=146)
  - 63 = 2^6-1 IS achievable (at n=8, transient gap at n=7)
  - 127 = 2^7-1 IS achievable (at n=7)
  - The k-nacci trace sequence [1,3,7,11,21,39,71,...] hits BOTH forbidden
    and achievable values — the trace connection is COINCIDENTAL
  - The Mersenne connection is WRONG: 3,15,31,63,127 are all achievable

WHAT IS TRUE about why 7 and 21 are forbidden:
  - THM-029: H=7 requires (α₁=3, i₂=0), but 3 pairwise-conflicting cycles
    always force additional cycles, making α₁≥4
  - THM-079: H=21 blocked by OCF decomposition analysis + tournament forcing
  - HYP-1231: 7·3^k forbidden for k=0,1; achievable for k≥2. "Nilpotency 2"
  - HYP-1317: 7=Φ₃(2), 21=Φ₃(4). Third cyclotomic polynomial at even args.
    But Φ₃(6)=43 IS achievable, so this is necessary but not sufficient.
  - HYP-1315: Both I-polynomials factor through I(K₃,x)=(1+3x)

FALSE CLAIMS TO CORRECT:
  1. casual-writeup.md line 122: "These Mersenne numbers control which H
     values are forbidden" — FALSE, 31 is achievable
  2. formal-writeup.md line 337: "connects forbidden H values (7=2^3-1,
     31=2^5-1, 127=2^7-1) to Mersenne primes" — FALSE, 31 and 127 achievable
  3. substack-hooks.md Hook N: "connects the forbidden H-spectrum values
     (7=2^3-1, 31=2^5-1)" — FALSE
  4. HYP-1600: "Connects forbidden H values (7=2^3-1, 31=2^5-1, 127=2^7-1)
     to Mersenne primes" — FALSE
  5. riemann-zeta-tournament.md: "connects the entire forbidden spectrum to
     Mersenne numbers: 7=2^3-1, 31=2^5-1, 127=2^7-1" — FALSE
  6. Various results files calling 31 "forbidden"
"""

# Verification
print("=" * 70)
print("FORBIDDEN VALUE AUDIT")
print("=" * 70)

print("""
PROVED FORBIDDEN (permanent gaps):
  H = 7   (THM-029, all n)
  H = 21  (THM-079, all n)

PROVED ACHIEVABLE (despite being k-nacci traces or Mersenne numbers):
  H = 1   (transitive tournament, any n)
  H = 3   (3-cycle C_3, n=3)
  H = 11  (various tournaments, n=5)
  H = 15  (regular tournament, n=5)     [Mersenne 2^4-1]
  H = 31  (bits=146, n=6)               [Mersenne 2^5-1, tribonacci p_5... wait]
  H = 39  (n=7)                         [tribonacci p_6]
  H = 63  (n=8, transient gap at n=7)   [Mersenne 2^6-1]
  H = 71  (n=7)                         [tribonacci p_7]
  H = 127 (n=7)                         [Mersenne 2^7-1]
  H = 131 (n=7)                         [tribonacci p_8]

WAIT — 31 is ALSO tribonacci p_5!
  Tribonacci: [1, 3, 7, 11, 21, 39, 71, 131, 241, ...]
  p_1=1 ✓, p_2=3 ✓, p_3=7 ✗, p_4=11 ✓, p_5=21 ✗, p_6=39 ✓, p_7=71 ✓

  So 7 = p_3 (forbidden) and 21 = p_5 (forbidden) are at ODD INDICES 3,5.
  But p_7 = 71 (achievable) is also odd index. So odd-index alone doesn't work.

  Actually: p_3 = 7, p_5 = 21. Both are forbidden. But 31 is NOT p_5!
  31 = Tr(M_5^5) = 2^5-1 (pentanacci trace, not tribonacci).
  In the tribonacci sequence: p_5 = 21, not 31.

  Let me recheck: tribonacci recurrence p_n = p_{n-1} + p_{n-2} + p_{n-3}
  p_1 = 1, p_2 = 3, p_3 = 7
  p_4 = 7+3+1 = 11
  p_5 = 11+7+3 = 21  ← THIS IS 21, not 31!
  p_6 = 21+11+7 = 39
  p_7 = 39+21+11 = 71

  And 31 = Tr(M_5^5) is the PENTANACCI trace, not tribonacci.
  31 is achievable. It's NOT in the tribonacci sequence at all.

CORRECT NARRATIVE:
  The tribonacci traces (k=3) are: 1, 3, 7, 11, 21, 39, 71, 131, ...
  Of these, exactly 7 and 21 are forbidden. The rest are achievable.
  This is NOT a coincidence but also NOT explained by the trace connection:
  - 7 is forbidden because of cycle-forcing (THM-029)
  - 21 is forbidden because of OCF decomposition blocking (THM-079)
  - The k-nacci trace merely HITS these values; it doesn't CAUSE forbiddenness

WHAT THM-227 ACTUALLY SHOWS:
  Tr(M_k^n) = 2^n - 1 for 1 ≤ n ≤ k. This is about k-nacci matrices.
  It does NOT say that 2^k-1 values are forbidden in the H-spectrum.
  The connection "forbidden H values are Mersenne numbers" was a false
  extrapolation from the single data point H=7 = 2^3-1 = Tr(M_3^3).

BEST CHARACTERIZATION OF FORBIDDEN VALUES:
  1. HYP-1231: {7, 21} = {7·3⁰, 7·3¹}. The 7-obstruction has nilpotency 2.
  2. HYP-1317: {7, 21} = {Φ₃(2), Φ₃(4)}. Third cyclotomic at even args.
  3. HYP-1315: Both have I-polynomials factoring through (1+3x).
  4. The KEY mechanism: H=7 requires 3 pairwise-sharing cycles (impossible);
     H=21 requires specific (α₁, i₂) decompositions (all blocked by forcing).
""")
