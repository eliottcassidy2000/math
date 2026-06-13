"""
Search for closed-form Omega_d for Paley tournaments.

Known:
  Omega_0 = 1
  Omega_1 = m
  Omega_2 = m(m-1)
  Omega_3 = m(m-1)(2m-3)/2

For P_7 (m=3): Omega = [1, 3, 6, 9, 9, 6, 3]
For P_11 (m=5): Omega = [1, 5, 20, 70, 205, 460, 700, 690, 450, 180, 30]

Question: is there a pattern?

Try: Omega_d as a polynomial in m, or as a product formula.
Also: the Omega_d are dimensions of per-eigenspace chain groups.
Total chain groups: p * Omega_d (except Omega_0^{total} = p).
"""
from math import comb, factorial
from fractions import Fraction

# Known data
omega_7 = [1, 3, 6, 9, 9, 6, 3]      # P_7, m=3
omega_11 = [1, 5, 20, 70, 205, 460, 700, 690, 450, 180, 30]  # P_11, m=5

print("=== OMEGA DIMENSION ANALYSIS ===\n")

# Check: are these related to falling factorials?
print("P_7 (m=3): Omega vs (m)_d / d! = C(m,d) and vs m^d / d!")
for d, o in enumerate(omega_7):
    cd = comb(3, d) if d <= 3 else 0
    fd = 3**d // factorial(d) if factorial(d) > 0 else 0
    print(f"  d={d}: Omega={o}, C(3,d)={cd}, m^d/d!={fd}")

print(f"\nP_11 (m=5): Omega / m")
for d, o in enumerate(omega_11):
    print(f"  d={d}: Omega={o}, Omega/m={Fraction(o,5)}")

# Check ratios Omega_d / Omega_{d-1}
print(f"\nRatios Omega_d / Omega_{{d-1}}:")
print("P_7:", [Fraction(omega_7[d], omega_7[d-1]) for d in range(1, len(omega_7))])
print("P_11:", [Fraction(omega_11[d], omega_11[d-1]) for d in range(1, len(omega_11))])

# Check: Omega_d = |A_d| - rank(C_d)
# |A_d| = number of d-diff-seqs = (number of d-tuples from QR with distinct partial sums)
# This is a derangement-type count.

# For P_7: |A_d| = [1, 3, 9, 21, 39, 45, 27]
A_7 = [1, 3, 9, 21, 39, 45, 27]
A_11 = [1, 5, 25, 110, 430, 1430, 3970, 8735, 14395, 15745, 8645]

print(f"\nRank(C_d) = |A_d| - Omega_d:")
print("P_7:", [A_7[d] - omega_7[d] for d in range(len(omega_7))])
print("P_11:", [A_11[d] - omega_11[d] for d in range(len(omega_11))])

# Alternating partial sums (boundary rank recursion)
print(f"\n=== BOUNDARY RANK RECURSION ===")
print("R_d = alternating sum of Omega: R_{d+1} = Omega_d - R_d for acyclic degrees")

for name, O, m_val in [("P_7", omega_7, 3), ("P_11", omega_11, 5)]:
    R = [0, 0]  # R_0, R_1
    for d in range(1, m_val):  # d=1,...,m-1
        R.append(O[d] - R[-1])
    print(f"\n{name} (m={m_val}): R from bottom = {R}")

    R_top = [0]  # R_{2m+1} = 0
    for d in range(2*m_val, m_val + 1, -1):  # d=2m,...,m+2
        R_top.insert(0, O[d] - R_top[0])
    print(f"  R from top (d={m_val+2}...{2*m_val+1}): {R_top}")

    R_m = R[-1]  # R_m from bottom
    R_m2 = R_top[0]  # R_{m+2} from top
    print(f"  R_m = {R_m}")
    print(f"  R_{{m+2}} = {R_m2}")
    print(f"  O_m - R_m = {O[m_val] - R_m} (= R_{{m+1}} + β_m)")
    print(f"  O_{{m+1}} - R_{{m+2}} = {O[m_val+1] - R_m2} (= R_{{m+1}} + β_{{m+1}})")
    # Since β_m = β_{m+1}:
    #  O_m - R_m = R_{m+1} + β
    #  O_{m+1} - R_{m+2} = R_{m+1} + β
    # These should be equal!
    val1 = O[m_val] - R_m
    val2 = O[m_val + 1] - R_m2
    print(f"  {val1} == {val2}? {'YES ✓' if val1 == val2 else 'NO ✗'}")

    # The total "budget" for R_{m+1} + β:
    budget = val1
    # And β_m = m(m-3)/2
    beta = m_val * (m_val - 3) // 2
    R_m1 = budget - beta
    print(f"  β_m = m(m-3)/2 = {beta}")
    print(f"  R_{{m+1}} = {budget} - {beta} = {R_m1}")

# Now let me look at the Omega generating function
print(f"\n=== GENERATING FUNCTION ANALYSIS ===")
print("Omega(x) = sum Omega_d * x^d")

for name, O, m_val in [("P_7", omega_7, 3), ("P_11", omega_11, 5)]:
    p = 2*m_val + 1
    print(f"\n{name} (m={m_val}, p={p}):")

    # Omega(1) = sum of all Omega_d
    omega_sum = sum(O)
    print(f"  Omega(1) = {omega_sum}")
    print(f"  Omega(1) / m = {Fraction(omega_sum, m_val)}")

    # Omega(-1) = alternating sum = chi
    chi = sum((-1)**d * o for d, o in enumerate(O))
    print(f"  Omega(-1) = chi = {chi}")

    # Check if Omega is related to (1+x)^m * something
    print(f"  Omega as polynomial coefficients:")
    for d, o in enumerate(O):
        # Try to write as sum of binomial coefficients
        # Omega_d = c_0 * C(m,d) + c_1 * C(m,d-1) + ...?
        pass

    # Check palindromic symmetry: Omega_d vs Omega_{2m-d}
    print(f"  Palindrome check (Omega_d vs Omega_{{2m-d}}):")
    for d in range(m_val + 1):
        print(f"    d={d}: Omega={O[d]}, Omega_{{2m-d}}={O[2*m_val-d]}, ratio={Fraction(O[d], O[2*m_val-d])}")

# Interesting: for P_7, Omega IS palindromic: [1,3,6,9,9,6,3]
# For P_11, Omega is NOT palindromic: [1,5,20,70,205,460,700,690,450,180,30]

print(f"\n=== |A_d| ANALYSIS ===")
# |A_d|/p counts d-diff-seqs divided by p
# Actually |A_d| already counts diff-seqs (not paths)
# |A_d| = number of d-tuples (s_1,...,s_d) in QR^d with partial sums distinct mod p

# For d small: |A_d| = m * (m-1) * ... with some correction
# A_0 = 1 (empty tuple)
# A_1 = m (choose s_1 from QR)
# A_2: choose s_1, s_2 from QR with 0, s_1, s_1+s_2 all distinct
#     = m * #{s_2 in QR : s_1+s_2 ≠ 0 and s_1+s_2 ≠ s_1}
#     = m * #{s_2 in QR : s_2 ≠ -s_1 and s_2 ≠ 0}
#   But s_2 ∈ QR so s_2 ≠ 0 always. And s_2 = -s_1 means s_2 = p-s_1.
#   If p-s_1 ∈ QR: m-1 choices. If p-s_1 ∉ QR: m choices.
#   For p ≡ 3 mod 4: -1 is NQR, so -s_1 = (-1)·s_1. If s_1 ∈ QR, then -s_1 ∈ NQR.
#   So p-s_1 ∉ QR for all s_1 ∈ QR. Hence A_2 = m * m = m².

print("For p ≡ 3 mod 4: -1 is NQR, so -QR = NQR.")
print("A_2 = m * m = m² (since s_2 = -s_1 is NQR, always allowed)")
for name, A, m_val in [("P_7", A_7, 3), ("P_11", A_11, 5)]:
    print(f"  {name}: A_2 = {A[2]}, m² = {m_val**2}. Match: {A[2] == m_val**2}")

# A_3: need 0, s_1, s_1+s_2, s_1+s_2+s_3 all distinct, all s_i in QR
# Forbidden: s_3 = -(s_1+s_2) [makes PS=0], s_3 = -s_2 [makes PS=s_1],
#            s_3 = 0 [impossible since QR doesn't contain 0]
# s_3 = -(s_1+s_2): in QR iff (s_1+s_2) has QR negation. Since -1 is NQR: -(s_1+s_2) ∈ NQR iff s_1+s_2 ∈ QR.
# s_3 = -s_2: always NQR (since s_2 ∈ QR and -1 is NQR).
# So: if s_1+s_2 ∈ QR, then -(s_1+s_2) ∈ NQR, no forbidden QR values from first condition.
#     if s_1+s_2 ∈ NQR, then -(s_1+s_2) ∈ QR, one forbidden QR value.
# Number of valid s_3: m (if s_1+s_2 ∈ QR) or m-1 (if s_1+s_2 ∈ NQR).
# But s_1+s_2 could also be 0: then -(s_1+s_2) = 0 ∉ QR, so m valid choices.
# Actually s_1+s_2 = 0 iff s_2 = -s_1, which means s_2 ∈ NQR (since -1 is NQR). But s_2 ∈ QR. Contradiction.
# So s_1+s_2 ≠ 0 always.

# Among A_2 = m² pairs (s_1,s_2): how many have s_1+s_2 ∈ QR?
# This is N_QR(QR+QR) = number of (a,b) ∈ QR×QR with a+b ∈ QR.
# By Jacobi sums: this equals m(m-1)/2 - 1/4 * (some Jacobi sum)
# Actually the exact count of (a,b) ∈ QR² with a+b = r depends on r.
# N(r ∈ QR) = (m-1)/2, N(r ∈ NQR) = (m+1)/2 for p ≡ 3 mod 4.
# Wait those are wrong signs. Let me use the exact Jacobi sum.
# For r ≠ 0: #{(a,b) ∈ QR² : a+b = r} = (p - 4χ(r) + 1)/4 where χ is Legendre symbol
# For p ≡ 3 mod 4: if r ∈ QR, N(r) = (p-3)/4 = (m-1)/2; if r ∈ NQR, N(r) = (p+1)/4 = (m+1)/2.

# Total (a,b) with a+b ∈ QR = sum_{r ∈ QR} N(r) = m * (m-1)/2
# Total (a,b) with a+b ∈ NQR = sum_{r ∈ NQR} N(r) = m * (m+1)/2

# So A_3 = m * (m-1)/2 * m + m * (m+1)/2 * (m-1) = m²(m-1)/2 + m(m+1)(m-1)/2
#        = m(m-1)/2 * (m + m + 1) = m(m-1)(2m+1)/2

# Wait let me redo this. Among A_2 = m² pairs:
# - m(m-1)/2 have s_1+s_2 ∈ QR, giving m choices for s_3 (none forbidden in QR)
#   Wait but we also need s_1+s_2+s_3 ≠ s_1, i.e., s_2+s_3 ≠ 0, i.e., s_3 ≠ -s_2.
#   Since -s_2 ∈ NQR, this is auto-satisfied.
# - m(m+1)/2 have s_1+s_2 ∈ NQR, giving m-1 choices for s_3
#   (s_3 = -(s_1+s_2) ∈ QR is forbidden since it would make PS_3 = 0)
#   Also need s_3 ≠ -s_2 (auto-satisfied).

# Hold on, we also need s_1+s_2+s_3 ≠ s_1, which means s_2+s_3 ≠ 0.
# But s_3 ≠ -s_2 is always true since -s_2 ∈ NQR.
# And s_1+s_2+s_3 ≠ 0 means s_3 ≠ -(s_1+s_2).
# If s_1+s_2 ∈ QR: -(s_1+s_2) ∈ NQR, so s_3 ≠ -(s_1+s_2) is auto-satisfied.
# If s_1+s_2 ∈ NQR: -(s_1+s_2) ∈ QR, forbidden, losing 1 choice.

# BUT WAIT: we also need the new partial sum s_1+s_2+s_3 to differ from
# all PREVIOUS partial sums: 0, s_1, s_1+s_2.
# s_1+s_2+s_3 ≠ 0: s_3 ≠ -(s_1+s_2) [handled above]
# s_1+s_2+s_3 ≠ s_1: s_3 ≠ -s_2 [handled above]
# s_1+s_2+s_3 ≠ s_1+s_2: s_3 ≠ 0 [always true]

# So the count is correct:
# A_3 = m*(m-1)/2 * m + m*(m+1)/2 * (m-1)

# Wait, the count of pairs with sum in QR vs NQR needs to exclude those where
# the PAIR itself fails (s_1 = s_2 mod p is impossible since QR allows repeats in diff-seqs).
# Actually diff-seqs allow repeated entries! Only partial sums must be distinct.

# Hmm, I had A_2 = m². But some of these m² pairs might have repeated s values.
# That's fine — we allow (s_1, s_2) with s_1 = s_2 as long as partial sums are distinct.
# s_1 = s_2 means PS = 0, s_1, 2s_1. All distinct iff s_1 ≠ 0 (always) and 2s_1 ≠ 0 (true for p>2).

# So the count of A_2 pairs with s_1+s_2 ∈ QR:
# = #{(s_1,s_2) ∈ QR² : s_1+s_2 ∈ QR}
# = sum_{r ∈ QR} #{(s_1,s_2) ∈ QR² : s_1+s_2 = r}
# = sum_{r ∈ QR} (m-1)/2 = m * (m-1)/2

# Count with s_1+s_2 ∈ NQR:
# = m² - m*(m-1)/2 = m(m+1)/2

# So A_3 = m*(m-1)/2 * m + m*(m+1)/2 * (m-1)
#        = m²(m-1)/2 + m(m-1)(m+1)/2 = m(m-1)(m + m + 1)/2 = m(m-1)(2m+1)/2

print("\nA_3 formula: m(m-1)(2m+1)/2")
for name, A, m_val in [("P_7", A_7, 3), ("P_11", A_11, 5)]:
    pred = m_val * (m_val - 1) * (2*m_val + 1) // 2
    print(f"  {name}: A_3 = {A[3]}, predicted = {pred}. Match: {A[3] == pred}")

# Omega_3 = A_3 - rank(C_3). rank(C_3) = 2m(m-1).
# Omega_3 = m(m-1)(2m+1)/2 - 2m(m-1) = m(m-1)(2m+1-4)/2 = m(m-1)(2m-3)/2. ✓

# Let me compute A_4 formula.
# At degree 4: need (s_1,s_2,s_3,s_4) ∈ QR^4 with 0, PS_1,...,PS_4 distinct.
# This gets complicated. Let me just verify numerically.

def count_A(p, d):
    """Count d-diff-seqs for Paley P_p."""
    QR = set(pow(x, 2, p) for x in range(1, p))
    QR_list = sorted(QR)
    count = [0]
    if d == 0:
        return 1
    def backtrack(seq, ps_list, ps_set, depth):
        if depth == d:
            count[0] += 1
            return
        for s in QR_list:
            new_ps = (ps_list[-1] + s) % p
            if new_ps in ps_set:
                continue
            ps_list.append(new_ps)
            ps_set.add(new_ps)
            backtrack(seq + [s], ps_list, ps_set, depth + 1)
            ps_list.pop()
            ps_set.remove(new_ps)
    backtrack([], [0], {0}, 0)
    return count[0]

# Verify A_d formulas
print("\n=== |A_d| FORMULAS ===")
for p in [7, 11]:
    m = (p - 1) // 2
    print(f"\nP_{p} (m={m}):")
    for d in range(min(5, 2*m+1)):
        Ad = count_A(p, d)
        # Known formulas:
        if d == 0:
            f = 1
        elif d == 1:
            f = m
        elif d == 2:
            f = m * m
        elif d == 3:
            f = m * (m-1) * (2*m+1) // 2
        else:
            f = None
        match = f"= {f} {'✓' if f == Ad else '✗'}" if f is not None else ""
        print(f"  d={d}: |A_d| = {Ad} {match}")

# Now try to find A_4 formula
# At d=4: for each (s_1,s_2,s_3) ∈ A_3:
#   Need s_4 ∈ QR with PS_4 = PS_3 + s_4 distinct from 0, PS_1, PS_2, PS_3
#   i.e., s_4 ≠ -PS_3, s_4 ≠ PS_1 - PS_3, s_4 ≠ PS_2 - PS_3, s_4 ≠ 0
#   s_4 ∈ QR, s_4 ≠ 0 always.
#   Need to count how many of the 3 forbidden values are in QR.

# The forbidden values are: -PS_3 = -(s_1+s_2+s_3), s_1-(s_1+s_2+s_3)=-(s_2+s_3), (s_1+s_2)-(s_1+s_2+s_3)=-s_3
# = {-(s_1+s_2+s_3), -(s_2+s_3), -s_3}
# Since -1 ∈ NQR: -x ∈ QR iff x ∈ NQR.
# So forbidden QR values = #{x ∈ {s_1+s_2+s_3, s_2+s_3, s_3} : x ∈ NQR}

# For PS_3 = s_1+s_2+s_3: this is the last partial sum
# For s_2+s_3 = PS_3 - s_1: this is PS_3 minus PS_1
# For s_3: this is in QR always

# So: s_3 ∈ QR means -s_3 ∈ NQR → not in QR → NOT a forbidden value.
# Number of forbidden QR values from {-(s_1+s_2+s_3), -(s_2+s_3)}:
#   Each of these is in QR iff the corresponding value (s_1+s_2+s_3 or s_2+s_3) is in NQR.

# Let α = s_1+s_2+s_3, β_val = s_2+s_3.
# Need: α ∈ QR or NQR, β_val ∈ QR or NQR.
# Four cases:
# α ∈ QR, β_val ∈ QR: 0 forbidden → m valid
# α ∈ QR, β_val ∈ NQR: 1 forbidden → m-1
# α ∈ NQR, β_val ∈ QR: 1 forbidden → m-1
# α ∈ NQR, β_val ∈ NQR: 2 forbidden → m-2

# Count the four cases among A_3 = m(m-1)(2m+1)/2 triples:
# This requires character sum analysis...

print(f"\n=== DEGREE 4: CASE ANALYSIS ===")
for p in [7, 11]:
    m = (p - 1) // 2
    QR = set(pow(x, 2, p) for x in range(1, p))
    QR_list = sorted(QR)

    cases = [[0]*4 for _ in range(1)]  # [QQ, QN, NQ, NN] where Q=QR, N=NQR
    # First letter: α = s_1+s_2+s_3, Second: β = s_2+s_3
    case_count = {'QQ': 0, 'QN': 0, 'NQ': 0, 'NN': 0}

    def backtrack3(seq, ps_list, ps_set, depth):
        if depth == 3:
            alpha = sum(seq) % p
            beta_val = (seq[1] + seq[2]) % p
            a_qr = 'Q' if alpha in QR else 'N'
            b_qr = 'Q' if beta_val in QR else 'N'
            case_count[a_qr + b_qr] += 1
            return
        for s in QR_list:
            new_ps = (ps_list[-1] + s) % p
            if new_ps in ps_set:
                continue
            ps_list.append(new_ps)
            ps_set.add(new_ps)
            backtrack3(seq + [s], ps_list, ps_set, depth + 1)
            ps_list.pop()
            ps_set.remove(new_ps)

    backtrack3([], [0], {0}, 0)
    print(f"\nP_{p} (m={m}): A_3 case split:")
    total_A4 = 0
    for case, cnt in case_count.items():
        forbidden = case.count('N')
        valid = m - forbidden
        total_A4 += cnt * valid
        print(f"  {case}: {cnt} triples, {forbidden} forbidden, {valid} valid → {cnt * valid}")
    print(f"  Total A_4 = {total_A4}")
    print(f"  Actual A_4 = {count_A(p, 4)}")
    print(f"  Match: {total_A4 == count_A(p, 4)}")

    # Also: A_3 = QQ + QN + NQ + NN
    A3 = sum(case_count.values())
    print(f"  A_3 = {A3} = {case_count}")

    # Check formulas for each case using Jacobi sums
    # QQ = #{triples : PS_3 ∈ QR, PS_2 ∈ QR}  (where PS_2 = s_2+s_3, oops PS_2 = s_1+s_2)
    # Hmm, I defined β_val = s_2+s_3, which is NOT a partial sum. Let me reconsider.
    # α = PS_3 = s_1+s_2+s_3 (partial sum 3)
    # β_val = s_2+s_3 = PS_3 - s_1 = PS_3 - PS_1
    # PS_2 = s_1+s_2 (this is a partial sum that appears in the distinct constraint)

    # The QR/NQR status of α and β_val depends on s_1, PS_2, PS_3.
    # α = PS_3 is a partial sum (already constrained to be distinct from 0,PS_1,PS_2)
    # β_val = PS_3 - PS_1 is the difference of two partial sums
    # Note: PS_3 could be ∈ QR or NQR (no constraint on which)
    # β_val = PS_3 - PS_1 is a difference of partial sums

    # This is getting complex. Let me just check the formula A_4 = ?
    # A_4 = m*QQ + (m-1)*QN + (m-1)*NQ + (m-2)*NN
    # = m*A_3 - QN - NQ - 2*NN
    # = m*A_3 - (A_3 - QQ) - NN
    # Hmm, not simplifying nicely.

print("\n=== A_4 DIRECT FORMULA SEARCH ===")
for p in [7, 11, 19, 23]:
    m = (p - 1) // 2
    if p > 19:
        continue
    A4 = count_A(p, 4)
    A3 = count_A(p, 3)
    r1 = Fraction(A4, m)
    r2 = Fraction(A4, m*(m-1))
    r3 = Fraction(A4, A3)
    print(f"  P_{p} (m={m}): A_4 = {A4}, A_4/m = {r1}, A_4/m(m-1) = {r2}, A_4/A_3 = {r3}")

# Also compute A_d for P_19 as far as feasible
print("\n=== P_19 (m=9) diff-seq counts ===")
p = 19
m = 9
for d in range(8):
    Ad = count_A(p, d)
    print(f"  d={d}: |A_d| = {Ad}", flush=True)
