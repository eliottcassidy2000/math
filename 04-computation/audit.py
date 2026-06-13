"""
audit.py -- kind-pasteur-2026-03-14-S108
SYSTEMATIC AUDIT of all claims from S105-S107

Go through every major claim and verify it computationally.
Flag anything suspicious, incorrect, or that needs caveats.
"""

import sys, math
import numpy as np
from fractions import Fraction
from itertools import combinations
from collections import Counter

sys.stdout.reconfigure(encoding='utf-8')

def count_ham_paths(adj, n):
    dp = {}
    for v in range(n):
        dp[(1 << v, v)] = 1
    for mask in range(1, 1 << n):
        for v in range(n):
            if not (mask & (1 << v)):
                continue
            if (mask, v) not in dp:
                continue
            for u in range(n):
                if mask & (1 << u):
                    continue
                if adj[v][u]:
                    key = (mask | (1 << u), u)
                    dp[key] = dp.get(key, 0) + dp[(mask, v)]
    full = (1 << n) - 1
    return sum(dp.get((full, v), 0) for v in range(n))

def main():
    print("=" * 70)
    print("SYSTEMATIC AUDIT — VERIFY ALL CLAIMS")
    print("kind-pasteur-2026-03-14-S108")
    print("=" * 70)

    errors = []
    warnings = []
    confirmed = []

    # ============================================================
    print(f"\n{'='*70}")
    print("CLAIM 1: Grand Energy Formula E_2k/E_0 = 2*(n-2k)^k/P(n,2k)")
    print(f"{'='*70}")

    # Verified at n=3..7 exactly. Let me re-verify independently.
    for n in [3, 4, 5, 6]:
        m = n*(n-1)//2
        h_vals = []
        for bits in range(1 << m):
            adj = [[0]*n for _ in range(n)]
            idx = 0
            for i in range(n):
                for j in range(i+1, n):
                    if bits & (1 << idx):
                        adj[i][j] = 1
                    else:
                        adj[j][i] = 1
                    idx += 1
            h_vals.append(count_ham_paths(adj, n))

        N = len(h_vals)
        mean = Fraction(sum(h_vals), N)
        var = Fraction(sum(h*h for h in h_vals), N) - mean*mean
        ratio_exact = var / (mean * mean)

        # Grand formula prediction
        pred = Fraction(0)
        for k in range(1, 100):
            if n - 2*k <= 0:
                break
            pred += Fraction(2*(n-2*k)**k, math.perm(n, 2*k))

        match = (ratio_exact == pred)
        status = "PASS" if match else "FAIL"
        if match:
            confirmed.append(f"Grand formula at n={n}: {ratio_exact} = {pred}")
        else:
            errors.append(f"Grand formula FAILS at n={n}: exact={ratio_exact}, pred={pred}")
        print(f"  n={n}: exact Var/Mean^2 = {ratio_exact}, formula = {pred}, {status}")

    # n=7 from stored computation
    n7_exact = Fraction(131, 504)
    n7_pred = Fraction(0)
    for k in range(1, 100):
        if 7 - 2*k <= 0:
            break
        n7_pred += Fraction(2*(7-2*k)**k, math.perm(7, 2*k))
    match7 = (n7_exact == n7_pred)
    status7 = "PASS" if match7 else "FAIL"
    print(f"  n=7: stored exact = {n7_exact}, formula = {n7_pred}, {status7}")
    if match7:
        confirmed.append("Grand formula at n=7: 131/504")
    else:
        errors.append("Grand formula FAILS at n=7")

    # ============================================================
    print(f"\n{'='*70}")
    print("CLAIM 2: Newton power sums S_3=7, S_5=21, S_8=131")
    print(f"{'='*70}")

    # The Newton sums of x^3-x^2-x-1 satisfy S_k = S_{k-1}+S_{k-2}+S_{k-3}
    # with S_0=3, S_1=1, S_2=3.
    S = [3, 1, 3]
    for i in range(3, 20):
        S.append(S[-1] + S[-2] + S[-3])

    checks = [(3, 7), (5, 21), (8, 131)]
    for k, expected in checks:
        actual = S[k]
        match = (actual == expected)
        status = "PASS" if match else "FAIL"
        print(f"  S_{k} = {actual}, expected {expected}, {status}")
        if match:
            confirmed.append(f"S_{k} = {expected}")
        else:
            errors.append(f"S_{k} = {actual} != {expected}")

    # Verify S_0, S_1, S_2 from Vieta's formulas for x^3-x^2-x-1
    # Roots: tau, sigma, sigma_bar
    # Sum of roots = 1 (coeff of x^2)
    # Sum of products of pairs = -1 (coeff of x)
    # Product of roots = 1 (constant term, with sign)
    print(f"\n  Vieta check for x^3 - x^2 - x - 1:")
    print(f"    S_0 = 3 (degree of polynomial)")
    print(f"    S_1 = sum of roots = 1 (= coeff of x^2)")
    print(f"    S_2 = S_1^2 - 2*(sum of products of pairs)")
    print(f"         = 1^2 - 2*(-1) = 1 + 2 = 3")
    s2_check = 1 - 2*(-1)
    if s2_check == 3:
        confirmed.append("S_2 = 3 from Vieta")
        print(f"    PASS: S_2 = {s2_check} = 3")
    else:
        errors.append(f"S_2 Vieta check failed: got {s2_check}")
        print(f"    FAIL")

    # ============================================================
    print(f"\n{'='*70}")
    print("CLAIM 3: S_k = 3*T_k + 4*T_{k-1} + T_{k-2}")
    print(f"{'='*70}")

    T = [0, 0, 1]
    for i in range(3, 20):
        T.append(T[-1] + T[-2] + T[-3])

    all_match = True
    for k in range(2, 15):
        pred = 3*T[k] + 4*T[k-1] + T[k-2]
        actual = S[k]
        if pred != actual:
            all_match = False
            errors.append(f"S_k = 3*T_k+4*T_(k-1)+T_(k-2) fails at k={k}: {pred} != {actual}")
            print(f"  k={k}: FAIL, pred={pred}, actual={actual}")

    if all_match:
        confirmed.append("S_k = 3*T_k + 4*T_{k-1} + T_{k-2} for k=2..14")
        print(f"  ALL PASS for k=2..14")

    # ============================================================
    print(f"\n{'='*70}")
    print("CLAIM 4: tau^3 = Phi_3(tau)")
    print(f"{'='*70}")

    tau = 1.8392867552141612
    phi3_tau = tau**2 + tau + 1
    tau_cubed = tau**3
    diff = abs(phi3_tau - tau_cubed)
    if diff < 1e-10:
        confirmed.append("tau^3 = Phi_3(tau)")
        print(f"  PASS: Phi_3(tau) = {phi3_tau:.10f}, tau^3 = {tau_cubed:.10f}, diff = {diff:.2e}")
    else:
        errors.append(f"tau^3 != Phi_3(tau), diff = {diff}")
        print(f"  FAIL: diff = {diff}")

    # ============================================================
    print(f"\n{'='*70}")
    print("CLAIM 5: E_2/Var -> 1 as n -> inf")
    print(f"{'='*70}")

    # E_2/Var should approach (n-2)/n? Let me check.
    # From the formula: E_2/E_0 = 2(n-2)/(n(n-1))
    # Var/E_0 = sum_k terms. The ratio E_2/Var = E_2/E_0 / (Var/E_0).
    # Let me compute this exactly.
    print(f"  E_2/Var for various n:")
    for n in [3, 5, 7, 10, 20, 50, 100]:
        e2 = Fraction(2*(n-2), n*(n-1))
        total = Fraction(0)
        for k in range(1, 100):
            if n - 2*k <= 0:
                break
            total += Fraction(2*(n-2*k)**k, math.perm(n, 2*k))
        ratio = e2 / total if total > 0 else Fraction(1)
        print(f"    n={n:4d}: E_2/Var = {float(ratio):.6f}")

    # WARNING: I claimed E_2/Var = (n-2)/n in the earlier script,
    # but the table showed it's NOT exactly (n-2)/n.
    # Let me check more carefully.
    print(f"\n  Is E_2/Var = (n-2)/n exactly?")
    for n in [5, 7, 10]:
        e2 = Fraction(2*(n-2), n*(n-1))
        total = Fraction(0)
        for k in range(1, 100):
            if n - 2*k <= 0:
                break
            total += Fraction(2*(n-2*k)**k, math.perm(n, 2*k))
        ratio = e2 / total
        pred = Fraction(n-2, n)
        print(f"    n={n}: E_2/Var = {ratio} = {float(ratio):.8f}, (n-2)/n = {pred} = {float(pred):.8f}, exact match: {ratio == pred}")

    # The earlier output showed they DON'T match!
    # E_2/Var at n=5 was 0.947368, but (n-2)/n = 3/5 = 0.6. Very different!
    warnings.append("CLAIM E_2/Var = (n-2)/n was stated in missing_insight.py but is WRONG. E_2/Var -> 1 is true but != (n-2)/n.")
    print(f"\n  WARNING: E_2/Var is NOT exactly (n-2)/n!")
    print(f"  The claim in missing_insight.py was INCORRECTLY stated.")
    print(f"  E_2/Var -> 1 is TRUE (monotonic from numerical evidence)")
    print(f"  But E_2/Var != (n-2)/n exactly.")

    # ============================================================
    print(f"\n{'='*70}")
    print("CLAIM 6: 360 = T(12) - F(12)")
    print(f"{'='*70}")

    fib = [0, 1]
    for i in range(2, 15):
        fib.append(fib[-1] + fib[-2])

    t12 = T[12]  # should be 504? Wait, indexing.
    # T uses 0-indexing: T[0]=0, T[1]=0, T[2]=1, T[3]=1, T[4]=2, T[5]=4,
    # T[6]=7, T[7]=13, T[8]=24, T[9]=44, T[10]=81, T[11]=149, T[12]=274, T[13]=504
    # Hmm, the earlier script said T(12) = 504. Let me check indexing.
    print(f"  T[12] = {T[12]}, T[13] = {T[13]}")
    print(f"  F[12] = {fib[12]}")

    # With 0-indexing: T[13]=504, F[12]=144. T[13]-F[12] = 504-144 = 360.
    # But the claim said T(12)-F(12) = 360. This depends on indexing convention!
    # In the earlier script, tribonacci used 1-indexing where T(1)=T(2)=0, T(3)=1.
    # So T(12) in 1-indexing = T[11] in 0-indexing = 149.
    # T(12)_1indexed = 149, F(12)_1indexed = 144. 149-144 = 5. NOT 360!

    # Let me check what the earlier script actually computed.
    # It said: TRIB = [0, 1, 1, 2, 4, 7, 13, 24, 44, 81, 149, 274, 504, 927, 1705]
    # This is 0-indexed with TRIB[0]=0. So TRIB[12] = 504.
    TRIB = [0, 1, 1, 2, 4, 7, 13, 24, 44, 81, 149, 274, 504, 927, 1705]
    FIB = [0, 1, 1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377]
    print(f"  TRIB[12] = {TRIB[12]}, FIB[12] = {FIB[12]}")
    print(f"  TRIB[12] - FIB[12] = {TRIB[12] - FIB[12]}")
    if TRIB[12] - FIB[12] == 360:
        confirmed.append("TRIB[12] - FIB[12] = 504 - 144 = 360")
        print(f"  PASS")
    else:
        errors.append(f"TRIB[12] - FIB[12] = {TRIB[12] - FIB[12]} != 360")
        print(f"  FAIL")

    # ============================================================
    print(f"\n{'='*70}")
    print("CLAIM 7: Forbidden values are base-6 repdigits (7='11', 21='33')")
    print(f"{'='*70}")

    def to_base(n, b):
        if n == 0: return '0'
        digits = []
        while n > 0:
            digits.append(str(n % b))
            n //= b
        return ''.join(reversed(digits))

    b6_7 = to_base(7, 6)
    b6_21 = to_base(21, 6)
    print(f"  7 in base 6: {b6_7}")
    print(f"  21 in base 6: {b6_21}")
    if b6_7 == '11' and b6_21 == '33':
        confirmed.append("7='11'_6, 21='33'_6 (repdigits)")
        print(f"  PASS")
    else:
        errors.append(f"Base-6: 7={b6_7}, 21={b6_21}")
        print(f"  FAIL")

    # ============================================================
    print(f"\n{'='*70}")
    print("CLAIM 8: Spectator freedom proof sketch")
    print(f"{'='*70}")

    print(f"""
  The "spectator freedom" argument says:
  E_2k/E_0 = 2*(n-2k)^k/P(n,2k) because fixing 2k arcs creates
  k interaction points each with (n-2k) free spectator choices.

  CAVEAT: This is a HEURISTIC, not a rigorous proof.
  The notion of "interaction point" and "spectator" has not been
  formally defined. The formula is VERIFIED but not PROVED.

  The claim that there are exactly k "interaction points" from 2k arcs
  needs justification. Why k and not 2k or k-1?

  POSSIBLE ISSUE: At level 4, n=6, the 90 full-coverage coefficients
  have |H_hat| = 1/4, but the 360 inherited ones have |H_hat| = 1/8.
  The spectator argument would need to explain this heterogeneity.
  The formula works for the TOTAL energy but individual coefficients
  vary in magnitude. The "spectator freedom" gives the AVERAGE, not
  individual values.

  STATUS: The proof sketch is SUGGESTIVE but not rigorous.
  The formula is empirically confirmed but lacks formal proof.""")
    warnings.append("Spectator freedom argument is heuristic, not a proof. Grand formula needs rigorous verification.")

    # ============================================================
    print(f"\n{'='*70}")
    print("CLAIM 9: Var(log H) grows (from missing_insight.py)")
    print(f"{'='*70}")

    for n in [3, 4, 5, 6]:
        m = n*(n-1)//2
        h_vals = []
        for bits in range(1 << m):
            adj = [[0]*n for _ in range(n)]
            idx = 0
            for i in range(n):
                for j in range(i+1, n):
                    if bits & (1 << idx):
                        adj[i][j] = 1
                    else:
                        adj[j][i] = 1
                    idx += 1
            h_vals.append(count_ham_paths(adj, n))
        h_arr = np.array(h_vals, dtype=float)
        var_log = np.var(np.log(h_arr))
        print(f"  n={n}: Var(log H) = {var_log:.6f}")

    # Var(log H) at n=3: 0.226, n=4: 0.502, n=5: 0.644, n=6: 0.640
    # It GROWS from n=3 to n=5, then DECREASES slightly at n=6!
    # So the claim "Var(log H) grows" is not clearly true for all n.
    warnings.append("Var(log H) grows n=3 to n=5 but decreases at n=6. Claim needs qualification.")
    print(f"  WARNING: Var(log H) DECREASES from n=5 to n=6!")
    print(f"  The claim 'Var(log H) grows' is only true for n=3..5.")

    # ============================================================
    print(f"\n{'='*70}")
    print("CLAIM 10: 189/360 = 21/40")
    print(f"{'='*70}")

    ratio = Fraction(189, 360)
    print(f"  189/360 = {ratio} = {float(ratio):.6f}")
    if ratio == Fraction(21, 40):
        confirmed.append("189/360 = 21/40")
        print(f"  PASS")
    else:
        errors.append(f"189/360 = {ratio} != 21/40")
        print(f"  FAIL")

    # ============================================================
    print(f"\n{'='*70}")
    print("CLAIM 11: 1+2+3 = 1*2*3 is unique")
    print(f"{'='*70}")

    found = []
    for a in range(1, 20):
        for b in range(a, 20):
            for c in range(b, 20):
                if a + b + c == a * b * c:
                    found.append((a, b, c))
    print(f"  Triples with sum=product: {found}")
    if found == [(1, 2, 3)]:
        confirmed.append("(1,2,3) is the unique triple with sum=product")
        print(f"  PASS: unique")
    else:
        warnings.append(f"Multiple triples found: {found}")

    # ============================================================
    print(f"\n{'='*70}")
    print("CLAIM 12: 504 is a tribonacci number")
    print(f"{'='*70}")

    if 504 in TRIB:
        confirmed.append("504 is tribonacci (TRIB[12])")
        print(f"  PASS: 504 = TRIB[{TRIB.index(504)}]")
    else:
        errors.append("504 is NOT in the tribonacci sequence")
        print(f"  FAIL")

    # ============================================================
    print(f"\n{'='*70}")
    print("CLAIM 13: Phi_3(5) = Phi_5(2) = 31")
    print(f"{'='*70}")

    phi3_5 = 5**2 + 5 + 1
    phi5_2 = 2**4 + 2**3 + 2**2 + 2 + 1
    print(f"  Phi_3(5) = {phi3_5}")
    print(f"  Phi_5(2) = {phi5_2}")
    if phi3_5 == 31 and phi5_2 == 31:
        confirmed.append("Phi_3(5) = Phi_5(2) = 31")
        print(f"  PASS")
    else:
        errors.append(f"Phi_3(5)={phi3_5}, Phi_5(2)={phi5_2}")
        print(f"  FAIL")

    # ============================================================
    print(f"\n{'='*70}")
    print("SUMMARY")
    print(f"{'='*70}")

    print(f"\n  CONFIRMED ({len(confirmed)}):")
    for c in confirmed:
        print(f"    + {c}")

    print(f"\n  WARNINGS ({len(warnings)}):")
    for w in warnings:
        print(f"    ! {w}")

    print(f"\n  ERRORS ({len(errors)}):")
    for e in errors:
        print(f"    X {e}")

    if not errors:
        print(f"\n  NO ERRORS FOUND. All numerical claims verified.")
    else:
        print(f"\n  {len(errors)} ERROR(S) FOUND. Investigate immediately.")

    print(f"\n{'='*70}")
    print("DONE — AUDIT COMPLETE")
    print(f"{'='*70}")

if __name__ == '__main__':
    main()
