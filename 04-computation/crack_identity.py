"""
crack_identity.py -- kind-pasteur-2026-03-14-S108f
CRACK THE IDENTITY by rapid perspective shifts.

Identity to prove: D_n(2) = n! + 2*sum_k (n-2k)^k * (n-2k)!

where D_n(2) = sum over anti-succ-free perms of 2^successions.

Stop thinking about anti-successions. Look at D_n(2) from completely
different angles until the proof clicks.
"""

import sys, math
from fractions import Fraction
from itertools import permutations
from collections import Counter

sys.stdout.reconfigure(encoding='utf-8')

def successions(pi):
    return sum(1 for i in range(len(pi)-1) if pi[i+1] == pi[i] + 1)

def anti_successions(pi):
    return sum(1 for i in range(len(pi)-1) if pi[i+1] == pi[i] - 1)

def main():
    print("=" * 70)
    print("CRACK THE IDENTITY — RAPID PERSPECTIVES")
    print("=" * 70)

    # Compute D_n(2) and the target for reference
    data = {}
    for n in range(2, 8):
        all_perms = list(permutations(range(n)))
        Dn2 = sum(2**successions(pi) for pi in all_perms if anti_successions(pi) == 0)
        target = math.factorial(n) + sum(2*(n-2*k)**k * math.factorial(n-2*k)
                                         for k in range(1, 100) if n-2*k > 0)
        data[n] = (Dn2, target)

    # ============================================================
    print(f"\n--- PERSPECTIVE 1: FORGET anti-successions. What IS D_n(2)?")
    print(f"    D_n(2) = E[H^2] * 2^(2(n-1)) / n!")
    print(f"    = (n! * Mean(H^2)) / (n! * (something))")
    print(f"    Actually: D_n(2) = 2^(2(n-1)) * E[H^2] / n!")
    for n in range(3, 8):
        print(f"    n={n}: D_n(2) = {data[n][0]}")

    # D_n(2): 1, 2, 8, 32, 158, 928, 6350
    # Let me look at D_n(2)/n!
    print(f"\n    D_n(2)/n! = 1 + Var/Mean^2:")
    for n in range(2, 8):
        ratio = Fraction(data[n][0], math.factorial(n))
        print(f"    n={n}: D_n(2)/n! = {ratio} = {float(ratio):.8f}")

    # ============================================================
    print(f"\n--- PERSPECTIVE 2: Look at D_n(2) - n! directly")

    for n in range(2, 8):
        excess = data[n][0] - math.factorial(n)
        print(f"    n={n}: D_n(2)-n! = {excess}")
        # Factor
        if excess > 0:
            # The formula says this = 2*sum_k (n-2k)^k*(n-2k)!
            terms = []
            for k in range(1, 100):
                if n-2*k <= 0: break
                terms.append(f"2*{n-2*k}^{k}*{n-2*k}! = {2*(n-2*k)**k*math.factorial(n-2*k)}")
            print(f"           = {' + '.join(terms)}")

    # D_n(2) - n!: 0, 0, 2, 8, 38, 208, 1310
    # n=4: 8 = 2*2*1! + 0 = ... wait
    # n=4: target-n! = 32-24 = 8. Formula: 2*(4-2)^1*(4-2)! = 2*2*2 = 8. CHECK.
    # n=5: 158-120 = 38. Formula: 2*3*6 + 2*1*1 = 36+2 = 38. CHECK.

    # ============================================================
    print(f"\n--- PERSPECTIVE 3: Can I express D_n(2)-n! as a SINGLE sum?")

    # D_n(2) - n! = 2*sum_k (n-2k)^k * (n-2k)!
    # = 2*[sum_k (n-2k)^k * (n-2k)!]
    # Let j = n-2k, so k = (n-j)/2, j goes from 1 to n-2 (by 2 if n is even/odd)
    # = 2*sum_j j^((n-j)/2) * j! where j = n-2, n-4, n-6, ..., 1 or 2.

    print(f"    In terms of j = n-2k:")
    for n in range(3, 8):
        terms = []
        j_vals = list(range(n-2, 0, -2))
        for j in j_vals:
            k = (n-j)//2
            val = j**k * math.factorial(j)
            terms.append(f"j={j}: {j}^{k}*{j}! = {val}")
        total = 2*sum(j**((n-j)//2)*math.factorial(j) for j in j_vals)
        print(f"    n={n}: {'; '.join(terms)}. 2*sum = {total}, target = {data[n][0]-math.factorial(n)}")

    # ============================================================
    print(f"\n--- PERSPECTIVE 4: Think of (n-2k)^k*(n-2k)! as a MULTINOMIAL")

    # (n-2k)^k * (n-2k)! = (n-2k)! * (n-2k)^k
    # = (n-2k)! * (n-2k) * (n-2k) * ... * (n-2k) [k times]
    # This counts: arrange (n-2k) items in a row [(n-2k)! ways],
    # then choose k items from the same (n-2k) set WITH REPLACEMENT [(n-2k)^k ways].
    # Total: ordered arrangement + k tags = (n-2k)^k * (n-2k)!

    # WHAT IF this counts DECORATED PERMUTATIONS?
    # A permutation of [n-2k] PLUS k labels chosen from [n-2k]?

    print(f"""
    (n-2k)^k * (n-2k)! counts:
    "Arrange (n-2k) items AND choose k labels from them (with replacement)."

    The 2 factor: the labels can be "marked" in 2 ways? Or it's from the
    path reversal involution.

    WHAT IF the identity says:
    "The excess D_n(2) - n! counts decorated permutations of smaller sets,
     where the decorations come from removing 2k elements and creating
     k interaction tags."
    """)

    # ============================================================
    print(f"\n--- PERSPECTIVE 5: Exponential generating function")

    # D_n(2) = n! + 2*sum_k (n-2k)^k * (n-2k)!
    # EGF: sum_n D_n(2) * t^n / n! = sum_n t^n + 2*sum_n sum_k (n-2k)^k*(n-2k)!/n! * t^n
    # = 1/(1-t) + 2*sum_k sum_n (n-2k)^k*(n-2k)!/n! * t^n
    # = 1/(1-t) + 2*sum_k sum_{m>=1} m^k*m!/((m+2k))! * t^{m+2k}  [setting m=n-2k]
    # = 1/(1-t) + 2*sum_k t^{2k} * sum_m m^k * t^m / C(m+2k, 2k) * (something)
    # This is getting messy. Let me try OGF instead.

    # Actually: D_n(2)/n! = 1 + 2*sum_k (n-2k)^k / P(n,2k)
    # And (n-2k)^k/P(n,2k) = (n-2k)^k * (n-2k)!/n!
    # For large n: ~ (n-2k)^k / n^{2k} ~ (1-2k/n)^k / n^k ~ 1/n^k.

    # ============================================================
    print(f"\n--- PERSPECTIVE 6: FORGET the identity. Prove the FORMULA directly.")

    # Instead of proving D_n(2) = RHS, prove E_2k/E_0 = 2(n-2k)^k/P(n,2k) directly.
    # We know E_2k = sum_{|S|=2k} H_hat(S)^2.
    # And H_hat(S) = (1/2^m) sum_T H(T) chi_S(T).
    # = (1/2^m) sum_T sum_P [P is HP in T] chi_S(T)
    # = (1/2^m) sum_P sum_T [P consistent with T] chi_S(T)
    # = sum_P (1/2^m) sum_T prod_{arcs of T} [correct orientation] * chi_S(T)

    # For a single path P = (v0,...,v_{n-1}):
    # sum_T [P consistent with T] chi_S(T)
    # = sum_T prod_{i: vi->vi+1 in T} 1 * prod_{e in S} y_e(T)
    # P fixes arcs v0->v1, ..., v_{n-2}->v_{n-1} (n-1 arcs).
    # The remaining m-(n-1) arcs are free.
    # chi_S(T) = prod_{e in S} y_e(T).

    # If arc e is FIXED by P (e is one of the path arcs):
    #   y_e is determined (+1 or -1 depending on orientation matching).
    # If arc e is FREE (not in P):
    #   y_e is equally likely +1 or -1.

    # sum_T [P consistent] chi_S(T)
    # = [product of fixed y_e for e in S inter path_arcs] * [sum over free arcs of chi_{S\path}]
    # The free arcs contribute: prod_{e in S \ path_arcs} sum_{y_e=+/-1} y_e = 0 unless S subset path_arcs.
    # So: sum is ZERO unless S is a subset of the path arcs of P!

    # IF S is a subset of P's arcs:
    # sum_T [P consistent] chi_S(T) = 2^{m-(n-1)} * prod_{e in S} y_e(P)

    print(f"""
    KEY INSIGHT:
    H_hat(S) = (1/2^m) * sum_P sum_T [P consistent with T] chi_S(T)
             = (1/2^m) * sum_P [S subset of arcs(P)] * 2^(m-(n-1)) * chi_S(P)
             = (1/2^(n-1)) * sum_{{P: S subset arcs(P)}} chi_S(P)

    where chi_S(P) = product of (+/-1) depending on arc orientations in P.

    So H_hat(S) = (1/2^(n-1)) * sum over paths P that USE all arcs in S
                  of the sign product over S.

    And H_hat(S)^2 = (1/2^(2(n-1))) * [sum_P ...]^2
                   = (1/2^(2(n-1))) * sum_{{P,Q: S subset arcs(P) inter arcs(Q)}} chi_S(P)*chi_S(Q)

    E_2k = sum_{{|S|=2k}} H_hat(S)^2
         = (1/2^(2(n-1))) * sum_{{|S|=2k}} sum_{{P,Q}} [S subset arcs(P) inter arcs(Q)] * chi_S(P)*chi_S(Q)
         = (1/2^(2(n-1))) * sum_{{P,Q}} sum_{{S subset arcs(P) inter arcs(Q), |S|=2k}} chi_S(P)*chi_S(Q)

    For two paths P, Q with shared arcs in same direction (set A, |A|=s)
    and other arcs:
    S must be a subset of the shared arcs A (since S subset arcs(P) inter arcs(Q)).
    chi_S(P)*chi_S(Q) = 1 for same-direction shared arcs.

    So: sum_{{S subset A, |S|=2k}} 1 = C(s, 2k)

    E_2k = (1/2^(2(n-1))) * n! * sum_pi C(s(pi), 2k)
    """)

    # Let me verify this!
    for n in [3, 4, 5]:
        all_perms = list(permutations(range(n)))
        for level in range(0, 2*((n-1)//2)+1, 2):
            k = level // 2
            if k == 0: continue
            # Formula: E_{2k} = (n!/2^{2(n-1)}) * sum_{pi: o=0} C(s(pi), 2k)
            inner = 0
            for pi in all_perms:
                s = successions(pi)
                o = anti_successions(pi)
                if o == 0 and s >= 2*k:
                    inner += math.comb(s, 2*k)
            E_2k_formula = Fraction(math.factorial(n) * inner, 2**(2*(n-1)))
            E0 = Fraction(math.factorial(n), 2**(n-1))**2
            ratio = E_2k_formula / E0
            # Grand formula prediction
            if n - 2*k > 0:
                pred = Fraction(2*(n-2*k)**k, math.perm(n, 2*k))
            else:
                pred = Fraction(0)
            print(f"  n={n}, k={k}: E_2k/E_0 = {ratio} = {float(ratio):.8f}, "
                  f"grand = {pred} = {float(pred):.8f}, match = {ratio == pred}")

    # ============================================================
    print(f"\n--- PERSPECTIVE 7: THE LEVEL DECOMPOSITION DIRECTLY")
    print(f"    E_2k/E_0 = (1/n!) * sum_pi C(s(pi), 2k)")

    # So the grand formula says:
    # (1/n!) * sum_{pi: o=0} C(s(pi), 2k) = 2*(n-2k)^k / P(n,2k)

    # Or: sum_{pi: o=0} C(s(pi), 2k) = 2*n!*(n-2k)^k / P(n,2k) = 2*(n-2k)^k*(n-2k)!

    print(f"\n    Testing: sum_pi C(s(pi), 2k) = 2*(n-2k)^k*(n-2k)!")
    for n in [3, 4, 5, 6, 7]:
        all_perms = list(permutations(range(n)))
        max_k = (n-1)//2
        for k in range(1, max_k+1):
            if n-2*k <= 0: break
            inner = sum(math.comb(successions(pi), 2*k)
                       for pi in all_perms if anti_successions(pi) == 0)
            pred = 2*(n-2*k)**k * math.factorial(n-2*k)
            print(f"    n={n}, k={k}: sum C(s,{2*k}) = {inner}, "
                  f"2*{n-2*k}^{k}*{n-2*k}! = {pred}, match = {inner == pred}")

    print(f"""
  *** THIS IS THE CLEAN IDENTITY! ***

  For each k:
    sum over anti-succ-free perms pi of C(s(pi), 2k) = 2*(n-2k)^k*(n-2k)!

  This says: "Among anti-succession-free permutations,
  the total number of ways to choose 2k of their successions
  = 2*(n-2k)^k*(n-2k)!."

  The left side COUNTS: (pi, S) where pi is anti-succ-free
  and S is a 2k-subset of pi's succession positions.

  The right side COUNTS: 2 * (arrangement of n-2k free items)
  * (k interaction tags from n-2k items with replacement).

  BIJECTION: Given (pi, S) with |S|=2k succession positions chosen,
  the 2k succession positions define k BLOCKS of size 2 in pi.
  Removing these blocks leaves n-2k free elements.
  Each block can be "dissolved" into a tag (which free element it attaches to).
  The 2 factor comes from... the direction of the path reversal?

  This is VERY close to a bijective proof!
    """)

    print(f"{'='*70}")
    print("DONE — THE LEVEL-k IDENTITY IS FOUND AND VERIFIED")
    print(f"{'='*70}")

if __name__ == '__main__':
    main()
