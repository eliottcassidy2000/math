#!/usr/bin/env python3
"""
Algebraic proof of sigma_adj = (n-3)! * [2*t_3 - C(n,3)/2]

sigma_adj = sigma({0,1}) = sum_P (A[p_0,p_1] - 1/2)(A[p_1,p_2] - 1/2)

Expanding:
sigma_adj = sum_P [A[p_0,p_1]*A[p_1,p_2] - (1/2)*A[p_0,p_1] - (1/2)*A[p_1,p_2] + 1/4]

= sum_P A[p_0,p_1]*A[p_1,p_2] - (1/2)*sum_P A[p_0,p_1] - (1/2)*sum_P A[p_1,p_2] + n!/4

For any position i: sum_P A[p_i,p_{i+1}] = n!/2
  (because summing over all permutations, each ordered pair (a,b) appears (n-2)! times
   at positions (i,i+1), and sum_{a!=b} A[a,b] = C(n,2), so total = (n-2)!*C(n,2) = n!/2)

So: sigma_adj = sum_P A[p_0,p_1]*A[p_1,p_2] - n!/4 - n!/4 + n!/4
              = sum_P A[p_0,p_1]*A[p_1,p_2] - n!/4

Now: sum_P A[p_0,p_1]*A[p_1,p_2]
For each ordered triple (a,b,c) of distinct vertices with A[a,b]=A[b,c]=1:
  how many permutations P have (p_0,p_1,p_2) = (a,b,c)? Answer: (n-3)!

So: sum_P A[p_0,p_1]*A[p_1,p_2] = (n-3)! * D
where D = #{ordered triples (a,b,c): A[a,b]=1 and A[b,c]=1, all distinct}
       = sum_b in_deg(b) * out_deg(b)   [since in-neighbors of b choose a, out-neighbors choose c]
       Wait: A[a,b]=1 means a->b, so a is in in-neighbors. But we need a!=c.
       Since a->b and b->c, and a,b,c distinct (but a could equal c? No, they're distinct in the permutation).

D = sum_b d_b^in * d_b^out (since a ranges over in-neighbors of b, c over out-neighbors, and a!=c automatically because in the tournament a->b means b doesn't point to a, while b->c, so c!=a unless there's a loop, impossible).

Wait, a and c can be equal... no, in a permutation positions 0,1,2 must be distinct.
And d_b^in = n-1-s_b, d_b^out = s_b where s_b is the score (out-degree) of b.

D = sum_b (n-1-s_b) * s_b = (n-1)*sum s_b - sum s_b^2

Now: sum s_b = C(n,2) (total arcs)
And: sum s_b^2 = ? (depends on tournament)

There's a known identity: t_3 = C(n,3) - (1/2)*sum s_b*(n-1-s_b)/1
Hmm, let me use the exact relation.

For a tournament: # 3-cycles = C(n,3) - (1/2)*sum_b C(s_b, 2)... no.
Actually: t_3 = C(n,3) - sum_b C(s_b, 2)? Let me check.

# transitive triples = sum_b C(s_b, 2) (choose 2 out-neighbors of b, they form a transitive triple with b as source)
Wait, that overcounts. Let's be careful.

For each triple {a,b,c}, it's either a 3-cycle or a transitive triple.
# transitive triples = C(n,3) - t_3
# transitive triples = sum_b C(s_b, 2) = (1/2)*sum s_b*(s_b-1)

So: t_3 = C(n,3) - (1/2)*sum s_b*(s_b-1) = C(n,3) - (1/2)*(sum s_b^2 - sum s_b)
       = C(n,3) - (1/2)*sum s_b^2 + (1/2)*C(n,2)

=> sum s_b^2 = 2*C(n,3) - 2*t_3 + C(n,2)

Let me verify this at n=5.
kind-pasteur-2026-03-06-S25g
"""
from itertools import permutations, combinations
from math import factorial, comb

print("ALGEBRAIC PROOF OF SIGMA_ADJ")
print("=" * 50)

for n in [5, 7]:
    print(f"\n--- n={n} ---")

    # Verify sum s_b^2 = 2*C(n,3) - 2*t3 + C(n,2) on a few tournaments
    for trial in range(5):
        import random
        random.seed(42 + trial)
        A = [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(i+1, n):
                if random.random() < 0.5:
                    A[i][j] = 1
                else:
                    A[j][i] = 1

        scores = [sum(A[v]) for v in range(n)]
        sum_s2 = sum(s*s for s in scores)

        t3 = sum(1 for v in combinations(range(n), 3)
                 for p in permutations(v)
                 if all(A[p[i]][p[(i+1)%3]] == 1 for i in range(3))) // 3

        pred_sum_s2 = 2*comb(n,3) - 2*t3 + comb(n,2)

        if trial < 3:
            print(f"  trial={trial}: t3={t3}, sum_s2={sum_s2}, pred={pred_sum_s2}, match={'YES' if sum_s2==pred_sum_s2 else 'NO'}")

    # Now compute sigma_adj algebraically
    # D = (n-1)*C(n,2) - sum_s2 = (n-1)*C(n,2) - (2*C(n,3) - 2*t3 + C(n,2))
    # sigma_adj = (n-3)!*D - n!/4

    print(f"\n  Algebraic formula:")
    print(f"  D = (n-1)*C(n,2) - sum_s2")
    print(f"    = (n-1)*C(n,2) - 2*C(n,3) + 2*t3 - C(n,2)")
    print(f"    = ((n-1)-1)*C(n,2) - 2*C(n,3) + 2*t3")
    print(f"    = (n-2)*C(n,2) - 2*C(n,3) + 2*t3")

    D_const = (n-2)*comb(n,2) - 2*comb(n,3)
    print(f"    = {D_const} + 2*t3")

    print(f"  sigma_adj = (n-3)!*D - n!/4")
    print(f"            = (n-3)!*({D_const} + 2*t3) - {factorial(n)}/4")
    print(f"            = {factorial(n-3)}*({D_const} + 2*t3) - {factorial(n)//4}")
    print(f"            = {factorial(n-3)*D_const - factorial(n)//4} + {2*factorial(n-3)}*t3")

    # Expected: sigma_adj = (n-3)! * [2*t3 - C(n,3)/2]
    #         = 2*(n-3)!*t3 - (n-3)!*C(n,3)/2
    pred_coeff = 2*factorial(n-3)
    pred_const = -factorial(n-3)*comb(n,3)//2
    print(f"  Expected: {pred_const} + {pred_coeff}*t3")

    actual_const = factorial(n-3)*D_const - factorial(n)//4
    actual_coeff = 2*factorial(n-3)
    print(f"  Computed:  {actual_const} + {actual_coeff}*t3")
    print(f"  Match: {'YES' if actual_const == pred_const and actual_coeff == pred_coeff else 'NO'}")

    if actual_const != pred_const:
        print(f"  MISMATCH in constant: {actual_const} vs {pred_const}")
        print(f"  D_const = (n-2)*C(n,2) - 2*C(n,3) = {(n-2)*comb(n,2)} - {2*comb(n,3)} = {D_const}")
        print(f"  (n-3)!*D_const = {factorial(n-3)*D_const}")
        print(f"  n!/4 = {factorial(n)//4}")
        print(f"  Difference = {factorial(n-3)*D_const - factorial(n)//4}")
        print(f"  Expected const = -(n-3)!*C(n,3)/2 = {-factorial(n-3)*comb(n,3)//2}")

        # Let me verify: is (n-2)*C(n,2) - 2*C(n,3) = C(n,3)?
        lhs = (n-2)*comb(n,2) - 2*comb(n,3)
        rhs = comb(n,3)
        print(f"  (n-2)*C(n,2) - 2*C(n,3) = {lhs}, C(n,3) = {rhs}, equal: {lhs == rhs}")

        # So D_const = C(n,3). Then:
        # sigma_adj = (n-3)!*C(n,3) + 2*(n-3)!*t3 - n!/4
        # Expected: -(n-3)!*C(n,3)/2 + 2*(n-3)!*t3
        # So: (n-3)!*C(n,3) - n!/4 should equal -(n-3)!*C(n,3)/2
        # (n-3)!*C(n,3) - n!/4 = (n-3)!*C(n,3) - (n*(n-1)*(n-2)*(n-3)!)/4
        # = (n-3)! * [C(n,3) - n(n-1)(n-2)/4]
        # = (n-3)! * [n(n-1)(n-2)/6 - n(n-1)(n-2)/4]
        # = (n-3)! * n(n-1)(n-2) * [1/6 - 1/4]
        # = (n-3)! * n(n-1)(n-2) * (-1/12)
        # = -(n-3)! * C(n,3) * 2/12 * 3
        # Hmm, let me just compute:
        val = factorial(n-3) * (comb(n,3) - n*(n-1)*(n-2)//4)
        print(f"  (n-3)! * [C(n,3) - n(n-1)(n-2)/4] = {val}")
        # C(n,3) = n(n-1)(n-2)/6
        # C(n,3) - n(n-1)(n-2)/4 = n(n-1)(n-2)*(1/6 - 1/4) = n(n-1)(n-2)*(-1/12)
        # = -C(n,3)/2
        print(f"  -C(n,3)/2 = {-comb(n,3)/2}")
        # So (n-3)! * [C(n,3) - n(n-1)(n-2)/4] = -(n-3)!*C(n,3)/2
        # This should match. But n(n-1)(n-2)/4 might not be integer...
        print(f"  n(n-1)(n-2)/4 = {n*(n-1)*(n-2)/4}")

print("\n" + "="*50)
print("COMPLETE ALGEBRAIC PROOF")
print("="*50)
print("""
sigma_adj = (n-3)! * D - n!/4
  where D = (n-2)*C(n,2) - 2*C(n,3) + 2*t_3 = C(n,3) + 2*t_3

So: sigma_adj = (n-3)! * (C(n,3) + 2*t_3) - n!/4
             = (n-3)!*C(n,3) + 2*(n-3)!*t_3 - n!/4

Now: n!/4 = n*(n-1)*(n-2)*(n-3)!/4 = C(n,3)*6*(n-3)!/4 = C(n,3)*(n-3)!*3/2

So: sigma_adj = (n-3)!*C(n,3) + 2*(n-3)!*t_3 - (n-3)!*C(n,3)*3/2
             = (n-3)! * [C(n,3) - 3*C(n,3)/2 + 2*t_3]
             = (n-3)! * [-C(n,3)/2 + 2*t_3]
             = (n-3)! * [2*t_3 - C(n,3)/2]

Therefore: w_{n-3} = (n-2) * sigma_adj = (n-2)*(n-3)! * [2*t_3 - C(n,3)/2]
                   = (n-2)! * [2*t_3 - C(n,3)/2]

QED!

The proof uses:
1. sigma_adj = sum_P s_{p_0,p_1} * s_{p_1,p_2} where s = A - 1/2
2. Expanding: uses sum_P A[p_i,p_{i+1}] = n!/2 and the 2-path count D
3. D = sum_b d_b^in * d_b^out = (n-1)*C(n,2) - sum s_b^2
4. sum s_b^2 = 2*C(n,3) - 2*t_3 + C(n,2) (from t_3 = C(n,3) - sum C(s_b,2))
5. Non-adjacent sigma = 0 (independent zero-mean factors)
""")
