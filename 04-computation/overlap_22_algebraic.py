#!/usr/bin/env python3
"""
ALGEBRAIC DERIVATION OF THE (2,2) CONTRIBUTION TO tr(c_{n-5})

THEOREM: At n=7, the (2,2) contribution to tr(c_2) = sum_P e_4(s_P) is:
  24*bc - 24*t_3 + 105

where bc = total complementary-cyclic-triple pairs across all 6-vertex subsets.

PROOF:
=====

1. A (2,2) position pattern is {i, i+1, j, j+1} with j >= i+3.
   At n=7: 3 such patterns ({0,1,3,4}, {0,1,4,5}, {1,2,4,5}).

2. For each pattern, sum_P s_i*s_{i+1}*s_j*s_{j+1} factorizes.
   The 2-block {i,i+1} involves 3 vertices: p_i, p_{i+1}, p_{i+2}.
   The 2-block {j,j+1} involves 3 vertices: p_j, p_{j+1}, p_{j+2}.
   Since j >= i+3, these 6 vertices are all distinct.

3. Define F(S) := sum_{orderings of 3-set S} (A[a,b]-1/2)*(A[b,c]-1/2)
   for S = {a,b,c}, summing over all 3! orderings.

   F(S) = H(S) - 3/2, where H(S) = #{directed Ham paths on S}.
   - Transitive triple: H=1, F = -1/2.
   - Cyclic triple: H=3, F = 3/2.

   Proof: sum s_0*s_1 = sum T_0*T_1 - sum T_0/2 - sum T_1/2 + 6/4.
   sum T_0*T_1 = H(S). sum T_0 = sum T_1 = 3 (each directed arc appears once).
   So F = H - 3/2 - 3/2 + 3/2 = H - 3/2.

4. The (2,2) sum for one pattern:
   sum_P s_i*s_{i+1}*s_j*s_{j+1}
   = (n-6)! * sum_{6-vertex subsets X} sum_{ordered (S1,S2) partition of X}
     F(S1)*F(S2)

   where S1 at positions {i,i+1,i+2} and S2 at {j,j+1,j+2}.

5. For each 6-vertex subset X, there are C(6,3)=20 ordered partitions
   into (S1,S2). Since F(S1)*F(S2) = F(S2)*F(S1), this equals
   2 * sum_{10 unordered partitions} F(T)*F(T').

6. Write F(T) = 2*c(T) - 1/2 where c(T)=1 if T is cyclic, 0 if transitive.
   Then F(T)*F(T') = 4*c(T)*c(T') - c(T) - c(T') + 1/4.

   Sum over 10 unordered partitions {T,T'} of X:
   - sum c(T)*c(T') = bc(X) [#partitions where both are cyclic]
   - sum (c(T)+c(T')) = t_3(X) [#cyclic triples in X, each appears in 1 partition]
   - constant = 10 * 1/4 = 5/2

   So: sum_{10} F*F' = 4*bc(X) - t_3(X) + 5/2.

7. One pattern contribution:
   = (n-6)! * sum_{X: |X|=6} 2*(4*bc(X) - t_3(X) + 5/2)
   = 2*(n-6)! * [4*sum_X bc(X) - sum_X t_3(X) + C(n,6)*5/2]

   At n=7: sum_X bc(X) = BC (= bc at n=7).
   sum_X t_3(X) = C(n-3, 6-3)*t_3 = C(4,3)*t_3 = 4*t_3.
   C(7,6) = 7.

   One pattern = 2*1*[4*BC - 4*t_3 + 7*5/2] = 8*BC - 8*t_3 + 35.

8. Total (2,2) at n=7: 3 patterns * (8*BC - 8*t_3 + 35) = 24*BC - 24*t_3 + 105.

QED.

Combined with the (4,) contribution (126 - 36*t_3 + 12*t_5) from OCF@n=5:
tr(c_2) = (126 - 36*t_3 + 12*t_5) + (24*BC - 24*t_3 + 105) = 24*BC - 60*t_3 + 12*t_5 + 231.

KEY OBSERVATION: The (2,2) derivation is PURELY ALGEBRAIC.
It does NOT use OCF at n=7 or any unproved result. It only uses:
- The definition of s_i
- The factorization of disjoint products
- The tournament identity A[a,b] + A[b,a] = 1
- Basic combinatorics

GENERAL FORMULA at any odd n >= 7:

(2,2) count = C(n-4, 2) position patterns.
Each contributes: (n-6)! * 2 * [4*sum_X bc(X) - sum_X t_3(X) + C(n,6)*5/2].

At general n:
sum_X t_3(X) = C(n-3, 3) * t_3  [each 3-cycle in C(n-3,3) six-element subsets]

So one pattern = (n-6)! * 2 * [4*BC - C(n-3,3)*t_3 + 5*C(n,6)/2].

Total (2,2) = C(n-4,2) * (n-6)! * 2 * [4*BC - C(n-3,3)*t_3 + 5*C(n,6)/2]
            = C(n-4,2) * (n-6)! * [8*BC - 2*C(n-3,3)*t_3 + 5*C(n,6)]

At n=7: C(3,2)*1*[8*BC - 2*4*t_3 + 5*7] = 3*[8BC-8t_3+35] = 24BC-24t_3+105 ✓

opus-2026-03-06-S11b (continued^5)
"""
print("Algebraic proof recorded. See derivation in docstring.")
print("DONE")
