#!/usr/bin/env python3
"""
Algebraic proof attempt for w_{n-3} = (n-2)! * [2*t_3 - C(n,3)/2]

W(r) = sum_P prod_{i=0}^{n-2} (r + s_{p_i, p_{i+1}})
     = sum_P prod (r + A[p_i,p_{i+1}] - 1/2)

Expanding the product:
  prod (r + s_e) = sum_{S subset [n-2]} r^{n-1-|S|} prod_{i in S} s_{p_i,p_{i+1}}

So W(r) = sum_k r^{n-1-k} * sigma_k
where sigma_k = sum_P sum_{|S|=k} prod_{i in S} s_{p_i,p_{i+1}}
             = sum_{|S|=k} sum_P prod_{i in S} s_{p_i,p_{i+1}}
             = sum_{|S|=k} sigma(S)

But W has only even-powered terms (odd n), so sigma_k = 0 for k even.

w_{n-3} is the coefficient of r^{n-1-2} = r^{n-3}, which corresponds to k=2:
  w_{n-3} = sum_{|S|=2} sigma(S)

S has 2 elements from {0, 1, ..., n-2}. These elements can be:
(a) Adjacent: S = {i, i+1} for some i. There are n-2 such sets.
(b) Non-adjacent: S = {i, j} with |i-j| >= 2.

For adjacent S: sigma({i,i+1}) = sum_P s_{p_i,p_{i+1}} * s_{p_{i+1},p_{i+2}}
  = sum_P (A[p_i,p_{i+1}]-1/2)(A[p_{i+1},p_{i+2}]-1/2)

  By position-translation invariance, this doesn't depend on i.
  = sigma_adj (one value for all adjacent pairs)

For non-adjacent S: sigma({i,j}) = sum_P s_{p_i,p_{i+1}} * s_{p_j,p_{j+1}}

  By translation invariance, this depends only on gap = j - i.
  For gap >= 2, the two factors involve DISJOINT position pairs.
  = sigma_{gap}

Now sigma_1 (singleton set, gap 1) = sum_P (A[p_i,p_{i+1}] - 1/2) = 0
  because sum_P A[p_i,p_{i+1}] = n!/2 (by symmetry, half the time the arc goes forward).

So the sigma(S) for |S|=1 is 0. Good, that's w_{n-2} = 0 (odd index).

For |S| = 2:
  Adjacent (gap 1): sigma({i,i+1}) for i=0,...,n-3. Count: n-2.
  Non-adjacent (gap g, g>=2): sigma({i,i+g}) for various i. Count: depends on g.

Total: w_{n-3} = (n-2) * sigma_adj + sum_{g>=2} (n-1-g) * sigma_g

The key is: what is sigma_adj and what is sigma_g?

sigma_adj = sum_P (A[p_0,p_1]-1/2)(A[p_1,p_2]-1/2)
  = sum_P A[p_0,p_1]*A[p_1,p_2] - (1/2)*sum_P A[p_0,p_1] - (1/2)*sum_P A[p_1,p_2] + (1/4)*n!
  = sigma({0,1}, A*A) - (1/2)*n!/2 - (1/2)*n!/2 + n!/4
  = sigma({0,1}, A*A) - n!/4

where sigma({0,1}, A*A) = sum_P A[p_0,p_1]*A[p_1,p_2] counts length-2 directed paths
starting at position (p_0, p_1, p_2) summed over all n! permutations.

This sum counts: for each ordered triple (a,b,c) of distinct vertices with A[a,b]=1 and A[b,c]=1,
how many permutations have (p_0,p_1,p_2) = (a,b,c)?
Answer: (n-3)! (the remaining n-3 positions can be any permutation).

So sigma({0,1}, A*A) = (n-3)! * sum_{a,b,c distinct} A[a,b]*A[b,c]
                     = (n-3)! * sum_b sum_{a!=b} A[a,b] * sum_{c!=b,a} A[b,c]

Hmm, this gets complicated. Let me use a different approach.

Actually: sum_{a,b,c distinct} A[a,b]*A[b,c] counts ordered directed 2-paths a->b->c.
For a tournament: each vertex b has in-degree d_in(b) and out-degree d_out(b) = n-1-d_in(b).
The number of 2-paths through b is d_in(b) * (d_out(b) - [if a=c, but a!=c guaranteed]).

Hmm, actually for ORDERED (a,b,c) with all distinct:
sum = sum_b (number of a with A[a,b]=1) * (number of c!=a,b with A[b,c]=1)
But c can equal a in theory... no, a,b,c are all distinct.

Let d_b^- = in-degree of b, d_b^+ = out-degree of b.
sum = sum_b d_b^- * d_b^+  ... wait, this overcounts because a!=c.

sum = sum_b [d_b^- * d_b^+ - (# of a with A[a,b]=1 and A[b,a]=1)]
But A[a,b]=1 and A[b,a]=1 can't both hold in a tournament! So no correction needed.

sum_{(a,b,c) ordered distinct} A[a,b]*A[b,c] = sum_b d_b^- * d_b^+

And d_b^- + d_b^+ = n-1, so d_b^- * d_b^+ = d_b^-(n-1-d_b^-).
Let S = score sequence. Then sum_b d_b^-(n-1-d_b^-) = (n-1)*sum d_b^- - sum d_b^-^2
= (n-1)*C(n,2) - sum s_b^2 where s_b = d_b^+ (out-degree).

Actually d_b^- = n-1-s_b, so d_b^-*d_b^+ = (n-1-s_b)*s_b.
sum = sum_b (n-1-s_b)*s_b = (n-1)*sum s_b - sum s_b^2
    = (n-1)*C(n,2) - sum s_b^2

And sum s_b = C(n,2) (total arcs), sum s_b^2 = C(n,2) + 4*C(n,3) - 4*t_3
(known identity: the sum of squared scores in a tournament).

Wait, let me recall the identity. For a tournament:
sum s_b^2 = sum s_b(s_b-1) + sum s_b = (number of 2-paths) + C(n,2)

Hmm, but # of 2-paths a->b->c is sum s_b choose 2 * 2! = sum s_b(s_b-1)... no.
# of ordered 2-paths b->c1, b->c2 (distinct) = sum s_b(s_b-1).

Actually this is getting circular. Let me just compute numerically.
kind-pasteur-2026-03-06-S25g
"""
from itertools import permutations
from math import factorial, comb

# Verify the algebraic decomposition at n=5
n = 5

# Compute sigma_adj and sigma_g for various gaps
for bits in [0, 7, 15, 31]:
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1

    t3 = sum(1 for i in range(n) for j in range(i+1,n) for k in range(j+1,n)
             if A[i][j]*A[j][k]*A[k][i] + A[i][k]*A[k][j]*A[j][i])

    # Compute sigma(S) for each S of size 2
    for gap in range(1, n-1):
        S_list = [(i, i+gap) for i in range(n-1-gap)]
        total = 0
        for S in S_list:
            for p in permutations(range(n)):
                prod = 1.0
                for pos in S:
                    prod *= (A[p[pos]][p[pos+1]] - 0.5)
                total += prod
        count = len(S_list)
        per_s = total / count if count > 0 else 0
        print(f"  bits={bits:2d} t3={t3}: gap={gap}, count={count}, sigma/count={per_s:.4f}")

    # Total w2 = w_{n-3}
    w2 = 0.0
    for S in [(i,j) for i in range(n-1) for j in range(i+1, n-1)]:
        for p in permutations(range(n)):
            prod = 1.0
            for pos in S:
                prod *= (A[p[pos]][p[pos+1]] - 0.5)
            w2 += prod

    pred = 2*factorial(n-2) * t3 - comb(n,3)*factorial(n-2)//2
    print(f"  w2={w2:.4f}, predicted={pred}, err={abs(w2-pred):.4f}")
    print()

# Now compute sigma_adj in terms of score sequence
print("\n" + "="*50)
print("SIGMA_ADJ DECOMPOSITION")
print("="*50)

for bits in range(1 << 10):
    A = [[0]*5 for _ in range(5)]
    idx = 0
    for i in range(5):
        for j in range(i+1, 5):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1

    t3 = sum(1 for i in range(5) for j in range(i+1,5) for k in range(j+1,5)
             if A[i][j]*A[j][k]*A[k][i] + A[i][k]*A[k][j]*A[j][i])

    scores = [sum(A[i]) for i in range(5)]
    sum_s2 = sum(s*s for s in scores)

    # sigma_adj: pick gap=1, one position set {0,1}
    sigma_adj = 0
    for p in permutations(range(5)):
        sigma_adj += (A[p[0]][p[1]] - 0.5) * (A[p[1]][p[2]] - 0.5)

    # sigma_gap2: gap=2, position set {0,2}
    sigma_gap2 = 0
    for p in permutations(range(5)):
        sigma_gap2 += (A[p[0]][p[1]] - 0.5) * (A[p[2]][p[3]] - 0.5)

    # sigma_gap3: gap=3, position set {0,3}
    sigma_gap3 = 0
    for p in permutations(range(5)):
        sigma_gap3 += (A[p[0]][p[1]] - 0.5) * (A[p[3]][p[4]] - 0.5)

    # w2 = 3*sigma_adj + 2*sigma_gap2 + 1*sigma_gap3
    w2 = 3*sigma_adj + 2*sigma_gap2 + 1*sigma_gap3
    pred = 12*t3 - 30

    if bits < 5:
        print(f"  bits={bits}: t3={t3}, sum_s2={sum_s2}")
        print(f"    sigma_adj={sigma_adj:.2f}, sigma_gap2={sigma_gap2:.2f}, sigma_gap3={sigma_gap3:.2f}")
        print(f"    w2={w2:.2f}, pred={pred}")

# Check if sigma_gap >= 2 always = 0
all_gap2_zero = all(True for _ in range(0))  # placeholder
gap2_vals = set()
gap3_vals = set()

for bits in range(1 << 10):
    A = [[0]*5 for _ in range(5)]
    idx = 0
    for i in range(5):
        for j in range(i+1, 5):
            if (bits >> idx) & 1:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1

    sigma_gap2 = sum(
        (A[p[0]][p[1]] - 0.5) * (A[p[2]][p[3]] - 0.5)
        for p in permutations(range(5))
    )
    sigma_gap3 = sum(
        (A[p[0]][p[1]] - 0.5) * (A[p[3]][p[4]] - 0.5)
        for p in permutations(range(5))
    )
    gap2_vals.add(round(sigma_gap2, 4))
    gap3_vals.add(round(sigma_gap3, 4))

print(f"\n  sigma_gap2 values: {sorted(gap2_vals)}")
print(f"  sigma_gap3 values: {sorted(gap3_vals)}")

if gap2_vals == {0.0} and gap3_vals == {0.0}:
    print("\n  BOTH non-adjacent sigmas are ZERO!")
    print("  => w_{n-3} = (n-2) * sigma_adj")
    print("  => sigma_adj = w_{n-3} / (n-2)")
    print("  => sigma_adj = (n-3)! * [2*t3 - C(n,3)/2]")
elif gap3_vals == {0.0}:
    print("\n  sigma_gap3 = 0 always. sigma_gap2 varies.")
else:
    print("\n  Both non-adjacent sigmas are non-zero.")

print("\nDONE")
