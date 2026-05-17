"""
Check if H values from all-0 staircase are Markov numbers.
Also search for the sequence in various databases.
"""

def enumerate_markov(limit=10**8):
    """BFS on the Markov tree to enumerate all Markov numbers up to limit."""
    markov_set = set()
    # BFS queue stores triples (a, b, c) with a <= b <= c
    from collections import deque
    queue = deque()
    seen_triples = set()

    # Seed: (1, 1, 1), (1, 1, 2), (1, 2, 5)
    initial = (1, 1, 2)
    queue.append(initial)
    seen_triples.add(initial)
    markov_set.update(initial)

    while queue:
        triple = queue.popleft()
        a, b, c = triple

        # Generate three children (one for each element replaced)
        for i in range(3):
            lst = list(triple)
            x, y, z = lst[(i+1)%3], lst[(i+2)%3], lst[i]
            # Replace z with z' = 3xy - z
            new_z = 3*x*y - z
            new_triple = tuple(sorted([x, y, new_z]))
            if new_z > 0 and new_z <= limit and new_triple not in seen_triples:
                seen_triples.add(new_triple)
                queue.append(new_triple)
                markov_set.add(new_z)

    return markov_set

print("Enumerating Markov numbers up to 10^8...")
markov_set = enumerate_markov(10**8)
print(f"Found {len(markov_set)} Markov numbers up to 10^8")
print(f"Sorted Markov numbers up to 100000:")
sorted_markov = sorted([m for m in markov_set if m <= 100000])
print(sorted_markov)

# H values from all-0 staircase
H_values = [5, 29, 233, 2489, 33773]
print(f"\nH values: {H_values}")
for H in H_values:
    in_markov = H in markov_set
    print(f"  H={H}: is Markov? {in_markov}")

# Look for Markov triples containing consecutive H values
print("\nSearching for Markov triples among H values:")
for i in range(len(H_values)):
    for j in range(i+1, len(H_values)):
        a, b = H_values[i], H_values[j]
        # Check if (a, b, c) could be a Markov triple for any c
        # a^2 + b^2 + c^2 = 3abc => c^2 - 3ab*c + (a^2+b^2) = 0
        discriminant = (3*a*b)**2 - 4*(a**2 + b**2)
        if discriminant >= 0:
            sqrt_d_approx = int(discriminant**0.5)
            for sd in range(max(0, sqrt_d_approx-2), sqrt_d_approx+3):
                if sd*sd == discriminant:
                    c1 = (3*a*b + sd) // 2
                    c2 = (3*a*b - sd) // 2
                    if (3*a*b + sd) % 2 == 0:
                        # Verify
                        if a**2 + b**2 + c1**2 == 3*a*b*c1:
                            print(f"  Markov triple ({a}, {b}, {c1}) ✓")
                        if c2 > 0 and a**2 + b**2 + c2**2 == 3*a*b*c2:
                            print(f"  Markov triple ({a}, {b}, {c2}) ✓")

# Sequence analysis
print("\n--- Sequence analysis for 5, 29, 233, 2489, 33773 ---")
H = H_values
print(f"Ratios H[k+1]/H[k]: {[H[i+1]/H[i] for i in range(len(H)-1)]}")
print(f"Differences: {[H[i+1]-H[i] for i in range(len(H)-1)]}")
print(f"Second differences: {[(H[i+2]-H[i+1])-(H[i+1]-H[i]) for i in range(len(H)-2)]}")

# Check for linear recurrence a(k) = p*a(k-1) + q*a(k-2) + ...
print("\nTesting linear recurrences:")
# Order 2: a[k] = p*a[k-1] + q*a[k-2]
# From H[2]=233, H[3]=2489, H[4]=33773: need H[1]=29, H[2]=233, H[3]=2489, H[4]=33773
# 2489 = p*233 + q*29
# 33773 = p*2489 + q*233
import numpy as np
A = np.array([[233, 29], [2489, 233]])
b = np.array([2489, 33773])
try:
    sol = np.linalg.solve(A, b)
    p, q = sol
    print(f"Order 2 recurrence: a(k) = {p:.4f}*a(k-1) + {q:.4f}*a(k-2)")
    # Verify
    check = p*H[1] + q*H[0]
    print(f"  Predicts H[2] = {check:.2f}, actual = {H[2]}")
    check2 = p*H[3] + q*H[2]
    print(f"  Predicts H[4] = {check2:.2f}, actual = {H[4]}")
except:
    print("Order 2 recurrence: no solution")

# Order 2 with integers
# Solve p*233 + q*29 = 2489 and p*2489 + q*233 = 33773
from math import gcd
# 2489 = 233p + 29q
# 33773 = 2489p + 233q
# Eliminate: 33773*29 - 2489*233 = 233*29*(p - p_prev) + ... actually:
# From row reduction: 233 * row2 - 2489 * row1:
# 233*33773 - 2489*2489 = 233*233*q - 2489*29*q
# (233^2 - 2489*29)*q = 233*33773 - 2489^2
lhs_q_coeff = 233**2 - 2489*29
rhs_q = 233*33773 - 2489**2
print(f"\nFor order-2 integer recurrence:")
print(f"Coefficient of q: {lhs_q_coeff}")
print(f"RHS for q: {rhs_q}")
if rhs_q % lhs_q_coeff == 0:
    q_int = rhs_q // lhs_q_coeff
    print(f"q = {q_int}")
    # Find p from 233p = 2489 - 29*q
    if (2489 - 29*q_int) % 233 == 0:
        p_int = (2489 - 29*q_int) // 233
        print(f"p = {p_int}")
        print(f"Recurrence: a(k) = {p_int}*a(k-1) + {q_int}*a(k-2)")
        # Verify all terms
        for i in range(2, len(H)):
            pred = p_int*H[i-1] + q_int*H[i-2]
            print(f"  H[{i}] predicted: {pred}, actual: {H[i]}, match: {pred==H[i]}")
    else:
        print(f"p not integer: {(2489 - 29*q_int)/233}")
else:
    print(f"q not integer: {rhs_q}/{lhs_q_coeff} = {rhs_q/lhs_q_coeff:.6f}")

# Check if this sequence matches any well-known sequences
print("\n--- Known sequence checks ---")
# Fibonacci-related
fib = [1, 1]
while fib[-1] < 10**6:
    fib.append(fib[-1] + fib[-2])
print(f"Fibonacci: {fib[:20]}")
for H_val in H_values:
    print(f"  H={H_val} is Fibonacci? {H_val in fib}")

# Pell numbers
pell = [0, 1]
while pell[-1] < 10**6:
    pell.append(2*pell[-1] + pell[-2])
print(f"Pell: {pell[:20]}")

# Catalan numbers
catalan = [1]
for n in range(1, 15):
    catalan.append(catalan[-1] * 2*(2*n-1) // (n+1))
print(f"Catalan: {catalan}")

# Alpha coefficients structure
print("\n--- Alpha coefficient structure ---")
alpha_by_k = {
    2: [1, 2],
    3: [1, 12, 1],
    4: [1, 68, 24],
    5: [1, 530, 317, 20],
    6: [1, 5750, 4244, 642, 10],
}

print("Cycle counts by length:")
cycle_counts = {
    2: {3: 2},
    3: {3: 6, 5: 6},
    4: {3: 12, 5: 28, 7: 28},
    5: {3: 20, 5: 80, 7: 220, 9: 210},
    6: {3: 30, 5: 180, 7: 906, 9: 2480, 11: 2154},
}

print("  3-cycles: 2, 6, 12, 20, 30 = k(k-1) for k=2..6 ✓")
print(f"  Check: {[k*(k-1) for k in range(2,7)]}")

print("\n  5-cycle counts by k:")
for k in range(3, 7):
    cnt = cycle_counts[k].get(5, 0)
    print(f"    k={k}: {cnt}")
# 6, 28, 80, 180
# Look for pattern: 6=6, 28=?, 80=?, 180=?
five_cyc = [6, 28, 80, 180]
print(f"  Ratios: {[five_cyc[i+1]/five_cyc[i] for i in range(len(five_cyc)-1)]}")
print(f"  Differences: {[five_cyc[i+1]-five_cyc[i] for i in range(len(five_cyc)-1)]}")
# 22, 52, 100 - differences of differences: 30, 48 - ratios: 100/52, 52/22...
# Actually: 6 = C(4,2), 28 = C(8,2), 80 = ?, 180 = ?
# Or: 6 = 2*3, 28 = 4*7, 80 = 5*16, 180 = 6*30?
# 3, 7, 16, 30 - differences: 4, 9, 14 - going by 5?
# Or: 6, 28, 80, 180 - look at ratios with k: 6/(3*2)=1, 28/(4*7)=1, 80/(5*16)=1, 180/(6*30)=1?
# 3*2, 7*4, 16*5, 30*6
print(f"  5-cycles/3-cycles: {[five_cyc[i]/cycle_counts[i+3][3] for i in range(4)]}")

print("\n  Independence number (max IS) per k:")
for k, alpha in alpha_by_k.items():
    max_is = len(alpha) - 1
    max_count = alpha[-1]
    print(f"    k={k}: max_IS={max_is}, count={max_count}")
