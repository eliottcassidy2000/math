#!/usr/bin/env python3
"""
CONNECTION: Key Identity (THM-030) and OCF / Independence Polynomial.

THM-030: B_b(W) + (-1)^m E_b(W) = 2r * cs_W(b)

OCF: H(T) = I(Omega(T), 2)

Can the Key Identity be used to prove or refine OCF?

APPROACH: The Key Identity gives a recurrence relating Hamiltonian path
counts to the transfer matrix column sums. If we can express the column
sums in terms of the independence polynomial, we get a new proof of OCF.

OBSERVATION: For VT tournaments at odd n:
  cs(b) = sum_{a != b} M[a,b] = 0 (since M = (H/n)*I, off-diag = 0)
  So B_b + (-1)^n E_b = 0, i.e., B_b = (-1)^{n+1} E_b

For odd n (VT): B_b = E_b (since (-1)^{n+1} = 1 when n is odd)
This says: every VT vertex has equal start-count and end-count!

For the KEY IDENTITY at the FULL vertex set W = V:
  B_b(V) + (-1)^n E_b(V) = 2r * cs_V(b) = cs_V(b) (at c=1)

Note: B_b(V) = sum_{a != b} H(a -> b) = total paths ENDING at b (not starting)
Wait, B_b(V) starts at b. So B_b(V) = sum over paths starting at b = E_b(V^op).
Actually no: B_b(V) is the total weight of Ham paths through V starting at b.
E_b(V) is the total weight of Ham paths through V ending at b.

So B_b(V) = #{Ham paths starting at b} and E_b(V) = #{Ham paths ending at b}.
(In the 0/1 tournament case.)

The identity B_b + (-1)^n E_b = cs(b) gives:
  n odd: B_b + E_b = cs(b)  [since (-1)^n = -1 for odd n... wait]

Actually (-1)^m where m = |W| = n. So:
  B_b + (-1)^n E_b = cs(b)

For n odd: (-1)^n = -1, so B_b - E_b = cs(b)
For n even: (-1)^n = 1, so B_b + E_b = cs(b)

And cs(b) = sum_{a != b} M[a,b].

For VT at odd n: cs(b) = 0, so B_b = E_b (each vertex starts and ends equally many paths).
For VT at even n: cs(b) = (n-1) * 2H/(n(n-1)) = 2H/n, so B_b + E_b = 2H/n.
Since VT implies B_b = E_b = H/n (by symmetry), this checks: 2H/n = 2H/n.

kind-pasteur-2026-03-06-S25 (continuation)
"""

from itertools import permutations
import numpy as np

def make_circulant(n, S):
    T = {}
    for i in range(n):
        for j in range(n):
            if i != j:
                T[(i,j)] = 1 if ((j - i) % n) in S else 0
    return T

def count_endpoint_stats(T, n):
    """Compute B_b and E_b for each vertex b."""
    B = [0] * n  # B[b] = paths starting at b
    E = [0] * n  # E[b] = paths ending at b
    for perm in permutations(range(n)):
        prod = 1
        for k in range(n-1):
            prod *= T.get((perm[k], perm[k+1]), 0)
        if prod > 0:
            B[perm[0]] += 1
            E[perm[-1]] += 1
    return B, E


# ============================================================
# Verify Key Identity at W = V for various tournaments
# ============================================================
print("=" * 70)
print("Key Identity at W = V (full vertex set)")
print("B_b(V) + (-1)^n E_b(V) = cs_V(b)")
print("=" * 70)

for n, name, S_set in [(3, "T_3", {1}), (5, "Paley T_5", {1,2}),
                         (5, "circ {1,3}", {1,3})]:
    T = make_circulant(n, S_set)
    B, E = count_endpoint_stats(T, n)
    H = sum(B)

    print(f"\n  {name} (n={n}): H = {H}")
    for b in range(min(n, 3)):
        lhs = B[b] + (-1)**n * E[b]
        print(f"    v={b}: B={B[b]}, E={E[b]}, B+(-1)^n*E = {lhs}")

    # For VT: B_b = E_b = H/n
    print(f"    H/n = {H/n}")
    print(f"    B = E = H/n: {all(B[i] == H//n and E[i] == H//n for i in range(n))}")


# ============================================================
# Non-VT tournaments: B_b != E_b
# ============================================================
print("\n" + "=" * 70)
print("Non-VT tournaments: B_b, E_b, and cs(b)")
print("=" * 70)

# Transitive tournament
for n in [3, 5, 7]:
    T = {}
    for i in range(n):
        for j in range(n):
            if i != j:
                T[(i,j)] = 1 if i < j else 0

    B, E = count_endpoint_stats(T, n)
    H = sum(B)
    print(f"\n  Transitive T_{n}: H = {H}")
    for b in range(n):
        lhs = B[b] + (-1)**n * E[b]
        print(f"    v={b}: B={B[b]}, E={E[b]}, B+(-1)^n*E = {lhs}")


# ============================================================
# CONNECTION TO OCF: can we relate cs(b) to Omega(T)?
# ============================================================
print("\n" + "=" * 70)
print("cs(b) and the independence polynomial")
print("=" * 70)

print("""
OCF: H(T) = I(Omega(T), 2) = 1 + 2*alpha_1 + 4*alpha_2 + ...

Key Identity: B_b + (-1)^n E_b = cs(b)

If we sum over b:
  sum_b B_b + (-1)^n sum_b E_b = sum_b cs(b)
  H + (-1)^n H = sum_b cs(b)

For odd n: 0 = sum_b cs(b) = Sigma (total off-diag sum) -- checks with THM-030.
For even n: 2H = sum_b cs(b) = Sigma -- also checks.

Now, the KEY QUESTION: can we express cs(b) in terms of Omega(T)?

Recall: cs(b) = sum_{a != b} M[a,b]
and M[a,b] = sum_S (-1)^|S| E_a(S+a) B_b(R+b)

The independence polynomial I(Omega(T), x) = sum_k alpha_k x^k
where alpha_k = #{independent sets of size k in Omega(T)}.

Can cs(b) be written as a "vertex-localized" independence polynomial?
For example, define I_b(x) = sum_k alpha_k(b) x^k
where alpha_k(b) counts independent sets of size k in Omega(T)
that "involve" vertex b (i.e., some cycle in the set contains b).

Then H(T) = sum_b delta_b where delta_b = (H(T) - H(T-b))/2 by OCF.
And delta_b = sum_{C containing b} mu(C) = I_b(2) - I_b(0) (approximately).

Is cs(b) related to delta_b?
""")

# Compute delta_b = (H(T) - H(T-b))/2 for various tournaments
from itertools import combinations

def count_odd_cycles(T, n):
    """Count directed odd cycles of all lengths."""
    cycles = []
    for length in range(3, n+1, 2):
        for vertices in combinations(range(n), length):
            # Check all cyclic orderings
            from itertools import permutations as perms
            for perm in perms(vertices):
                if perm[0] != min(vertices): continue  # canonical
                is_cycle = True
                for k in range(length):
                    if T.get((perm[k], perm[(k+1) % length]), 0) != 1:
                        is_cycle = False
                        break
                if is_cycle:
                    cycles.append(set(vertices))
                    break  # only one canonical rep per vertex set
    return cycles


for n, name, S_set in [(5, "Paley T_5", {1,2}), (5, "circ {1,3}", {1,3})]:
    T = make_circulant(n, S_set)
    B, E = count_endpoint_stats(T, n)
    H = sum(B)

    print(f"\n  {name} (n={n}): H = {H}")

    # Compute H(T-b) for each b
    for b in range(min(n, 3)):
        # Delete vertex b
        T_minus_b = {}
        remaining = [v for v in range(n) if v != b]
        for i in remaining:
            for j in remaining:
                if i != j:
                    T_minus_b[(i,j)] = T.get((i,j), 0)

        H_minus_b = 0
        for perm in permutations(remaining):
            prod = 1
            for k in range(len(perm)-1):
                prod *= T_minus_b.get((perm[k], perm[k+1]), 0)
            H_minus_b += prod

        delta_b = (H - H_minus_b)
        cs_b = B[b] + (-1)**n * E[b]
        print(f"    v={b}: H(T-b)={H_minus_b}, delta=H-H(T-b)={delta_b}, delta/2={delta_b/2}")
        print(f"           B={B[b]}, E={E[b]}, cs(b)=B+(-1)^n*E={cs_b}")


# ============================================================
# Vertex deletion and Key Identity
# ============================================================
print("\n" + "=" * 70)
print("Vertex deletion: H(T) - H(T-b) vs cs(b)")
print("=" * 70)

print("""
By OCF (Claim A): H(T) - H(T-b) = 2 * sum_{C containing b} mu(C)

And by Key Identity: cs(b) = B_b + (-1)^n E_b

For VT at odd n: cs(b) = 0 for all b.
So B_b = E_b (each vertex starts and ends equally many paths).

But H(T) - H(T-b) = 2*delta_b is generally NONZERO.
So cs(b) and delta_b are DIFFERENT quantities.

HOWEVER: sum_b cs(b) = 0 (odd n) and sum_b delta_b = H*(n-1) - sum_b H(T-b).
These are measuring different things.

cs(b) measures the NET endpoint imbalance of vertex b in Ham paths.
delta_b measures the contribution of cycles through b to H(T).
""")


print("\n" + "=" * 70)
print("DONE")
print("=" * 70)
