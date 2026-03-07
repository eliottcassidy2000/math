#!/usr/bin/env python3
"""
Can we prove OCF (H(T) = I(Omega(T), 2)) by induction on n?

Induction step: Assume OCF holds for n-1, prove for n.

H(T) = H(T-v) + 2*sum_{C through v} mu(C)  [Claim A]

If Claim A holds AND H(T-v) = I(Omega(T-v), 2) (induction hypothesis), then:

H(T) = I(Omega(T-v), 2) + 2*sum_{C through v} mu(C)

We need: I(Omega(T), 2) = I(Omega(T-v), 2) + 2*sum_{C through v} mu(C)

This is EXACTLY Claim B! Which is already proved (THM-003).

But wait — this circular: we need Claim A to do the induction, and
Claim A is what we're trying to prove.

ALTERNATIVE: Can we prove Claim A directly by induction?

Let delta_H(T,v) = H(T) - H(T-v).
Let delta_I(T,v) = I(Omega(T), 2) - I(Omega(T-v), 2).

By Claim B: delta_I(T,v) = 2*sum_{C through v} mu(C).

So OCF (H = I) is equivalent to: delta_H(T,v) = delta_I(T,v) for all T, v.

If we can show delta_H = delta_I for all T,v, then by induction
(base case n=1: H=1=I=1), OCF follows.

But delta_H = delta_I IS Claim A (since delta_I = 2*sum mu by Claim B).

So proving Claim A DIRECTLY would give an independent proof of OCF!

The question is: can we prove Claim A without using OCF?

From THM-070, we showed: IF OCF holds, THEN Claim A follows (via through-v clique).
But can we go the other direction?

Actually, delta_H has a DIRECT combinatorial interpretation:
  delta_H(T,v) = sum_{P' Ham path of T-v} insact(v, P')
where insact counts valid insertion positions.

And delta_I has:
  delta_I(T,v) = 2*sum_{C through v} mu(C)   [Claim B]

So Claim A = "sum of insertion counts = 2 * sum of mu values"

We already know inshat(v,P') is always odd [Redei].
And (inshat-1)/2 = #{TypeII positions} = #{3-cycle embeddings at P'} [THM-004/005].

But insact != inshat in general...

Let me compute the relationship between insact and the OCF deletion
formula more carefully.

opus-2026-03-07-S36
"""
from itertools import permutations, combinations

def tournament_from_bits(n, bits_int):
    m = n*(n-1)//2
    b = [(bits_int >> k) & 1 for k in range(m)]
    A = [[0]*n for _ in range(n)]
    idx = 0
    for i in range(n):
        for j in range(i+1, n):
            if b[idx]:
                A[i][j] = 1
            else:
                A[j][i] = 1
            idx += 1
    return A

def H(A, n, vertices=None):
    """Count Hamiltonian paths on given vertices."""
    if vertices is None:
        vertices = list(range(n))
    nv = len(vertices)
    if nv <= 1:
        return 1
    count = 0
    for perm in permutations(vertices):
        if all(A[perm[i]][perm[i+1]] for i in range(nv-1)):
            count += 1
    return count

def insertion_count(A, n, v, path):
    """Count valid positions to insert v into path to make a Ham path of T."""
    nv = len(path)
    count = 0
    # Position 0: v at start
    if A[v][path[0]]:
        count += 1
    # Position nv: v at end
    if A[path[-1]][v]:
        count += 1
    # Position j (1 <= j <= nv-1): v between path[j-1] and path[j]
    for j in range(1, nv):
        if A[path[j-1]][v] and A[v][path[j]]:
            count += 1
    return count

def main():
    n = 5
    m = n*(n-1)//2
    print(f"=== Insertion Count Analysis at n={n} ===\n")

    # For each tournament, compute delta_H via direct counting
    # and check against insertion counts
    total = 0
    match = 0

    for bits_int in range(2**m):
        A = tournament_from_bits(n, bits_int)
        H_T = H(A, n)

        for v in range(n):
            rem = [u for u in range(n) if u != v]
            H_Tv = H(A, n, rem)
            delta_H = H_T - H_Tv

            # Sum of insertion counts
            ins_sum = 0
            for perm in permutations(rem):
                if all(A[perm[i]][perm[i+1]] for i in range(len(rem)-1)):
                    ins_sum += insertion_count(A, n, v, list(perm))

            total += 1

            if delta_H == ins_sum:
                match += 1
            elif total <= 10:
                print(f"  bits={bits_int}, v={v}: delta_H={delta_H}, ins_sum={ins_sum}")

    print(f"\nTotal: {total}, Match: {match}")
    if match == total:
        print("delta_H = sum insact for ALL (T,v) at n=5!")
        print("\nThis confirms: H(T) = H(T-v) + sum_P' insact(v, P')")
        print("(trivial identity: each Ham path of T corresponds to")
        print(" inserting v into a Ham path of T-v at some position)")
    print()

    # Now investigate: insact vs inshat
    print("=== insact vs inshat comparison ===\n")

    mismatch_count = 0
    for bits_int in [0, 100, 200, 500, 920]:
        if bits_int >= 2**m:
            continue
        A = tournament_from_bits(n, bits_int)
        for v in range(n):
            rem = [u for u in range(n) if u != v]
            for perm in permutations(rem):
                if not all(A[perm[i]][perm[i+1]] for i in range(len(rem)-1)):
                    continue
                path = list(perm)
                insact_val = insertion_count(A, n, v, path)

                # inshat = boundary + typeI + typeII
                sig = [1 if A[v][path[j]] else 0 for j in range(len(path))]
                boundary = sig[0] + (1 - sig[-1])
                typeI = sum(1 for j in range(len(sig)-1) if sig[j]==0 and sig[j+1]==1)
                typeII = sum(1 for j in range(len(sig)-1) if sig[j]==1 and sig[j+1]==0)
                inshat_val = boundary + typeI + typeII

                if insact_val != inshat_val:
                    mismatch_count += 1
                    if mismatch_count <= 5:
                        print(f"  bits={bits_int}, v={v}, path={path}")
                        print(f"    insact={insact_val}, inshat={inshat_val}")
                        print(f"    sig={sig}, boundary={boundary}, tI={typeI}, tII={typeII}")

    if mismatch_count == 0:
        print("insact = inshat for ALL tested cases at n=5!")
        print("(This is expected: inshat IS the insertion count formula)")
    else:
        print(f"{mismatch_count} mismatches found")

if __name__ == "__main__":
    main()
