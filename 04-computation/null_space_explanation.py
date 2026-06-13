#!/usr/bin/env python3
"""
WHY the null space exists: g(t) depends only on S = sum(len_i - 1).

The correction function for an independent set of cycle type (l_1,...,l_p) is:
  g(t) = A_{f+1}(t) * (t-1)^{d-f}
where f = d - S, S = sum(l_i - 1), d = n-1.

So ANY two cycle types with the same S get the SAME correction function g(t).

At n=7 (d=6):
  S=2: (3,) only. f=4, g=A_5(t)*(t-1)^2
  S=4: (5,) AND (3,3). f=2, g=A_3(t)*(t-1)^4
  S=6: (7,) only. f=0, g=(t-1)^6

The S=4 degeneracy explains the null vector (0,-2,0,1):
  2*(5-cycle contribution) = (3,3)-pair contribution in the G_T formula.

At n=9 (d=8):
  S=2: (3,)
  S=4: (5,), (3,3)           <-- same g(t)
  S=6: (7,), (3,5)           <-- same g(t)!
  S=8: (9,), (3,7), (5,5), (3,3,3)  <-- all same g(t)!

This predicts the null space dimension at each n!

opus-2026-03-07-S34
"""
from math import comb

def eulerian_number(n, k):
    return sum((-1)**j * comb(n+1, j) * (k+1-j)**n for j in range(k+1))

def eulerian_poly_eval(n, t):
    return sum(eulerian_number(n, k) * t**k for k in range(n))

def cycle_types_by_S(n):
    """Find all possible independent set cycle-type partitions grouped by S."""
    d = n - 1
    max_S = d

    S_groups = {}
    for S in range(2, max_S + 1, 2):
        types = []
        find_types(S, 3, n, [], types)
        if types:
            S_groups[S] = types
    return S_groups

def find_types(S, min_k, max_verts, current, result):
    """Find all ways to partition S into (l_i - 1) values where l_i >= 3 odd."""
    if S == 0:
        result.append(tuple(sorted(current, reverse=True)))
        return
    if S < 2:
        return
    total_verts = sum(current)  # total vertices used by disjoint cycles
    for k in range(min_k, min(S + 1, max_verts - total_verts) + 1, 2):
        if k - 1 <= S:
            new_verts = total_verts + k
            if new_verts <= max_verts:
                find_types(S - (k - 1), k, max_verts, current + [k], result)

def main():
    print("=" * 70)
    print("NULL SPACE STRUCTURE: WHY CYCLE TYPES ARE INVISIBLE TO G_T")
    print("=" * 70)

    for n in [5, 7, 9, 11, 13]:
        d = n - 1
        print(f"\n{'='*50}")
        print(f"n = {n}, d = {d}")
        print(f"{'='*50}")

        S_groups = cycle_types_by_S(n)

        # Count invariants
        total_invariants = sum(len(types) for types in S_groups.values())
        num_S_values = len(S_groups)

        print(f"\nCycle types grouped by S = sum(l_i - 1):")
        for S in sorted(S_groups.keys()):
            types = S_groups[S]
            type_strs = [str(t) for t in types]
            f = d - S
            print(f"  S={S} (f={f}): {', '.join(type_strs)}")

        print(f"\nTotal invariant types: {total_invariants}")
        print(f"Distinct S values: {num_S_values}")
        print(f"Null space dimension: {total_invariants - num_S_values}")
        print(f"  (these many cycle-type combinations are invisible to G_T)")

        # At x=2, the deformed Eulerian polynomial a_k has independent equations
        # equal to floor((n-1)/2) (palindromy halves them)
        num_ak = n - 1  # coefficients a_0 through a_{n-2}
        independent_ak = (n - 1) // 2 + 1  # palindromy: a_k = a_{n-1-k}
        # But a_0 = 1 always, so effectively
        # Wait: a_k = A(n,k) + corrections. The corrections are determined by
        # the cycle invariants. Palindromy reduces the independent equations.

        print(f"\nForward-edge distribution: {num_ak} coefficients a_0,...,a_{n-2}")
        print(f"After palindromy: {independent_ak} independent")
        print(f"Each S-group gives 1 equation -> {num_S_values} equations")
        print(f"But each S-group contributes to multiple a_k values")

    print(f"\n\n{'='*70}")
    print(f"KEY THEORETICAL RESULT")
    print(f"{'='*70}")
    print()
    print("G_T(t,x) can only 'see' cycle-type information through")
    print("the aggregate S = sum(len_i - 1) = total excess.")
    print()
    print("This means G_T(t,x) is STRICTLY COARSER than GS's U_T.")
    print("U_T separates all cycle types; G_T only separates by S.")
    print()
    print("However, G_T(t,x) adds EULERIAN polynomial structure (t variable)")
    print("that U_T's power-sum expansion doesn't directly capture.")
    print()
    print("The null space of the a_k -> invariant map at each n is:")
    for n in [5, 7, 9, 11, 13]:
        S_groups = cycle_types_by_S(n)
        total = sum(len(t) for t in S_groups.values())
        distinct = len(S_groups)
        print(f"  n={n}: null dim = {total - distinct}")

    print()
    print("OPEN QUESTION: Can we define a 'typed G_T' that separates all")
    print("cycle lengths AND captures descent structure?")
    print("  G_T^typed(t; y_3, y_5, ...) uses individual y_k weights")
    print("  but the correction function g(t) is STILL determined by S alone.")
    print("  So typed G_T can't distinguish cycles with same S in t-coefficients.")
    print("  The only way to separate is via DIFFERENT x-weights (y_k != y_l).")

if __name__ == "__main__":
    main()
