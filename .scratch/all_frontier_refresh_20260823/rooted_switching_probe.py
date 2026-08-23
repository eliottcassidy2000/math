#!/usr/bin/env python3
"""Exact rooted-switching hostile for the all-frontier refresh.

Find the smallest labelled tournament switching class containing two strong
representatives with the same unrooted pair (H, disc) but different fixed-root
odd Pfaffian energies.  This is a bounded diagnostic for the shared
"quotient needs a moving root/owner" pattern; it proves no tournament, LRC,
or factorial conjecture.
"""

from fractions import Fraction
from itertools import permutations


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def det_bareiss(matrix):
    """Fraction-free exact determinant of an integer matrix."""
    n = len(matrix)
    if n == 0:
        return 1
    a = [row[:] for row in matrix]
    sign = 1
    previous = 1
    for k in range(n - 1):
        pivot_row = next((i for i in range(k, n) if a[i][k] != 0), None)
        if pivot_row is None:
            return 0
        if pivot_row != k:
            a[k], a[pivot_row] = a[pivot_row], a[k]
            sign = -sign
        pivot = a[k][k]
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                numerator = a[i][j] * pivot - a[i][k] * a[k][j]
                require(numerator % previous == 0, "Bareiss division failed")
                a[i][j] = numerator // previous
        previous = pivot
        for i in range(k + 1, n):
            a[i][k] = 0
    return sign * a[n - 1][n - 1]


def identity_plus(k):
    n = len(k)
    return [[(1 if i == j else 0) + k[i][j] for j in range(n)] for i in range(n)]


def tournament_from_bits(n, bits):
    k = [[0] * n for _ in range(n)]
    e = 0
    for i in range(n):
        for j in range(i + 1, n):
            value = 1 if (bits >> e) & 1 else -1
            k[i][j] = value
            k[j][i] = -value
            e += 1
    return k


def upper_signature(k):
    return tuple(k[i][j] for i in range(len(k)) for j in range(i + 1, len(k)))


def switch(k, signs):
    n = len(k)
    return [[signs[i] * k[i][j] * signs[j] for j in range(n)] for i in range(n)]


def switching_key(k):
    n = len(k)
    signatures = []
    for mask in range(1 << (n - 1)):
        signs = [1] + [(-1 if (mask >> (i - 1)) & 1 else 1) for i in range(1, n)]
        signatures.append(upper_signature(switch(k, signs)))
    return min(signatures)


def is_strong(k):
    n = len(k)
    for start in range(n):
        seen = {start}
        stack = [start]
        while stack:
            i = stack.pop()
            for j in range(n):
                if k[i][j] == 1 and j not in seen:
                    seen.add(j)
                    stack.append(j)
        if len(seen) != n:
            return False
    return True


def hamiltonian_paths(k):
    n = len(k)
    return sum(
        all(k[p[i]][p[i + 1]] == 1 for i in range(n - 1))
        for p in permutations(range(n))
    )


def disc(k):
    n = len(k)
    return Fraction(det_bareiss(identity_plus(k)), 1 << (n - 1))


def rooted_odd_energy(k, root):
    n = len(k)
    khat = [k[i][:] + [root[i]] for i in range(n)]
    khat.append([-root[i] for i in range(n)] + [0])
    raw = det_bareiss(identity_plus(khat)) - det_bareiss(identity_plus(k))
    return Fraction(raw, 1 << (n - 1))


def pfaffian(matrix):
    """Independent recursive Pfaffian for the final small hostile."""
    n = len(matrix)
    require(n % 2 == 0, "Pfaffian matrix must have even order")
    if n == 0:
        return 1
    total = 0
    for j in range(1, n):
        minor = [
            [matrix[r][c] for c in range(n) if c not in (0, j)]
            for r in range(n)
            if r not in (0, j)
        ]
        total += (1 if j % 2 == 1 else -1) * matrix[0][j] * pfaffian(minor)
    return total


def rooted_odd_energy_pfaffian(k, root):
    n = len(k)
    total = 0
    for mask in range(1 << n):
        subset = [i for i in range(n) if (mask >> i) & 1]
        if len(subset) % 2 == 0:
            continue
        augmented = [[k[i][j] for j in subset] + [root[i]] for i in subset]
        augmented.append([-root[j] for j in subset] + [0])
        total += pfaffian(augmented) ** 2
    return Fraction(total, 1 << (n - 1))


def arcs(k):
    return [(i + 1, j + 1) for i in range(len(k)) for j in range(i + 1, len(k)) if k[i][j] == 1] + [
        (j + 1, i + 1) for i in range(len(k)) for j in range(i + 1, len(k)) if k[i][j] == -1
    ]


def find_hostile():
    fallback = None
    for n in range(3, 7):
        classes = {}
        edge_count = n * (n - 1) // 2
        for bits in range(1 << edge_count):
            k = tournament_from_bits(n, bits)
            if not is_strong(k):
                continue
            item = {
                "bits": bits,
                "k": k,
                "H": hamiltonian_paths(k),
                "disc": disc(k),
                "E": rooted_odd_energy(k, [1] * n),
            }
            classes.setdefault(switching_key(k), []).append(item)
        for items in classes.values():
            by_unrooted = {}
            for item in items:
                by_unrooted.setdefault((item["H"], item["disc"]), []).append(item)
            for same in by_unrooted.values():
                values = {}
                for item in same:
                    values.setdefault(item["E"], item)
                if len(values) >= 2:
                    pair = list(values.values())[:2]
                    return n, pair[0], pair[1], "same H and disc"
            if fallback is None:
                energies = {}
                for item in items:
                    energies.setdefault(item["E"], item)
                if len(energies) >= 2:
                    pair = list(energies.values())[:2]
                    fallback = (n, pair[0], pair[1], "same switching class and disc")
    if fallback is not None:
        return fallback
    raise RuntimeError("No rooted switching hostile found through n=6")


def switching_signs_between(k1, k2):
    n = len(k1)
    for mask in range(1 << (n - 1)):
        signs = [1] + [(-1 if (mask >> (i - 1)) & 1 else 1) for i in range(1, n)]
        if switch(k1, signs) == k2:
            return signs
    raise RuntimeError("Representatives were not switching equivalent")


def main():
    n, first, second, strength = find_hostile()
    signs = switching_signs_between(first["k"], second["k"])
    covariance = rooted_odd_energy(first["k"], signs)
    require(covariance == second["E"], "Rooted switching covariance failed")
    require(first["disc"] == second["disc"], "Switching changed discriminant")
    require(first["E"] != second["E"], "Fixed-root energies did not separate")
    require(is_strong(first["k"]) and is_strong(second["k"]), "Hostile must remain in strong locus")
    require(rooted_odd_energy_pfaffian(first["k"], [1] * n) == first["E"], "first Pfaffian audit failed")
    require(rooted_odd_energy_pfaffian(second["k"], [1] * n) == second["E"], "second Pfaffian audit failed")

    print("ROOTED SWITCHING / OWNER HOSTILE")
    print(f"minimal searched order: n={n}")
    print(f"separation strength: {strength}")
    print(f"switch signs D: {tuple(signs)}")
    print(f"shared discriminant: {first['disc']}")
    print(f"first:  H={first['H']}, E_odd(1)={first['E']}, arcs={sorted(arcs(first['k']))}")
    print(f"second: H={second['H']}, E_odd(1)={second['E']}, arcs={sorted(arcs(second['k']))}")
    print(f"covariant check E_odd(K,D1)=E_odd(DKD,1): {covariance}")
    print("independent odd-subset Pfaffian-square reconstruction: pass")
    print("VERDICT: the unrooted switching data do not determine the fixed-root response.")
    print("BOUNDARY: this is a finite exact observer hostile, not a proof of H>=disc or any LRC/FC claim.")


if __name__ == "__main__":
    main()
