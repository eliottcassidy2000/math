"""Exact floor-sum / nilpotent / zero-divisor controls; standard library only.

Run from any directory; the JSON certificate is written beside this source.
All theorem claims are proved in the companion Markdown, not inferred from
these finite checks. Checks remain active under Python -O.
"""
from fractions import Fraction
from itertools import product
from math import gcd, prod
from pathlib import Path
import hashlib
import json


def check(condition, message):
    if not condition:
        raise RuntimeError(message)


def factor(n):
    factors = {}
    p = 2
    while p*p <= n:
        while n % p == 0:
            factors[p] = factors.get(p, 0)+1
            n //= p
        p += 1
    if n > 1:
        factors[n] = 1
    return factors


def zero_power_count(n, d):
    return prod(p**(a-(a+d-1)//d) for p, a in factor(n).items())


def zero_pair_count(n):
    total = prod(p**(a-1)*((a+1)*p-a) for p, a in factor(n).items())
    return total-2*n+1


def nth_root_floor(x, d, upper):
    """Certified root by integer binary search, with a verified upper bound."""
    lo, hi = 0, upper
    check(hi**d > x, "root upper bound")
    while lo+1 < hi:
        mid = (lo+hi)//2
        if mid**d <= x:
            lo = mid
        else:
            hi = mid
    check(lo**d <= x < (lo+1)**d, "root boundary")
    return lo


def main():
    certificate = {
        "status": "FINITE-EXACT controls for independently proved all-modulus identities",
        "source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "universes": {},
    }
    odd_power_pairs = 0
    odd_power_residues = 0
    for n in range(2, 501):
        squarefree = all(a == 1 for a in factor(n).values())
        for d in (1, 3, 5, 7, 9, 11):
            residues = [pow(k, d, n) for k in range(n)]
            count = zero_power_count(n, d)
            check(residues.count(0) == count, "power-zero CRT formula")
            floor_sum = sum(k**d//n for k in range(1, n))
            baseline = Fraction(sum(k**d for k in range(1, n)), n)-Fraction(n-1, 2)
            check(floor_sum == baseline+Fraction(count-1, 2), "odd power defect")
            check(sum(residues) == n*(n-count)//2, "antipodal residue sum")
            if d >= 3:
                check((floor_sum == baseline) == squarefree, "squarefree iff")
            else:
                check(floor_sum == 0 and count == 1, "degree-one hostile")
            odd_power_pairs += 1
            odd_power_residues += n
    certificate["universes"]["odd_power"] = {
        "moduli": "2<=n<=500", "degrees": [1, 3, 5, 7, 9, 11],
        "parameter_pairs": odd_power_pairs, "residue_evaluations": odd_power_residues,
    }

    inverse_pairs = 0
    direct_root_terms = 0
    # Directly enumerate root terms; no reciprocity formula is used here.
    for d, max_n in ((1, 30), (2, 30), (3, 80), (5, 12), (7, 8), (9, 6)):
        for n in range(2, max_n+1):
            ceiling = (n-1)**d//n
            root_sum = sum(nth_root_floor(n*k, d, n-1) for k in range(1, ceiling+1))
            power_sum = sum(k**d//n for k in range(1, n))
            count = zero_power_count(n, d)
            check(root_sum+power_sum == (n-1)*ceiling+count-1, "lattice reciprocity")
            if d % 2:
                check(ceiling == ((n-1)**d-(n-1))//n, "odd inverse rectangle")
            if d == 3:
                check(ceiling == (n-1)*(n-2), "cube rectangle")
                baseline = Fraction((3*n-5)*(n-1)*(n-2), 4)
                check(root_sum == baseline+Fraction(count-1, 2), "cube-root defect")
            inverse_pairs += 1
            direct_root_terms += ceiling
    certificate["universes"]["inverse_roots"] = {
        "degree_and_max_modulus": [[1, 30], [2, 30], [3, 80], [5, 12], [7, 8], [9, 6]],
        "parameter_pairs": inverse_pairs, "direct_root_terms": direct_root_terms,
        "method": "integer binary search for every individual root term",
    }

    pair_entries = 0
    for n in range(2, 151):
        actual = 0
        zero_pairs = 0
        for i in range(1, n):
            row = sum(i*j//n for j in range(1, n))
            check(row == Fraction((i-1)*(n-1)+gcd(i, n)-1, 2), "row floor formula")
            actual += row
            zero_pairs += sum(i*j % n == 0 for j in range(1, n))
        predicted = zero_pair_count(n)
        check(zero_pairs == predicted == sum(gcd(i, n)-1 for i in range(1, n)),
              "zero-pair independent counts")
        baseline = Fraction((n-2)*(n-1)**2, 4)
        check(actual == baseline+Fraction(zero_pairs, 2), "product floor defect")
        prime = factor(n) == {n: 1}
        check((actual == baseline) == prime, "primality iff")
        pair_entries += (n-1)**2
    certificate["universes"]["product_pairs"] = {
        "moduli": "2<=n<=150", "ordered_pair_entries": pair_entries,
        "methods": ["direct row floor sums", "direct modular zero pairs", "gcd sum", "CRT factor product"],
    }

    higher_entries = 0
    for arity, max_n in ((3, 25), (4, 12)):
        for n in range(2, max_n+1):
            actual = zeros = 0
            for coordinates in product(range(1, n), repeat=arity):
                value = prod(coordinates)
                actual += value//n
                zeros += value % n == 0
                higher_entries += 1
            baseline = Fraction((n-1)**arity*(n**(arity-1)-2**(arity-1)), 2**arity)
            check(actual == baseline+Fraction(zeros, 2), "higher product involution")
    certificate["universes"]["higher_products"] = {
        "arity_and_max_modulus": [[3, 25], [4, 12]], "tuple_entries": higher_entries,
    }

    examples = []
    for n in (2, 3, 4, 6, 7, 8, 9, 12, 16, 18, 27, 30, 36, 60, 64, 72):
        root_count, zero_pairs = zero_power_count(n, 3), zero_pair_count(n)
        cube_base = Fraction((n-2)*(n-1)*(n+1), 4)
        inverse_base = Fraction((3*n-5)*(n-2)*(n-1), 4)
        product_base = Fraction((n-2)*(n-1)**2, 4)
        examples.append({
            "n": n, "factorization": factor(n), "cube_zero_roots_including_zero": root_count,
            "nonzero_zero_product_pairs": zero_pairs,
            "cube_floor_sum": int(cube_base+Fraction(root_count-1, 2)),
            "cube_root_floor_sum": int(inverse_base+Fraction(root_count-1, 2)),
            "product_floor_sum": int(product_base+Fraction(zero_pairs, 2)),
            "cube_and_root_defect": str(Fraction(root_count-1, 2)),
            "product_defect": str(Fraction(zero_pairs, 2)),
        })
    certificate["examples"] = examples
    check(sum(k*k//3 for k in range(1, 3)) != Fraction(1+4, 3)-1, "even-power hostile")
    cube_residues_7 = [pow(k, 3, 7) for k in range(1, 7)]
    check(len(set(cube_residues_7)) == 2, "cube map is not a permutation hostile")
    certificate["hostiles"] = {
        "literal_cube_p3": {"without_floor": "3", "claimed_rhs": "2"},
        "literal_product_p3": {"without_floor": "3", "claimed_rhs": "1"},
        "literal_root_p3": "cuberoot(3)+cuberoot(6)>2=claimed_rhs",
        "squarefree_composite_separator": 6,
        "first_composite_defect": 4,
        "even_power_antipodal_failure": {"n": 3, "d": 2},
        "nonpermutation_prime_cube_residues": {"n": 7, "residues": cube_residues_7},
        "degree_one_squarefree_detection_failure": "degree one has zero defect at every modulus",
    }
    squarefree = lambda n: all(a == 1 for a in factor(n).values())
    def accelerated_collatz(n):
        image = 3*n+1
        return image//(image & -image)
    loss = next((n, accelerated_collatz(n)) for n in range(3, 501, 2)
                if squarefree(n) and not squarefree(accelerated_collatz(n)))
    gain = next((n, accelerated_collatz(n)) for n in range(3, 501, 2)
                if not squarefree(n) and squarefree(accelerated_collatz(n)))
    check(loss == (33, 25) and gain == (9, 7), "squarefree status is not monotone")
    certificate["hostiles"]["first_odd_collatz_loss_of_squarefree"] = loss
    certificate["hostiles"]["first_odd_collatz_gain_of_squarefree"] = gain
    certificate["result"] = "PASS"
    path = Path(__file__).with_suffix(".json")
    path.write_text(json.dumps(certificate, indent=2, sort_keys=True)+"\n", encoding="utf-8")
    print(f"odd-power identities: {odd_power_pairs} pairs, {odd_power_residues} residues; PASS")
    print(f"inverse lattice identities: {inverse_pairs} pairs, {direct_root_terms} direct root terms; PASS")
    print(f"product floor identities: {pair_entries} ordered pair entries; PASS")
    print(f"higher product identities: {higher_entries} tuple entries; PASS")
    print("hostiles: literal sums, n=4, n=6, even powers, nonpermuting cubes modulo 7; PASS")
    print(f"certificate: {path.name}")


if __name__ == "__main__":
    main()
