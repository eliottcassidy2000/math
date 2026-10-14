"""Marked Hamiltonian paths, prime necklace defects, and typed Burnside counts.

No spectrum-completeness or new primality theorem is inferred. Literal graph
APIs have explicit small caps; the proofs in the companion note are all-order.
Run normally or with -O. All checks remain active. Stdout is deterministic JSON.
"""
from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations, permutations
from math import factorial, gcd, isqrt
import json

CHECKS = 0


def require(test, message):
    global CHECKS
    CHECKS += 1
    if not test:
        raise ValueError(message)


def integer(x, lo, hi=None):
    if type(x) is not int or x < lo or (hi is not None and x > hi):
        raise ValueError("exact integer outside declared range")
    return x


def prime(p, cap=100000):
    integer(p, 3, cap)
    if p % 2 == 0 or any(p % d == 0 for d in range(3, isqrt(p)+1, 2)):
        raise ValueError("odd prime required")
    return p


def endpoint_rows(n):
    integer(n, 3, 64)
    rows = [sum(1 << j for j in range(i+1, n)) for i in range(n)]
    rows[0] ^= 1 << (n-1)
    rows[-1] |= 1
    return tuple(rows)


def tournament(rows, cap=64):
    if type(rows) is not tuple:
        raise ValueError("tuple of exact adjacency bit rows required")
    n = integer(len(rows), 1, cap)
    for i, row in enumerate(rows):
        integer(row, 0, (1 << n)-1)
        if row >> i & 1:
            raise ValueError("no loops")
        for j in range(i):
            if ((row >> j) & 1) + ((rows[j] >> i) & 1) != 1:
                raise ValueError("exactly one orientation per pair")
    return rows


def path_valid(rows, path):
    tournament(rows)
    n = len(rows)
    if type(path) is not tuple or len(path) != n:
        return False
    if any(type(v) is not int for v in path) or set(path) != set(range(n)):
        return False
    return all(rows[a] >> b & 1 for a, b in zip(path, path[1:]))


@dataclass(frozen=True)
class Address:
    order: int
    subset: int | None  # None is the unique path avoiding the reversed edge.


def address(a):
    if type(a) is not Address:
        raise ValueError("marked endpoint-family address required")
    integer(a.order, 3, 64)
    if a.subset is not None:
        integer(a.subset, 0, (1 << (a.order-2))-1)
    return a


def decode(a):
    address(a)
    n, mask = a.order, a.subset
    if mask is None:
        return tuple(range(n))
    left = tuple(i+1 for i in range(n-2) if mask >> i & 1)
    right = tuple(i+1 for i in range(n-2) if not (mask >> i & 1))
    return left + (n-1, 0) + right


def encode(n, path):
    rows = endpoint_rows(n)
    if not path_valid(rows, path):
        raise ValueError("not a native Hamiltonian path")
    if path == tuple(range(n)):
        return Address(n, None)
    cut = path.index(n-1)
    if cut+1 >= n or path[cut+1] != 0:
        raise ValueError("missing distinguished reversed edge")
    out = Address(n, sum(1 << (v-1) for v in path[:cut]))
    if decode(out) != path:
        raise ValueError("ordered path sidecar mismatch")
    return out


def rotate(a, shift=1):
    address(a)
    integer(shift, 0)
    if a.subset is None:
        return a
    m = a.order-2
    shift %= m
    mask = ((a.subset << shift) | (a.subset >> (m-shift))) & ((1 << m)-1)
    return Address(a.order, mask)


def transport(rows, path, shift=1):
    tournament(rows)
    if len(rows) < 3 or rows != endpoint_rows(len(rows)):
        raise ValueError("the original endpoint-family source is required")
    return decode(rotate(encode(len(rows), path), shift))


def count_paths(rows):
    tournament(rows, 12)
    n = len(rows)
    full = (1 << n)-1
    dp = [[0]*n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, full+1):
        for v, value in enumerate(dp[mask]):
            if value:
                choices = rows[v] & (full ^ mask)
                while choices:
                    bit = choices & -choices
                    choices -= bit
                    dp[mask | bit][bit.bit_length()-1] += value
    return sum(dp[full])


def graph_from_bits(n, mask):
    integer(n, 1, 12)
    integer(mask, 0, (1 << (n*(n-1)//2))-1)
    rows = [0]*n
    for k, (i, j) in enumerate(combinations(range(n), 2)):
        if mask >> k & 1:
            rows[i] |= 1 << j
        else:
            rows[j] |= 1 << i
    return tuple(rows)


def automorphism(rows, perm):
    tournament(rows)
    n = len(rows)
    if type(perm) is not tuple or len(perm) != n:
        return False
    if any(type(x) is not int for x in perm) or set(perm) != set(range(n)):
        return False
    return all(((rows[i] >> j) & 1) == ((rows[perm[i]] >> perm[j]) & 1)
               for i in range(n) for j in range(n))


def necklace_receipt(p):
    prime(p)
    q = ((1 << (p-1))-1)//p
    return {"prime": p, "tournament_order": p+2,
            "fixed_path_addresses": 3, "free_path_orbits": 2*q,
            "fermat_quotient_mod_p": q % p,
            "free_orbits_mod_p": (2*q) % p,
            "wieferich": q % p == 0}


def odd_partitions(n, lower=1):
    if n == 0:
        yield ()
    for part in range(lower, n+1, 2):
        for rest in odd_partitions(n-part, part):
            yield (part,) + rest


def permutation_term(n, parts):
    multiplicities = Counter(parts)
    z = 1
    for size, amount in multiplicities.items():
        z *= size**amount * factorial(amount)
    exponent = sum((a-1)//2 for a in parts)
    exponent += sum(gcd(a, b) for a, b in combinations(parts, 2))
    return factorial(n)//z * (1 << exponent)


def unlabelled_count(n):
    integer(n, 1, 23)
    value = Fraction(sum(permutation_term(n, parts) for parts in odd_partitions(n)),
                     factorial(n))
    if value.denominator != 1:
        raise ValueError("Burnside integrality failure")
    return int(value)


def fixed_orientation_count(perm):
    """Independent pair-orbit reader, including a possible orientation reversal."""
    n = len(perm)
    seen, cycles = set(), 0
    for pair in combinations(range(n), 2):
        if pair in seen:
            continue
        a, b = pair
        while True:
            seen.add(tuple(sorted((a, b))))
            a, b = perm[a], perm[b]
            if (a, b) == pair:
                break
            if (b, a) == pair:
                return 0
        cycles += 1
    return 1 << cycles


def burnside_defect(p):
    prime(p, 23)
    chi = 1 if pow(2, (p-1)//2, p) == 1 else -1
    f = ((1 << ((p-1)//2))-chi)//p
    w = (factorial(p-1)+1)//p
    rest = sum(permutation_term(p, parts) for parts in odd_partitions(p)
               if parts != (p,) and parts != (1,)*p)
    if rest % p:
        raise ValueError("remaining permutation classes must be p-divisible")
    return {"prime": p, "unlabelled_count": unlabelled_count(p),
            "legendre_2": chi, "euler_quotient_mod_p": f % p,
            "wilson_quotient_mod_p": w % p, "remaining_sum_mod_p": (rest//p) % p,
            "formula_mod_p": (f-chi*w-rest//p) % p}


def rejected(fn, *args):
    try:
        fn(*args)
    except (ValueError, TypeError):
        return True
    return False


def main():
    # Independent brute listings check every address through order eight.
    family = []
    for n in range(3, 13):
        rows = endpoint_rows(n)
        addresses = [Address(n, None)] + [Address(n, s) for s in range(1 << (n-2))]
        paths = {decode(a) for a in addresses}
        require(len(paths) == (1 << (n-2))+1, "path bijection")
        require(count_paths(rows) == len(paths), "independent subset DP")
        for a in addresses:
            path = decode(a)
            require(path_valid(rows, path), "native path")
            require(encode(n, path) == a, "source address round trip")
            require(transport(rows, path, n-2) == path, "cyclic action returns")
            require(rotate(rotate(a, 2), 3) == rotate(a, 5), "action law")
        if n <= 8:
            literal = {p for p in permutations(range(n)) if path_valid(rows, p)}
            require(literal == paths, "complete permutation enumeration")
        if n <= 7:
            aut = sum(automorphism(rows, p) for p in permutations(range(n)))
            require(aut == (3 if n == 3 else 1), "graph symmetry boundary")
        family.append({"order": n, "hamiltonian_paths": len(paths)})

    necklaces = []
    for p in (3, 5, 7, 11, 13):
        seen, sizes = set(), Counter()
        for mask in range(1 << p):
            if mask in seen:
                continue
            a = Address(p+2, mask)
            orbit = {rotate(a, j).subset for j in range(p)}
            seen |= orbit
            sizes[len(orbit)] += 1
        receipt = necklace_receipt(p)
        require(sizes[1] == 2 and set(sizes) == {1, p}, "prime subset orbit sizes")
        require(sizes[p] == receipt["free_path_orbits"], "Fermat quotient orbit count")
        necklaces.append({"prime": p, "subset_orbits_by_size": dict(sizes)})

    # An artificial action on paths is not the corresponding vertex permutation.
    a = Address(5, 1)
    source = decode(a)
    sigma = (0, 2, 3, 1, 4)
    literal = tuple(sigma[v] for v in source)
    require(not automorphism(endpoint_rows(5), sigma), "interior rotation is not an automorphism")
    require(not path_valid(endpoint_rows(5), literal), "unsorted literal path fails")
    require(path_valid(endpoint_rows(5), transport(endpoint_rows(5), source)), "resorted path action succeeds")
    require(literal != transport(endpoint_rows(5), source), "vertex and address actions differ")
    # At composite length four a nonconstant address has period two.
    a4 = Address(6, 5)
    require(rotate(a4, 2) == a4 and rotate(a4, 1) != a4, "composite-length hostile")

    # Complete labelled small universe; no inference beyond inherited all-order theorem.
    spectra = []
    for n in range(1, 7):
        histogram = Counter(count_paths(graph_from_bits(n, b))
                            for b in range(1 << (n*(n-1)//2)))
        require(7 not in histogram and 21 not in histogram, "small forbidden controls")
        require(all(h % 2 for h in histogram), "small Redei parity")
        spectra.append({"order": n, "labelled_tournaments": sum(histogram.values()),
                        "attained_counts": sorted(histogram)})
    h63 = (238, 220, 216, 16, 225, 206, 8, 72)
    require(count_paths(h63) == 63, "explicit seven-divisible count")
    padded = h63 + (255,)
    require(count_paths(padded) == 63, "non-strong H63 refutes gap tower")

    fermat = []
    for j in range(6):
        n = (1 << j)+2
        h = (1 << (n-2))+1
        residue = pow(2, h-1, h*h)
        require(residue % h == 1, "Fermat-number base-two congruence")
        quotient_residue = (residue-1)//h
        require(quotient_residue == pow(2, -j, h), "all Fermat-number quotient residue")
        fermat.append({"index": j, "order": n, "paths": h,
                       "base_two_quotient_mod_value": quotient_residue})
        require(h == (1 << (1 << j))+1, "Fermat subsequence")
    require(fermat[-1]["paths"] == 641*6700417, "composite Fermat value")
    arithmetic = [necklace_receipt(p) for p in (3, 5, 7, 17, 257, 1093, 3511, 65537)]
    for row in arithmetic:
        p = row["prime"]
        require(row["wieferich"] == (pow(2, p-1, p*p) == 1), "independent modular Wieferich test")
        require(row["free_orbits_mod_p"] == 2*row["fermat_quotient_mod_p"] % p, "correct residue observable")
        require(((1 << p)+1) % p == 3 % p, "path count itself has different residue")
        row.pop("free_path_orbits")  # Compact output; the exact API retains it.

    known = [1, 1, 2, 4, 12, 56, 456, 6880]
    require([unlabelled_count(n) for n in range(1, 9)] == known, "Burnside count controls")
    for n in range(1, 7):
        total = sum(fixed_orientation_count(p) for p in permutations(range(n)))
        require(total == factorial(n)*unlabelled_count(n), "independent pair-orbit Burnside reader")
    burnside = [burnside_defect(p) for p in (3, 5, 7, 11, 13, 17, 19, 23)]
    for row in burnside:
        p = row["prime"]
        require(row["formula_mod_p"] == row["unlabelled_count"] % p, "THM307 full residual formula")
        q = necklace_receipt(p)["fermat_quotient_mod_p"]
        require(2*row["legendre_2"]*row["euler_quotient_mod_p"] % p == q,
                "Euler-Fermat quotient relation")
    require(burnside_defect(5)["remaining_sum_mod_p"] == 4, "remaining Burnside term cannot be dropped")
    require(unlabelled_count(7) == 456 and count_paths(endpoint_rows(7)) == 33,
            "class count versus path count")

    malformed = [
        (endpoint_rows, (True,)), (endpoint_rows, (2,)), (endpoint_rows, (65,)),
        (endpoint_rows, (5.0,)),
        (decode, (Address(5, True),)), (decode, (Address(5, 8),)),
        (decode, (Address(5, 1.0),)),
        (rotate, (Address(5, 1), -1)), (encode, (5, (0, 1, 2, 3, True))),
        (transport, (graph_from_bits(5, 0), source)),
        (necklace_receipt, (True,)), (necklace_receipt, (9,)),
        (necklace_receipt, (2,)), (burnside_defect, (29,)),
        (count_paths, ((False,),)), (tournament, ([0],)),
    ]
    for fn, args in malformed:
        require(rejected(fn, *args), "malformed typed receipt rejected")
    result = {"status": "PROVED scoped elementary bridges; inherited H!=7,21; full odd spectrum OPEN",
              "checks": CHECKS, "endpoint_family": family, "prime_necklace_controls": necklaces,
              "finite_labelled_spectra": spectra, "H63_rows": h63,
              "Fermat_subsequence": fermat, "prime_path_defects": arithmetic,
              "Burnside_residuals": burnside,
              "nonautomorphism_hostile": {"path": source, "literal_vertex_image": literal,
                                           "actual_address_image": transport(endpoint_rows(5), source)},
              "scope": "no aliquot graph is converted to a tournament; no new spectrum completeness"}
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
