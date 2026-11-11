"""Exact small controls; run from any working directory. No dependencies."""
from hashlib import sha256
import json
from math import gcd, isqrt
from pathlib import Path


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def prime(n):
    return n >= 2 and all(n % d for d in range(2, isqrt(n) + 1))


def sieve(limit):
    flags = bytearray(b"\x01") * (limit + 1)
    flags[0:2] = b"\x00\x00"
    for d in range(2, isqrt(limit) + 1):
        if flags[d]:
            flags[d*d:limit+1:d] = b"\x00" * ((limit-d*d)//d+1)
    return flags


def main():
    P = 30030
    flags = sieve(P + 289)
    for d in range(1, 290):
        for sign in (-1, 1):
            value = P + sign*d
            require(prime(value) == bool(flags[value]), "primality paths differ")
    F = next(d for d in range(2, 290) if flags[P+d])
    L = next(d for d in range(2, 290) if flags[P-d])
    require(F == L == 17, "first-hit offset")
    require(flags[30029] and not flags[30031], "offset-one boundary")
    require(30031 == 59*509 and 30049 == 151*199, "endpoint factors")
    require(flags[30319] and not prime(289), "strict-square hostile")
    require(flags[30013] and 30013 == 122**2+123**2, "Gaussian norm")
    require(2*30013 == 245**2+1, "twice-prime identity")
    require(173**2 < 30012 < 174**2, "not a Landau value")
    require((123+122, 123-122) == (245, 1), "Gaussian multiplication")
    require(6**2+1 == 37 and prime(37), "positive control")
    require(30**2+1 == 17*53, "small quadratic composite")
    require(P**2+1 == 73*12353437, "large quadratic composite")
    require(prime(73) and prime(12353437), "factor primality")
    require(gcd(P**2+1, P) == 1, "wheel survival")

    quadratic = [n*n+1 for n in range(1, 175) if prime(n*n+1)]
    quadratic2 = [n*n+1 for n in range(1, 175) if flags[n*n+1]]
    require(quadratic == quadratic2, "quadratic primality paths differ")
    require(quadratic[:7] == [2, 5, 17, 37, 101, 197, 257], "quadratic prefix")
    roots = {}
    for q in range(2, 98):
        if prime(q):
            roots[q] = sum((r*r+1) % q == 0 for r in range(q))
            expected = 1 if q == 2 else (2 if q % 4 == 1 else 0)
            require(roots[q] == expected, "root-count formula")

    cases = 0
    admissible = 0
    for W in range(1, 31):
        for a in range(1, W+1):
            A = a*a+1
            for t in range(1, 6):
                n = a+t*W*A
                B = 1+2*a*t*W+t*t*W*W*A
                require(n % W == a % W and n*n+1 == A*B, "congruence identity")
                require(A > 1 and B > 1, "proper factors")
                if gcd(A, W) == 1:
                    require(gcd(n*n+1, W) == 1, "admissible class lost")
                    admissible += 1
                cases += 1
    result = {
        "status": "FINITE-EXACT controls passed; infinite claims use the written proof",
        "script_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
        "primorial": P, "F_6": F, "L_6": L,
        "nearest_below": 30029, "nearest_above": 30047, "crossing_gap": 18,
        "primality_crosschecks": 578,
        "quadratic_n_range": [1, 174], "quadratic_prime_values": quadratic,
        "quadratic_root_counts": roots,
        "composite_construction_cases": cases,
        "locally_admissible_composite_cases": admissible,
        "factorization": [901800901, 73, 12353437],
    }
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
