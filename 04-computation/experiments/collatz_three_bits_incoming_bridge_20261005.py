#!/usr/bin/env python3
"""Exact integration of e3b54e1376 with the critical-flow construction."""
from fractions import Fraction as F
from hashlib import sha256
from math import isqrt
from pathlib import Path
import json

checks = 0


def require(condition, witness=None):
    global checks
    checks += 1
    if not condition:
        raise RuntimeError(witness)


def v2(n):
    return (n & -n).bit_length() - 1


def U(n):
    z = 3 * n + 1
    return z >> v2(z)


def nu(n):
    return F(8, 3 * 4 ** n.bit_length())


def mu(n):
    return F(2, 4 ** ((n + 1) // 2).bit_length())


def obs(n):
    m = n + 1
    b = (3 * m - isqrt(5 * m * m) - 1) // 2
    fs = [1, 2]
    while fs[-1] <= n:
        fs.append(fs[-1] + fs[-2])
    rem, lowest = n, 0
    for f in reversed(fs):
        d = int(f <= rem)
        rem -= d * f
        if f == 1:
            lowest = d
    return n.bit_length(), v2(n + 1), n % 2, b % 2, lowest


def main():
    for d in range(1, 14):
        images = set()
        for n in range(1, 2 ** d, 2):
            a = d - n.bit_length() + 1
            y = 2 ** a * n - 1
            m = (y + 1) // 2
            images.add(m)
            require(v2(y + 1) == a and (y + 1) >> a == n, (d, n))
            require(mu(y) == F(2, 4 ** d) == nu(n) * F(3, 4 ** a), (d, n))
            require(m.bit_length() == d, (d, n))
            if a > 1:
                require(U(y) == 2 ** (a - 1) * (3 * n) - 1, (d, n))
        require(images == set(range(2 ** (d - 1), 2 ** d)), d)
    for n in range(1, 1 << 14, 2):
        h = v2(n + 1)
        t = (n + 1) >> h
        require(mu(n) == F(3, 4 ** h) * nu(t), n)
    path = [107, 161, 121, 91]
    require([U(n) for n in path[:-1]] == path[1:])
    require(obs(path[0]) == obs(path[-1]) == (7, 2, 1, 1, 0))
    require(U(429) == 161 and 429 not in path)
    require(nu(107) == nu(91) and mu(107) == mu(91))
    require(U(37) == 7)

    # Independently validate the supplied 88-base graph and reconstruct the
    # equality weights for rho=1. Only this finite supported union is covered.
    kernel_path = Path(__file__).resolve().parents[2] / "05-knowledge/results/collatz_three_bit_sibling_flow_20261005.json"
    payload = json.loads(kernel_path.read_text())
    entries = payload["bases"]
    by_base = {row["b"]: row for row in entries}
    require(len(by_base) == len(entries) == 88)
    require(by_base[1] == {"b": 1, "parent": None, "depth": 0, "valuation": 2})
    weights = {1: F(1)}
    lengths, depths = {1: 0}, {1: 0}
    for b in sorted(by_base):
        row = by_base[b]
        require(type(b) is int and b > 0 and b % 2 == 1)
        require(type(row["depth"]) is int and row["depth"] >= 0)
        require(type(row["valuation"]) is int and v2(3 * b + 1) == row["valuation"])
        require(row["valuation"] in (1, 2))
        if b == 1:
            continue
        parent, k = row["parent"], row["depth"]
        require(type(parent) is int and parent in by_base)
        require(U(b) == 4 ** k * parent + (4 ** k - 1) // 3, b)
    for b in sorted(by_base):
        path = []
        c = b
        while c not in weights:
            require(c not in path, (b, c))
            path.append(c)
            c = by_base[c]["parent"]
        for c in reversed(path):
            row = by_base[c]
            p, k = row["parent"], row["depth"]
            lengths[c] = lengths[p] + 1
            depths[c] = depths[p] + k
            weights[c] = F(15, 16) ** lengths[c] * F(1, 16) ** depths[c]
            require(weights[c] == F(15, 16) * F(1, 16) ** k * weights[p], c)
    total = F(16, 15) * sum(weights.values(), F())
    for b, row in by_base.items():
        if b > 1:
            incoming = weights[b] * F(16, 15)
            target_mass = F(1, 16) ** row["depth"] * weights[row["parent"]]
            require(incoming == target_mass, b)
    capacity_rows = []
    for r in (F(1, 16), F(1, 4), F(1, 2)):
        denominator = 1 + r + r * r
        kappa = 1 - r * r / denominator
        for forbidden in range(3):
            exact = 1 - r ** forbidden / denominator
            partial = sum(((1-r) * r ** k for k in range(36) if k % 3 != forbidden), F())
            require(partial == exact * (1 - r ** 36))
            require(exact <= kappa < 1)
        capacity_rows.append({"r": str(r), "kappa": str(kappa),
                              "rooted_base_mass_bound": str(1 / (1 - kappa))})
    for c in range(1, 1024, 2):
        if v2(3*c+1) not in (1, 2):
            continue
        for k in range(13):
            y = 4 ** k * c + (4 ** k - 1) // 3
            require(y % 3 == (c + k) % 3)
            if y % 3 == 0:
                continue
            b = ((2*y-1)//3) if y % 6 == 5 else ((4*y-1)//3)
            require(U(b) == y and v2(3*b+1) in (1, 2))
            require((b == 1) == (c == 1 and k == 0))
    require(sum(weights.values(), F()) <= 273)
    require(U(9) == 7 and F(1,9)/(F(3,4)*F(1,7)) == F(28,27))
    report = {
        "status": "FINITE-EXACT integration; universal coverage remains OPEN",
        "source_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
        "kernel_sha256": sha256(kernel_path.read_bytes()).hexdigest(),
        "coding_universe": "all mixture descriptions at heights 1..13; all gamma sources below 2^14",
        "coding_bijection": "(D,n) -> Y=2^(D-bit_length(n)+1)*n-1",
        "plateau": {"path": [107, 161, 121, 91], "observer": list(obs(107)),
                    "additional_predecessor": [429, 161], "nu_endpoint_atom": str(nu(107))},
        "critical_branch_capacity": capacity_rows,
        "size_carry_hostile": {"base_edge": [9,7], "r": "1/4", "rho": "1", "s": 1,
                                "source_to_allowed_weight_ratio": "28/27"},
        "critical_finite_kernel": {"bases": 88, "r": "1/16", "rho": "1",
            "root_weight_before_normalization": "1", "total_full_sibling_flow_mass": str(total),
            "maximum_base_path_length": max(lengths.values()),
            "weights": [{"b": b, "g": str(weights[b]), "L": lengths[b], "K": depths[b]}
                        for b in sorted(weights)],
            "scope": "same certified infinite sibling union; no new convergence coverage"},
        "checks_total": checks}
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
