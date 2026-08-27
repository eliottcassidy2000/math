#!/usr/bin/env python3
"""Exact scout for the 2-adic reflection-orbit tower of Sun 2-4-6-8 fibres."""

from __future__ import annotations

from hashlib import sha256
from math import comb


RANKS = (2, 4, 6, 8)
TARGET = 896_315_812_331_399


def factor(n: int) -> list[tuple[int, int]]:
    ans: list[tuple[int, int]] = []
    p = 2
    while p * p <= n:
        if n % p:
            p = 3 if p == 2 else p + 2
            continue
        a = 0
        while n % p == 0:
            n //= p
            a += 1
        ans.append((p, a))
        p = 3 if p == 2 else p + 2
    if n > 1:
        ans.append((n, 1))
    return ans


def floor_log_p(k: int, p: int) -> int:
    e = 0
    while k >= p:
        k //= p
        e += 1
    return e


def period(m: int, k: int) -> int:
    ans = 1
    for p, a in factor(m):
        ans *= p ** (a + floor_log_p(k, p))
    return ans


def b(n: int, k: int) -> int:
    return comb(n, k) if n >= k else 0


def circular_convolution(a: list[int], c: list[int]) -> list[int]:
    m = len(a)
    assert len(c) == m
    out = [0] * m
    aa = [(i, v) for i, v in enumerate(a) if v]
    cc = [(j, v) for j, v in enumerate(c) if v]
    for i, av in aa:
        for j, cv in cc:
            out[(i + j) % m] += av * cv
    return out


def role_data(m: int, k: int) -> tuple[int, list[int], list[int], list[int], list[int]]:
    """Return period, fixed indices, fixed, paired-orbit, and full histograms."""
    p = period(m, k)
    fixed_indices = [n for n in range(p) if (2 * n - (k - 1)) % p == 0]
    assert len(fixed_indices) == (1 if p % 2 else 0)
    fixed_hist = [0] * m
    g = [0] * m
    full = [0] * m
    seen = set()
    for n in range(p):
        full[b(n, k) % m] += 1
        mate = (k - 1 - n) % p
        assert b(n, k) % m == b(mate, k) % m
        if n == mate:
            fixed_hist[b(n, k) % m] += 1
            continue
        if n in seen:
            continue
        assert mate != n
        seen.add(n)
        seen.add(mate)
        g[b(n, k) % m] += 1
    assert len(seen) == p - len(fixed_indices)
    expected = [2 * v for v in g]
    expected = [x + y for x, y in zip(expected, fixed_hist)]
    assert expected == full
    return p, fixed_indices, fixed_hist, g, full


def tower(m: int) -> tuple[list[list[int]], list[int], list[int], list[tuple[int, list[int]]]]:
    assert m >= 2
    roles = [role_data(m, k) for k in RANKS]
    fixed_hists = [row[2] for row in roles]
    gs = [row[3] for row in roles]
    fulls = [row[4] for row in roles]

    beta = [[0] * m for _ in range(len(RANKS) + 1)]
    for mask in range(1 << len(RANKS)):
        poly = [0] * m
        poly[0] = 1
        for i in range(len(RANKS)):
            poly = circular_convolution(poly, gs[i] if (mask >> i) & 1 else fixed_hists[i])
        d = mask.bit_count()
        for residue, count in enumerate(poly):
            beta[d][residue] += count

    reconstructed = [sum((1 << d) * beta[d][t] for d in range(5)) for t in range(m)]
    direct = [0] * m
    direct[0] = 1
    for hist in fulls:
        direct = circular_convolution(direct, hist)
    assert reconstructed == direct
    fixed_values = [next((i for i, count in enumerate(row[2]) if count), -1) for row in roles]
    role_summary = [(row[0], row[1]) for row in roles]
    return beta, direct, fixed_values, role_summary


def v2(n: int) -> int:
    if n == 0:
        return -1
    return (n & -n).bit_length() - 1


def main() -> None:
    print("SUN_2468_REFLECTION_ORBIT_TOWER")
    print("ranks=2,4,6,8 target=896315812331399")
    audit = sha256()
    moduli = list(range(2, 101))
    odd_moduli = list(range(3, 100, 2))
    total_residues = 0
    target_zero_counts = [0] * 5
    target_zero_moduli: list[list[int]] = [[] for _ in range(5)]
    target_min_positive: list[int | None] = [None] * 5
    target_v2_hist: dict[int, int] = {}
    for m in moduli:
        beta, direct, fixed_values, role_summary = tower(m)
        total_residues += m
        for t in range(m):
            row = f"{m},{t}," + ",".join(str(beta[d][t]) for d in range(5)) + f",{direct[t]}\n"
            audit.update(row.encode("ascii"))
            if m % 2:
                assert direct[t] % 2 == (1 if t == sum(fixed_values) % m else 0)
            else:
                assert direct[t] % 16 == 0
        if m % 2:
            assert all(p % 2 and len(fixed) == 1 for p, fixed in role_summary)
        else:
            assert all(p % 2 == 0 and not fixed for p, fixed in role_summary)
            assert all(not beta[d][t] for d in range(4) for t in range(m))
        target_residue = TARGET % m
        if m % 2:
            for d in range(5):
                count = beta[d][target_residue]
                if count == 0:
                    target_zero_counts[d] += 1
                    target_zero_moduli[d].append(m)
                elif target_min_positive[d] is None or count < target_min_positive[d]:
                    target_min_positive[d] = count
            vv = v2(direct[target_residue])
            target_v2_hist[vv] = target_v2_hist.get(vv, 0) + 1
    print(f"audit_moduli=2..100 count={len(moduli)} residues={total_residues}")
    print(f"audit_odd_moduli=3..99 count={len(odd_moduli)}")
    print(f"audit_sha256={audit.hexdigest()}")
    print(f"target_zero_counts_by_dimension={target_zero_counts}")
    print(f"target_zero_moduli_by_dimension={target_zero_moduli}")
    print(f"target_min_positive_by_dimension={target_min_positive}")
    print(f"target_v2_hist={sorted(target_v2_hist.items())}")

    small = sha256()
    for m in (3, 5, 7, 9, 11, 13):
        beta, direct, _, _ = tower(m)
        for t in range(m):
            row = f"{m},{t}," + ",".join(str(beta[d][t]) for d in range(5)) + f",{direct[t]}\n"
            small.update(row.encode("ascii"))
    print(f"independent_range_sha256={small.hexdigest()}")

    for m in (2, 3, 4, 5, 7, 8, 11, 13, 16, 17, 23, 31, 499):
        beta, direct, fixed_values, role_summary = tower(m)
        t = TARGET % m
        bvals = [beta[d][t] for d in range(5)]
        print(
            f"m={m} target_residue={t} periods_fixed={role_summary} "
            f"fingerprint={(sum(fixed_values) % m) if m % 2 else 'none'} beta={bvals} "
            f"R={direct[t]} v2={v2(direct[t])}"
        )
    print("verdict=exact_boolean_orbit_tower_native_orbit_graphs_bipartite")


if __name__ == "__main__":
    main()
