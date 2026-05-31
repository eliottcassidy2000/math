#!/usr/bin/env python3
"""
lrc_bruhat_tits_tree_s399.py

codex-2026-05-31 S399

Read the n=16 Lonely Runner dyadic endpoint law as a finite shadow of the
Bruhat-Tits tree for PGL_2(Q_2).

This is an exploratory computation, not a proof.  It takes the exact THM-367
formula from S391 and reorganizes it as:

* owner u=2^k              -> a dyadic horosphere / tree depth k;
* protector p=2^j q        -> a step from depth k with drop L=k-j;
* odd q mod 16             -> a boundary residue direction;
* protected endpoint count -> a finite Hecke-kernel capacity.

The main surprise is a mass-conservation pattern: after normalizing by the
2u endpoints in a pure dyadic horosphere, every drop L has total active
odd-residue mass exactly one.
"""

from __future__ import annotations

from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from fractions import Fraction
from importlib.machinery import SourceFileLoader
from math import gcd
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
S391 = SourceFileLoader(
    "lrc_n16_dyadic_endpoint_formula_s391",
    str(ROOT / "04-computation" / "lrc_n16_dyadic_endpoint_formula_s391.py"),
).load_module()


N = 16
ODD_RESIDUES = S391.ODD_RESIDUES


@dataclass(frozen=True)
class DropKernel:
    drop: int
    active_odd_residues: tuple[int, ...]
    capacity_ratio: Fraction
    total_active_mass: Fraction


@dataclass(frozen=True)
class FanLedger:
    owner: int
    cover: tuple[int, ...]
    normalized_cover: tuple[int, ...]
    cover_gcd: int
    incidence_hist: tuple[tuple[int, int], ...]
    all_lower_indeg_hist: tuple[tuple[int, int], ...]
    global_private_endpoint_count: int
    all_lower_core: tuple[int, int]
    cylinder_signatures: tuple[tuple[tuple[int, int, int, int], int], ...]


def pow2_depth(value: int) -> int:
    if value <= 0 or value & (value - 1):
        raise ValueError(f"not a positive power of two: {value}")
    return value.bit_length() - 1


def odd_part_mod16(value: int) -> int:
    return (value >> S391.v2(value)) % N


def kernel_for_drop(k: int, drop: int) -> DropKernel:
    """Return the full residue-class kernel at owner 2^k and depth drop."""

    owner = 1 << k
    active: list[int] = []
    ratios: set[Fraction] = set()
    for odd in ODD_RESIDUES:
        protector = (1 << (k - drop)) * odd
        count = S391.pure_dyadic_formula(k, protector)
        if count:
            active.append(odd)
            ratios.add(Fraction(count, 2 * owner))
    if len(ratios) > 1:
        raise AssertionError(f"drop {drop} has nonuniform active capacities")
    ratio = next(iter(ratios), Fraction(0, 1))
    return DropKernel(
        drop=drop,
        active_odd_residues=tuple(active),
        capacity_ratio=ratio,
        total_active_mass=ratio * len(active),
    )


def full_tree_kernels(k: int) -> tuple[DropKernel, ...]:
    return tuple(kernel_for_drop(k, drop) for drop in range(1, k + 1))


def lower_truncation_profile(owner: int) -> tuple[tuple[int, tuple[int, ...], tuple[int, ...]], ...]:
    """Available lower protectors p<u, grouped by tree drop and odd q mod 16."""

    k = pow2_depth(owner)
    by_drop: dict[int, set[int]] = defaultdict(set)
    for protector in range(1, owner):
        drop = k - S391.v2(protector)
        by_drop[drop].add(odd_part_mod16(protector))

    rows: list[tuple[int, tuple[int, ...], tuple[int, ...]]] = []
    for kernel in full_tree_kernels(k):
        available = tuple(sorted(by_drop.get(kernel.drop, ())))
        active_available = tuple(q for q in available if q in kernel.active_odd_residues)
        rows.append((kernel.drop, available, active_available))
    return tuple(rows)


def endpoint_indegrees(owner: int, protectors: tuple[int, ...]) -> Counter[int]:
    indeg: Counter[tuple[int, int]] = Counter()
    for protector in protectors:
        for endpoint in S391.protected_endpoints(owner, protector):
            indeg[endpoint] += 1
    return Counter(indeg[endpoint] for endpoint in S391.endpoint_labels(owner))


def all_lower_private_endpoint_count(owner: int) -> int:
    lower = tuple(range(1, owner))
    covers = {protector: S391.protected_endpoints(owner, protector) for protector in lower}
    counts: Counter[tuple[int, int]] = Counter()
    for endpoints in covers.values():
        counts.update(endpoints)
    return sum(1 for endpoint in S391.endpoint_labels(owner) if counts[endpoint] == 1)


def peel_all_lower_core(owner: int) -> tuple[int, int]:
    """Leaf-peel the bipartite graph: lower protectors versus endpoints."""

    endpoints = set(S391.endpoint_labels(owner))
    protectors = {
        protector
        for protector in range(1, owner)
        if S391.protected_endpoints(owner, protector)
    }
    p_to_e = {
        protector: set(S391.protected_endpoints(owner, protector))
        for protector in protectors
    }
    e_to_p: dict[tuple[int, int], set[int]] = {endpoint: set() for endpoint in endpoints}
    for protector, covered in p_to_e.items():
        for endpoint in covered:
            e_to_p[endpoint].add(protector)

    queue: deque[tuple[str, object]] = deque()
    for endpoint in endpoints:
        if len(e_to_p[endpoint]) <= 1:
            queue.append(("e", endpoint))
    for protector in protectors:
        if len(p_to_e[protector]) <= 1:
            queue.append(("p", protector))

    live_e = set(endpoints)
    live_p = set(protectors)
    while queue:
        kind, node = queue.popleft()
        if kind == "e":
            endpoint = node  # type: ignore[assignment]
            if endpoint not in live_e:
                continue
            live_e.remove(endpoint)
            for protector in tuple(e_to_p[endpoint]):
                if protector not in live_p:
                    continue
                p_to_e[protector].discard(endpoint)
                if len(p_to_e[protector]) <= 1:
                    queue.append(("p", protector))
        else:
            protector = node  # type: ignore[assignment]
            if protector not in live_p:
                continue
            live_p.remove(protector)
            for endpoint in tuple(p_to_e[protector]):
                if endpoint not in live_e:
                    continue
                e_to_p[endpoint].discard(protector)
                if len(e_to_p[endpoint]) <= 1:
                    queue.append(("e", endpoint))

    return len(live_e), len(live_p)


def cylinder_signature(owner: int, protectors: tuple[int, ...], depth: int = 4) -> tuple[tuple[tuple[int, int, int, int], int], ...]:
    """Summarize fan load on 2-adic boundary cylinders (m mod 2^depth, sign)."""

    modulus = 1 << depth
    indeg: Counter[tuple[int, int]] = Counter()
    for protector in protectors:
        indeg.update(S391.protected_endpoints(owner, protector))

    bucket_loads: dict[tuple[int, int], list[int]] = defaultdict(list)
    for endpoint in S391.endpoint_labels(owner):
        m, sign = endpoint
        bucket_loads[(m % modulus, sign)].append(indeg[endpoint])

    signatures: Counter[tuple[int, int, int, int]] = Counter()
    for loads in bucket_loads.values():
        signatures[(len(loads), min(loads), max(loads), sum(loads))] += 1
    return tuple(sorted(signatures.items()))


def fan_ledger(owner: int) -> FanLedger:
    cover = S391.constructive_nine_cover(owner)
    if cover is None:
        raise ValueError(f"no constructive fan for owner {owner}")
    cover_gcd = 0
    for protector in cover:
        cover_gcd = gcd(cover_gcd, protector)
    normalized = tuple(protector // cover_gcd for protector in cover)
    all_lower = tuple(range(1, owner))
    return FanLedger(
        owner=owner,
        cover=cover,
        normalized_cover=normalized,
        cover_gcd=cover_gcd,
        incidence_hist=tuple(sorted(endpoint_indegrees(owner, cover).items())),
        all_lower_indeg_hist=tuple(sorted(endpoint_indegrees(owner, all_lower).items())),
        global_private_endpoint_count=all_lower_private_endpoint_count(owner),
        all_lower_core=peel_all_lower_core(owner),
        cylinder_signatures=cylinder_signature(owner, cover),
    )


def fmt_residues(values: tuple[int, ...]) -> str:
    return "{" + ",".join(str(value) for value in values) + "}"


def fmt_fraction(frac: Fraction) -> str:
    if frac.denominator == 1:
        return str(frac.numerator)
    return f"{frac.numerator}/{frac.denominator}"


def print_dictionary() -> None:
    print("BRUHAT-TITS DICTIONARY FOR THE n=16 DYADIC ROW")
    print("=" * 78)
    print("  owner u=2^k        -> dyadic horosphere / tree depth k")
    print("  endpoint (16m+-1)/(16u) -> two unit boundary fibers over m")
    print("  protector p=2^j q  -> tree step with drop L=k-j")
    print("  odd q mod 16       -> boundary residue direction")
    print("  protected count    -> finite Hecke-kernel capacity")
    print()


def print_mass_conservation() -> None:
    print("1. THM-367 AS A TREE-KERNEL MASS LAW")
    print("=" * 78)
    print("For every checked owner depth k, active odd residue classes at a fixed")
    print("drop L have equal capacity, and their normalized total mass is 1.")
    print()
    print("  k owner drop active q mod 16          each capacity  active mass")
    for k in range(4, 9):
        owner = 1 << k
        for kernel in full_tree_kernels(k):
            print(
                f"  {k:>1} {owner:>5} {kernel.drop:>4} "
                f"{fmt_residues(kernel.active_odd_residues):<24} "
                f"{fmt_fraction(kernel.capacity_ratio):>13}  "
                f"{fmt_fraction(kernel.total_active_mass):>11}"
            )
        print()
    print(
        "Read: L=1 has two active boundary residues, each covering half the "
        "horosphere; L=2 has four active residues, each covering a quarter; "
        "L>=3 has eight active residues, each covering an eighth."
    )
    print()


def print_lower_truncation() -> None:
    print("2. LOWER-PROTECTOR TRUNCATION")
    print("=" * 78)
    print("The full residue kernel is symmetric in q mod 16.  A maximal owner,")
    print("however, can only use lower protectors p<u.  This takes a one-sided")
    print("truncation of the Bruhat-Tits apartment.")
    print()
    for owner in (16, 32, 64, 128):
        print(f"  owner u={owner}:")
        for drop, available, active_available in lower_truncation_profile(owner):
            print(
                f"    drop {drop:>2}: available={fmt_residues(available):<24} "
                f"active lower={fmt_residues(active_available)}"
            )
        print()


def print_fan_ledgers() -> None:
    print("3. NINE-SPEED FAN AS A FINITE TREE STAR")
    print("=" * 78)
    print("For u>=32 the constructive cover is:")
    print("  u/2 plus (u/32)*{1,3,5,7,9,11,13,15}.")
    print("It is one radial L=1 direction plus the full eight-direction boundary")
    print("sphere at L=5.  Its raw incidence mass is always 3u, so average")
    print("endpoint indegree is 3/2.")
    print()
    print("  owner gcd normalized cover                         fan indeg hist")
    for owner in (16, 32, 64, 128, 256, 512):
        ledger = fan_ledger(owner)
        print(
            f"  {owner:>5} {ledger.cover_gcd:>3} "
            f"{str(ledger.normalized_cover):<40} {ledger.incidence_hist}"
        )
    print()
    print("All-lower local cores and private endpoints:")
    print("  owner all-lower indeg hist                 private endpoints  core(E,P)")
    for owner in (16, 32, 64, 128, 256):
        ledger = fan_ledger(owner)
        print(
            f"  {owner:>5} {str(ledger.all_lower_indeg_hist):<34} "
            f"{ledger.global_private_endpoint_count:>17}  {ledger.all_lower_core}"
        )
    print()
    print("Boundary cylinder signatures for the fan, using cylinders (m mod 16, sign):")
    print("  signature tuple = (bucket size, min indeg, max indeg, total indeg)")
    for owner in (16, 32, 64, 128, 256):
        ledger = fan_ledger(owner)
        print(f"  owner {owner:>3}: {ledger.cylinder_signatures}")
    print()


def print_proof_targets() -> None:
    print("4. PROOF TARGETS SUGGESTED BY THE TREE VIEW")
    print("=" * 78)
    print("A. The local dyadic endpoint theorem is a normalized Markov/Hecke kernel:")
    print("   every tree drop has total active residue mass exactly one.")
    print("B. The known nine-cover is a star current: one radial half-horosphere")
    print("   plus eight deep boundary eighth-horospheres.")
    print("C. The u=16 private-leaf proof is exceptional.  From u=32 onward, the")
    print("   all-lower local incidence graph has no private endpoints and has a")
    print("   nonempty 2-core, so the n=16 proof needs a global flow inequality.")
    print("D. A would-be disproof should look like a harmonic current on a truncated")
    print("   2-adic tree.  The conjectural obstruction is that primitive gcd-breaker")
    print("   speeds create positive divergence, visible as gap or endpoint leaves.")


def main() -> None:
    print_dictionary()
    print_mass_conservation()
    print_lower_truncation()
    print_fan_ledgers()
    print_proof_targets()


if __name__ == "__main__":
    main()
