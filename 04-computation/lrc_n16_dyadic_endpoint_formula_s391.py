#!/usr/bin/env python3
"""
lrc_n16_dyadic_endpoint_formula_s391.py

codex-2026-05-31 S391

Exact local endpoint-count algebra for the n=16 Lonely Runner row.

The previous S389/S390 passes isolated a dyadic endpoint-debt law.  This file
turns the cleanest part into a theorem-shaped finite statement:

* for a pure dyadic owner u=2^k, the number of endpoints protected by a
  residue p mod 16u has an explicit 2-adic formula;
* u=16 is the first owner whose lower protectors can cover all endpoints;
* the u=16 lower cover has an exact nine-residue private-endpoint certificate.

The last section also probes odd payloads u=16w.  Those scans are evidence,
not a theorem: the pure dyadic formula appears to tensor by the odd payload at
the level of maximum layer capacity.
"""

from __future__ import annotations

from dataclasses import dataclass


N = 16
SIGNS = (1, -1)
ODD_RESIDUES = (1, 3, 5, 7, 9, 11, 13, 15)


@dataclass(frozen=True)
class FormulaRow:
    k: int
    owner: int
    residues_checked: int
    nonzero_residues: int
    distinct_counts: tuple[int, ...]


def v2(value: int) -> int:
    if value == 0:
        return 10**9
    out = 0
    while value % 2 == 0:
        out += 1
        value //= 2
    return out


def circular_residue_distance(value: int, modulus: int) -> int:
    residue = value % modulus
    return min(residue, modulus - residue)


def endpoint_labels(owner: int) -> tuple[tuple[int, int], ...]:
    return tuple((m, sign) for m in range(owner) for sign in SIGNS)


def protects_endpoint(owner: int, protector: int, m: int, sign: int) -> bool:
    modulus = N * owner
    return circular_residue_distance(protector * (N * m + sign), modulus) < owner


def protected_endpoints(owner: int, protector: int) -> frozenset[tuple[int, int]]:
    return frozenset(
        (m, sign)
        for m, sign in endpoint_labels(owner)
        if protects_endpoint(owner, protector, m, sign)
    )


def pure_dyadic_formula(k: int, protector: int) -> int:
    """Predicted protected endpoint count for owner u=2^k at n=16."""
    owner = 1 << k
    modulus = N * owner
    residue = protector % modulus
    if residue == 0:
        return 2 * owner

    j = v2(residue)
    if j >= k:
        return 0

    odd_part = (residue >> j) % N
    layer_drop = k - j

    if layer_drop >= 3:
        return 1 << (k - 2)
    if layer_drop == 2:
        return (1 << (k - 1)) if odd_part in (1, 3, 13, 15) else 0
    if layer_drop == 1:
        return (1 << k) if odd_part in (1, 15) else 0
    raise AssertionError("unreachable dyadic layer")


def verify_pure_dyadic_formula(kmax: int = 8) -> tuple[FormulaRow, ...]:
    rows: list[FormulaRow] = []
    for k in range(kmax + 1):
        owner = 1 << k
        modulus = N * owner
        actual_counts: list[int] = []
        for protector in range(modulus):
            actual = len(protected_endpoints(owner, protector))
            predicted = pure_dyadic_formula(k, protector)
            if actual != predicted:
                raise AssertionError(
                    f"k={k}, owner={owner}, p={protector}: "
                    f"actual {actual}, predicted {predicted}"
                )
            actual_counts.append(actual)
        rows.append(
            FormulaRow(
                k=k,
                owner=owner,
                residues_checked=modulus,
                nonzero_residues=sum(1 for count in actual_counts if count),
                distinct_counts=tuple(sorted(set(actual_counts))),
            )
        )
    return tuple(rows)


def formula_description() -> str:
    return "\n".join(
        [
            "For owner u=2^k and residue p mod 16u:",
            "  p = 0 mod 16u                         -> 2u endpoints",
            "  p = 2^j q, j >= k, p not 0 mod 16u    -> 0 endpoints",
            "  L = k-j >= 3                           -> 2^(k-2) endpoints",
            "  L = 2 and q mod 16 in {1,3,13,15}     -> 2^(k-1) endpoints",
            "  L = 1 and q mod 16 in {1,15}          -> 2^k endpoints",
            "  the other L=1,2 odd classes            -> 0 endpoints",
        ]
    )


def lower_union(owner: int) -> frozenset[tuple[int, int]]:
    covered: set[tuple[int, int]] = set()
    for protector in range(1, owner):
        covered.update(protected_endpoints(owner, protector))
    return frozenset(covered)


def exact_min_lower_cover(owner: int) -> tuple[int, tuple[int, ...], int] | None:
    """Exact set-cover solver for the small pure-dyadic lower-owner cases."""
    universe = endpoint_labels(owner)
    index = {endpoint: i for i, endpoint in enumerate(universe)}
    full_mask = (1 << len(universe)) - 1

    sets: list[tuple[int, int, int]] = []
    for protector in range(1, owner):
        mask = 0
        for endpoint in protected_endpoints(owner, protector):
            mask |= 1 << index[endpoint]
        if mask:
            sets.append((protector, mask, mask.bit_count()))

    union_mask = 0
    for _protector, mask, _size in sets:
        union_mask |= mask
    if union_mask != full_mask:
        return None

    candidates_by_endpoint: list[list[int]] = [[] for _ in universe]
    for set_index, (_protector, mask, _size) in enumerate(sets):
        work = mask
        while work:
            bit = work & -work
            candidates_by_endpoint[bit.bit_length() - 1].append(set_index)
            work -= bit

    # The constructive nine-cover gives a strong initial upper bound once it
    # exists.  For the small exact cases this makes the lower-bound search fast.
    initial = constructive_nine_cover(owner)
    if initial is not None:
        initial_mask = 0
        by_protector = {protector: mask for protector, mask, _size in sets}
        for protector in initial:
            initial_mask |= by_protector.get(protector, 0)
        if initial_mask == full_mask:
            best_size = len(initial)
            best_cover = tuple(initial)
        else:
            best_size = len(sets) + 1
            best_cover = ()
    else:
        best_size = len(sets) + 1
        best_cover = ()

    max_set_size = max(size for _protector, _mask, size in sets)
    seen: dict[int, int] = {}
    calls = 0

    def dfs(covered: int, chosen: tuple[int, ...]) -> None:
        nonlocal best_size, best_cover, calls
        calls += 1
        if len(chosen) >= best_size:
            return
        if covered == full_mask:
            best_size = len(chosen)
            best_cover = chosen
            return
        remaining = (full_mask & ~covered).bit_count()
        lower_bound = (remaining + max_set_size - 1) // max_set_size
        if len(chosen) + lower_bound >= best_size:
            return
        if seen.get(covered, 10**9) <= len(chosen):
            return
        seen[covered] = len(chosen)

        uncovered = full_mask & ~covered
        best_options: list[int] | None = None
        work = uncovered
        while work:
            bit = work & -work
            endpoint_index = bit.bit_length() - 1
            work -= bit
            options = [
                set_index
                for set_index in candidates_by_endpoint[endpoint_index]
                if sets[set_index][1] & ~covered
            ]
            if best_options is None or len(options) < len(best_options):
                best_options = options
                if len(options) == 1:
                    break

        assert best_options is not None
        ordered_options = sorted(
            best_options,
            key=lambda set_index: (sets[set_index][1] & ~covered).bit_count(),
            reverse=True,
        )
        for set_index in ordered_options:
            protector, mask, _size = sets[set_index]
            dfs(covered | mask, chosen + (protector,))

    dfs(0, ())
    return best_size, best_cover, calls


def constructive_nine_cover(owner: int) -> tuple[int, ...] | None:
    if owner < 16 or owner & (owner - 1):
        return None
    k = owner.bit_length() - 1
    scale = 1 << max(0, k - 5)
    cover = (owner // 2,) + tuple(scale * residue for residue in ODD_RESIDUES)
    if all(0 < protector < owner for protector in cover):
        return cover
    return None


def private_endpoint_certificate(owner: int) -> tuple[tuple[int, int, tuple[tuple[int, int], ...]], ...]:
    lower = tuple(range(1, owner))
    covers = {protector: protected_endpoints(owner, protector) for protector in lower}
    rows: list[tuple[int, int, tuple[tuple[int, int], ...]]] = []
    for protector in lower:
        others: set[tuple[int, int]] = set()
        for other in lower:
            if other != protector:
                others.update(covers[other])
        private = tuple(sorted(covers[protector] - others))
        if private:
            rows.append((protector, len(covers[protector]), private))
    return tuple(rows)


def odd_payload_scan(odd_payloads: tuple[int, ...] = (1, 3, 5, 7, 9, 15)) -> tuple[str, ...]:
    lines: list[str] = []
    pure_max = {0: 4, 1: 4, 2: 8, 3: 16}
    for odd_payload in odd_payloads:
        owner = 16 * odd_payload
        modulus = N * owner
        pieces: list[str] = []
        for layer in range(4):
            counts = sorted(
                {
                    len(protected_endpoints(owner, protector))
                    for protector in range(1, modulus)
                    if v2(protector) == layer
                }
            )
            pieces.append(
                f"v2={layer}: counts={counts}, max={max(counts)}, "
                f"w*pure={odd_payload * pure_max[layer]}"
            )
        lines.append(f"u={owner} (w={odd_payload}): " + "; ".join(pieces))
    return tuple(lines)


def print_formula_verification() -> None:
    print("N=16 PURE DYADIC ENDPOINT-COUNT FORMULA")
    print("=" * 78)
    print(formula_description())
    print()
    rows = verify_pure_dyadic_formula(kmax=8)
    print("Verified against direct endpoint protection:")
    print("  k  owner  residues  nonzero  distinct protected-counts")
    for row in rows:
        print(
            f"  {row.k:>1} {row.owner:>6} {row.residues_checked:>9} "
            f"{row.nonzero_residues:>8}  {row.distinct_counts}"
        )
    print(
        f"Total residues checked: {sum(row.residues_checked for row in rows)}"
    )
    print()


def print_lower_cover_atlas() -> None:
    print("LOWER-PROTECTOR COVER ATLAS")
    print("=" * 78)
    print("A maximal owner cannot use a super-gate.  Its endpoints must be")
    print("covered by lower protectors p<u, so lower-cover size is the local")
    print("debt bill for that owner.")
    print()
    print("  owner  lower union  exact min lower cover")
    for owner in (2, 4, 8, 16, 32, 64):
        union_size = len(lower_union(owner))
        total = 2 * owner
        exact = exact_min_lower_cover(owner)
        exact_text = "no cover" if exact is None else f"{exact[0]} {exact[1]} (calls={exact[2]})"
        print(f"  {owner:>5}  {union_size:>4}/{total:<4}  {exact_text}")
    print()

    print("Constructive nine-cover self-similarity:")
    print("  owner  nine lower protectors                         covered?")
    for owner in (16, 32, 64, 128, 256, 512):
        cover = constructive_nine_cover(owner)
        assert cover is not None
        covered: set[tuple[int, int]] = set()
        for protector in cover:
            covered.update(protected_endpoints(owner, protector))
        print(f"  {owner:>5}  {cover!s:<47} {len(covered) == 2 * owner}")
    print()


def print_u16_private_certificate() -> None:
    print("u=16 PRIVATE-ENDPOINT CERTIFICATE")
    print("=" * 78)
    print("Among lower protectors p=1,...,15, these and only these have")
    print("endpoints no other lower protector touches.  Hence every lower cover")
    print("must contain all nine; the nine also cover all 32 endpoints.")
    print()
    print("  p   protected  private endpoints (m, sign)")
    necessary: list[int] = []
    for protector, protected_count, private in private_endpoint_certificate(16):
        necessary.append(protector)
        print(f"  {protector:>2}  {protected_count:>9}  {private}")
    covered: set[tuple[int, int]] = set()
    for protector in necessary:
        covered.update(protected_endpoints(16, protector))
    print()
    print(f"Necessary set: {tuple(necessary)}")
    print(f"Covers every u=16 endpoint: {len(covered) == 32}")
    print()


def print_odd_payload_scan() -> None:
    print("ODD-PAYLOAD STRESS SCAN")
    print("=" * 78)
    print("For u=16w with odd w, the maximum capacity in each v2 layer")
    print("matches w times the pure u=16 capacity in these exact scans.")
    print("The lower counts vary once the odd payload is present, so this is")
    print("recorded as proof guidance rather than a theorem.")
    print()
    for line in odd_payload_scan():
        print("  " + line)
    print()


def print_proof_synthesis() -> None:
    print("PROOF SYNTHESIS")
    print("=" * 78)
    print(
        "1. THM-366 already closes the no-16-gate branch: without a speed "
        "divisible by 16, the odd unit endpoints a/16 are boundary witnesses."
    )
    print(
        "2. A maximal pure owner u=16 cannot invoke a super-gate.  The local "
        "endpoint bill is exactly nine lower residues, forced by private "
        "endpoints."
    )
    print(
        "3. Owners 2,4,8 cannot be closed by lower protectors at all.  Owners "
        "16,32,64 have exact lower-cover number 9 in the bounded solver, and "
        "there is a uniform nine-cover construction for all pure dyadic "
        "owners u>=16."
    )
    print(
        "4. This is the recursive wrinkle: the first dyadic gate has private "
        "leaves, but higher dyadic owners can self-similarly close their "
        "local endpoint layer.  A full n=16 proof therefore needs a global "
        "debt-flow inequality using the speed budget and maximality, not only "
        "a pointwise private-leaf argument."
    )
    print()


def main() -> None:
    print_formula_verification()
    print_lower_cover_atlas()
    print_u16_private_certificate()
    print_odd_payload_scan()
    print_proof_synthesis()


if __name__ == "__main__":
    main()
