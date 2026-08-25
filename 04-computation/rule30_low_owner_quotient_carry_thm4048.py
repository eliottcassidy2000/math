#!/usr/bin/env python3
"""Finite exact low-owner and quotient-carry controls for THM-4048."""

from __future__ import annotations

from collections import defaultdict
import hashlib


HORIZON = 1 << 18
EXPECTED_PREFIX = "75be8407c265798fa046baa3ba0f9336e4cbe27bff0be21aeba3e87a7681bea4"
KS = (1, 2, 3, 4, 5, 6, 8)


def require(ok: bool, label: object) -> None:
    if not ok:
        raise RuntimeError(label)


def phi(x: int) -> int:
    return x ^ ((x << 1) | (x << 2))


def center_and_power_states(stop: int):
    word = bytearray(stop)
    powers: dict[int, int] = {}
    row = 1
    for t in range(stop):
        word[t] = (row >> t) & 1
        if t and t & (t - 1) == 0:
            powers[t] = row
        row = phi(row)
    return word, powers


def valuation(x: int) -> int:
    return (x & -x).bit_length() - 1


def low_owners(q: int, v: int, k: int) -> list[int]:
    modulus = 1 << (v + k)
    mask = modulus - 1
    rows = [0] * q
    row = 1
    for t in range(2 * q):
        if t < q:
            rows[t] = row
        else:
            rows[t - q] = ((row - rows[t - q]) >> v) & ((1 << k) - 1)
        row = phi(row) & mask
    require(all(x & 1 for x in rows), ("owners not odd", q, v, k))
    return rows


def owner_window(q: int, v: int, precision: int, count: int) -> list[int]:
    """Return U_m(t) modulo 2^precision for 0 <= t < count."""
    modulus = 1 << (v + precision)
    mask = modulus - 1
    base = [0] * count
    row = 1
    for t in range(count + q):
        if t < count:
            base[t] = row
        if t >= q and t - q < count:
            base[t - q] = ((row - base[t - q]) >> v) & ((1 << precision) - 1)
        row = phi(row) & mask
    require(all(x & 1 for x in base), ("owner window not odd", q, v, precision))
    return base


def first_scale_closure_hostile(powers: dict[int, int], k: int):
    """Find physical failure of K-bit sibling data to determine the Haar child."""
    low_mask = (1 << k) - 1
    for m in range(1, 12):
        q = 1 << m
        v = valuation(powers[q] - 1)
        v_next = valuation(powers[2 * q] - 1)
        d = v_next - v
        owners = owner_window(q, v, k + d, 2 * q)
        seen_pair: dict[tuple[int, int], tuple[int, int]] = {}
        seen_left: dict[int, tuple[int, int]] = {}
        for r in range(q):
            left = owners[r]
            right = owners[r + q]
            child = ((left + right) >> d) & low_mask
            pair_key = (left & low_mask, right & low_mask)
            left_key = left & low_mask
            old = seen_pair.get(pair_key)
            if old is not None and old[1] != child:
                return ("sibling", k, m, q, v, d, pair_key, old, (r, child))
            seen_pair.setdefault(pair_key, (r, child))
            old = seen_left.get(left_key)
            if old is not None and old[1] != child:
                # Retain the first one-sided hostile as a fallback, but keep
                # searching for the stronger sibling-pair collision.
                one_sided = ("left", k, m, q, v, d, left_key, old, (r, child))
            seen_left.setdefault(left_key, (r, child))
        if "one_sided" in locals():
            fallback = one_sided
    return fallback if "fallback" in locals() else None


def minimal_k2_witness(powers: dict[int, int]):
    q = 4
    v = valuation(powers[q] - 1)
    v_next = valuation(powers[2 * q] - 1)
    d = v_next - v
    owners = owner_window(q, v, 4, 2 * q)
    records = tuple(
        (
            r,
            owners[r],
            owners[r + q],
            (owners[r] + owners[r + q]) >> d,
        )
        for r in (1, 2)
    )
    require((v, d, records) == (4, 2, ((1, 15, 5, 5), (2, 15, 13, 7))), "K=2 witness")
    return q, v, d, records


def target_blind_k2_witness(powers: dict[int, int]):
    """First K=2 collision with every target digit at or above cutoff K."""
    q = 8
    v = valuation(powers[q] - 1)
    v_next = valuation(powers[2 * q] - 1)
    d = v_next - v
    owners = owner_window(q, v, 3, 2 * q)
    records = tuple(
        (
            r,
            owners[r],
            owners[r + q],
            (owners[r] + owners[r + q]) >> d,
        )
        for r in (0, 3)
    )
    require(
        (v, d, q - v, records) == (6, 1, 2, ((0, 7, 7, 7), (3, 7, 3, 5))),
        "target-blind K=2 witness",
    )
    return q, v, d, q - v, records


def shell_base_carry_audit(stop: int = 512):
    """Audit c_(q+r)=bit_(q+r-v)(B+U) on every complete shell."""
    rows = [1]
    for _ in range(stop):
        rows.append(phi(rows[-1]))

    first_carry = None
    shell_counts = []
    m = 0
    while 2 * (1 << m) - 1 <= stop:
        q = 1 << m
        v = valuation(rows[q] - 1)
        carry_count = 0
        for r in range(q):
            n = q + r
            j_base = (rows[r] - 1) >> 1
            owner = (rows[n] - rows[r]) >> v
            actual = (rows[n] >> n) & 1
            if n < v:
                predicted = 0
            else:
                target = n - v
                base = j_base >> (v - 1)
                require(((base >> target) & 1) == 0, ("base target support", m, r))
                predicted = ((base + owner) >> target) & 1
                if target:
                    mask = (1 << target) - 1
                    carry = int((base & mask) + (owner & mask) >= 1 << target)
                else:
                    carry = 0
                direct = (owner >> target) & 1
                require(predicted == direct ^ carry, ("base carry split", m, r))
                if carry:
                    carry_count += 1
                    if first_carry is None:
                        first_carry = (
                            m,
                            q,
                            v,
                            r,
                            n,
                            target,
                            j_base,
                            base,
                            owner,
                            direct,
                            carry,
                            actual,
                        )
            require(predicted == actual, ("repaired shell readout", m, r))
        shell_counts.append((m, q, v, carry_count))
        m += 1

    require(
        first_carry == (2, 4, 4, 2, 6, 2, 12, 1, 399, 1, 1, 0),
        "minimal base-addition carry",
    )
    return stop, tuple(shell_counts), first_carry


def fibre_variation(labels: bytearray, owners: list[int]):
    delta: dict[int, int] = defaultdict(int)
    maximum = 0
    final = 0
    for bit, owner in zip(labels, owners):
        delta[owner] += 1 if bit else -1
        final = sum(abs(x) for x in delta.values())
        maximum = max(maximum, final)
    return maximum, final, len(delta)


def consecutive_pairing(labels: bytearray, owners: list[int]):
    pending: dict[int, tuple[int, int]] = {}
    bad_pairs = 0
    good_pairs = 0
    max_fixed_unmatched = 0
    fixed_unmatched = 0
    for index, (bit, owner) in enumerate(zip(labels, owners)):
        old = pending.pop(owner, None)
        if old is None:
            pending[owner] = (index, bit)
            fixed_unmatched += 1
        else:
            fixed_unmatched -= 1
            if old[1] == bit:
                bad_pairs += 1
                fixed_unmatched += 2
            else:
                good_pairs += 1
        max_fixed_unmatched = max(max_fixed_unmatched, fixed_unmatched)
    return good_pairs, bad_pairs, len(pending), max_fixed_unmatched, fixed_unmatched


def toggle_pairing(labels: bytearray, owners: list[int], j: int):
    q = len(labels)
    bit = 1 << j
    good = bad = owner_fail = 0
    for r in range(q):
        s = r ^ bit
        if r > s:
            continue
        if owners[r] != owners[s]:
            owner_fail += 1
        elif labels[r] != labels[s]:
            good += 1
        else:
            bad += 1
    unmatched = 2 * (bad + owner_fail)
    return j, good, bad, owner_fail, unmatched


def main() -> None:
    center, powers = center_and_power_states(HORIZON)
    center_hash = hashlib.sha256(center).hexdigest()
    require(center_hash == EXPECTED_PREFIX, "frozen center-prefix hash")
    print("RULE30_LOW_OWNER_QUOTIENT_CARRY_THM4048")
    print(f"universe=center_0..{HORIZON-1};ks={KS};center_sha256={center_hash}")
    print(
        "scale_closure_hostiles="
        f"{tuple(first_scale_closure_hostile(powers, k) for k in (1, 2, 3, 4, 6, 8))}"
    )
    print(f"minimal_K2_full_residue_witness={minimal_k2_witness(powers)}")
    print(f"target_blind_K2_full_residue_witness={target_blind_k2_witness(powers)}")
    print(f"shell_base_carry_audit={shell_base_carry_audit()}")
    for m in range(12, HORIZON.bit_length() - 1):
        q = 1 << m
        v = valuation(powers[q] - 1)
        labels = center[q : 2 * q]
        signed = sum(1 if x else -1 for x in labels)
        rows = []
        for k in KS:
            owners = low_owners(q, v, k)
            variation = fibre_variation(labels, owners)
            consecutive = consecutive_pairing(labels, owners)
            toggles = tuple(
                toggle_pairing(labels, owners, j)
                for j in sorted(set((0, 1, max(0, m // 2), m - 2, m - 1)))
            )
            rows.append(
                (
                    k,
                    variation,
                    consecutive[3:],
                    min(toggle[-1] for toggle in toggles),
                )
            )
        print(
            f"shell={(m,q,v,signed)};"
            "row=(K,(max_fibre_variation,final_fibre_variation,fibres),"
            "(max_successive_unmatched,final_successive_unmatched),"
            f"best_tested_toggle_unmatched);rows={rows}"
        )
    print("scope=FINITE-EXACT controls; universal reductions are proved in THM-4048;no asymptotic inference")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
