"""Exact tests for the equal-rank upset skeleton and neutral boundary bills."""
from __future__ import annotations
import argparse
import itertools as it
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
COUNTS: dict[str, int] = {}


def need(value: bool, tag: str) -> None:
    COUNTS[tag] = COUNTS.get(tag, 0) + 1
    if not value:
        raise RuntimeError(tag)


def bit_ids(mask: int):
    while mask:
        bit = mask & -mask
        yield bit.bit_length() - 1
        mask ^= bit


def closures(n: int):
    """All orders with the displayed order as a linear extension, once each."""
    edges = tuple(it.combinations(range(n), 2))
    seen = set()
    for code in range(1 << len(edges)):
        succ = [1 << i for i in range(n)]
        for p, (i, j) in enumerate(edges):
            if code >> p & 1:
                succ[i] |= 1 << j
        for i in range(n - 1, -1, -1):
            for j in range(i + 1, n):
                if succ[i] >> j & 1:
                    succ[i] |= succ[j]
        row = tuple(succ)
        if row not in seen:
            seen.add(row)
            yield row


def upset(mask: int, succ: tuple[int, ...]) -> bool:
    return all(not (succ[i] & ~mask) for i in bit_ids(mask))


def route(a: int, b: int, succ: tuple[int, ...]) -> list[int]:
    if a.bit_count() != b.bit_count() or not upset(a, succ) or not upset(b, succ):
        raise ValueError("Route requires equally ranked upsets")
    states = [a]
    while a != b:
        remove = next(i for i in bit_ids(a & ~b)
                      if not any(j != i and succ[j] >> i & 1 for j in bit_ids(a)))
        add = next(j for j in bit_ids(b & ~a)
                   if not ((succ[j] & ~(1 << j)) & (b & ~a)))
        a = (a ^ (1 << remove)) | (1 << add)
        states.append(a)
    return states


def mass(mask: int, weights: tuple[int, ...]) -> int:
    return sum(weights[i] for i in bit_ids(mask))


def verify(nmax: int = 5) -> dict:
    poset_counts = {}
    total_routes = 0
    max_length = 0
    # Labels respecting a linear extension represent every finite poset up to
    # relabelling. We also explicitly reverse labels below to test covariance.
    for n in range(nmax + 1):
        count = 0
        for succ in closures(n):
            count += 1
            upsets = [x for x in range(1 << n) if upset(x, succ)]
            for a in upsets:
                for b in upsets:
                    if a.bit_count() != b.bit_count():
                        continue
                    states = route(a, b, succ)
                    total_routes += 1
                    max_length = max(max_length, len(states) - 1)
                    need(states[0] == a and states[-1] == b, "endpoints")
                    need(len(states) - 1 == (a & ~b).bit_count(), "optimal_exchange_count")
                    changed = 0
                    for x, y in zip(states, states[1:]):
                        need(upset(y, succ), "literal_upset")
                        need(y.bit_count() == a.bit_count(), "rank")
                        need((x ^ y).bit_count() == 2, "two_components")
                        need((y ^ b).bit_count() + 2 == (x ^ b).bit_count(), "monotone_distance")
                        need(not (changed & (x ^ y)), "no_core_coordinate_repeated")
                        changed |= x ^ y
                    need(changed == a ^ b, "exact_core_difference")
            # Covariance is checked by a non-order-preserving labelling.
            def rev(mask):
                return sum(1 << (n - 1 - i) for i in bit_ids(mask))
            revsucc = tuple(rev(succ[n - 1 - i]) for i in range(n))
            for a, b in it.combinations(upsets, 2):
                if a.bit_count() == b.bit_count():
                    rr = route(rev(a), rev(b), revsucc)
                    need(all(upset(x, revsucc) for x in rr), "relabelled_route")
                    need(rr[-1] == rev(b), "relabelled_endpoint")
            # Full interval minimizers on all posets through size four.
            if n <= 4:
                intervals = [(lo, hi) for lo in upsets for hi in upsets if not (lo & ~hi)]
                weights = tuple((i % 3) + 1 for i in range(n))
                for (la, ra), (lb, rb) in it.product(intervals, repeat=2):
                    fa = [x for x in upsets if not (la & ~x) and not (x & ~ra)]
                    fb = [x for x in upsets if not (lb & ~x) and not (x & ~rb)]
                    exact = min(mass(x ^ y, weights) for x in fa for y in fb)
                    wa, wb = la | (lb & ra), lb | (la & rb)
                    formula = mass(la & ~rb, weights) + mass(lb & ~ra, weights)
                    need(formula == exact, "boundary_exact_enumeration")
                    need(wa in fa and wb in fb, "boundary_attainment_membership")
                    need(mass(wa ^ wb, weights) == exact, "boundary_attainment_cost")
        poset_counts[str(n)] = count

    # Independent literal guard model: two visible choices v0,v1 of rank one,
    # one neutral shared component z of weight7, and a parity-paying helper h.
    # h=v0, z=v0; total charge v0+h vanishes in Z/2. Two accepted full states.
    states = []
    for v0, v1, h, z in it.product((0, 1), repeat=4):
        if v0 + v1 == 1 and h == v0 and z == v0 and (v0 + h) % 2 == 0:
            states.append((v0, v1, h, z))
    weights = (2, 3, 5, 7)
    exact = sum(w for x, y, w in zip(*states, weights) if x != y)
    need(len(states) == 2 and exact == 17, "guarded_bridge_17")
    need(2 + 3 + 5 < exact, "ignored_neutral_bill_is_false")
    # Unequal exact margins require a three-component exchange.
    unequal = [mask for mask in range(8) if mass(mask, (1, 1, 2)) == 2]
    need(unequal == [3, 4], "unequal_margin_hostile_universe")
    need((unequal[0] ^ unequal[1]).bit_count() == 3, "unequal_margin_hostile_support")
    # Wrong chain direction: {1,2} is an upset, removing2 is illegal.
    chain = (7, 6, 4)
    need(upset(6, chain) and not upset(2, chain), "orientation_hostile")
    return {"status": "FINITE-EXACT", "max_order": nmax,
            "natural_poset_counts": poset_counts, "all_equal_rank_ordered_pairs": total_routes,
            "largest_exchange_count": max_length, "guarded_bridge_cost": exact,
            "checks": COUNTS, "total_checks": sum(COUNTS.values())}


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--write", action="store_true")
    args = p.parse_args()
    result = verify()
    target = HERE / "receipt.json"
    if args.write:
        target.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    elif json.loads(target.read_text(encoding="utf-8")) != result:
        raise RuntimeError("Frozen receipt mismatch")
    print(json.dumps({"status": "PASS", "checks": result["total_checks"],
                      "posets": result["natural_poset_counts"], "pairs": result["all_equal_rank_ordered_pairs"]}))


if __name__ == "__main__":
    main()
