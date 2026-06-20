#!/usr/bin/env python3
"""S57 local tile-address recurrence verifier.

The full tiling recurrence, the even half-tiling recurrence, and the odd
half-tiling recurrence become local once each tile is assigned two clocks:

  beta = upper endpoint a          (full-staircase birth strip)
  tau  = a+b-1                     (mirror/fixed-line crossing time)

The pair (beta,tau) recovers the tile (a,b), so the full recurrence and the
half recurrences are complementary address projections, not competing count
identities.
"""

from __future__ import annotations

from collections import Counter
from math import comb


Tile = tuple[int, int]
Address = tuple[int, int]


def full_count(n: int) -> int:
    return comb(n - 1, 2) if n >= 1 else 0


def half_count(n: int) -> int:
    return ((n - 1) * (n - 1)) // 4 if n >= 1 else 0


def fixed_count(n: int) -> int:
    return (n - 1) // 2 if n >= 1 else 0


def tiles(n: int) -> list[Tile]:
    return [(a, b) for a in range(3, n + 1) for b in range(1, a - 1)]


def addr(tile: Tile) -> Address:
    a, b = tile
    return (a, a + b - 1)


def tile_from_addr(address: Address) -> Tile:
    beta, tau = address
    return (beta, tau - beta + 1)


def valid_addr(address: Address) -> bool:
    beta, tau = address
    return beta >= 3 and beta <= tau <= 2 * beta - 3


def gap_from_addr(address: Address) -> int:
    beta, tau = address
    return 2 * beta - tau - 2


def lower_from_addr(address: Address) -> int:
    beta, tau = address
    return tau - beta + 1


def reflect_tile(n: int, tile: Tile) -> Tile:
    a, b = tile
    return (n + 1 - b, n + 1 - a)


def reflect_addr(n: int, address: Address) -> Address:
    beta, tau = address
    return (n + beta - tau, 2 * n - tau)


def side(n: int, address: Address) -> str:
    beta, tau = address
    if beta > n:
        return "absent"
    if tau < n:
        return "kept"
    if tau == n:
        return "fixed"
    return "discarded"


def canonical_rep_addr(n: int, address: Address) -> Address:
    s = side(n, address)
    if s in {"kept", "fixed"}:
        return address
    if s == "discarded":
        return reflect_addr(n, address)
    raise ValueError(f"tile with address {address} absent at n={n}")


def full_ie_word(n: int, tile: Tile) -> str:
    """Membership in the three size-(n-1) sub-staircases.

    Pin coordinates are r=a-b-1, c=b, with the full triangle r+c<=n-1.
    The three shifted subtriangles are:
      A: r+c<=n-2, B: r>=2, C: c>=2.
    """
    a, b = tile
    r = a - b - 1
    c = b
    bits = []
    if r + c <= n - 2:
        bits.append("A")
    if r >= 2:
        bits.append("B")
    if c >= 2:
        bits.append("C")
    return "".join(bits)


def full_ie_weight(word: str) -> int:
    k = len(word)
    if k == 0:
        return 0
    # Sum singles - pair overlaps + triple overlap for this local membership set.
    return k - comb(k, 2) + comb(k, 3)


def one_flip_facts(address: Address) -> tuple[int, int, int]:
    """Return (lower endpoint, c3 one-flip, H one-flip) from THM-513."""
    lower = lower_from_addr(address)
    gap = gap_from_addr(address)
    return lower, gap, 1 + 2**gap


def verify_addresses(nmax: int = 16) -> None:
    for n in range(2, nmax + 1):
        ts = tiles(n)
        addresses = [addr(t) for t in ts]
        assert len(ts) == full_count(n)
        assert len(set(addresses)) == len(addresses)
        assert all(valid_addr(a) for a in addresses)
        assert all(tile_from_addr(addr(t)) == t for t in ts)

        kept = [t for t in ts if side(n, addr(t)) in {"kept", "fixed"}]
        fixed = [t for t in ts if side(n, addr(t)) == "fixed"]
        assert len(kept) == half_count(n)
        assert len(fixed) == fixed_count(n)

        for t in ts:
            a = addr(t)
            rt = reflect_tile(n, t)
            ra = reflect_addr(n, a)
            assert addr(rt) == ra
            assert valid_addr(ra)
            assert reflect_addr(n, ra) == a
            rep = canonical_rep_addr(n, a)
            assert side(n, rep) in {"kept", "fixed"}
            if n >= 4:
                assert full_ie_weight(full_ie_word(n, t)) == 1


def crossing_layer_count(tau: int) -> int:
    return fixed_count(tau)


def print_increment_ledger() -> None:
    print("INCREMENT LEDGER")
    print(" n full_new fixed_crossing half_new full_bits half_bits")
    for n in range(2, 14):
        full_new = n - 2
        half_new = fixed_count(n)
        assert full_count(n) - full_count(n - 1) == full_new
        assert half_count(n) - half_count(n - 1) == half_new
        print(
            f"{n:2d} {full_new:8d} {fixed_count(n):14d} "
            f"{half_new:8d} {full_count(n):9d} {half_count(n):9d}"
        )


def print_layer_identities(nmax: int = 14) -> None:
    print()
    print("LAYER IDENTITIES")
    for n in range(2, nmax + 1):
        by_birth = Counter(addr(t)[0] for t in tiles(n))
        by_crossing = Counter(a[1] for a in (addr(t) for t in tiles(n)) if a[1] <= n)
        assert sum(by_birth.values()) == full_count(n)
        assert sum(by_crossing.values()) == half_count(n)
        assert all(by_birth[b] == b - 2 for b in by_birth)
        assert all(by_crossing[t] == crossing_layer_count(t) for t in by_crossing)
    print(f"  verified n=2..{nmax}:")
    print("  full_count(n)=sum_{beta<=n}(beta-2)")
    print("  half_count(n)=sum_{tau<=n}floor((tau-1)/2)")


def print_local_examples() -> None:
    examples = [(7, 3), (9, 1), (9, 6)]
    print()
    print("LOCAL TILE TIMELINES")
    for t in examples:
        a = addr(t)
        lower, c3, hpaths = one_flip_facts(a)
        print(
            f"  tile={t} addr=(beta={a[0]},tau={a[1]}) "
            f"lower={lower} gap/c3={c3} oneflip_H={hpaths}"
        )
        start = a[0]
        stop = min(a[1] + 2, 14)
        for n in range(start, stop + 1):
            s = side(n, a)
            rep = canonical_rep_addr(n, a)
            print(
                f"    n={n:2d} parity={'even' if n % 2 == 0 else 'odd ':>4} "
                f"side={s:9s} rep={tile_from_addr(rep)} rep_addr={rep} "
                f"full_word={full_ie_word(n, t)}"
            )


def print_nine_snapshot() -> None:
    n = 9
    print()
    print("N=9 ADDRESS SNAPSHOT")
    words = Counter(full_ie_word(n, t) for t in tiles(n))
    sides = Counter(side(n, addr(t)) for t in tiles(n))
    crossings = Counter(addr(t)[1] for t in tiles(n) if addr(t)[1] <= n)
    print(f"  full IE word histogram={dict(sorted(words.items()))}")
    print(f"  mirror side histogram={dict(sorted(sides.items()))}")
    print(f"  crossing layers through n={n}: {dict(sorted(crossings.items()))}")
    print("  first half-carrier addresses:")
    rows = []
    for t in tiles(n):
        a = addr(t)
        if side(n, a) in {"kept", "fixed"}:
            rows.append((a[1], a[0], t, gap_from_addr(a), full_ie_word(n, t), side(n, a)))
    for tau, beta, t, gap, word, s in rows[:12]:
        print(f"    tile={t} beta={beta} tau={tau} gap={gap} side={s:5s} full_word={word}")


def tournament_analysis() -> None:
    print()
    print("TOURNAMENT ANALYSIS: address layers, not runners")
    vertices = {
        "birth_strip": {"full_recursion", "new_vertex", "upper_endpoint"},
        "crossing_spine": {"half_recursion", "fixed_line", "mirror_event"},
        "mirror_orbit": {"half_recursion", "complement_pair", "pure_coordinate"},
        "gap_root": {"interval_root", "oneflip_c3", "oneflip_H"},
        "endpoint_score": {"score_defect", "lower_endpoint", "upper_endpoint"},
        "cycle_packet": {"interval_root", "cycle_space", "oneflip_c3"},
        "complement_even_dp": {"complement_pair", "half_recursion", "cycle_space"},
    }
    names = list(vertices)
    adj = {name: set() for name in names}
    scores = Counter()
    for i, a in enumerate(names):
        for b in names[i + 1 :]:
            key_a = (len(vertices[a] & vertices[b]), len(vertices[a]), -i)
            key_b = (len(vertices[a] & vertices[b]), len(vertices[b]), -names.index(b))
            if key_a >= key_b:
                adj[a].add(b)
                scores[a] += 1
            else:
                adj[b].add(a)
                scores[b] += 1
    path = sorted(names, key=lambda x: (scores[x], len(vertices[x]), x), reverse=True)
    cycles3 = 0
    for i, a in enumerate(names):
        for j in range(i + 1, len(names)):
            b = names[j]
            for k in range(j + 1, len(names)):
                c = names[k]
                if b in adj[a] and c in adj[b] and a in adj[c]:
                    cycles3 += 1
                if c in adj[a] and b in adj[c] and a in adj[b]:
                    cycles3 += 1
    print("  pairwise observable=(shared predicates, predicate count, declaration order)")
    print("  switch/gauge=larger observable orients the edge")
    print(f"  scores={dict(scores)}")
    print(f"  directed_3cycles={cycles3}")
    print("  tie Hamiltonian path=" + " > ".join(path))


def main() -> None:
    print("S57 local tile-address recurrence verifier")
    print("address(tile=(a,b))=(beta=a, tau=a+b-1)")
    print("half side uses KPS convention: kept iff a+b<=n+1, equivalently tau<=n")
    verify_addresses()
    print("VERIFIED: address bijection, mirror formula, and canonical reps n=2..16")
    print("VERIFIED: local full IE weights n=4..16, where the shifted-triangle cover is nondegenerate")
    print_increment_ledger()
    print_layer_identities()
    print_local_examples()
    print_nine_snapshot()
    tournament_analysis()
    print()
    print("SYNTHESIS")
    print("  full recurrence clock: beta=a, with beta-layer size beta-2.")
    print("  half recurrence clock: tau=a+b-1, with crossing-layer size floor((tau-1)/2).")
    print("  combined address (beta,tau) recovers the tile and its one-flip root facts:")
    print("    b=tau-beta+1, gap=c3=2*beta-tau-2, H_oneflip=1+2^gap.")
    print("  dynamic payoff: n->n+1 adds n-1 full strip bits but only floor(n/2)")
    print("    new independent half-quotient orbit bits for complement-even invariants.")


if __name__ == "__main__":
    main()
