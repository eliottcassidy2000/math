#!/usr/bin/env python3
"""S643: Goldbach/Lemoine same-pair bridge.

Even Goldbach shadow:

    E = p + q

Odd Lemoine shadow:

    O = a + 2*b

When the same odd-prime pair is used, the odd shadow is ordered because one
prime is doubled.  If O=a+2*b and E=a+b, then

    b = O - E,       a = 2E - O.

So an even target plus one of its odd companions reconstructs the ordered prime
pair.  Duplicates p=q are the branch locus: E=2p, O=3p.  The lone even prime 2
is an apex/boundary channel because pairs {2,q} do not sum to an even target.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from functools import lru_cache
from itertools import combinations_with_replacement


LIMIT = 220


def primes_upto(n: int) -> list[int]:
    sieve = [True] * (n + 1)
    sieve[0:2] = [False, False]
    for p in range(2, int(n**0.5) + 1):
        if sieve[p]:
            step = p
            start = p * p
            sieve[start : n + 1 : step] = [False] * (((n - start) // step) + 1)
    return [i for i, is_prime in enumerate(sieve) if is_prime]


PRIMES = primes_upto(LIMIT)
PRIME_SET = set(PRIMES)
ODD_PRIMES = [p for p in PRIMES if p % 2 == 1]


def bridge_pairs(max_even: int):
    rows = []
    by_even = defaultdict(list)
    by_odd = defaultdict(list)
    for p, q in combinations_with_replacement(ODD_PRIMES, 2):
        E = p + q
        if E > max_even:
            continue
        companions = sorted({E + p, E + q})
        duplicate = p == q
        row = (p, q, E, companions, duplicate)
        rows.append(row)
        by_even[E].append(row)
        for O in companions:
            by_odd[O].append((E, p, q, duplicate))
    return rows, dict(by_even), dict(by_odd)


def lemoine_reps(max_odd: int):
    reps = defaultdict(list)
    for O in range(7, max_odd + 1, 2):
        for b in PRIMES:
            a = O - 2 * b
            if a < 3:
                continue
            if a in PRIME_SET and a % 2 == 1:
                reps[O].append((a, b))
    return dict(reps)


def valid_same_pair(E: int, O: int) -> tuple[bool, tuple[int, int] | None]:
    b = O - E
    a = 2 * E - O
    ok = E % 2 == 0 and O % 2 == 1 and a in PRIME_SET and b in PRIME_SET and a % 2 == 1 and b % 2 == 1
    return ok, (a, b) if ok else None


@dataclass(frozen=True)
class Route:
    name: str
    exactness: int
    bridge_value: int
    repo_fit: int
    computable: int
    risk: int


ROUTES = [
    Route("linear_pair_vieta", 5, 5, 5, 5, 1),
    Route("duplicate_branch_locus", 5, 4, 5, 5, 1),
    Route("bridge_graph_even_odd", 4, 5, 5, 5, 1),
    Route("apex_prime_two_boundary", 4, 4, 5, 5, 1),
    Route("mod6_residue_wheel", 4, 4, 4, 5, 2),
    Route("side_channel_jackknife", 3, 5, 5, 4, 2),
    Route("lemoine_to_goldbach_transfer", 3, 4, 4, 4, 3),
    Route("raw_representation_counts", 2, 2, 3, 5, 3),
    Route("scalar_parity_slogan", 1, 1, 2, 5, 4),
]

CRITERIA = ["exactness", "bridge_value", "repo_fit", "computable"]


def route_tournament():
    idx = {r.name: i for i, r in enumerate(ROUTES)}
    edges: dict[tuple[int, int], int] = {}
    wins = Counter({r.name: 0 for r in ROUTES})
    for a_i, b_i in combinations_with_replacement(range(len(ROUTES)), 2):
        if a_i == b_i:
            continue
        a = ROUTES[a_i]
        b = ROUTES[b_i]
        av = 0
        bv = 0
        for criterion in CRITERIA:
            if getattr(a, criterion) > getattr(b, criterion):
                av += 1
            elif getattr(b, criterion) > getattr(a, criterion):
                bv += 1
        if a.risk < b.risk:
            av += 1
        elif b.risk < a.risk:
            bv += 1
        winner = a if av >= bv else b
        i, j = sorted((a_i, b_i))
        edges[(i, j)] = 1 if winner is ROUTES[i] else 0
        wins[winner.name] += 1
    return edges, wins


def beats(edges: dict[tuple[int, int], int], a: int, b: int) -> bool:
    if a < b:
        return edges[(a, b)] == 1
    return edges[(b, a)] == 0


def hamiltonian_count(edges: dict[tuple[int, int], int], n: int) -> int:
    @lru_cache(maxsize=None)
    def dp(mask: int, last: int) -> int:
        if mask == (1 << last):
            return 1
        prev_mask = mask ^ (1 << last)
        total = 0
        for prev in range(n):
            if ((prev_mask >> prev) & 1) and beats(edges, prev, last):
                total += dp(prev_mask, prev)
        return total

    full = (1 << n) - 1
    return sum(dp(full, last) for last in range(n))


def directed_3cycles(edges: dict[tuple[int, int], int], n: int) -> int:
    total = 0
    for a in range(n):
        for b in range(a + 1, n):
            for c in range(b + 1, n):
                out = Counter()
                for i, j in [(a, b), (a, c), (b, c)]:
                    winner = i if beats(edges, i, j) else j
                    out[winner] += 1
                if sorted(out.values()) == [1, 1, 1]:
                    total += 1
    return total


def main() -> str:
    rows, by_even, by_odd = bridge_pairs(120)
    lemoine = lemoine_reps(121)
    lines: list[str] = []
    lines.append("S643 Goldbach/Lemoine Same-Pair Bridge")
    lines.append("======================================")
    lines.append("")
    lines.append("Linear carrier")
    lines.append("--------------")
    lines.append("Even shadow: E = p + q, unordered for odd-prime pairs.")
    lines.append("Odd shadow:  O = a + 2*b, ordered because b is doubled.")
    lines.append("Same-pair reconstruction from (E,O): doubled prime b=O-E; single prime a=2E-O.")
    lines.append("Validity window: E < O < 2E, with a,b odd primes.")
    lines.append("For an unordered pair {p,q}, the odd companions are O=E+p and O=E+q.")
    lines.append("Duplicate branch p=q: E=2p, O=3p; the two odd orientations fold together.")
    lines.append("Boundary: q=2 Lemoine reps are apex reps and do not share an even Goldbach pair.")
    lines.append("")
    lines.append("First duplicate/branch pairs")
    lines.append("----------------------------")
    dupes = [(p, 2 * p, 3 * p) for p in ODD_PRIMES if 3 * p <= 160]
    lines.append("p -> (even=2p, odd=3p): " + ", ".join(f"{p}->({E},{O})" for p, E, O in dupes[:18]))
    lines.append("The prime 2 is exceptional: 2->(4,6), but 6 is even, so it exits the odd bridge.")
    lines.append("")
    lines.append("Sample even nodes and odd companions")
    lines.append("------------------------------------")
    for E in range(6, 62, 2):
        if E not in by_even:
            continue
        chunks = []
        for p, q, _E, companions, duplicate in by_even[E]:
            tag = "dup" if duplicate else "two"
            chunks.append(f"({p},{q})->{companions}:{tag}")
        lines.append(f"E={E:<3} " + "; ".join(chunks))
    lines.append("")
    lines.append("Odd nodes with many even parents")
    lines.append("--------------------------------")
    dense_odds = sorted(by_odd.items(), key=lambda item: (len(item[1]), item[0]), reverse=True)
    for O, parents in dense_odds[:12]:
        parent_bits = ", ".join(f"E={E} via ({p},{q})" for E, p, q, _dup in parents[:6])
        lines.append(f"O={O:<3} parents={len(parents):<2} {parent_bits}")
    lines.append("")
    lines.append("Lemoine reps up to 121: bridge vs apex q=2")
    lines.append("-------------------------------------------")
    apex_only = []
    bridge_missing = []
    for O in range(7, 122, 2):
        reps = lemoine.get(O, [])
        apex = [rep for rep in reps if rep[1] == 2]
        bridge = [rep for rep in reps if rep[1] % 2 == 1]
        if apex and not bridge:
            apex_only.append(O)
        if not bridge:
            bridge_missing.append(O)
    lines.append(f"odd targets with only apex q=2 reps: {apex_only}")
    lines.append(f"odd targets with no odd-prime bridge rep: {bridge_missing}")
    for O in range(7, 42, 2):
        reps = lemoine.get(O, [])
        apex = [rep for rep in reps if rep[1] == 2]
        bridge = [rep for rep in reps if rep[1] % 2 == 1]
        lines.append(f"O={O:<3} reps={reps} bridge={bridge} apex={apex}")
    lines.append("")
    lines.append("Pair recovery checks")
    lines.append("--------------------")
    checks = [(10, 15), (10, 17), (16, 21), (22, 33), (34, 51), (50, 75)]
    for E, O in checks:
        ok, pair = valid_same_pair(E, O)
        lines.append(f"(E,O)=({E},{O}) -> {pair if ok else 'invalid'}")
    lines.append("")
    lines.append("Tournament Analysis")
    lines.append("-------------------")
    lines.append("Assumption challenge: I considered vertices as evens, odds, primes, unordered pairs, ordered pairs, gaps, residues, bridge edges, and proof lenses.  The chosen vertices are carrier lenses because they preserve the same-pair predicate and expose what is destroyed by forgetting orientation or the q=2 boundary.")
    edges, wins = route_tournament()
    n = len(ROUTES)
    order = sorted(range(n), key=lambda i: wins[ROUTES[i].name], reverse=True)
    score_hist = Counter(wins.values())
    lines.append(f"vertices={n}")
    lines.append(f"score_hist={dict(sorted(score_hist.items()))}")
    lines.append(f"directed_3cycles={directed_3cycles(edges, n)}")
    lines.append(f"Hamiltonian_paths={hamiltonian_count(edges, n)}")
    lines.append("ranking:")
    for rank, idx in enumerate(order, start=1):
        route = ROUTES[idx]
        lines.append(f"  {rank}. {route.name} (score={wins[route.name]})")
    lines.append("")
    lines.append("New applications")
    lines.append("----------------")
    lines.append("1. Treat (E,O) as a linear Vieta carrier: same-pair Goldbach/Lemoine data reconstructs the ordered prime pair exactly.")
    lines.append("2. Use the duplicate line (2p,3p) as the ramification/branch locus, analogous to discriminant zero in the pi/e carrier.")
    lines.append("3. Track q=2 Lemoine reps as apex boundary terms; they are odd representations without an even same-pair companion.")
    lines.append("4. Build an even-odd companion graph and jackknife duplicate, twin-gap, and q=2 boundary channels separately.")
    return "\n".join(lines)


if __name__ == "__main__":
    print(main())
