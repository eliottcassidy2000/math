#!/usr/bin/env python3
"""S662: pi as a sum/product/fraction representation carrier.

This script is deliberately small.  The goal is not to compete with numerical
pi algorithms, but to make three classical infinite representations behave as
three different side-channel carriers:

* sums retain additive packets and cancellation/moment order;
* products retain local factors, zeros, norms, and sieve channels;
* continued fractions retain recursive boundary state and convergent owners.

Tournament Analysis is over representation faces, not over numerical values.
The pairwise observable is: which side channel survives this face better?
"""

from __future__ import annotations

from collections import Counter, deque
from dataclasses import dataclass
from decimal import Decimal, getcontext


getcontext().prec = 90

PI = Decimal(
    "3.14159265358979323846264338327950288419716939937510582097494459230781640628620899"
)


def dec(n: int) -> Decimal:
    return Decimal(n)


def leibniz_pi(terms: int) -> Decimal:
    total = Decimal(0)
    sign = Decimal(1)
    for k in range(terms):
        total += sign / dec(2 * k + 1)
        sign = -sign
    return Decimal(4) * total


def atan_inv(q: int, terms: int) -> Decimal:
    """Return arctan(1/q) from its alternating power series."""

    inv = Decimal(1) / dec(q)
    xpow = inv
    x2 = inv * inv
    total = Decimal(0)
    sign = Decimal(1)
    for k in range(terms):
        total += sign * xpow / dec(2 * k + 1)
        xpow *= x2
        sign = -sign
    return total


def machin_pi(terms: int) -> Decimal:
    # pi/4 = 4 arctan(1/5) - arctan(1/239)
    return Decimal(4) * (Decimal(4) * atan_inv(5, terms) - atan_inv(239, terms))


def basel_pi(terms: int) -> Decimal:
    total = Decimal(0)
    for n in range(1, terms + 1):
        total += Decimal(1) / dec(n * n)
    return (Decimal(6) * total).sqrt()


def wallis_pi(factors: int) -> Decimal:
    product = Decimal(1)
    for n in range(1, factors + 1):
        product *= dec(4 * n * n) / dec(4 * n * n - 1)
    return Decimal(2) * product


def brouncker_pi(depth: int) -> Decimal:
    """Finite truncation of Brouncker's continued fraction for 4/pi."""

    if depth < 1:
        raise ValueError("depth must be positive")
    tail = Decimal(2)
    for k in range(depth, 1, -1):
        odd = 2 * k - 1
        tail = Decimal(2) + dec(odd * odd) / tail
    four_over_pi = Decimal(1) + Decimal(1) / tail
    return Decimal(4) / four_over_pi


def cot_partial_fraction_pi(terms: int) -> Decimal:
    """Use pi*cot(pi*z)=1/z+sum_n(1/(z+n)+1/(z-n)) at z=1/4."""

    z = Decimal(1) / Decimal(4)
    total = Decimal(1) / z
    for n in range(1, terms + 1):
        nd = dec(n)
        total += Decimal(1) / (z + nd) + Decimal(1) / (z - nd)
    return total


def sci(x: Decimal, digits: int = 4) -> str:
    return f"{x:.{digits}E}"


def short(x: Decimal, places: int = 18) -> str:
    return f"{x:.{places}f}"


@dataclass(frozen=True)
class Face:
    label: str
    additive_packet: int
    local_factor: int
    recursive_boundary: int
    payload: str
    repo_transfer: str


FACES = [
    Face(
        "sum",
        additive_packet=3,
        local_factor=1,
        recursive_boundary=2,
        payload="terms, moments, cancellation order",
        repo_transfer="Basel power sums; LRC pair sums; OCF coefficient sums",
    ),
    Face(
        "product",
        additive_packet=2,
        local_factor=3,
        recursive_boundary=1,
        payload="local factors, zeros, norm/sieve channels",
        repo_transfer="sine/Wallis factors; Euler products; p-adic obstruction ledgers",
    ),
    Face(
        "fraction",
        additive_packet=1,
        local_factor=2,
        recursive_boundary=3,
        payload="convergents, continuants, boundary owner state",
        repo_transfer="continued-fraction clocks; carry/owner derivatives; Euclidean descent",
    ),
    Face(
        "raw_decimal",
        additive_packet=0,
        local_factor=0,
        recursive_boundary=0,
        payload="scalar digits only",
        repo_transfer="diagnostic checksum; not a proof carrier",
    ),
]


def face_scores(face: Face) -> tuple[int, int, int]:
    return (face.additive_packet, face.local_factor, face.recursive_boundary)


def face_beats(a: Face, b: Face) -> bool:
    wins = sum(x > y for x, y in zip(face_scores(a), face_scores(b)))
    losses = sum(x < y for x, y in zip(face_scores(a), face_scores(b)))
    if wins != losses:
        return wins > losses
    # Tie Hamiltonian path for the quotient: sum -> product -> fraction -> raw.
    tie_order = {face.label: i for i, face in enumerate(FACES)}
    return tie_order[a.label] < tie_order[b.label]


def tournament_fingerprint(faces: list[Face]) -> dict[str, object]:
    n = len(faces)
    out = [0] * n
    edges: dict[tuple[int, int], int] = {}
    for i in range(n):
        for j in range(i + 1, n):
            winner = i if face_beats(faces[i], faces[j]) else j
            out[winner] += 1
            edges[(i, j)] = winner

    def arc(i: int, j: int) -> bool:
        if i < j:
            return edges[(i, j)] == i
        return edges[(j, i)] == i

    c3 = 0
    cycle_triples: list[tuple[str, str, str]] = []
    for a in range(n):
        for b in range(a + 1, n):
            for c in range(b + 1, n):
                wins = Counter()
                for x, y in ((a, b), (a, c), (b, c)):
                    wins[edges[(min(x, y), max(x, y))]] += 1
                if sorted(wins.values()) == [1, 1, 1]:
                    c3 += 1
                    cycle_triples.append((faces[a].label, faces[b].label, faces[c].label))

    # Hamiltonian path count by DP.
    dp: dict[tuple[int, int], int] = {}
    for last in range(n):
        dp[(1 << last, last)] = 1
    for mask in range(1 << n):
        for last in range(n):
            val = dp.get((mask, last), 0)
            if not val:
                continue
            for nxt in range(n):
                if (mask >> nxt) & 1:
                    continue
                if arc(last, nxt):
                    dp[(mask | (1 << nxt), nxt)] = dp.get((mask | (1 << nxt), nxt), 0) + val
    full = (1 << n) - 1
    hamiltonian_paths = sum(dp.get((full, last), 0) for last in range(n))

    # Strongly connected components.
    adj = [[j for j in range(n) if i != j and arc(i, j)] for i in range(n)]
    radj = [[i for i in range(n) if i != j and arc(i, j)] for j in range(n)]
    seen = [False] * n
    order: list[int] = []
    for start in range(n):
        if seen[start]:
            continue
        stack = [(start, False)]
        while stack:
            v, done = stack.pop()
            if done:
                order.append(v)
                continue
            if seen[v]:
                continue
            seen[v] = True
            stack.append((v, True))
            for w in adj[v]:
                if not seen[w]:
                    stack.append((w, False))
    seen = [False] * n
    sccs: list[list[str]] = []
    for start in reversed(order):
        if seen[start]:
            continue
        comp: list[str] = []
        q = deque([start])
        seen[start] = True
        while q:
            v = q.popleft()
            comp.append(faces[v].label)
            for w in radj[v]:
                if not seen[w]:
                    seen[w] = True
                    q.append(w)
        sccs.append(sorted(comp))

    edge_lines = []
    for i in range(n):
        for j in range(i + 1, n):
            winner = edges[(i, j)]
            loser = j if winner == i else i
            edge_lines.append(f"{faces[winner].label} -> {faces[loser].label}")

    return {
        "score_hist": dict(sorted(Counter(out).items())),
        "outscores": {faces[i].label: out[i] for i in range(n)},
        "directed_3cycles": c3,
        "cycle_triples": cycle_triples,
        "sccs": sccs,
        "hamiltonian_paths": hamiltonian_paths,
        "edges": edge_lines,
    }


def approximation_rows() -> list[tuple[str, str, Decimal]]:
    rows: list[tuple[str, str, Decimal]] = []
    for terms in (10, 100, 1000, 10000):
        rows.append(("sum: Leibniz", f"N={terms}", leibniz_pi(terms)))
    for terms in (3, 5, 10, 20):
        rows.append(("sum: Machin/arctan", f"terms={terms}", machin_pi(terms)))
    for terms in (10, 100, 1000, 10000):
        rows.append(("sum: Basel sqrt", f"N={terms}", basel_pi(terms)))
    for factors in (10, 100, 1000, 10000):
        rows.append(("product: Wallis", f"N={factors}", wallis_pi(factors)))
    for depth in (3, 5, 10, 20, 50, 100):
        rows.append(("fraction: Brouncker", f"depth={depth}", brouncker_pi(depth)))
    for terms in (10, 100, 1000, 10000):
        rows.append(("fraction: cot poles", f"N={terms}", cot_partial_fraction_pi(terms)))
    return rows


def main() -> None:
    print("=" * 78)
    print("S662 pi trinity: infinite sum, infinite product, infinite fraction")
    print("=" * 78)
    print()
    print("Core formulas")
    print("  sum       pi/4 = sum_{k>=0} (-1)^k/(2k+1)")
    print("  sum-fast  pi/4 = 4*atan(1/5) - atan(1/239), with atan by power series")
    print("  sum       pi^2/6 = sum_{n>=1} 1/n^2")
    print("  product   pi/2 = product_{n>=1} (2n)^2/((2n-1)(2n+1))")
    print("  fraction  4/pi = 1 + 1^2/(2 + 3^2/(2 + 5^2/(2 + ...)))")
    print("  fraction  pi*cot(pi*z)=1/z+sum_n(1/(z+n)+1/(z-n)); set z=1/4")
    print()

    print("Approximation ledger")
    print(f"{'face':<24} {'cutoff':<14} {'pi estimate':>24} {'abs error':>14}")
    for face, cutoff, value in approximation_rows():
        err = abs(value - PI)
        print(f"{face:<24} {cutoff:<14} {short(value):>24} {sci(err):>14}")

    print()
    print("Representation side-channel ledger")
    print(f"{'face':<12} {'metric triple':<16} {'payload':<48} transfer")
    for face in FACES:
        metric = face_scores(face)
        print(f"{face.label:<12} {str(metric):<16} {face.payload:<48} {face.repo_transfer}")

    fp = tournament_fingerprint(FACES)
    print()
    print("Tournament Analysis over representation faces")
    print("  vertices=sum, product, fraction, raw_decimal")
    print("  pairwise observable=which side channel is preserved better")
    print("  switch=gauge majority over additive_packet, local_factor, recursive_boundary")
    print("  tie Hamiltonian path=sum -> product -> fraction -> raw_decimal")
    print(f"  outscores={fp['outscores']}")
    print(f"  score_hist={fp['score_hist']}")
    print(f"  directed_3cycles={fp['directed_3cycles']} {fp['cycle_triples']}")
    print(f"  sccs={fp['sccs']}")
    print(f"  hamiltonian_paths={fp['hamiltonian_paths']}")
    print("  edges:")
    for edge in fp["edges"]:
        print(f"    {edge}")

    print()
    print("Assumption challenge")
    print("  Vertex candidates considered: values, formulas, terms, factors, convergents,")
    print("  proof lenses, side-channel types, LRC residues, OCF packets, and owner states.")
    print("  This script uses representation faces because the predicate is not numerical")
    print("  equality; it is which information survives scalar collapse.")

    print()
    print("Hypothesis update")
    print("  A scalar value should be stored as a trinity record when possible:")
    print("    sum packets + product factors + fraction boundary state.")
    print("  HYP-2229 already supplies the product-to-sum bridge through the sine product;")
    print("  S662 adds the missing fraction face as a recursive owner/convergent channel.")
    print("  For LRC14 specifically, read the same trinity as:")
    print("    odd-wall sums + C=27 product/gcd card + paired carry-owner fraction branch.")
    print("  The trinity is globally cyclic, but the carry-seam task is fraction-first.")


if __name__ == "__main__":
    main()
