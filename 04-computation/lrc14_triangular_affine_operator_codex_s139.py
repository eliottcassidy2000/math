"""S139: triangular/perfect-number affine-operator carriers for LRC14.

The prompt proposes the equality

    x*(2*x - 1) = 2*log2(x) + 1 - x,

which is equivalent to x^2 = 1/2 + log2(x).  This script records the
correction that the equality has no positive real solution, then keeps the
useful part: triangular/perfect-number and affine-composition data form a
labelled product/depth carrier, not a replacement for exact LRC14 M/Farey data.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
import math
from typing import Callable, Iterable


def log2(x: float) -> float:
    return math.log(x, 2.0)


def log_quadratic_gap(x: float) -> float:
    """Return x^2 - 1/2 - log2(x)."""

    return x * x - 0.5 - log2(x)


def cubic_barrier(x: Fraction) -> Fraction:
    """The cubic carrier named by the prompt."""

    return ((x - 1) / 2) * (x / 2) * (x - Fraction(1, 2))


def cubic_gap(x: Fraction) -> Fraction:
    return x * x - Fraction(1, 2) - cubic_barrier(x)


def bisection_root(f: Callable[[float], float], lo: float, hi: float) -> float:
    flo = f(lo)
    for _ in range(80):
        mid = (lo + hi) / 2
        fmid = f(mid)
        if flo * fmid <= 0:
            hi = mid
        else:
            lo = mid
            flo = fmid
    return (lo + hi) / 2


def affine_profile(word: str) -> tuple[Fraction, list[int], int]:
    """Left-to-right action of a(x)=x/2 and b(x)=x+1.

    If a b occurs before a later a, that +1 is halved.  The depths list records
    how many future halvings each b increment suffers.
    """

    depths: list[int] = []
    future_halves = 0
    for letter in reversed(word):
        if letter == "a":
            future_halves += 1
        elif letter == "b":
            depths.append(future_halves)
        else:
            raise ValueError(f"unexpected letter {letter!r}")
    depths.reverse()

    beta = sum(Fraction(1, 2**depth) for depth in depths)
    return beta, depths, sum(depths)


def triangular(n: int) -> int:
    return n * (n + 1) // 2


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    d = 3
    while d * d <= n:
        if n % d == 0:
            return False
        d += 2
    return True


@dataclass(frozen=True)
class Carrier:
    key: str
    name: str
    score: tuple[int, int, int, int, int, int]
    role: str


CARRIERS = [
    Carrier(
        "M",
        "exact M/Farey branch",
        (6, 6, 6, 5, 5, 6),
        "keeps p/q, q, and excess e=14p-q",
    ),
    Carrier(
        "C27",
        "C27 two-swap splice",
        (5, 6, 5, 5, 5, 6),
        "unit petal plus 12-branch shell transfer",
    ),
    Carrier(
        "K33",
        "Kpq/K33 owner packet",
        (5, 5, 4, 4, 6, 5),
        "p>=3 incidence and forbidden-H state-lift address",
    ),
    Carrier(
        "aff",
        "affine depth profile",
        (4, 6, 5, 5, 4, 6),
        "keeps a/b composition order and dyadic halving depths",
    ),
    Carrier(
        "prod",
        "product triangular/perfect lane",
        (4, 4, 3, 5, 3, 5),
        "p*(Np-1), perfect-number triangular side channel",
    ),
    Carrier(
        "cubic",
        "cubic barrier carrier",
        (3, 4, 4, 4, 3, 5),
        "three-wall signed factor packet",
    ),
    Carrier(
        "log",
        "log/quadratic curvature guardrail",
        (3, 3, 3, 5, 2, 6),
        "proves quadratic/product growth dominates log depth",
    ),
    Carrier(
        "raw",
        "raw scalar equality numerology",
        (1, 1, 1, 2, 1, 1),
        "false equality if used without labels",
    ),
]


def tournament_edges(carriers: Iterable[Carrier]) -> dict[tuple[str, str], str]:
    items = list(carriers)
    edges: dict[tuple[str, str], str] = {}
    order = {carrier.key: i for i, carrier in enumerate(items)}
    for i, left in enumerate(items):
        for right in items[i + 1 :]:
            if left.score > right.score:
                winner = left.key
            elif right.score > left.score:
                winner = right.key
            else:
                winner = left.key if order[left.key] < order[right.key] else right.key
            edges[(left.key, right.key)] = winner
    return edges


def count_directed_triangles(carriers: list[Carrier], edges: dict[tuple[str, str], str]) -> int:
    def beats(a: str, b: str) -> bool:
        key = (a, b) if (a, b) in edges else (b, a)
        return edges[key] == a

    total = 0
    keys = [c.key for c in carriers]
    for i in range(len(keys)):
        for j in range(i + 1, len(keys)):
            for k in range(j + 1, len(keys)):
                a, b, c = keys[i], keys[j], keys[k]
                if (beats(a, b) and beats(b, c) and beats(c, a)) or (
                    beats(a, c) and beats(c, b) and beats(b, a)
                ):
                    total += 1
    return total


def score_hist(carriers: list[Carrier], edges: dict[tuple[str, str], str]) -> dict[int, int]:
    scores = {c.key: 0 for c in carriers}
    for winner in edges.values():
        scores[winner] += 1
    hist: dict[int, int] = {}
    for score in scores.values():
        hist[score] = hist.get(score, 0) + 1
    return dict(sorted(hist.items()))


def main() -> None:
    print("S139 LRC14 TRIANGULAR / AFFINE-OPERATOR CARRIER")
    print("=" * 78)

    print("[0] Assumption challenge / quotient declaration")
    print("  considered vertices:")
    print("    runners, gaps, affine words in a,b, dyadic depth profiles,")
    print("    triangular/product lanes, perfect-number indices, C27 shells,")
    print("    Kpq owners, cubic factor packets, logarithmic depth, and proof obligations.")
    print("  chosen quotient:")
    print("    proof carriers plus affine-word depth profiles; exact M/Farey data stays attached.")
    print("  preserved predicate:")
    print("    branch label, unit/composition order, pair-owner packet, and state-lift fit.")
    print("  destroyed predicate:")
    print("    exact safe-time geometry and row-specific witnesses.")
    print("  challenged assumption:")
    print("    the scalar equality or a raw triangular number can replace q or M(S).")

    print("\n[1] Log/quadratic equality audit")
    print("  x*(2x-1)=2*log2(x)+1-x is exactly equivalent to")
    print("  x^2 = 1/2 + log2(x).")
    min_x = 1.0 / math.sqrt(2.0 * math.log(2.0))
    min_gap = log_quadratic_gap(min_x)
    print(f"  gap g(x)=x^2-1/2-log2(x) has its minimum at x={min_x:.15f}.")
    print(f"  minimum gap = {min_gap:.15f} > 0, so there is no positive real solution.")
    print("  sample gaps:")
    for x in [0.5, 0.75, 1.0, 2.0, 4.0]:
        print(f"    x={x:4.2f}: g(x)={log_quadratic_gap(x): .12f}")
    print("  readout: this is a curvature guardrail, not an identity.")

    print("\n[2] Cubic factor packet from the prompt")
    print("  C(x)=((x-1)/2)*(x/2)*(x-1/2) has roots 0, 1/2, 1.")
    print("  It is a signed three-wall packet, analogous to depth-3/Bonferroni labels.")
    print("  exact comparison with x^2-1/2:")
    for x in [Fraction(1, 2), Fraction(3, 4), Fraction(1), Fraction(2), Fraction(4)]:
        print(
            "    x={:<4}: C(x)={:<8} x^2-1/2-C(x)={}".format(
                str(x), str(cubic_barrier(x)), str(cubic_gap(x))
            )
        )
    roots = [
        bisection_root(lambda y: float(cubic_gap(Fraction.from_float(y))), 0.5, 0.8),
        bisection_root(lambda y: float(cubic_gap(Fraction.from_float(y))), 5.0, 5.8),
    ]
    print(f"  intersections with x^2-1/2 occur near {roots[0]:.12f} and {roots[1]:.12f}.")
    print("  readout: useful only as a typed three-factor carrier; it is not log2(x).")

    print("\n[3] Triangular and perfect-number product lane")
    print("  x*(2x-1)=T_{2x-1}.")
    print("  If x=2^(r-1) and 2^r-1 is prime, this is the even perfect number")
    print("  2^(r-1)*(2^r-1).")
    for r in range(2, 14):
        x = 2 ** (r - 1)
        mersenne = 2**r - 1
        value = x * (2 * x - 1)
        tag = "perfect" if is_prime(mersenne) else "composite-index"
        print(f"    r={r:2d}: x={x:5d}  2^r-1={mersenne:5d}  T={value:8d}  {tag}")
    print("  LRC14 product lane comparison on p/(14p-1):")
    print("      p   p*(2p-1)   p*(14p-1)   difference from 7*p*(2p-1)")
    for p in range(1, 8):
        tri_lane = p * (2 * p - 1)
        lrc_lane = p * (14 * p - 1)
        print(f"    {p:3d} {tri_lane:10d} {lrc_lane:12d} {lrc_lane - 7 * tri_lane:12d}")
    print("  readout: perfect/triangular arithmetic belongs to the product side channel,")
    print("  the same family as p*q, not to the binding denominator q.")

    print("\n[4] Affine operators a(x)=x/2 and b(x)=x+1")
    print("  left-to-right words produce x -> x/2^h + beta,")
    print("  beta=sum 2^{-future_halves_at_each_b}.")
    print("  noncommutator: ab(x)-ba(x)=1/2.")
    print("  word families:")
    for n in range(1, 7):
        staircase = "ba" * n
        block = "b" * n + "a" * n
        tail = "a" * n + "b" * n
        for label, word in [("stair", staircase), ("block", block), ("tail", tail)]:
            beta, depths, depth_sum = affine_profile(word)
            if label == "stair":
                expect = f"T_{n}={triangular(n)}"
            elif label == "block":
                expect = f"{n}^2={n*n}"
            else:
                expect = "0"
            print(
                f"    n={n} {label:<5} beta={str(beta):>8} "
                f"depths={depths!s:<22} depth_sum={depth_sum:<3} expected={expect}"
            )
    print("  readout: triangular numbers arise from the depth profile of a staircase")
    print("  composition, while squares arise from block order.  Same letter counts,")
    print("  different order labels.  That is exactly the anti-scalarization lesson.")

    print("\n[5] Tournament Analysis on proof carriers")
    print("  pairwise observable:")
    print("    (branch retention, composition-order retention, unit preservation,")
    print("     finite checkability, state-lift fit, scalar-decoy resistance).")
    print("  switch/gauge: lexicographically larger role score wins; listed order breaks ties.")
    edges = tournament_edges(CARRIERS)
    for carrier in CARRIERS:
        print(
            f"    {carrier.name:<36} score={carrier.score} role={carrier.role}"
        )
    print(
        "  fingerprint: score_hist={} c3={}".format(
            score_hist(CARRIERS, edges), count_directed_triangles(CARRIERS, edges)
        )
    )
    print("  Hamiltonian carrier order:")
    print("    " + " > ".join(carrier.name for carrier in CARRIERS))

    print("\n[6] Proof readout")
    print("  The triangular/perfect-number idea adds a product-depth carrier:")
    print("    p*(Np-1) generalizes x*(2x-1), and the affine a/b words explain")
    print("    why order/depth labels matter as much as scalar size.")
    print("  Apply it to LRC14 only after exact M/Farey branch and C27/K33 labels:")
    print("    p=2 branch -> C27 unit-petal/two-swap splice;")
    print("    p>=3 branch -> K33 owner packet plus octahedral/Clebsch side labels.")
    print("  Candidate lemma refinement:")
    print("    a low-gap non-AP/GW residual should carry an affine-depth packet whose")
    print("    unit-visible entries force the C27 petal splice, or whose nonunit")
    print("    depth entries expose the K33/octahedral/Clebsch state-lift packet.")
    print("  This is not a proof of LRC14.  It is a new labelled packet interface")
    print("  and a warning that the false log/quadratic equality must not be used")
    print("  as a scalar theorem.")


if __name__ == "__main__":
    main()
