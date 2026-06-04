#!/usr/bin/env python3
"""
sequence_shadow_lab_s633.py

Small shadow-sequence lab for hard isomorphism counts.

The motivating hard row is the self-converse tournament sequence

    1, 1, 2, 2, 8, 12, 88, ...

Instead of only asking for the next expensive term, this scout records nearby
sequences that are still nontrivial but carry more recursion:

* total tournament classes T(n)=A000568(n);
* fixed layer SC(n)=A002785(n);
* merged and nonfixed layers (T +/- SC)/2;
* q-deformed odd-cycle Burnside shadows A(m,q);
* even/odd SC bisection factors from HYP-2074;
* round/LRC self-converse shadow counts;
* LRC transporter shell orbit counts under <2,-1>;
* a proof-lens Tournament Analysis over the sequence-shadow methods.
"""

from __future__ import annotations

from collections import Counter, deque
from fractions import Fraction
from math import factorial, gcd
import os


ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULT_PATH = os.path.join(
    ROOT, "05-knowledge", "results", "sequence_shadow_lab_s633.out"
)


def partitions_nonincreasing(n: int, max_part: int | None = None):
    if max_part is None:
        max_part = n

    def rec(rem: int, mx: int, cur: list[int]):
        if rem == 0:
            yield tuple(cur)
            return
        for p in range(min(rem, mx), 0, -1):
            cur.append(p)
            yield from rec(rem - p, p, cur)
            cur.pop()

    yield from rec(n, max_part, [])


def odd_partitions(n: int):
    for lam in partitions_nonincreasing(n):
        if all(x % 2 == 1 for x in lam):
            yield lam


def z_factor(lam: tuple[int, ...]) -> int:
    out = 1
    for part, mult in Counter(lam).items():
        out *= (part**mult) * factorial(mult)
    return out


def pair_orbit_exponent(lam: tuple[int, ...]) -> int:
    total = sum(part // 2 for part in lam)
    for i in range(len(lam)):
        for j in range(i + 1, len(lam)):
            total += gcd(lam[i], lam[j])
    return total


def A_q(n: int, q: int) -> int:
    """Odd-cycle Burnside q-shadow.  A_q(n,2)=A000568(n)."""
    total = Fraction(0)
    for lam in odd_partitions(n):
        total += Fraction(q ** pair_orbit_exponent(lam), z_factor(lam))
    assert total.denominator == 1
    return total.numerator


def A000568(n: int) -> int:
    return A_q(n, 2)


def sc_partitions(n: int):
    extra = 1 if n % 2 else 0
    target = n - extra
    for lam in partitions_nonincreasing(target):
        if all(x % 4 == 2 for x in lam):
            yield lam + ((1,) if extra else ())


def SC(n: int) -> int:
    if n <= 1:
        return 1
    total = Fraction(0)
    for lam in sc_partitions(n):
        total += Fraction(2 ** pair_orbit_exponent(lam), z_factor(lam))
    assert total.denominator == 1
    return total.numerator


def G_between(lam: tuple[int, ...]) -> int:
    total = 0
    for i in range(len(lam)):
        for j in range(i + 1, len(lam)):
            total += gcd(lam[i], lam[j])
    return total


def c2_odd(lam: tuple[int, ...]) -> int:
    return sum((x - 1) // 2 for x in lam) + G_between(lam)


def sc_even_bisection(m: int) -> int:
    return A_q(m, 4)


def sc_odd_bisection(m: int) -> int:
    """SC(2m+1)=sum 2^ell 4^c2/z over odd partitions of m."""
    total = Fraction(0)
    for lam in odd_partitions(m):
        total += Fraction((2 ** len(lam)) * (4 ** c2_odd(lam)), z_factor(lam))
    assert total.denominator == 1
    return total.numerator


def dominant_T(n: int) -> Fraction:
    """Identity-permutation term in A000568."""
    return Fraction(2 ** (n * (n - 1) // 2), factorial(n))


def dominant_SC_even(m: int) -> Fraction:
    """All-one odd partition term in SC(2m)=A(m,4)."""
    return Fraction(4 ** (m * (m - 1) // 2), factorial(m))


def dominant_SC_odd(m: int) -> Fraction:
    """All-one odd partition term in SC(2m+1)."""
    return Fraction((2**m) * (4 ** (m * (m - 1) // 2)), factorial(m))


def generated_group(mod: int, gens: tuple[int, ...]) -> list[int]:
    seen = {1 % mod}
    q = deque([1 % mod])
    while q:
        a = q.popleft()
        for gen in gens:
            b = (a * gen) % mod
            if b not in seen:
                seen.add(b)
                q.append(b)
    return sorted(seen)


def fold_residue(x: int, mod: int) -> int:
    x %= mod
    return min(x, mod - x)


def lrc_shell_orbits(n: int):
    mod = 2 * n - 1
    group = generated_group(mod, (2, -1 % mod))
    shells = list(range(1, n))
    seen: set[int] = set()
    orbits: list[list[int]] = []
    for shell in shells:
        if shell in seen:
            continue
        orbit = sorted({fold_residue(g * shell, mod) for g in group})
        seen.update(orbit)
        orbits.append(orbit)
    return mod, group, orbits


def round_self_converse_count(m: int) -> int:
    """Repo HYP-2086/HYP-2097 skinny round shadow for m runner vertices."""
    if m <= 1:
        return 1
    return 2 ** ((m - 1) // 2)


def centered_hex(r: int) -> int:
    return 3 * r * (r + 1) + 1


def hamiltonian_count(edges: dict[tuple[int, int], int], n: int) -> int:
    # edges[(i,j)] = 1 iff i -> j for i<j convention, else 0 means j -> i.
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1 << n):
        for last in range(n):
            ways = dp[mask][last]
            if not ways:
                continue
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                a, b = sorted((last, nxt))
                bit = edges[(a, b)]
                last_beats = bit == 1 if last == a else bit == 0
                if last_beats:
                    dp[mask | (1 << nxt)][nxt] += ways
    return sum(dp[-1])


def tournament_fingerprint(orderings: list[list[str]]):
    vertices = sorted({x for order in orderings for x in order})
    idx = {v: i for i, v in enumerate(vertices)}
    wins = Counter()
    edges: dict[tuple[int, int], int] = {}
    for i, a in enumerate(vertices):
        for j, b in enumerate(vertices):
            if i >= j:
                continue
            a_wins = 0
            b_wins = 0
            for order in orderings:
                pos = {name: k for k, name in enumerate(order)}
                if pos[a] < pos[b]:
                    a_wins += 1
                else:
                    b_wins += 1
            winner = a if a_wins >= b_wins else b
            wins[winner] += 1
            edges[(idx[a], idx[b])] = 1 if winner == a else 0

    directed_3 = 0
    for a in range(len(vertices)):
        for b in range(a + 1, len(vertices)):
            for c in range(b + 1, len(vertices)):
                pairs = [(a, b), (a, c), (b, c)]
                out = Counter()
                for x, y in pairs:
                    bit = edges[(x, y)]
                    winner = x if bit else y
                    out[winner] += 1
                if sorted(out.values()) == [1, 1, 1]:
                    directed_3 += 1

    h = hamiltonian_count(edges, len(vertices))
    ranked = sorted(vertices, key=lambda v: (-wins[v], v))
    return {
        "scores": {v: wins[v] for v in vertices},
        "score_hist": dict(sorted(Counter(wins.values()).items())),
        "directed_3cycles": directed_3,
        "H": h,
        "tie_path": ranked,
    }


def fmt_int(x: int) -> str:
    s = str(x)
    if len(s) <= 18:
        return s
    return s[:9] + "..." + s[-6:]


def main() -> str:
    lines: list[str] = []
    lines.append("S633 Sequence Shadow Lab")
    lines.append("========================")
    lines.append("")
    lines.append("Principle")
    lines.append("---------")
    lines.append(
        "When a count is hard, compute adjacent shadows: fixed layer, "
        "merged layer, nonfixed-pair layer, q-deformation, skinny quotient, "
        "and transporter quotient."
    )
    lines.append("")

    lines.append("1. Tournament / self-converse mirror table")
    lines.append("------------------------------------------")
    lines.append(
        "n  T=A000568     SC=fixed      merged=(T+SC)/2  nonfixed=(T-SC)/2  SC/T"
    )
    for n in range(1, 15):
        t = A000568(n)
        sc = SC(n)
        merged = (t + sc) // 2
        nonfixed = (t - sc) // 2
        frac = Fraction(sc, t)
        lines.append(
            f"{n:2d} {fmt_int(t):>14} {fmt_int(sc):>14} "
            f"{fmt_int(merged):>18} {fmt_int(nonfixed):>19} "
            f"{frac.numerator}/{frac.denominator}"
        )
    sc_seq = [SC(n) for n in range(1, 11)]
    lines.append("")
    lines.append("SC(1..10) = " + ", ".join(str(x) for x in sc_seq))
    lines.append(
        "Reading: SC is the fixed layer of complement/converse; merged and "
        "nonfixed layers are easier companions when SC itself is expensive."
    )
    lines.append("")

    lines.append("2. q-deformed odd-cycle Burnside shadows")
    lines.append("----------------------------------------")
    lines.append(
        "m  A(m,2)=T(m)   A(m,4)=SC(2m)  A(m,6) shadow  SC(2m+1)  odd/even"
    )
    for m in range(1, 13):
        a2 = A_q(m, 2)
        a4 = A_q(m, 4)
        a6 = A_q(m, 6)
        odd = sc_odd_bisection(m)
        ratio = Fraction(odd, a4)
        lines.append(
            f"{m:2d} {fmt_int(a2):>14} {fmt_int(a4):>15} "
            f"{fmt_int(a6):>15} {fmt_int(odd):>15} "
            f"{ratio.numerator}/{ratio.denominator}"
        )
    lines.append(
        "Reading: even SC is exactly A(m,4); odd SC is the same sum with an "
        "extra 2^(#parts) fixed-vertex factor.  A(m,6) is not a theorem target, "
        "but it is a nearby pressure count with the same Burnside skeleton."
    )
    lines.append("")

    lines.append("3. Dominant-term recursion feel")
    lines.append("-------------------------------")
    lines.append(
        "m  SC(2m+2)/SC(2m)  dominant 4^m/(m+1)  "
        "SC(2m+3)/SC(2m+1)  dominant 2*4^m/(m+1)"
    )
    for m in range(1, 10):
        even_ratio = Fraction(SC(2 * m + 2), SC(2 * m))
        odd_ratio = Fraction(SC(2 * m + 3), SC(2 * m + 1))
        even_dom = Fraction(4**m, m + 1)
        odd_dom = Fraction(2 * (4**m), m + 1)
        lines.append(
            f"{m:2d} {even_ratio.numerator}/{even_ratio.denominator:<12} "
            f"{even_dom.numerator}/{even_dom.denominator:<10} "
            f"{odd_ratio.numerator}/{odd_ratio.denominator:<12} "
            f"{odd_dom.numerator}/{odd_dom.denominator}"
        )
    lines.append(
        "Dominant partition heuristic: all singleton cycles give the recursion "
        "4^m/(m+1) on the even bisection and twice that on the odd bisection. "
        "It is crude early but points at the right growth currency."
    )
    lines.append("")

    lines.append("4. Skinny and transporter shadows")
    lines.append("---------------------------------")
    lines.append(
        "n  round-SC(n)  C=2n-1  |<2,-1>|  folded shell orbits  gcd reps  hex?"
    )
    hex_values = {centered_hex(r): r for r in range(0, 8)}
    for n in range(3, 25):
        mod, group, orbits = lrc_shell_orbits(n)
        gcds = [gcd(orb[0], mod) for orb in orbits]
        hex_mark = f"C_hex({hex_values[n]})" if n in hex_values else "-"
        lines.append(
            f"{n:2d} {round_self_converse_count(n):11d} {mod:7d} "
            f"{len(group):8d} {len(orbits):20d} {str(gcds):>14} {hex_mark}"
        )
    lines.append(
        "Reading: the LRC round/self-converse shadow is exponentially thinner "
        "than A000568, while transporter shell-orbits expose modulus strata. "
        "Prime C often collapses to one shell orbit; composite C creates the "
        "gcd side channels seen at C=27."
    )
    lines.append("")

    lines.append("5. How to extend a hard sequence without only extending it")
    lines.append("---------------------------------------------------------")
    lines.append("- Fixed layer: count structures with Anti(X) nonempty.")
    lines.append("- Merged layer: quotient X by J, giving (total+fixed)/2 when J is an involution.")
    lines.append("- Nonfixed-pair layer: (total-fixed)/2, often the generic sea.")
    lines.append("- q-shadow: keep the Burnside skeleton but vary the color/pressure parameter.")
    lines.append("- Bisection: split by parity/fixed-point factor before asking for a recurrence.")
    lines.append("- Transporter quotient: count Aut(X), Anti(X), and what the quotient forgets.")
    lines.append("- Skinny quotient: restrict to round/circular/LRC-realizable carriers.")
    lines.append("")

    orderings = [
        [
            "bisection factor",
            "fixed layer",
            "transporter quotient",
            "merged layer",
            "q-shadow",
            "skinny quotient",
            "raw next term",
        ],
        [
            "transporter quotient",
            "fixed layer",
            "bisection factor",
            "skinny quotient",
            "merged layer",
            "q-shadow",
            "raw next term",
        ],
        [
            "skinny quotient",
            "transporter quotient",
            "bisection factor",
            "q-shadow",
            "fixed layer",
            "merged layer",
            "raw next term",
        ],
    ]
    fp = tournament_fingerprint(orderings)
    lines.append("6. Tournament Analysis over sequence-extension methods")
    lines.append("------------------------------------------------------")
    lines.append(
        "Vertices are methods, not tournament vertices. Pairwise observable: "
        "which method gives a related hard sequence while preserving proof "
        "side channels. Gauges: Burnside recursion, transporter data, "
        "LRC/unit-distance portability."
    )
    lines.append(
        f"score_hist={fp['score_hist']} directed_3cycles={fp['directed_3cycles']} "
        f"H={fp['H']}"
    )
    for name in fp["tie_path"]:
        lines.append(f"  score={fp['scores'][name]}: {name}")
    lines.append("")
    lines.append("Assumption Challenge")
    lines.append("--------------------")
    lines.append(
        "Alternate vertices considered: raw values, partitions, automorphism "
        "groups, anti-cosets, LRC shells, q-colors, unit-distance shells, "
        "and proof obligations. Chosen vertices: extension methods."
    )
    lines.append(
        "Preserved predicate: each hard value is surrounded by fixed, merged, "
        "bisection, q-shadow, skinny, or transporter companions. Destroyed "
        "data: individual isomorphism representatives and exact LRC/unit-distance "
        "embeddings."
    )
    lines.append(
        "Challenged assumption: progress means the next term of the original "
        "sequence. Sometimes the better next term is a neighboring sequence "
        "whose recursion exposes the missing side channel."
    )
    lines.append("")
    lines.append("Interpretation")
    lines.append("--------------")
    lines.append(
        "- A000568 and SC already form a mirror pair: total versus fixed under converse."
    )
    lines.append(
        "- SC's even bisection is A(m,4); the odd bisection is the same base "
        "with a fixed-vertex 2^(#parts) tax."
    )
    lines.append(
        "- The merged and nonfixed layers are useful generic companions when "
        "the fixed layer is sparse."
    )
    lines.append(
        "- LRC round counts and shell transporter orbits are much thinner "
        "shadows that preserve proof-relevant side channels."
    )
    lines.append(
        "- This suggests a general recipe for famous hard spaces: do not only "
        "ask for A(n+1); ask for fixed(A), merged(A), transporter(A), and "
        "skinny(A) as separate recursive probes."
    )
    return "\n".join(lines)


if __name__ == "__main__":
    text = main()
    print(text)
    with open(RESULT_PATH, "w", encoding="utf-8") as f:
        f.write(text)
        f.write("\n")
