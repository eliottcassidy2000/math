#!/usr/bin/env python3
"""Beatty-Pell decomposition of the triangular tower crossover word.

codex-2026-06-13

This is a sharpening of the HYP-2453/HYP-2454 observation that the overlap
word between the square-shell additive tower A and the triangular square tower
B is "Sturmian/Beatty-like".

The useful correction is:

* The shell address itself is exactly an inhomogeneous Beatty sequence.
* The address increments are the genuine Sturmian binary clock.
* The visible two-side token word is richer: a six-interval rotation coding
  decorated by two zero-density Pell boundary atoms.

The output is meant to be a proof scaffold and a transfer pattern for LRC14
address ledgers, not a finished theorem about LRC.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from itertools import combinations
from math import isqrt, sqrt


@dataclass(frozen=True)
class Interval:
    lo: int
    hi: int

    @property
    def size(self) -> int:
        return self.hi - self.lo + 1

    def contains(self, other: "Interval") -> bool:
        return self.lo <= other.lo and other.hi <= self.hi

    def __str__(self) -> str:
        return f"[{self.lo},{self.hi}]"


def T(n: int) -> int:
    return n * (n + 1) // 2


def A_shell(n: int) -> Interval:
    return Interval(n * n, n * n + 2 * n)


def A_L(n: int) -> Interval:
    return Interval(n * n, n * n + n)


def A_R(n: int) -> Interval:
    return Interval(n * n + n + 1, n * n + 2 * n)


def B_L(m: int) -> Interval:
    a = T(2 * m)
    return Interval(a, a + m)


def B_R(m: int) -> Interval:
    a = T(2 * m)
    return Interval(a + m + 1, a + 2 * m)


def square_shell_index(x: int) -> int:
    return isqrt(x)


def direct_side_state(interval: Interval) -> tuple[str, int]:
    """Direct floor-sqrt classifier from HYP-2453.

    State codes:
      L = contained in A_n.L
      R = contained in A_n.R
      M = crosses the A_n left/right midpoint
      S = crosses the square-shell boundary after A_n
    """

    n = square_shell_index(interval.lo)
    shell = A_shell(n)
    if not shell.contains(interval):
        return ("S", n)
    if A_L(n).contains(interval):
        return ("L", n)
    if A_R(n).contains(interval):
        return ("R", n)
    return ("M", n)


def beatty_address(m: int) -> tuple[int, int, int]:
    """Return (n,d,r) for B_L(m)'s A-shell address.

    n = floor(sqrt(2m^2+m))
      = floor(sqrt(2)*(m+1/4)).

    d = n-m is the inhomogeneous Beatty address
      floor((sqrt(2)-1)m + sqrt(2)/4).

    r = m^2 - 2md - d^2 is the Pell/carry remainder.
    """

    n = isqrt(2 * m * m + m)
    d = n - m
    r = m * m - 2 * m * d - d * d
    return n, d, r


def verify_beatty_identity(limit: int) -> tuple[bool, int | None]:
    """Verify n=floor(sqrt(2)*(m+1/4)) using integer inequalities.

    Avoid floating point: X=sqrt(2)*(4m+1)/4, so
    floor(X)=n iff 16*n^2 <= 2*(4m+1)^2 < 16*(n+1)^2.
    """

    for m in range(1, limit + 1):
        n, _, _ = beatty_address(m)
        x_sq_num = 2 * (4 * m + 1) * (4 * m + 1)
        if not (16 * n * n <= x_sq_num < 16 * (n + 1) * (n + 1)):
            return False, m
        if n != isqrt(2 * m * m + m):
            return False, m
    return True, None


def formula_states(m: int) -> tuple[str, str, int, int, int, int]:
    """Exact token from the Beatty/Pell normal form.

    Returns (left_state, right_state, n, d, r, epsilon), where epsilon records
    whether B_R starts in A_n (0) or A_{n+1} (1).
    """

    n, d, r = beatty_address(m)

    if r <= d - m:
        left = "L"
    elif r <= d:
        left = "M"
    elif r <= 2 * d:
        left = "R"
    else:
        left = "S"

    epsilon = 1 if r >= 2 * d else 0
    if epsilon:
        right = "L" if r <= 3 * d + 2 else "M"
    else:
        if r <= d - 2 * m:
            right = "L"
        elif r < d - m:
            right = "M"
        elif r <= 2 * d - m:
            right = "R"
        else:
            right = "S"

    return left, right, n, d, r, epsilon


def direct_token(m: int) -> tuple[str, int, int]:
    left, n_left = direct_side_state(B_L(m))
    right, n_right = direct_side_state(B_R(m))
    return left + right, n_left, n_right


def verify_formula(limit: int) -> tuple[bool, tuple[int, str, str] | None]:
    for m in range(1, limit + 1):
        left, right, n, _, _, epsilon = formula_states(m)
        token, n_left, n_right = direct_token(m)
        expected = left + right
        if token != expected or n != n_left or n + epsilon != n_right:
            return False, (m, token, expected)
    return True, None


def address_rows(count: int = 24) -> list[str]:
    rows: list[str] = []
    prev_d = 0
    for m in range(1, count + 1):
        left, right, n, d, r, eps = formula_states(m)
        inc = d - prev_d
        prev_d = d
        rows.append(
            f"m={m:2d}: n={n:2d}, d={d:2d}, delta_d={inc}, r={r:3d}, "
            f"eps={eps}, token={left}{right}, B_L={B_L(m)}, B_R={B_R(m)}"
        )
    return rows


def token_counts(limit: int) -> Counter[str]:
    counts: Counter[str] = Counter()
    for m in range(1, limit + 1):
        left, right, *_ = formula_states(m)
        counts[left + right] += 1
    return counts


def expected_token_lengths() -> dict[str, float]:
    rt2 = sqrt(2.0)
    a = (2.0 - rt2) / 4.0
    b = (2.0 - rt2) / 2.0
    c = 0.5
    d = 1.0 - 1.0 / (2.0 * rt2)
    e = (3.0 * rt2 - 2.0) / (2.0 * rt2)
    return {
        "LM": a,
        "MR": b - a,
        "MS": c - b,
        "RS": d - c,
        "SL": e - d,
        "SM": 1.0 - e,
        "LR": 0.0,
        "RL": 0.0,
    }


def interval_model_lines() -> list[str]:
    rt2 = sqrt(2.0)
    pts = [
        ("0", 0.0),
        ("a=(2-sqrt(2))/4", (2.0 - rt2) / 4.0),
        ("b=(2-sqrt(2))/2", (2.0 - rt2) / 2.0),
        ("1/2", 0.5),
        ("c=1-1/(2sqrt(2))", 1.0 - 1.0 / (2.0 * rt2)),
        ("e=(3sqrt(2)-2)/(2sqrt(2))", (3.0 * rt2 - 2.0) / (2.0 * rt2)),
        ("1", 1.0),
    ]
    tokens = ["LM", "MR", "MS", "RS", "SL", "SM"]
    rows: list[str] = []
    for i, token in enumerate(tokens):
        left_name, left_val = pts[i]
        right_name, right_val = pts[i + 1]
        rows.append(
            f"  eta in [{left_name}, {right_name}) -> {token} "
            f"(length {right_val-left_val:.9f})"
        )
    return rows


def rare_tokens(limit: int = 1_000_000, want: int = 8) -> dict[str, list[tuple[int, int, int, int]]]:
    hits: dict[str, list[tuple[int, int, int, int]]] = defaultdict(list)
    for m in range(1, limit + 1):
        left, right, n, d, r, _ = formula_states(m)
        token = left + right
        if token in {"LR", "RL"} and len(hits[token]) < want:
            hits[token].append((m, n, d, r))
        if all(len(hits[token]) >= want for token in ("LR", "RL")):
            break
    return hits


def recurrence_residuals(seq: list[int]) -> list[int]:
    return [seq[i + 2] - 6 * seq[i + 1] + seq[i] - 2 for i in range(len(seq) - 2)]


def pell_wall_explanation(token: str, hit: tuple[int, int, int, int]) -> str:
    m, n, d, r = hit
    if token == "LR":
        x = 2 * n + 1
        y = 2 * m + 1
        return (
            f"m={m}, n={n}, d={d}, r={r}: r=d-m; "
            f"(2n+1)^2-2(2m+1)^2={x*x-2*y*y}"
        )
    x = 2 * (n + 1)
    y = 2 * m + 1
    return (
        f"m={m}, n={n}, d={d}, r={r}: r=2d; "
        f"2(n+1)^2-(2m+1)^2={x*x//2-y*y if x*x % 2 == 0 else 'n/a'}, "
        f"(2(n+1))^2-2(2m+1)^2={x*x-2*y*y}"
    )


def factor_complexity(word: list[str], k_max: int) -> list[int]:
    return [
        len({tuple(word[i : i + k]) for i in range(0, len(word) - k + 1)})
        for k in range(1, k_max + 1)
    ]


def word_complexity_sample(limit: int = 10_000, k_max: int = 20) -> tuple[list[int], list[int], list[int]]:
    token_word: list[str] = []
    carry_word: list[str] = []
    inc_word: list[str] = []
    prev_d = 0
    for m in range(1, limit + 1):
        left, right, _, d, _, eps = formula_states(m)
        token_word.append(left + right)
        carry_word.append(str(eps))
        inc_word.append(str(d - prev_d))
        prev_d = d
    return (
        factor_complexity(inc_word, k_max),
        factor_complexity(carry_word, k_max),
        factor_complexity(token_word, k_max),
    )


def count_directed_3cycles(vertices: list[str], edges: dict[tuple[str, str], str]) -> int:
    cycles = 0
    for a, b, c in combinations(vertices, 3):
        winners = {edges[(a, b)], edges[(a, c)], edges[(b, c)]}
        scores = Counter(winners)
        if sorted(scores.values()) == [1, 1, 1]:
            cycles += 1
    return cycles


def scc_sizes(vertices: list[str], edges: dict[tuple[str, str], str]) -> list[int]:
    graph: dict[str, set[str]] = {v: set() for v in vertices}
    rev: dict[str, set[str]] = {v: set() for v in vertices}
    for (a, b), winner in edges.items():
        loser = b if winner == a else a
        graph[winner].add(loser)
        rev[loser].add(winner)

    seen: set[str] = set()
    order: list[str] = []

    def dfs(v: str) -> None:
        seen.add(v)
        for w in graph[v]:
            if w not in seen:
                dfs(w)
        order.append(v)

    for v in vertices:
        if v not in seen:
            dfs(v)

    seen.clear()
    sizes: list[int] = []

    def rdfs(v: str, bucket: list[str]) -> None:
        seen.add(v)
        bucket.append(v)
        for w in rev[v]:
            if w not in seen:
                rdfs(w, bucket)

    for v in reversed(order):
        if v not in seen:
            bucket: list[str] = []
            rdfs(v, bucket)
            sizes.append(len(bucket))
    return sorted(sizes, reverse=True)


def count_hamiltonian_paths(vertices: list[str], edges: dict[tuple[str, str], str]) -> int:
    index = {v: i for i, v in enumerate(vertices)}
    n = len(vertices)
    dp = [[0] * n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for mask in range(1 << n):
        for last in range(n):
            if not dp[mask][last]:
                continue
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                a, b = vertices[last], vertices[nxt]
                winner = edges[(a, b) if (a, b) in edges else (b, a)]
                if winner == a:
                    dp[mask | (1 << nxt)][nxt] += dp[mask][last]
    return sum(dp[(1 << n) - 1])


def tournament_analysis() -> list[str]:
    vertices = [
        "beatty_shell_address",
        "sturmian_increment_clock",
        "pell_wall_atoms",
        "six_interval_rotation",
        "exact_carry_normal_form",
        "token_complexity_warning",
        "lrc14_address_ledger",
        "higher_moment_fractional_address",
    ]
    features = {
        "beatty_shell_address": (5, 5, 4, 4, 5, 4),
        "sturmian_increment_clock": (4, 5, 3, 4, 3, 3),
        "pell_wall_atoms": (5, 4, 5, 3, 4, 5),
        "six_interval_rotation": (4, 4, 4, 5, 3, 3),
        "exact_carry_normal_form": (5, 5, 5, 5, 4, 4),
        "token_complexity_warning": (3, 4, 4, 4, 3, 2),
        "lrc14_address_ledger": (2, 5, 5, 3, 5, 5),
        "higher_moment_fractional_address": (3, 4, 5, 3, 4, 5),
    }
    edges: dict[tuple[str, str], str] = {}
    scores: Counter[str] = Counter()
    for i, a in enumerate(vertices):
        for j, b in enumerate(vertices):
            if i >= j:
                continue
            wins_a = sum(1 for x, y in zip(features[a], features[b]) if x > y)
            wins_b = sum(1 for x, y in zip(features[a], features[b]) if y > x)
            winner = a if wins_a >= wins_b else b
            edges[(a, b)] = winner
            scores[winner] += 1

    hist = Counter(scores[v] for v in vertices)
    leader = max(vertices, key=lambda v: (scores[v], -vertices.index(v)))
    return [
        "Tournament vertices: "
        + ", ".join(vertices),
        "Observable: majority comparison over "
        "(exactness, clock clarity, hidden-address value, computability, "
        "LRC transfer, proof potential).",
        "Switch/gauge: orient toward the carrier retaining more hidden-address data; "
        "tie path is the listed vertex order.",
        f"score_hist = {dict(sorted(hist.items()))}",
        f"directed_3cycles = {count_directed_3cycles(vertices, edges)}",
        f"scc_sizes = {scc_sizes(vertices, edges)}",
        f"hamiltonian_paths = {count_hamiltonian_paths(vertices, edges)}",
        f"leader = {leader}",
    ]


def main() -> None:
    print("TRIANGULAR TOWER BEATTY-PELL DECOMPOSITION")
    print("==========================================")
    print("Codex/HYP-2456 addendum to HYP-2453/HYP-2454 and concrete HYP-2455 boundary-lift instance.")
    print()

    print("EXACT BEATTY ADDRESS")
    print("====================")
    print("For B_L(m)=[T_{2m},T_{2m}+m], set n=floor(sqrt(T_{2m})).")
    print("Then exactly:")
    print("  n_m = floor(sqrt(2m^2+m)) = floor(sqrt(2)*(m+1/4))")
    print("  d_m = n_m-m = floor((sqrt(2)-1)m + sqrt(2)/4)")
    print("  r_m = m^2 - 2m*d_m - d_m^2")
    ok, bad = verify_beatty_identity(250_000)
    print(f"Integer-inequality verification to m=250000: {ok}; first_bad={bad}")
    print()
    for row in address_rows(28):
        print(f"  {row}")
    print()

    print("EXACT NORMAL FORM FOR THE TWO-SIDE TOKEN")
    print("========================================")
    print("Let token be B_L-state followed by B_R-state, with")
    print("  L=in A_n.L, R=in A_n.R, M=crosses A midpoint, S=crosses square boundary.")
    print("B_L:")
    print("  L iff r <= d-m")
    print("  M iff d-m < r <= d")
    print("  R iff d < r <= 2d")
    print("  S iff r > 2d")
    print("B_R uses epsilon=1 iff r>=2d, i.e. B_R starts in A_{n+1}.")
    print("  if epsilon=0:")
    print("    L iff r <= d-2m")
    print("    M iff d-2m < r < d-m")
    print("    R iff d-m <= r <= 2d-m")
    print("    S iff 2d-m < r < 2d")
    print("  if epsilon=1:")
    print("    L iff r <= 3d+2")
    print("    M iff r > 3d+2")
    ok, bad = verify_formula(250_000)
    print(f"Formula-vs-floor-sqrt verification to m=250000: {ok}; first_bad={bad}")
    print()

    print("ASYMPTOTIC SIX-INTERVAL ROTATION MODEL")
    print("======================================")
    print("Let eta_m = {sqrt(2)*(m+1/4)}.  Ignoring zero-density equality walls,")
    print("the exact carry inequalities collapse to this interval coding:")
    for row in interval_model_lines():
        print(row)
    print()
    counts = token_counts(500_000)
    expected = expected_token_lengths()
    print("Token frequencies through m=500000:")
    for token, count in sorted(counts.items(), key=lambda item: (-item[1], item[0])):
        freq = count / 500_000
        print(
            f"  {token}: {count:7d}, freq={freq:.9f}, "
            f"expected_length={expected.get(token, 0.0):.9f}, "
            f"diff={freq-expected.get(token, 0.0):+.9f}"
        )
    print()

    print("PELL BOUNDARY ATOMS")
    print("===================")
    hits = rare_tokens()
    for token in ("LR", "RL"):
        print(f"{token} first {len(hits[token])} hits:")
        for hit in hits[token]:
            print("  " + pell_wall_explanation(token, hit))
        seq = [m for m, _, _, _ in hits[token]]
        print(f"  m-sequence: {seq}")
        print(f"  recurrence residuals m[k+2]-6m[k+1]+m[k]-2: {recurrence_residuals(seq)}")
        print()

    print("WORD COMPLEXITY")
    print("===============")
    inc_complexity, carry_complexity, token_complexity = word_complexity_sample()
    print("Factor complexity p(k), k=1..20:")
    print(f"  delta_d_m Sturmian clock: {inc_complexity}")
    print(f"  epsilon carry word: {carry_complexity}")
    print(f"  two-side token word: {token_complexity}")
    print("Conclusion: the Beatty address increment is genuinely Sturmian (p(k)=k+1 in the sample),")
    print("but the visible crossover word is a finite-interval/carry decoration, not a Sturmian word.")
    print()

    print("TRANSFER NOTE")
    print("=============")
    print("This is the smallest exact model of a recurring repo pattern:")
    print("  scalar/product collision -> attach address/carry coordinate -> unroll into a word.")
    print("For LRC14, the analogous move is to stop treating a blocked denominator/fiber as")
    print("a scalar wall and instead record its Beatty-like shell address, owner, carry, and")
    print("rare boundary atoms.  The Pell walls here are not noise; they are the exact cases")
    print("where a coarse containment statement becomes a rigid equality/endpoint event.")
    print()

    print("TOURNAMENT ANALYSIS")
    print("===================")
    for row in tournament_analysis():
        print(row)


if __name__ == "__main__":
    main()
