#!/usr/bin/env python3
"""Grothendieck-Katz p-curvature support-gate toy atlas.

codex-2026-06-12

This script is a structured merge note, not a proof of the p-curvature
conjecture.  It keeps one exact lesson from Grothendieck-Katz in a form that
can be tested against the repo's LRC14 support-gate work:

    the mod-p obstruction is an operator, not a scalar residue.

Part 1 computes rank-one p-curvature on local p-jets R = F_p[z]/(z^p) for
small connections d/dz + a(z).  The p-curvature on this p-jet is (D+M_a)^p,
since D^p=0 on R.  The deliberately lossy scalar comparison is multiplication
by a(z)^p.  The two ranks can disagree in both directions:

  * a=m/(1-z): scalar a^p is full rank, but operator p-curvature is zero;
  * a=z/(1-z): scalar a^p is zero on p-jets, but operator p-curvature is full.

That is exactly the warning for LRC14: raw "q is blocked" is not the proof
object.  The missing operator/carry coordinate is the marked support ledger
(Q27/Pisano class/divisor fiber/Bprime owner).

Part 2 prints the matching LRC14 denominator table for AP and the recent
single-stranger residual rows.  Part 3 records Tournament Analysis over proof
routes.
"""
from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from itertools import combinations, permutations
from math import gcd


N = 14
CORE = [7 * k for k in range(1, 13)]
PRIMES = [3, 5, 7, 11, 13]


def add_poly(a: list[int], b: list[int], p: int, n: int) -> list[int]:
    return [((a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)) % p for i in range(n)]


def mul_poly(a: list[int], b: list[int], p: int, n: int) -> list[int]:
    out = [0] * n
    for i, ai in enumerate(a[:n]):
        if ai % p == 0:
            continue
        for j, bj in enumerate(b[: n - i]):
            if bj % p:
                out[i + j] = (out[i + j] + ai * bj) % p
    return out


def inv_series(q: list[int], p: int, n: int) -> list[int]:
    q = [(x % p) for x in q] + [0] * max(0, n - len(q))
    if q[0] % p == 0:
        raise ValueError("denominator is not invertible at this prime")
    inv0 = pow(q[0], -1, p)
    out = [0] * n
    out[0] = inv0
    for k in range(1, n):
        s = 0
        for i in range(1, k + 1):
            s = (s + q[i] * out[k - i]) % p
        out[k] = (-inv0 * s) % p
    return out


def rational_series(num: list[int], den: list[int], p: int, n: int) -> list[int]:
    return mul_poly([(x % p) for x in num], inv_series(den, p, n), p, n)


def matmul(A: list[list[int]], B: list[list[int]], p: int) -> list[list[int]]:
    n = len(A)
    return [
        [sum(A[i][k] * B[k][j] for k in range(n)) % p for j in range(n)]
        for i in range(n)
    ]


def matpow(A: list[list[int]], e: int, p: int) -> list[list[int]]:
    n = len(A)
    out = [[1 if i == j else 0 for j in range(n)] for i in range(n)]
    base = A
    while e:
        if e & 1:
            out = matmul(out, base, p)
        base = matmul(base, base, p)
        e //= 2
    return out


def rank_mod_p(A: list[list[int]], p: int) -> int:
    A = [row[:] for row in A]
    m, n = len(A), len(A[0])
    r = 0
    for c in range(n):
        pivot = next((i for i in range(r, m) if A[i][c] % p), None)
        if pivot is None:
            continue
        A[r], A[pivot] = A[pivot], A[r]
        inv = pow(A[r][c] % p, -1, p)
        A[r] = [(x * inv) % p for x in A[r]]
        for i in range(m):
            if i != r and A[i][c] % p:
                f = A[i][c] % p
                A[i] = [(A[i][j] - f * A[r][j]) % p for j in range(n)]
        r += 1
        if r == m:
            break
    return r


def multiplication_matrix(a: list[int], p: int, n: int) -> list[list[int]]:
    M = [[0] * n for _ in range(n)]
    for col in range(n):
        for i, ai in enumerate(a[: n - col]):
            M[i + col][col] = (M[i + col][col] + ai) % p
    return M


def derivative_matrix(p: int, n: int) -> list[list[int]]:
    D = [[0] * n for _ in range(n)]
    for col in range(1, n):
        D[col - 1][col] = col % p
    return D


def pth_power_poly(a: list[int], p: int, n: int) -> list[int]:
    out = [1] + [0] * (n - 1)
    for _ in range(p):
        out = mul_poly(out, a, p, n)
    return out


@dataclass(frozen=True)
class ConnectionCase:
    name: str
    numerator: tuple[int, ...]
    denominator: tuple[int, ...]
    reading: str


CONNECTIONS = [
    ConnectionCase("zero", (0,), (1,), "trivial finite monodromy"),
    ConnectionCase("constant_1", (1,), (1,), "exp(z), non-algebraic"),
    ConnectionCase("one_over_1_minus_z", (1,), (1, -1), "(1-z)^-1, rational"),
    ConnectionCase("two_over_1_minus_z", (2,), (1, -1), "(1-z)^-2, rational"),
    ConnectionCase("one_over_1_minus_z_squared", (1,), (1, -2, 1), "exp(1/(1-z)), non-algebraic"),
    ConnectionCase("z_over_1_minus_z", (0, 1), (1, -1), "exp(-z)/(1-z), non-algebraic"),
]


def connection_ranks(case: ConnectionCase, p: int) -> tuple[int, int]:
    n = p
    a = rational_series(list(case.numerator), list(case.denominator), p, n)
    D = derivative_matrix(p, n)
    M = multiplication_matrix(a, p, n)
    L = add_matrix(D, M, p)
    pcurv = matpow(L, p, p)
    naive = multiplication_matrix(pth_power_poly(a, p, n), p, n)
    return rank_mod_p(pcurv, p), rank_mod_p(naive, p)


def add_matrix(A: list[list[int]], B: list[list[int]], p: int) -> list[list[int]]:
    n = len(A)
    return [[(A[i][j] + B[i][j]) % p for j in range(n)] for i in range(n)]


def witness_at_q(S: list[int], q: int, n: int = N) -> int | None:
    band = q // n
    residues = [v % q for v in S]
    for a in range(1, q):
        if gcd(a, q) != 1:
            continue
        if all(min((a * r) % q, q - ((a * r) % q)) > band for r in residues):
            return a
    return None


def first_plain_witness(S: list[int], qmax: int) -> tuple[int, int] | None:
    for q in range(2, qmax + 1):
        a = witness_at_q(S, q)
        if a is not None:
            return q, a
    return None


def q_lattice(max_m: int) -> list[int]:
    return sorted({d * m for d in (1, 2, 7, 14) for m in range(1, max_m + 1)} - {1})


def first_fiber_witness(S: list[int], max_m: int) -> tuple[int, int] | None:
    for q in q_lattice(max_m):
        a = witness_at_q(S, q)
        if a is not None:
            return q, a
    return None


def fmt_witness(w: tuple[int, int] | None) -> str:
    if w is None:
        return "none"
    q, a = w
    return f"{a}/{q}"


@dataclass(frozen=True)
class LrcRow:
    name: str
    speeds: tuple[int, ...]
    reading: str


LRC_ROWS = [
    LrcRow("AP", tuple(range(1, 14)), "wall atom: finite-monodromy-like floor row"),
    LrcRow("Vstar", tuple([1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 24]), "wall atom with apex carry"),
    LrcRow("S611", tuple(sorted(CORE + [611])), "plain shell blocked; Q27 rescued by q=91"),
    LrcRow("S702", tuple(sorted(CORE + [702])), "plain shell blocked; Q27 rescued by q=91"),
    LrcRow("S793", tuple(sorted(CORE + [793])), "plain shell blocked; Q27 rescued by q=40"),
    LrcRow("S1053", tuple(sorted(CORE + [1053])), "plain shell blocked; Q27 rescued by q=40"),
]


@dataclass(frozen=True)
class Route:
    name: str
    scalar_shadow: str
    retained_channel: str
    operator: int
    exactness: int
    lrc_leverage: int
    computable: int
    risk: int

    def score_tuple(self) -> tuple[int, int, int, int, int]:
        return (self.operator, self.exactness, self.lrc_leverage, self.computable, -self.risk)


ROUTES = [
    Route(
        "p_curvature_operator",
        "A(z) mod p or A(z)^p",
        "(D+M_A)^p; derivative/carry terms retained",
        5,
        5,
        4,
        4,
        1,
    ),
    Route(
        "lrc14_Q27_curvature_ledger",
        "denominator q is blocked",
        "Pisano class, 13-clock, divisor fiber, Bprime owner",
        5,
        4,
        5,
        5,
        2,
    ),
    Route(
        "marked_blocker_hypergraph",
        "all unit twists at q blocked",
        "tau_q, canonical cover load, universal blockers",
        4,
        4,
        5,
        5,
        2,
    ),
    Route(
        "partial_frobenius_support_gate",
        "supersingular scalar reductions",
        "diagonal forms on all partial Frobenius twists",
        5,
        5,
        3,
        3,
        2,
    ),
    Route(
        "raw_q_blocking",
        "q blocked/not blocked",
        "none",
        1,
        3,
        2,
        5,
        4,
    ),
    Route(
        "naive_scalar_modp",
        "A(z)^p or scalar residue",
        "none; shown to lie both ways in toy cases",
        1,
        2,
        1,
        5,
        5,
    ),
]


def route_tournament(routes: list[Route]) -> dict[str, object]:
    wins: dict[int, set[int]] = {i: set() for i in range(len(routes))}
    flips_vs_scalar = 0
    for i, j in combinations(range(len(routes)), 2):
        a, b = routes[i], routes[j]
        if a.score_tuple() >= b.score_tuple():
            winner, loser = i, j
        else:
            winner, loser = j, i
        wins[winner].add(loser)
        scalar_pref = i if (a.exactness, a.computable) >= (b.exactness, b.computable) else j
        if scalar_pref != winner:
            flips_vs_scalar += 1

    score_hist = Counter(len(v) for v in wins.values())
    cycles = 0
    for i, j, k in combinations(range(len(routes)), 3):
        if j in wins[i] and k in wins[j] and i in wins[k]:
            cycles += 1
        if k in wins[i] and j in wins[k] and i in wins[j]:
            cycles += 1

    ham = 0
    first_path: tuple[str, ...] | None = None
    for perm in permutations(range(len(routes))):
        if all(perm[t + 1] in wins[perm[t]] for t in range(len(perm) - 1)):
            ham += 1
            if first_path is None:
                first_path = tuple(routes[i].name for i in perm)
    return {
        "score_hist": dict(sorted(score_hist.items())),
        "directed_3cycles": cycles,
        "hamiltonian_paths": ham,
        "first_path": first_path,
        "flips_vs_scalar": flips_vs_scalar,
        "scores": sorted(((len(wins[i]), routes[i].name) for i in range(len(routes))), reverse=True),
    }


def main() -> None:
    print("=" * 78)
    print("Grothendieck-Katz p-curvature x LRC14 support gate")
    print("=" * 78)
    print("HYP-2446 / T790")
    print("External anchors: arXiv:1610.05674, arXiv:1412.7875, arXiv:2501.13175")
    print()

    print("[1] Rank-one p-curvature on local p-jets R=F_p[z]/(z^p)")
    print("    pcurv_rank = rank((D+M_a)^p); naive_rank = rank(M_{a^p})")
    print("    The mismatch is the operator/carry lesson.")
    for case in CONNECTIONS:
        pairs = [connection_ranks(case, p) for p in PRIMES]
        pcurv = [x for x, _ in pairs]
        naive = [y for _, y in pairs]
        mismatch = sum(1 for x, y in pairs if x != y)
        print(f"  {case.name:28s} pcurv={pcurv!s:20s} naive={naive!s:20s} mismatch={mismatch}/5")
        print(f"    reading: {case.reading}")
    print()

    print("[2] LRC14 denominator-curvature table")
    print("    plain<=27 is the scalar shell horizon; Q27 is the fibered/operator ledger.")
    print("    row       plain<=27  plain<=41  Q27       Q41       reading")
    for row in LRC_ROWS:
        S = list(row.speeds)
        plain27 = first_plain_witness(S, 27)
        plain41 = first_plain_witness(S, 41)
        q27 = first_fiber_witness(S, 27)
        q41 = first_fiber_witness(S, 41)
        print(
            f"    {row.name:8s} {fmt_witness(plain27):10s} {fmt_witness(plain41):10s} "
            f"{fmt_witness(q27):9s} {fmt_witness(q41):9s} {row.reading}"
        )
    print()

    print("[3] Transfer dictionary")
    dictionary = [
        ("finite monodromy", "finite wall family / descending quotient atom"),
        ("almost all primes p", "persistent denominator/fiber ladder tests"),
        ("p-curvature vanishes", "blocked ledger is compatible under carry/Frobenius lifts"),
        ("nonzero p-curvature", "a denominator witness or Bprime/owner-private opening appears"),
        ("operator (D+A)^p", "marked support: Pisano class + 13-clock + divisor fiber + owner"),
        ("naive A^p", "raw q-blocked scalar, too lossy"),
    ]
    for left, right in dictionary:
        print(f"  {left:24s} -> {right}")
    print()

    print("[4] Proposed LRC14 p-curvature principle")
    print("  If a primitive LRC14 row blocks the fibered ladder for a long range,")
    print("  its blocker data should define compatible local sections under q -> p*q")
    print("  and q -> d*m fiber maps.  Compatibility is rare: it should force AP,")
    print("  Vstar, 2AP, or a lower resource via Bprime(any runner)/owner-private")
    print("  deletion.  Incompatibility is a nonzero curvature certificate, i.e. a")
    print("  finite denominator witness.")
    print()

    print("[5] Tournament Analysis over proof-route vertices")
    analysis = route_tournament(ROUTES)
    print("  observable = (operator_retained, exactness, LRC leverage, computability, -risk)")
    print("  score histogram:", analysis["score_hist"])
    print("  directed 3-cycles:", analysis["directed_3cycles"])
    print("  Hamiltonian paths:", analysis["hamiltonian_paths"])
    print("  first Hamiltonian path:", analysis["first_path"])
    print("  edge flips versus scalar-only ranking:", analysis["flips_vs_scalar"])
    print("  scores:")
    for score, name in analysis["scores"]:
        print(f"    {score}: {name}")
    print()

    print("[6] Assumption challenge")
    print("  Alternate vertices considered: primes, denominators, unit twists, runners,")
    print("  local p-jets, blocker supports, Pisano classes, Bprime targets, Frobenius")
    print("  twists, monodromy generators, and proof obligations.")
    print("  Chosen vertices: proof routes.  This preserves whether scalar local tests")
    print("  retain an operator/carry side channel, and destroys actual differential")
    print("  geometry plus arbitrary multi-stranger LRC interactions.")


if __name__ == "__main__":
    main()
