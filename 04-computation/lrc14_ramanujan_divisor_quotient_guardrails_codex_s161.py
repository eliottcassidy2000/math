#!/usr/bin/env python3
"""
LRC14 Ramanujan/divisor quotient guardrail audit.

This script stress-tests what a quotient of an LRC14 row is allowed to forget.
The deliberately sharp comparison is:

  scalar divisor data      tau, sigma, omega, unitary divisor counts
  exact-period data        q-cover profile, Ramanujan sums c_q(v)
  packet data              Ramanujan sums c_q(v_i +/- v_j), residue multisets
  proof predicate          first qdiv route and exact strict-safe measure bucket

The audit is finite and exact.  It does not prove LRC14.  Its purpose is to
produce admissibility constraints for future proofs: a quotient can be a proof
carrier only when the proof predicate is constant on its fibers, or when the
quotient explicitly carries a certificate for the labels it has forgotten.

Tournament Analysis:
  vertices: quotient candidates / proof carriers, not runners;
  pairwise observable: fewer mixed proof-route fibers wins; if tied, the more
  compact quotient wins;
  tie Hamiltonian path: lexical order after the metric tuple.
"""

from __future__ import annotations

from collections import Counter, defaultdict, deque
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from math import gcd


AP = tuple(range(1, 14))
UNITS14 = tuple(a for a in range(1, 14) if gcd(a, 14) == 1)
Q_MODES = (2, 3, 4, 6, 7, 14, 27, 41)


@dataclass(frozen=True)
class SeedRow:
    label: str
    speeds: tuple[int, ...]
    source: str
    note: str


@dataclass(frozen=True)
class AuditedRow:
    label: str
    speeds: tuple[int, ...]
    source: str
    note: str
    qdiv: int
    safe_measure: Fraction


_factor_cache: dict[int, dict[int, int]] = {}


def factor(n: int) -> dict[int, int]:
    if n in _factor_cache:
        return _factor_cache[n]
    original = n
    out: dict[int, int] = {}
    p = 2
    while p * p <= n:
        while n % p == 0:
            out[p] = out.get(p, 0) + 1
            n //= p
        p += 1 if p == 2 else 2
    if n > 1:
        out[n] = out.get(n, 0) + 1
    _factor_cache[original] = out
    return out


def divisors(n: int) -> list[int]:
    out = [1]
    for p, e in factor(n).items():
        old = out[:]
        pe = 1
        for _ in range(e):
            pe *= p
            out.extend(d * pe for d in old)
    return sorted(out)


def mobius(n: int) -> int:
    f = factor(n)
    if any(e > 1 for e in f.values()):
        return 0
    return -1 if len(f) % 2 else 1


def tau(n: int) -> int:
    out = 1
    for e in factor(n).values():
        out *= e + 1
    return out


def sigma(n: int, k: int = 1) -> int:
    out = 1
    for p, e in factor(n).items():
        if k == 0:
            out *= e + 1
        else:
            out *= (p ** ((e + 1) * k) - 1) // (p**k - 1)
    return out


def omega(n: int) -> int:
    return len(factor(n))


def bigomega(n: int) -> int:
    return sum(factor(n).values())


def rad(n: int) -> int:
    out = 1
    for p in factor(n):
        out *= p
    return out


def ramanujan_c(q: int, n: int) -> int:
    g = gcd(q, n)
    return sum(mobius(q // d) * d for d in divisors(g))


def qdiv(speeds: tuple[int, ...], cap: int = 240) -> int:
    for d in range(2, cap + 1):
        if all(v % d for v in speeds):
            return d
    return cap + 1


def strict_safe_measure(speeds: tuple[int, ...]) -> Fraction:
    """Haar measure of {t in [0,1]: all ||v t|| >= 1/14}.

    Danger intervals are open, but endpoints have measure zero, so closed
    interval union arithmetic gives the exact measure.
    """
    intervals: list[tuple[Fraction, Fraction]] = []
    for v in speeds:
        for m in range(v + 1):
            lo = Fraction(14 * m - 1, 14 * v)
            hi = Fraction(14 * m + 1, 14 * v)
            if hi <= 0 or lo >= 1:
                continue
            lo = max(lo, Fraction(0))
            hi = min(hi, Fraction(1))
            if lo < hi:
                intervals.append((lo, hi))

    intervals.sort()
    merged: list[list[Fraction]] = []
    for lo, hi in intervals:
        if not merged or lo > merged[-1][1]:
            merged.append([lo, hi])
        elif hi > merged[-1][1]:
            merged[-1][1] = hi
    danger = sum((hi - lo for lo, hi in merged), Fraction(0))
    return Fraction(1) - danger


def replace_one(drop: int, add: int) -> tuple[int, ...]:
    return tuple(sorted((set(AP) - {drop}) | {add}))


def replace_many(drops: tuple[int, ...], adds: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(sorted((set(AP) - set(drops)) | set(adds)))


def named_rows() -> list[SeedRow]:
    return [
        SeedRow("AP", AP, "named", "strict boundary equality atom"),
        SeedRow("GW 12->24", replace_one(12, 24), "named", "Goddyn-Wong boundary atom"),
        SeedRow("near 12->36", replace_one(12, 36), "named", "K33 near-miss"),
        SeedRow("petal 10->20", replace_one(10, 20), "named", "unit-petal row"),
        SeedRow("petal 13->26", replace_one(13, 26), "named", "unit-petal row"),
        SeedRow(
            "two-swap 10,12->20,24",
            replace_many((10, 12), (20, 24)),
            "named",
            "P10 plus GW splice",
        ),
        SeedRow(
            "two-swap 10,12->20,36",
            replace_many((10, 12), (20, 36)),
            "named",
            "P10 plus K33 splice",
        ),
        SeedRow("covering 6->98", replace_one(6, 98), "named", "small covering repair"),
        SeedRow("covering 12->84", replace_one(12, 84), "named", "lcm-tail covering repair"),
        SeedRow("covering 12->168", replace_one(12, 168), "named", "same lcm-tail packet"),
    ]


def one_swap_bank(max_add: int = 220) -> list[SeedRow]:
    rows: list[SeedRow] = []
    for drop in AP:
        for add in range(14, max_add + 1):
            if add in AP:
                continue
            speeds = replace_one(drop, add)
            if len(speeds) == 13 and gcd(*speeds) == 1:
                rows.append(
                    SeedRow(
                        f"swap {drop}->{add}",
                        speeds,
                        "one-swap AP bank",
                        "bounded quotient-stress row",
                    )
                )
    return rows


def audit_rows(seed_rows: list[SeedRow]) -> list[AuditedRow]:
    out: list[AuditedRow] = []
    seen: set[tuple[int, ...]] = set()
    for row in seed_rows:
        if row.speeds in seen:
            continue
        seen.add(row.speeds)
        out.append(
            AuditedRow(
                label=row.label,
                speeds=row.speeds,
                source=row.source,
                note=row.note,
                qdiv=qdiv(row.speeds),
                safe_measure=strict_safe_measure(row.speeds),
            )
        )
    return out


def safe_bucket(mu: Fraction) -> str:
    if mu == 0:
        return "zero"
    if mu <= Fraction(1, 1260):
        return "tiny<=1/1260"
    if mu <= Fraction(1, 100):
        return "small<=1/100"
    return "open"


def proof_route(row: AuditedRow) -> str:
    if row.safe_measure == 0:
        return "boundary-zero"
    q = row.qdiv
    if q > 14:
        return f"qdiv>14:{safe_bucket(row.safe_measure)}"
    return f"qdiv={q}:{safe_bucket(row.safe_measure)}"


def gcd14_counts(speeds: tuple[int, ...]) -> tuple[tuple[int, int], ...]:
    return tuple(sorted(Counter(gcd(v, 14) for v in speeds).items()))


def scalar_divisor_sig(row: AuditedRow) -> tuple:
    s = row.speeds
    return (
        len(s),
        sum(tau(v) for v in s),
        sum(sigma(v, 1) for v in s),
        sum(omega(v) for v in s),
        sum(bigomega(v) for v in s),
        gcd14_counts(s),
    )


def unitary_divisor_sig(row: AuditedRow) -> tuple:
    s = row.speeds
    return scalar_divisor_sig(row) + (
        sum(2 ** omega(v) for v in s),
        sum(rad(v) for v in s),
        tuple(sorted(Counter(rad(v) % 14 for v in s).items())),
    )


def qcover_sig(row: AuditedRow) -> tuple:
    s = row.speeds
    return (
        row.qdiv,
        tuple(any(v % d == 0 for v in s) for d in range(2, 15)),
        tuple(sum(1 for v in s if v % d == 0) for d in (2, 7, 14)),
    )


def ramanujan_speed_sig(row: AuditedRow) -> tuple:
    s = row.speeds
    return tuple(sum(ramanujan_c(q, v) for v in s) for q in Q_MODES)


def ramanujan_pair_sig(row: AuditedRow) -> tuple:
    s = row.speeds
    packet = []
    for q in Q_MODES:
        sum_packet = sum(ramanujan_c(q, a + b) for a, b in combinations(s, 2))
        diff_packet = sum(ramanujan_c(q, abs(a - b)) for a, b in combinations(s, 2))
        packet.append((sum_packet, diff_packet))
    return tuple(packet)


def exact_period_packet_sig(row: AuditedRow) -> tuple:
    residues = tuple(sorted(Counter(v % 14 for v in row.speeds).items()))
    return (
        qcover_sig(row),
        ramanujan_speed_sig(row),
        ramanujan_pair_sig(row),
        residues,
    )


def endpoint_measure_sig(row: AuditedRow) -> tuple:
    return (row.qdiv, safe_bucket(row.safe_measure), row.safe_measure)


def full_row_sig(row: AuditedRow) -> tuple[int, ...]:
    return row.speeds


QUOTIENTS = (
    ("scalar_divisor", scalar_divisor_sig, "tau/sigma/omega/gcd14 sums"),
    ("unitary_divisor", unitary_divisor_sig, "adds unitary divisor and radical data"),
    ("qcover", qcover_sig, "first missing divisor and q=2..14 cover profile"),
    ("ramanujan_speed", ramanujan_speed_sig, "primitive-root trace sums c_q(v)"),
    ("ramanujan_pair", ramanujan_pair_sig, "primitive traces of pair sums/differences"),
    ("exact_period_packet", exact_period_packet_sig, "qcover + Ramanujan + residues"),
    ("endpoint_measure", endpoint_measure_sig, "qdiv plus exact safe-measure bucket"),
    ("full_row", full_row_sig, "no quotient; row itself"),
)


@dataclass(frozen=True)
class QuotientAudit:
    name: str
    note: str
    classes: int
    bad_fibers: int
    bad_pair_collisions: int
    largest_fiber: int
    examples: tuple[tuple[str, ...], ...]


def pair_collision_count(route_counts: Counter[str]) -> int:
    total = sum(route_counts.values())
    same = sum(c * (c - 1) // 2 for c in route_counts.values())
    return total * (total - 1) // 2 - same


def audit_quotient(
    rows: list[AuditedRow], name: str, fn, note: str, example_limit: int = 3
) -> QuotientAudit:
    buckets: dict[tuple, list[AuditedRow]] = defaultdict(list)
    for row in rows:
        buckets[fn(row)].append(row)

    bad: list[tuple[int, Counter[str], list[AuditedRow]]] = []
    bad_pairs = 0
    largest = 0
    for bucket_rows in buckets.values():
        largest = max(largest, len(bucket_rows))
        routes = Counter(proof_route(row) for row in bucket_rows)
        if len(routes) > 1:
            bad.append((len(bucket_rows), routes, bucket_rows))
            bad_pairs += pair_collision_count(routes)

    bad.sort(key=lambda x: (-x[0], -pair_collision_count(x[1]), sorted(x[1])))
    examples = []
    for _, routes, bucket_rows in bad[:example_limit]:
        shown = []
        for row in bucket_rows[:6]:
            shown.append(
                f"{row.label} | route={proof_route(row)} | safe={fmt(row.safe_measure)}"
            )
        shown.append("routes=" + dict_repr(routes))
        examples.append(tuple(shown))

    return QuotientAudit(
        name=name,
        note=note,
        classes=len(buckets),
        bad_fibers=len(bad),
        bad_pair_collisions=bad_pairs,
        largest_fiber=largest,
        examples=tuple(examples),
    )


def metric(audit: QuotientAudit) -> tuple[int, int, int, str]:
    return (
        audit.bad_pair_collisions,
        audit.bad_fibers,
        audit.classes,
        audit.name,
    )


def tournament_fingerprint(audits: list[QuotientAudit]) -> dict[str, object]:
    n = len(audits)
    adj = [[False] * n for _ in range(n)]
    for i, ai in enumerate(audits):
        for j, aj in enumerate(audits):
            if i == j:
                continue
            adj[i][j] = metric(ai) < metric(aj)

    scores = [sum(1 for j in range(n) if adj[i][j]) for i in range(n)]
    c3 = 0
    for i, j, k in combinations(range(n), 3):
        if (adj[i][j] and adj[j][k] and adj[k][i]) or (
            adj[i][k] and adj[k][j] and adj[j][i]
        ):
            c3 += 1

    # Reachability-based SCCs is fine at this tiny n.
    reach = [[adj[i][j] or i == j for j in range(n)] for i in range(n)]
    for k in range(n):
        for i in range(n):
            if reach[i][k]:
                for j in range(n):
                    reach[i][j] = reach[i][j] or reach[k][j]
    seen: set[int] = set()
    scc_sizes = []
    for i in range(n):
        if i in seen:
            continue
        comp = [j for j in range(n) if reach[i][j] and reach[j][i]]
        seen.update(comp)
        scc_sizes.append(len(comp))

    # Hamiltonian paths by DP.
    dp = [[0] * n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for mask in range(1 << n):
        for last in range(n):
            count = dp[mask][last]
            if not count:
                continue
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                if adj[last][nxt]:
                    dp[mask | (1 << nxt)][nxt] += count
    hp = sum(dp[(1 << n) - 1])

    order = sorted(range(n), key=lambda i: metric(audits[i]))
    return {
        "score_hist": dict(Counter(scores)),
        "c3": c3,
        "scc_sizes": tuple(sorted(scc_sizes, reverse=True)),
        "hamiltonian_paths": hp,
        "tie_hamiltonian_path": " > ".join(audits[i].name for i in order),
    }


def fmt(fr: Fraction) -> str:
    if fr.denominator == 1:
        return str(fr.numerator)
    return f"{fr.numerator}/{fr.denominator}"


def dict_repr(counter: Counter[str]) -> str:
    items = sorted(counter.items())
    return "{" + ", ".join(f"{k}:{v}" for k, v in items) + "}"


def row_packet_summary(row: AuditedRow) -> str:
    s = row.speeds
    c14_speed = sum(ramanujan_c(14, v) for v in s)
    c14_sum = sum(ramanujan_c(14, a + b) for a, b in combinations(s, 2))
    c14_diff = sum(ramanujan_c(14, abs(a - b)) for a, b in combinations(s, 2))
    return (
        f"{row.label:26s} qdiv={row.qdiv:3d} safe={fmt(row.safe_measure):>13s} "
        f"tau={sum(tau(v) for v in s):3d} sigma={sum(sigma(v) for v in s):5d} "
        f"c14(v)={c14_speed:4d} c14(sum)={c14_sum:5d} c14(diff)={c14_diff:5d} "
        f"route={proof_route(row)}"
    )


def print_sequence_table(limit: int = 28) -> None:
    print("DIVISOR/RAMANUJAN SEQUENCE SNAPSHOT n<=28")
    print(" n  tau sigma omega bigO unitary c14  gcd14")
    for n in range(1, limit + 1):
        print(
            f"{n:2d} {tau(n):4d} {sigma(n):5d} {omega(n):5d} "
            f"{bigomega(n):4d} {2 ** omega(n):7d} {ramanujan_c(14, n):4d} "
            f"{gcd(14, n):5d}"
        )


def main() -> None:
    seed = named_rows() + one_swap_bank(max_add=220)
    rows = audit_rows(seed)
    audits = [
        audit_quotient(rows, name, fn, note)
        for name, fn, note in QUOTIENTS
    ]

    named = [row for row in rows if row.source == "named"]

    print("LRC14 RAMANUJAN/DIVISOR QUOTIENT GUARDRAIL AUDIT")
    print("=" * 72)
    print("External seed facts used:")
    print("- sigma_k is multiplicative and is Id * 1 under Dirichlet convolution.")
    print("- Ramanujan c_q(n) is the primitive q-th root trace and has a Mobius/gcd formula.")
    print("- sigma_k(n) has a Ramanujan-sum expansion, so scalar divisor data hides phase packets.")
    print()
    print(f"Rows audited: {len(rows)} ({len(named)} named + one-swap AP bank through add<=220)")
    route_counts = Counter(proof_route(row) for row in rows)
    print("Proof-route counts:", dict_repr(route_counts))
    print()

    print_sequence_table(28)
    print()

    print("NAMED ROW PACKETS")
    for row in named:
        print("  " + row_packet_summary(row))
    print()

    print("QUOTIENT FIBER AUDIT")
    print("name                    classes bad_fibers bad_pair_collisions largest_fiber")
    for audit in audits:
        print(
            f"{audit.name:24s} {audit.classes:7d} {audit.bad_fibers:10d} "
            f"{audit.bad_pair_collisions:19d} {audit.largest_fiber:13d}"
        )
    print()

    print("REPRESENTATIVE MIXED FIBERS")
    for audit in audits:
        print(f"\n[{audit.name}] {audit.note}")
        if not audit.examples:
            print("  no mixed proof-route fibers in this audit")
            continue
        for ex_i, example in enumerate(audit.examples, 1):
            print(f"  example {ex_i}:")
            for line in example:
                print("    " + line)
    print()

    print("TOURNAMENT ANALYSIS ON QUOTIENT CANDIDATES")
    fp = tournament_fingerprint(audits)
    print("  pairwise observable: lower (bad_pair_collisions, bad_fibers, classes)")
    print("  vertices: quotient/proof carriers, not runners")
    print(f"  score_hist={fp['score_hist']}")
    print(f"  directed_3_cycles={fp['c3']}")
    print(f"  SCC_sizes={fp['scc_sizes']}")
    print(f"  Hamiltonian_path_count={fp['hamiltonian_paths']}")
    print(f"  tie_Hamiltonian_path={fp['tie_hamiltonian_path']}")
    print()

    print("GUARDRAIL THEOREM TARGET")
    print(
        "  A quotient Q is admissible for an LRC14 proof step P only when P is "
        "constant on Q-fibers, or when Q carries a named certificate for every "
        "forgotten label needed to recover P."
    )
    print(
        "  In this audit, scalar divisor and unitary-divisor quotients mix rows "
        "with different qdiv/safe-route predicates.  They may be features, but "
        "they are not proof carriers for the LRC route predicate."
    )
    print(
        "  Ramanujan speed traces improve the arithmetic interface but still "
        "forget pair/end-point ownership.  Pair and exact-period packets are the "
        "first quotient types that approach admissibility without using the full row."
    )
    print()

    print("BUCKET OF INQUIRY")
    print("  1. Prove a fiber theorem: exact-period Ramanujan packets are route-homogeneous")
    print("     after adding endpoint owner labels, AP/GW equality labels, and K33 state-lift debt.")
    print("  2. Compare c_14(v_i+v_j) against HYP-2970 zero-credit endpoint sums")
    print("     K=14(rm-sn)+r+s; this should be the primitive-root trace of the same equality.")
    print("  3. Test shifted Carmichael/Ramanujan autocorrelation of danger multiplicity")
    print("     against HYP-2973 count moments and HYP-2974 Toeplitz PSD certificates.")
    print("  4. Keep multiplicative functions as irreducibility ledgers, not scalar endpoints:")
    print("     perfect/Farey product, unitary design, Pollock, tiling, and unit-distance analogies")
    print("     all require their incidence/unit labels before quotienting.")


if __name__ == "__main__":
    main()
