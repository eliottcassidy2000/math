#!/usr/bin/env python3
"""HYP-2676 / T916: signed Erdos-Turan packet estimate scout.

The previous packet-alignment scout (HYP-2674) found that the dangerous
finite far-peel rows have six positive one-missed-sector packets.  This file
keeps that exact packet decomposition, then overlays additive-doubling and
small-model diagnostics inspired by Ruzsa modeling and Plunnecke-Ruzsa growth.

It is a scout, not a proof.  It tests the sharper claim:

  same-sign packet bias is a low-rank/Freiman-model phenomenon; once additive
  growth is genuinely broad, the signed packet estimate gains cancellation or
  small mass before the LRC(14) cap margin is threatened.

All LRC quantities below use exact Fraction arithmetic.  The only floating
display is the already-proved THM-546 comparison constant
    kappa = 1.856901...
for |Delta_w| <= kappa*V(E')/(pi^2*w).
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from functools import reduce
from itertools import combinations
from math import gcd, pi
from typing import Iterable


F = Fraction
KAPPA = 1.856901
CAP9_MINUS_Q8 = F(129643, 980980)
ALL_INNER = 0b1111110


def lcm(a: int, b: int) -> int:
    return a // gcd(a, b) * b


def primitive(row: Iterable[int]) -> bool:
    nonzero = [abs(x) for x in row if x]
    return bool(nonzero) and reduce(gcd, nonzero) == 1


def normalized(row: Iterable[int]) -> tuple[int, ...]:
    return tuple(sorted(set(row)))


def wall_breakpoints(row: tuple[int, ...]) -> tuple[int, list[int]]:
    nonzero = [x for x in row if x]
    if not nonzero:
        return 1, [0, 1]
    ell = reduce(lcm, nonzero)
    den = 7 * ell
    bps = {0, den}
    for e in nonzero:
        step = den // (7 * e)
        for a in range(0, 7 * e + 1):
            bps.add(a * step)
    return den, sorted(bps)


def missed_distribution(row: tuple[int, ...]) -> tuple[F, ...]:
    nonzero = [x for x in row if x]
    if not nonzero:
        return tuple([F(0)] * 6 + [F(1)])
    den, bps = wall_breakpoints(row)
    ell = den // 7
    den2 = 2 * ell
    nums = [0] * 7
    for lo, hi in zip(bps, bps[1:]):
        if hi <= lo:
            continue
        midnum = lo + hi
        mask = 0
        for e in nonzero:
            mask |= 1 << ((e * midnum // den2) % 7)
        missed = 6 - (mask & ALL_INNER).bit_count()
        nums[missed] += hi - lo
    return tuple(F(num, den) for num in nums)


def phi(row: tuple[int, ...]) -> F:
    p = missed_distribution(row)
    return p[0] + p[1] / 7


def p0(row: tuple[int, ...]) -> F:
    return missed_distribution(row)[0]


def frac01(x: F) -> F:
    return x - (x.numerator // x.denominator)


def g0(y: F) -> F:
    y = frac01(y)
    if y < F(1, 7):
        return y * F(6, 7)
    return F(6, 49) - (y - F(1, 7)) * F(1, 7)


def one_missed_runs(row: tuple[int, ...]) -> tuple[int, list[tuple[int, int, int]]]:
    """Maximal intervals where row misses exactly one inner sector.

    The return denominator den means an interval (a,b,s) is [a/den,b/den].
    """

    nonzero = [x for x in row if x]
    if not nonzero:
        return 1, []
    den, bps = wall_breakpoints(row)
    ell = den // 7
    den2 = 2 * ell
    cells: list[tuple[int, int, int | None]] = []
    for lo, hi in zip(bps, bps[1:]):
        if hi <= lo:
            continue
        midnum = lo + hi
        mask = 0
        for e in nonzero:
            mask |= 1 << ((e * midnum // den2) % 7)
        missed = [j for j in range(1, 7) if not (mask & (1 << j))]
        cells.append((lo, hi, missed[0] if len(missed) == 1 else None))

    runs: list[tuple[int, int, int]] = []
    i = 0
    while i < len(cells):
        lo, hi, s = cells[i]
        if s is None:
            i += 1
            continue
        a = lo
        j = i
        while j + 1 < len(cells) and cells[j + 1][2] == s and cells[j + 1][0] == cells[j][1]:
            j += 1
        runs.append((a, cells[j][1], s))
        i = j + 1
    return den, runs


def packet_from_runs(
    den: int, runs: list[tuple[int, int, int]], w: int
) -> tuple[F, dict[int, F], dict[int, F], list[F]]:
    packets = {s: F(0) for s in range(1, 7)}
    packet_abs = {s: F(0) for s in range(1, 7)}
    terms: list[F] = []
    for lo, hi, s in runs:
        term = g0(F(w * hi, den) - F(s, 7)) - g0(F(w * lo, den) - F(s, 7))
        packets[s] += term
        packet_abs[s] += abs(term)
        terms.append(term)
    return sum(packets.values()), packets, packet_abs, terms


def packet_report(Ep: tuple[int, ...], w: int, check_direct: bool = True) -> tuple[F, dict[int, F], dict[int, F], list[F], int]:
    den, runs = one_missed_runs(Ep)
    wdelta, packets, packet_abs, terms = packet_from_runs(den, runs, w)
    delta = wdelta / w
    if check_direct:
        direct = p0(normalized(Ep + (w,))) - phi(Ep)
        if direct != delta:
            raise AssertionError((Ep, w, direct, delta))
    return delta, packets, packet_abs, terms, len(runs)


def sign_word(packets: dict[int, F]) -> str:
    out = []
    for s in range(1, 7):
        value = packets[s]
        out.append("+" if value > 0 else "-" if value < 0 else "0")
    return "".join(out)


def fmt(q: F) -> str:
    return f"{q} ({float(q):.9f})"


def sumset(row: tuple[int, ...], times: int = 2) -> set[int]:
    out = {0}
    for _ in range(times):
        out = {x + a for x in out for a in row}
    return out


def sumset_excess(row: tuple[int, ...]) -> int:
    return len(sumset(row, 2)) - (2 * len(row) - 1)


def additive_energy(row: tuple[int, ...]) -> int:
    counts = Counter(a + b for a in row for b in row)
    return sum(c * c for c in counts.values())


def longest_run(row: tuple[int, ...]) -> int:
    best = cur = 1
    for a, b in zip(row, row[1:]):
        if b == a + 1:
            cur += 1
            best = max(best, cur)
        else:
            cur = 1
    return best


def longest_ap(row: tuple[int, ...]) -> int:
    present = set(row)
    best = 1
    span = row[-1] - row[0]
    for a in row:
        for d in range(1, span + 1):
            length = 0
            x = a
            while x in present:
                length += 1
                x += d
            best = max(best, length)
    return best


def squarefree_kernel(n: int) -> int:
    if n == 0:
        return 0
    n = abs(n)
    out = 1
    p = 2
    while p * p <= n:
        parity = 0
        while n % p == 0:
            n //= p
            parity ^= 1
        if parity:
            out *= p
        p += 1 if p == 2 else 2
    if n > 1:
        out *= n
    return out


def squarefree_profile(row: tuple[int, ...]) -> tuple[tuple[int, int], ...]:
    return tuple(sorted(Counter(squarefree_kernel(x) for x in row if x).items()))


def additive_profile(row: tuple[int, ...]) -> dict[str, object]:
    k = len(row)
    k2 = F(len(sumset(row, 2)), k)
    k3 = F(len(sumset(row, 3)), k)
    return {
        "span": row[-1] - row[0],
        "density": F(k, row[-1] - row[0] + 1),
        "sumset_excess": sumset_excess(row),
        "K2": k2,
        "K3": k3,
        "K3_over_K2sq": F(k3.numerator * k2.denominator * k2.denominator, k3.denominator * k2.numerator * k2.numerator),
        "energy": additive_energy(row),
        "longest_run": longest_run(row),
        "longest_ap": longest_ap(row),
        "gap_count": sum(1 for a, b in zip(row, row[1:]) if b > a + 1),
        "sqfree": squarefree_profile(row),
    }


def excess_bucket(exc: int) -> str:
    if exc == 0:
        return "AP"
    if exc <= 4:
        return "near_AP"
    if exc <= 10:
        return "small"
    if exc <= 18:
        return "medium"
    return "large"


@dataclass(frozen=True)
class RowPacket:
    name: str
    Ep: tuple[int, ...]
    w: int
    delta: F
    wdelta: F
    sign: str
    sector_abs: F
    run_abs: F
    runs: int
    profile: dict[str, object]

    @property
    def thm546_bound(self) -> float:
        return KAPPA * self.runs / (pi * pi * self.w)

    @property
    def risk_ratio(self) -> float:
        bound = self.thm546_bound
        return abs(float(self.delta)) / bound if bound else 0.0


def build_row_packet(name: str, Ep: tuple[int, ...], w: int, check_direct: bool = True) -> RowPacket:
    Ep = normalized(Ep)
    delta, packets, packet_abs, terms, runs = packet_report(Ep, w, check_direct=check_direct)
    wdelta = delta * w
    return RowPacket(
        name=name,
        Ep=Ep,
        w=w,
        delta=delta,
        wdelta=wdelta,
        sign=sign_word(packets),
        sector_abs=sum(abs(v) for v in packets.values()),
        run_abs=sum(abs(v) for v in terms),
        runs=runs,
        profile=additive_profile(Ep),
    )


def print_row_packet(rep: RowPacket) -> None:
    print(rep.name)
    print(f"  E'={rep.Ep}, w={rep.w}, full={normalized(rep.Ep + (rep.w,))}")
    print(f"  Delta_w={fmt(rep.delta)}  wDelta={fmt(rep.wdelta)}  margin9-Delta={fmt(CAP9_MINUS_Q8 - rep.delta)}")
    print(f"  packet sign={rep.sign}  sector_abs={fmt(rep.sector_abs)}  run_abs={fmt(rep.run_abs)}")
    print(
        "  cancellation: "
        f"|sum|/sector_abs={fmt(abs(rep.wdelta) / rep.sector_abs) if rep.sector_abs else 'n/a'}; "
        f"|sum|/run_abs={fmt(abs(rep.wdelta) / rep.run_abs) if rep.run_abs else 'n/a'}"
    )
    print(
        "  THM546: "
        f"V={rep.runs}, bound={rep.thm546_bound:.9f}, "
        f"|Delta|/bound={rep.risk_ratio:.6f}"
    )
    prof = rep.profile
    print(
        "  Ruzsa/Plunnecke profile: "
        f"span={prof['span']} density={prof['density']} exc={prof['sumset_excess']} "
        f"K2={prof['K2']} K3={prof['K3']} K3/K2^2={prof['K3_over_K2sq']} "
        f"energy={prof['energy']} run={prof['longest_run']} ap={prof['longest_ap']} gaps={prof['gap_count']}"
    )
    print(f"  squarefree profile={prof['sqfree']}")
    print()


def named_rows() -> list[RowPacket]:
    specs = [
        ("consec8_w10", (0, 1, 2, 3, 4, 5, 6, 7), 10),
        ("B13_same_sign_packet", (0, 1, 2, 4, 6, 7, 8, 10), 12),
        ("HYP2671_dyadic_block", (0, 1, 2, 4, 8, 12, 16, 20), 24),
        ("HYP2672_doubled_odd_tail", (0, 1, 2, 4, 8, 14, 26, 34), 38),
        ("nonshell_warning", (0, 2, 3, 5, 6, 15), 18),
        ("HYP2675_boundary_leader_base", (0, 2, 4, 6, 8, 10, 12, 14), 15),
        ("HYP2675_true_wide_base", (0, 4, 6, 8, 10, 12, 14, 15), 16),
        ("KPS_third_pocket", (0, 3, 5, 16, 28, 30, 33), 35),
    ]
    return [build_row_packet(name, Ep, w) for name, Ep, w in specs]


def dyadic_family_scan(s_max: int = 160) -> list[RowPacket]:
    rows: list[RowPacket] = []
    for mult in (6, 10):
        for s in range(3, s_max + 1):
            Ep = (0, 1, 2, 4, 8, 3 * s, 4 * s, 5 * s)
            rows.append(build_row_packet(f"dyadic_m{s}_w{mult}s", Ep, mult * s, check_direct=False))
    return rows


def print_dyadic_summary(rows: list[RowPacket]) -> list[RowPacket]:
    print("Dyadic/Bertrand scale ledger")
    print("  family E_s={0,1,2,4,8,3s,4s,5s}, exact packet scan for s=3..160")
    leaders: list[RowPacket] = []
    for mult in (6, 10):
        subset = [r for r in rows if r.name.endswith(f"w{mult}s")]
        global_leader = max(subset, key=lambda r: r.delta)
        leaders.append(global_leader)
        print(f"  w={mult}s global: s={global_leader.w // mult}, Delta={fmt(global_leader.delta)}, sign={global_leader.sign}, exc={global_leader.profile['sumset_excess']}")
        print("  dyadic blocks:")
        lo = 3
        while lo <= 128:
            hi = min(2 * lo - 1, 160)
            block = [r for r in subset if lo <= r.w // mult <= hi]
            if block:
                leader = max(block, key=lambda r: r.delta)
                leaders.append(leader)
                print(
                    f"    s={lo:3d}..{hi:3d}: best s={leader.w // mult:3d}, "
                    f"Delta={fmt(leader.delta)}, sign={leader.sign}, "
                    f"|sum|/run_abs={fmt(abs(leader.wdelta) / leader.run_abs) if leader.run_abs else 'n/a'}"
                )
            lo *= 2
    print()
    return leaders


def scan_b14_bank() -> tuple[list[RowPacket], dict[str, dict[str, object]]]:
    """Near-speed B14 bank: E'={0}+7-subsets of [1,14], w=max+1..max+4."""

    top: list[RowPacket] = []
    buckets: dict[str, dict[str, object]] = {}
    for comb in combinations(range(1, 15), 7):
        Ep = (0,) + comb
        if not primitive(Ep):
            continue
        den, runs = one_missed_runs(Ep)
        profile = additive_profile(Ep)
        bucket = excess_bucket(int(profile["sumset_excess"]))
        if bucket not in buckets:
            buckets[bucket] = {
                "count": 0,
                "same_positive": 0,
                "max_delta": F(-10),
                "best": None,
                "max_eff": F(0),
                "eff_best": None,
            }
        for w in range(Ep[-1] + 1, Ep[-1] + 5):
            if not primitive(Ep + (w,)):
                continue
            wdelta, packets, _packet_abs, terms = packet_from_runs(den, runs, w)
            delta = wdelta / w
            sign = sign_word(packets)
            run_abs = sum(abs(v) for v in terms)
            sector_abs = sum(abs(v) for v in packets.values())
            rep = RowPacket(
                name=f"B14_near_Emax{Ep[-1]}_w{w}",
                Ep=Ep,
                w=w,
                delta=delta,
                wdelta=wdelta,
                sign=sign,
                sector_abs=sector_abs,
                run_abs=run_abs,
                runs=len(runs),
                profile=profile,
            )
            data = buckets[bucket]
            data["count"] = int(data["count"]) + 1
            if sign == "++++++":
                data["same_positive"] = int(data["same_positive"]) + 1
            if delta > data["max_delta"]:
                data["max_delta"] = delta
                data["best"] = rep
            eff = abs(wdelta) / run_abs if run_abs else F(0)
            if eff > data["max_eff"]:
                data["max_eff"] = eff
                data["eff_best"] = rep
            top.append(rep)
            top.sort(key=lambda r: (r.delta, -int(r.profile["sumset_excess"]), tuple(-x for x in r.Ep)), reverse=True)
            del top[12:]
    return top, buckets


def print_b14_summary(top: list[RowPacket], buckets: dict[str, dict[str, object]]) -> list[RowPacket]:
    print("B14 near-speed packet/Ruzsa bank")
    print("  rows: E'={0}+7-subsets of [1,14], w=max(E')+1..max(E')+4, primitive")
    print("  bucket summary by additive sumset excess:")
    leaders: list[RowPacket] = []
    for bucket in ["AP", "near_AP", "small", "medium", "large"]:
        if bucket not in buckets:
            continue
        data = buckets[bucket]
        best = data["best"]
        eff_best = data["eff_best"]
        if isinstance(best, RowPacket):
            leaders.append(best)
        if isinstance(eff_best, RowPacket):
            leaders.append(eff_best)
        print(
            f"    {bucket:8s}: count={data['count']:5d}, ++++++={data['same_positive']:4d}, "
            f"maxDelta={fmt(data['max_delta'])}, maxEff={fmt(data['max_eff'])}"
        )
        if isinstance(best, RowPacket):
            print(
                f"      best Delta: E'={best.Ep}, w={best.w}, sign={best.sign}, "
                f"exc={best.profile['sumset_excess']}, K2={best.profile['K2']}"
            )
        if isinstance(eff_best, RowPacket) and eff_best != best:
            print(
                f"      best signed efficiency: E'={eff_best.Ep}, w={eff_best.w}, sign={eff_best.sign}, "
                f"Delta={fmt(eff_best.delta)}"
            )
    print("  top positive Delta rows:")
    for i, rep in enumerate(top, 1):
        print(
            f"    {i:2d}. E'={rep.Ep}, w={rep.w}, Delta={fmt(rep.delta)}, "
            f"sign={rep.sign}, exc={rep.profile['sumset_excess']}, "
            f"K2={rep.profile['K2']}, eff={fmt(abs(rep.wdelta) / rep.run_abs) if rep.run_abs else 'n/a'}"
        )
    print()
    return leaders + top[:8]


def tournament_analysis(rows: list[RowPacket]) -> None:
    dedup: dict[tuple[tuple[int, ...], int], RowPacket] = {}
    for row in rows:
        key = (row.Ep, row.w)
        if key not in dedup or row.risk_ratio > dedup[key].risk_ratio:
            dedup[key] = row
    verts = sorted(dedup.values(), key=lambda r: (-r.risk_ratio, -float(r.delta), r.Ep, r.w))[:12]
    score = Counter()
    wins: set[tuple[int, int]] = set()
    for i, a in enumerate(verts):
        for j, b in enumerate(verts):
            if i >= j:
                continue
            # Pairwise observable: larger share of the THM-546 BV budget is the
            # sharper signed-packet proof obligation.  Ties use larger Delta and
            # then lexicographic row order as the Hamiltonian tie path.
            a_key = (a.risk_ratio, float(a.delta), tuple(-x for x in a.Ep), -a.w)
            b_key = (b.risk_ratio, float(b.delta), tuple(-x for x in b.Ep), -b.w)
            winner = i if a_key >= b_key else j
            loser = j if winner == i else i
            wins.add((winner, loser))
            score[winner] += 1
            score.setdefault(loser, score[loser])
    cycles = 0
    for i, j, k in combinations(range(len(verts)), 3):
        if (
            (i, j) in wins and (j, k) in wins and (k, i) in wins
        ) or (
            (i, k) in wins and (k, j) in wins and (j, i) in wins
        ):
            cycles += 1
    path = sorted(range(len(verts)), key=lambda i: (-score[i], -verts[i].risk_ratio, -float(verts[i].delta), verts[i].Ep))
    print("Tournament Analysis")
    print("  vertices: packet proof obligations / additive profiles, not runners or arcs")
    print("  pairwise observable: larger |Delta_w| share of the THM-546 BV budget wins")
    print("  switch/gauge: sign of the far packet; ties follow Delta then lexicographic row")
    print(f"  score_hist={dict(sorted(Counter(score.values()).items()))}, directed_3cycles={cycles}")
    print("  Hamiltonian path:")
    for rank, idx in enumerate(path, 1):
        rep = verts[idx]
        print(
            f"    {rank:2d}. {rep.name}: risk={rep.risk_ratio:.6f}, "
            f"Delta={fmt(rep.delta)}, sign={rep.sign}, exc={rep.profile['sumset_excess']}, E'={rep.Ep}, w={rep.w}"
        )
    print("  challenged assumption: the signed-packet vertices do not have to be runners.")
    print("    This quotient preserves the far-peel discrepancy and additive model class;")
    print("    it destroys the full endpoint state word, so it must feed back into HYP-2648/2675.")
    print()


def main() -> None:
    print("HYP-2676 LRC14 signed Erdos-Turan packet / Ruzsa-model scout")
    print("Arithmetic: exact Fractions for packet Delta; float display only for THM-546.")
    print(f"k=9 comparison margin cap_9-Q(8) = {fmt(CAP9_MINUS_Q8)}")
    print()

    named = named_rows()
    print("Named packet rows")
    for rep in named:
        print_row_packet(rep)

    family_rows = dyadic_family_scan()
    family_leaders = print_dyadic_summary(family_rows)

    top_b14, buckets = scan_b14_bank()
    bank_leaders = print_b14_summary(top_b14, buckets)

    tournament_analysis(named + family_leaders + bank_leaders)

    print("Proof-route reading")
    print("  1. Same-sign packet rows are not random ET errors; they live in compressed additive models.")
    print("  2. Ruzsa modeling should split low-K rows into finite cyclic/GAP models before scalar bounds.")
    print("  3. Plunnecke-Ruzsa growth is useful in the complementary branch: larger growth spreads endpoints,")
    print("     making the signed sector packet a cancellation problem rather than a same-sign pocket.")
    print("  4. Landau/Chebyshev/Bertrand analogy: local packet bias can persist on a finite scale,")
    print("     but the scale ledger asks for a reset in each dyadic block after the finite resonance.")
    print("  5. No LRC(14) proof is claimed here; this turns the next obligation into a finite-model")
    print("     same-sign classification plus a signed ET packet bound on high-growth rows.")
    print("PASS: signed packet/Ruzsa scout complete.")


if __name__ == "__main__":
    main()
