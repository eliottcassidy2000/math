"""LRC14 petal themes and marked-Farey proof attack.

This scout extends the AP/Goddyn-Wong petal story in three directions:

1. Flower denominators: compare the 22-family from 1/pi with the LRC petal
   denominators 23, 27, 41 and the pair-sum modulus 27.
2. Petal migration: exact max-min marks for AP single replacements
   {1,...,13} - {h} + {r}, using the marked Farey node M(S).
3. Proof attack: treat filters as vertices in a tournament-style proof carrier
   and isolate the terminal single-petal lemma still needed for LRC14.

The computation is exact for all LRC max-min values.  Floating point is only
used in the pi/flower readout.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import combinations
from math import gcd, pi
import sys


sys.stdout.reconfigure(line_buffering=True)

N = 14
THR = F(1, N)
AP = tuple(range(1, N))
GW = tuple(list(range(1, 12)) + [13, 24])
UNITS = (1, 3, 5, 9, 11, 13)
EVENS = (2, 4, 6, 8, 10, 12)
ODD_SKELETON = (1, 3, 5, 7, 9, 11, 13)
PAIR_SUM_MODULUS = 2 * N - 1
FILTER_LIMIT = 300
EXACT_LIMIT = 80


@dataclass(frozen=True)
class Record:
    S: tuple[int, ...]
    mark: F
    argmax: tuple[F, ...]
    qdiv: int
    status: str


def norm1(x: F) -> F:
    n = x.numerator % x.denominator
    d = x.denominator
    if 2 * n <= d:
        return F(n, d)
    return F(d - n, d)


def primitive(S: tuple[int, ...]) -> bool:
    g = 0
    for s in S:
        g = gcd(g, s)
    return g == 1


def candidate_taus(S: tuple[int, ...]) -> set[F]:
    S = tuple(sorted(set(S)))
    cands: set[F] = {F(1, 2)}
    for i, a in enumerate(S):
        k = 0
        while True:
            t = F(2 * k + 1, 2 * a)
            if t > F(1, 2):
                break
            cands.add(t)
            k += 1
        for b in S[i + 1 :]:
            for d in (a + b, b - a):
                if d <= 0:
                    continue
                k = 1
                while True:
                    t = F(k, d)
                    if t > F(1, 2):
                        break
                    cands.add(t)
                    k += 1
    return cands


_M_CACHE: dict[tuple[int, ...], tuple[F, tuple[F, ...]]] = {}


def M_exact(S: tuple[int, ...]) -> tuple[F, tuple[F, ...]]:
    S = tuple(sorted(set(S)))
    if S in _M_CACHE:
        return _M_CACHE[S]
    best = F(0)
    pts: list[F] = []
    for t in candidate_taus(S):
        val = min(norm1(s * t) for s in S)
        if val > best:
            best = val
            pts = [t]
        elif val == best:
            pts.append(t)
    out = (best, tuple(sorted(pts)))
    _M_CACHE[S] = out
    return out


def q_threshold(S: tuple[int, ...]) -> int:
    q = 2
    while any(s % q == 0 for s in S):
        q += 1
    return q


def residue_counts(S: tuple[int, ...]) -> tuple[int, ...]:
    c = [0] * N
    for s in S:
        c[s % N] += 1
    return tuple(c)


def q_punctured_cover(S: tuple[int, ...]) -> bool:
    return all(any(s % q == 0 for s in S) for q in range(2, N)) and all(s % N for s in S)


def unit_cover(S: tuple[int, ...]) -> bool:
    c = residue_counts(S)
    return all(c[u] >= 1 for u in UNITS)


def exact_unit_shell(S: tuple[int, ...]) -> bool:
    c = residue_counts(S)
    return c[0] == 0 and all(c[u] == 1 for u in UNITS)


def exact_odd_skeleton(S: tuple[int, ...]) -> bool:
    c = residue_counts(S)
    return c[0] == 0 and all(c[r] == 1 for r in ODD_SKELETON)


def cofinite_residue_set(S: tuple[int, ...]) -> bool:
    c = residue_counts(S)
    return c[0] == 0 and sum(1 for r in range(1, N) if c[r]) >= 12


def zsum_max(S: tuple[int, ...]) -> bool:
    return sum(1 for i, a in enumerate(S) for b in S[i + 1 :] if (a + b) % N == 0) == 6


def complement_sum14(S: tuple[int, ...]) -> bool:
    by_res: dict[int, list[int]] = {r: [] for r in range(N)}
    for s in S:
        by_res[s % N].append(s)
    for r in (1, 3, 5):
        a = by_res[r]
        b = by_res[N - r]
        if len(a) != 1 or len(b) != 1:
            return False
        if a[0] + b[0] != N:
            return False
    return True


def even_dipole_data(S: tuple[int, ...]) -> tuple[str, tuple[int, ...], tuple[int, ...]]:
    c = residue_counts(S)
    defect = [c[r] - (1 if r else 0) for r in range(N)]
    support = tuple(r for r, d in enumerate(defect) if d)
    missing = tuple(r for r, d in enumerate(defect) if d < 0)
    extra = tuple(r for r, d in enumerate(defect) if d > 0)
    if not support:
        return ("zero", missing, extra)
    if (
        set(support).issubset(set(EVENS))
        and len(missing) == 1
        and len(extra) == 1
        and sum(abs(defect[r]) for r in support) == 2
    ):
        return ("single_even_dipole", missing, extra)
    return ("other", missing, extra)


def even_dipole_or_zero(S: tuple[int, ...]) -> bool:
    return even_dipole_data(S)[0] in {"zero", "single_even_dipole"}


def one_petal_or_ap(S: tuple[int, ...]) -> bool:
    kind, missing, extra = even_dipole_data(S)
    if kind == "zero":
        return True
    if kind != "single_even_dipole":
        return False
    m = missing[0]
    e = extra[0]
    outside = sorted(s for s in S if s > 13 and s % N == e)
    return len(outside) == 1 and outside[0] == 2 * m


def minimal_one_petal(S: tuple[int, ...]) -> bool:
    kind, missing, extra = even_dipole_data(S)
    if kind == "zero":
        return set(S) == set(AP)
    if kind != "single_even_dipole":
        return False
    m = missing[0]
    e = extra[0]
    return set(S) == (set(AP) - {m}) | {2 * m} and (2 * m) % N == e


def top_petal_or_ap(S: tuple[int, ...]) -> bool:
    return set(S) == set(AP) or set(S) == set(GW)


FILTERS: tuple[tuple[str, object], ...] = (
    ("punctured_q_cover", q_punctured_cover),
    ("unit_cover", unit_cover),
    ("exact_unit_shell", exact_unit_shell),
    ("exact_odd_skeleton", exact_odd_skeleton),
    ("cofinite_residue_set", cofinite_residue_set),
    ("zsum_max", zsum_max),
    ("complement_sum14", complement_sum14),
    ("even_dipole_or_zero", even_dipole_or_zero),
    ("one_petal_or_ap", one_petal_or_ap),
    ("minimal_one_petal", minimal_one_petal),
    ("top_petal_or_ap", top_petal_or_ap),
)


def first_failed_filter(S: tuple[int, ...]) -> str:
    for name, pred in FILTERS:
        if not pred(S):  # type: ignore[operator]
            return name
    return "PASS"


def gw_window(v: int) -> tuple[int, int]:
    return (N - v, 2 * N - 1 - 2 * v)


def gw_window_coprime_blockers(v: int) -> tuple[int, ...]:
    lo, hi = gw_window(v)
    return tuple(x for x in range(lo, hi + 1) if gcd(x, v) == 1)


def gw_gate(v: int) -> bool:
    return len(gw_window_coprime_blockers(v)) == 0


def single_swap_rows(limit: int) -> list[tuple[int, ...]]:
    rows: set[tuple[int, ...]] = set()
    for h in AP:
        kept = set(AP) - {h}
        for r in range(1, limit + 1):
            S = tuple(sorted(kept | {r}))
            if len(S) == 13 and primitive(S):
                rows.add(S)
    return sorted(rows)


def row_name(S: tuple[int, ...]) -> str:
    missing = sorted(set(AP) - set(S))
    added = sorted(set(S) - set(AP))
    if not missing and not added:
        return "AP"
    if len(missing) == 1 and len(added) == 1:
        return f"{missing[0]}->{added[0]}"
    return str(S)


def mark_status(S: tuple[int, ...], mark: F, qdiv: int) -> str:
    if mark == THR:
        return "tight-apex"
    if any(s % N == 0 for s in S):
        return "apex-blocked"
    if qdiv < N:
        return "loose-up"
    excess = N * mark.numerator - mark.denominator
    if excess == 1:
        return "unit-child"
    return "nonunit-child"


def packet(mark: F) -> str:
    a = mark.numerator
    b = mark.denominator
    planar = "nonplanar" if min(a, b) >= 3 else "planar"
    det = b - N * a
    return f"{a}/{b}: sum={a + b}, prod={a * b}, det={det}, K_{{{a},{b}}} {planar}"


def make_records(rows: list[tuple[int, ...]]) -> list[Record]:
    records: list[Record] = []
    for S in rows:
        mark, pts = M_exact(S)
        qdiv = q_threshold(S)
        records.append(Record(S, mark, pts, qdiv, mark_status(S, mark, qdiv)))
    return records


def record_for(S: tuple[int, ...]) -> Record:
    mark, pts = M_exact(S)
    qdiv = q_threshold(S)
    return Record(S, mark, pts, qdiv, mark_status(S, mark, qdiv))


def print_flower_extension() -> None:
    theta = 1 / pi
    delta = 1 / (2 * pi * pi)
    print("[Flower denominators next to LRC petal marks]")
    print(f"  22/7-pi={22 / 7 - pi:+.12e}; 31^(1/3)-pi={31 ** (1 / 3) - pi:+.12e}")
    print(f"  1/pi={theta:.15f}; 7/22={7 / 22:.15f}; error={theta - 7 / 22:+.12e}")
    print("  k-petal returns at theta=1/pi radians:")
    print("    k  k/pi-nearest-integer   normalized-turns-nearest-integer")
    for k in (14, 20, 22, 23, 24, 27, 31, 41):
        rad = k * theta
        rad_err = rad - round(rad)
        turns = k * delta
        turn_err = turns - round(turns)
        print(f"    {k:2d} {rad_err:+.12e} {turn_err:+.12e}")
    print("  Readout: 22 is a radian-denominator family; 27 is the LRC14 pair-sum")
    print("  modulus 2n-1; 41 is the first nonplanar unit Farey child 3/(14*3-1).")
    print()


def print_doubling_petals() -> None:
    print("[Minimal doubling petals h -> 2h]")
    print(
        "   h  row       q  gate blockers        M        status        packet"
    )
    for h in range(7, N):
        r = 2 * h
        S = tuple(sorted((set(AP) - {h}) | {r}))
        rec = record_for(S)
        blockers = gw_window_coprime_blockers(h)
        print(
            f"  {h:2d}  {h}->{r:<3d} {rec.qdiv:2d} {str(gw_gate(h)):>5} "
            f"{str(blockers):<15} {str(rec.mark):<8} {rec.status:<13} {packet(rec.mark)}"
        )
    print("  Only 12->24 is both blocker-free in the GW window and marked at 1/14.")
    print()


def first_add2_tails(h: int, count: int = 3) -> list[int]:
    out: list[int] = []
    r = h + 2
    kept = set(AP) - {h}
    while len(out) < count:
        if r not in kept:
            out.append(r)
        r += 2
    return out


def print_add2_probe() -> None:
    print("[n+2 recursion probe: first valid same-parity tails]")
    print("   h  tails tested        best mark       best row    status")
    for h in AP:
        best: tuple[F, int, Record] | None = None
        tested = first_add2_tails(h)
        for r in tested:
            S = tuple(sorted((set(AP) - {h}) | {r}))
            rec = record_for(S)
            item = (rec.mark, r, rec)
            if best is None or item < best:
                best = item
        assert best is not None
        _, r_best, rec_best = best
        print(
            f"  {h:2d}  {str(tested):<18} {str(rec_best.mark):<14} "
            f"{h}->{r_best:<3d} {rec_best.status}"
        )
    print("  The +2 mode mostly falls into apex-blocked or loose-up residue lifts;")
    print("  it does not reproduce the mark-unital 12->24 doubling petal.")
    print()


def print_single_swap_atlas(records: list[Record], exact_limit: int) -> None:
    print(f"[AP single-replacement exact marked atlas, unique rows v<={exact_limit}]")
    by_status = Counter(rec.status for rec in records)
    by_mark = Counter(rec.mark for rec in records)
    tight = [rec for rec in records if rec.mark == THR]
    below = [rec for rec in records if rec.mark < THR]
    print(f"  rows={len(records)} status_counts={dict(sorted(by_status.items()))}")
    print(f"  tight={len(tight)} below_threshold={len(below)}")
    for rec in tight:
        print(
            f"    TIGHT {row_name(rec.S):<7} q={rec.qdiv} M={rec.mark} "
            f"denoms={sorted({t.denominator for t in rec.argmax})} first_fail={first_failed_filter(rec.S)}"
        )
    print("  most common marked Farey nodes:")
    for mark, count in by_mark.most_common(12):
        print(f"    {str(mark):<7} count={count:4d} {packet(mark)}")
    print()


def print_site_profiles(limit: int) -> None:
    print(f"[Single missing-site profiles for tails 14..{limit}]")
    print("   h  status counts                                      closest mark above 1/14")
    for h in AP:
        site_records = []
        for r in range(N, limit + 1):
            S = tuple(sorted((set(AP) - {h}) | {r}))
            if len(S) == 13 and primitive(S):
                site_records.append(record_for(S))
        counts = Counter(rec.status for rec in site_records)
        loose = [rec for rec in site_records if rec.mark > THR]
        closest = min(loose, key=lambda rec: (rec.mark - THR, max(rec.S))) if loose else None
        if closest is None:
            tail = "-"
        else:
            tail = (
                f"{row_name(closest.S)} M={closest.mark} "
                f"gap={closest.mark - THR} {closest.status}"
            )
        count_text = ", ".join(f"{k}:{counts[k]}" for k in sorted(counts))
        print(f"  {h:2d}  {count_text:<50} {tail}")
    print()


def filter_survivors(rows: list[tuple[int, ...]], through: str) -> list[tuple[int, ...]]:
    current = set(rows)
    for name, pred in FILTERS:
        current = {S for S in current if pred(S)}  # type: ignore[operator]
        if name == through:
            return sorted(current)
    raise ValueError(through)


def print_terminal_core(rows: list[tuple[int, ...]], limit: int) -> None:
    terminal = filter_survivors(rows, "minimal_one_petal")
    print(f"[Terminal core inside the minimal-one-petal premise, bank v<={limit}]")
    print(f"  terminal rows={len(terminal)}")
    for S in terminal:
        rec = record_for(S)
        kind, missing, extra = even_dipole_data(S)
        print(
            f"    {row_name(S):<7} S={S} M={str(rec.mark):<5} status={rec.status:<13} "
            f"defect={kind} missing={missing} extra={extra} first_fail={first_failed_filter(S)}"
        )
    print("  Exact terminal implication in this bank:")
    print("    minimal one-petal + mark-unital  <=>  AP or 12->24.")
    print("  The remaining proof target is to turn this bounded terminal fact into")
    print("  a lemma, or prove every tight survivor must enter this premise.")
    print()


def condition_tournament(rows: list[tuple[int, ...]]) -> None:
    gates: list[tuple[str, object]] = list(FILTERS)
    pass_sets: list[set[tuple[int, ...]]] = []
    for _, pred in gates:
        pass_sets.append({S for S in rows if pred(S)})  # type: ignore[operator]

    n = len(gates)
    adj = [[False] * n for _ in range(n)]
    incomparable = 0
    equal_pairs: list[tuple[str, str]] = []
    for i in range(n):
        for j in range(i + 1, n):
            Ai, Aj = pass_sets[i], pass_sets[j]
            if Ai == Aj:
                adj[i][j] = True
                equal_pairs.append((gates[i][0], gates[j][0]))
            elif Ai < Aj:
                adj[i][j] = True
            elif Aj < Ai:
                adj[j][i] = True
            else:
                incomparable += 1
                if (len(Ai), i) <= (len(Aj), j):
                    adj[i][j] = True
                else:
                    adj[j][i] = True
    cyc3 = 0
    for a, b, c in combinations(range(n), 3):
        if adj[a][b] and adj[b][c] and adj[c][a]:
            cyc3 += 1
        if adj[a][c] and adj[c][b] and adj[b][a]:
            cyc3 += 1
    print("[Tournament Analysis: filter-strength carrier]")
    print("  vertex pass counts:")
    for (name, _), passed in zip(gates, pass_sets):
        print(f"    {name:<22} {len(passed):5d}")
    print(f"  directed_3cycles={cyc3}; incomparable_pairs={incomparable}")
    print(f"  equal pass-set pairs={equal_pairs}")
    print("  The mark-unital check is then applied only to the terminal core: there")
    print("  the top-petal gate is exactly AP/GW, and exact M separates those from")
    print("  the loose terminal petals.")
    print()


def print_farey_exemplars() -> None:
    examples = [
        ("AP", AP),
        ("GW", GW),
        ("apex blocked 7->14", tuple(sorted((set(AP) - {7}) | {14}))),
        ("same-residue liar 12->26", tuple(sorted((set(AP) - {12}) | {26}))),
        ("nonunit petal 8->16", tuple(sorted((set(AP) - {8}) | {16}))),
        ("unit child 10->20", tuple(sorted((set(AP) - {10}) | {20}))),
        ("nonplanar unit child 12->36", tuple(sorted((set(AP) - {12}) | {36}))),
    ]
    print("[Farey/summand/multiplicand exemplars]")
    for label, S in examples:
        rec = record_for(S)
        print(
            f"  {label:<28} {row_name(S):<7} q={rec.qdiv:<2} M={str(rec.mark):<5} "
            f"{rec.status:<13} {packet(rec.mark)}"
        )
    print("  The first child with numerator 3 is 3/41, hence the first child whose")
    print("  complete-bipartite packet K_{3,41} contains the K_{3,3} obstruction.")
    print()


def print_proof_attack() -> None:
    print("[Proof attack distilled]")
    print("  1. Existing theorem-level gates: q=2..13 covered, q=14 omitted, and")
    print("     unit residues covered. Failure here is loose-up or apex-blocked.")
    print("  2. Existing AP-tail measure theorems THM-542/543/544 close the one-tail")
    print("     and two-tail below-second AP collar; the petal attack should not redo")
    print("     those exact interval certificates.")
    print("  3. Petal lemma target: once the proof enters the minimal one-petal branch,")
    print("     the only mark-unital terminal rows are AP and 12->24. The exact escape")
    print("     witnesses are 8->16 with M=2/23 and 10->20 with M=2/27.")
    print("  4. The GW window explains why 12 is special: [14-12,27-24]=[2,3] has no")
    print("     coprime blocker for 12, while every other h has a blocker.")
    print("  5. Remaining global work: prove every primitive tight survivor either")
    print("     enters that one-petal branch or is killed by a multi-dipole/far-tail")
    print("     theorem preserving the marked Farey node.")


def main() -> None:
    print("LRC14 PETAL THEMES AND PROOF ATTACK")
    print("=" * 78)
    print_flower_extension()
    print_doubling_petals()
    print_add2_probe()

    exact_rows = single_swap_rows(EXACT_LIMIT)
    records = make_records(exact_rows)
    print_single_swap_atlas(records, EXACT_LIMIT)
    print_site_profiles(EXACT_LIMIT)

    rows = single_swap_rows(FILTER_LIMIT)
    print_terminal_core(rows, FILTER_LIMIT)
    print_farey_exemplars()
    condition_tournament(rows)
    print_proof_attack()


if __name__ == "__main__":
    main()
