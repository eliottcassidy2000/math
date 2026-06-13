#!/usr/bin/env python3
"""
lrc_prime_power_zero_branch_core_s546b.py

Prime-power zero-branch cover cores for LRC.

This supplements the product-tree HYP-2036 scan.  For each prime power
q=p^d <= n, THM-369
says that if no speed is divisible by q then t=u/q, (u,q)=1, is a sieve
witness: every runner is at closed distance at least 1/q >= 1/n from the
observer.  Thus a possible open-cover counterexample must put at least one
speed in the zero branch 0 mod q.

The question here is what that zero branch buys.  It certainly kills the unit
point u/q, because a q-divisible runner lands at the observer.  But the local
danger intervals centered at those unit points are nested stars.  This script
checks, exactly, whether those zero-branch local covers have any nonpeeling
endpoint-protection core.  In all audited prime-power branches they peel to
empty.  The obstruction is therefore exported to descendant/event layers, not
stored in a local q-zero branch core.

Tournament Analysis declaration:
    vertices:
        prime-power branches q=p^d (plus an optional n-gate diagnostic)
    pairwise observable:
        zero-branch occupancy, local endpoint-core size, and branch radius
    switch/gauge:
        a branch with more zero carriers / larger local radius / larger core
        beats another; exact equality is a t(r)ienerment tie
    tie Hamiltonian path:
        increasing q
    fingerprints:
        tie counts, score histograms, strict directed triangles, SCC sizes,
        and Hamiltonian path counts after resolving ties by q-order.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations
from math import gcd


ZERO = Fraction(0)
ONE = Fraction(1)


@dataclass(frozen=True)
class Case:
    name: str
    n: int
    speeds: tuple[int, ...]


CASES = (
    Case("n6_AP_wall", 6, (1, 2, 3, 4, 5)),
    Case("n6_Z3_bridge_open", 6, (1, 3, 6, 9, 12)),
    Case("n8_AP_wall", 8, tuple(range(1, 8))),
    Case("n14_AP_wall", 14, tuple(range(1, 14))),
    Case("n18_AP_wall", 18, tuple(range(1, 18))),
    Case("n18_no_9_zero_branch", 18, (1, 2, 3, 4, 5, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16, 17, 19)),
    Case("n18_18_gate_skip8", 18, tuple(sorted(set(range(1, 18)) - {8} | {18}))),
    Case("n18_36_double_gate_skip8", 18, tuple(sorted(set(range(1, 18)) - {8} | {36}))),
    Case("n18_deep_3_zero_branch", 18, tuple([1] + [3 * k for k in range(1, 17)])),
)


def nstar(n: int) -> int:
    return n // 2 if n % 2 == 0 else n


def factor(n: int) -> dict[int, int]:
    out: dict[int, int] = {}
    d = 2
    while d * d <= n:
        while n % d == 0:
            out[d] = out.get(d, 0) + 1
            n //= d
        d += 1
    if n > 1:
        out[n] = out.get(n, 0) + 1
    return out


def prime_powers_upto(n: int) -> list[int]:
    powers = set()
    for p in range(2, n + 1):
        if factor(p) == {p: 1}:  # p is prime
            q = p
            while q <= n:
                powers.add(q)
                q *= p
    return sorted(powers)


def is_prime_power(q: int) -> bool:
    f = factor(q)
    return len(f) == 1


def p_adic_val(x: int, p: int, cap: int | None = None) -> int:
    if x == 0:
        return cap if cap is not None else 10**9
    v = 0
    while x % p == 0:
        v += 1
        x //= p
        if cap is not None and v >= cap:
            return cap
    return v


def frac(x: Fraction) -> Fraction:
    return x - Fraction(x.numerator // x.denominator)


def add_interval_mod(out: list[tuple[Fraction, Fraction, str]], left: Fraction, right: Fraction, owner: str) -> None:
    while left < ZERO:
        out.append((left + ONE, ONE, owner))
        left += ONE
        right += ONE
    while right > ONE:
        out.append((ZERO, right - ONE, owner))
        right -= ONE
        left -= ONE
    if left < right:
        out.append((left, right, owner))


def danger_intervals(speeds: tuple[int, ...], n: int) -> list[tuple[Fraction, Fraction, str]]:
    intervals: list[tuple[Fraction, Fraction, str]] = []
    threshold = Fraction(1, n)
    for s in speeds:
        for k in range(s):
            center = Fraction(k, s)
            add_interval_mod(intervals, center - threshold / s, center + threshold / s, f"v={s}")
    return intervals


def branch_intervals(speeds: tuple[int, ...], n: int, q: int) -> list[tuple[Fraction, Fraction, str]]:
    """Intervals centered at primitive unit points u/q from q-divisible speeds."""
    intervals: list[tuple[Fraction, Fraction, str]] = []
    zero_speeds = [s for s in speeds if s % q == 0]
    for s in zero_speeds:
        half = Fraction(1, n * s)
        for u in range(1, q):
            if gcd(u, q) != 1:
                continue
            center = Fraction(u, q)
            intervals.append((center - half, center + half, f"q={q},v={s},u={u}"))
    return intervals


def inside_open_interval(point: Fraction, interval: tuple[Fraction, Fraction, str]) -> bool:
    left, right, _ = interval
    return left < point < right


def endpoint_core(intervals: list[tuple[Fraction, Fraction, str]]) -> tuple[int, int, tuple[int, ...]]:
    """Peel intervals whose left or right endpoint is not strictly protected."""
    active = set(range(len(intervals)))
    layer_sizes = []
    while True:
        remove = set()
        for idx in active:
            left, right, _ = intervals[idx]
            left_protected = any(
                j != idx and j in active and inside_open_interval(left, intervals[j])
                for j in active
            )
            right_protected = any(
                j != idx and j in active and inside_open_interval(right, intervals[j])
                for j in active
            )
            if not (left_protected and right_protected):
                remove.add(idx)
        if not remove:
            break
        layer_sizes.append(len(remove))
        active -= remove
    return len(active), len(layer_sizes), tuple(layer_sizes)


def nz_flow_count_mod(weights: tuple[int, ...], k: int) -> int:
    dist = {0: 1}
    for w in weights:
        nd: defaultdict[int, int] = defaultdict(int)
        for residue, count in dist.items():
            for m in range(1, k):
                nd[(residue + m * w) % k] += count
        dist = dict(nd)
    return dist.get(0, 0)


def branch_radius(speeds: tuple[int, ...], n: int, q: int) -> Fraction:
    zero_speeds = [s for s in speeds if s % q == 0]
    if not zero_speeds:
        return ZERO
    return Fraction(1, n * min(zero_speeds))


def branch_record(speeds: tuple[int, ...], n: int, q: int) -> dict[str, object]:
    zspeeds = tuple(s for s in speeds if s % q == 0)
    intervals = branch_intervals(speeds, n, q)
    core, rounds, layers = endpoint_core(intervals)
    radius = branch_radius(speeds, n, q)
    return {
        "q": q,
        "prime_power": is_prime_power(q),
        "phi": sum(1 for u in range(1, q) if gcd(u, q) == 1),
        "zero_count": len(zspeeds),
        "zero_speeds": zspeeds,
        "closed_witness_if_empty": len(zspeeds) == 0 and q <= n,
        "open_witness_if_empty": len(zspeeds) == 0 and q < n,
        "local_intervals": len(intervals),
        "local_core": core,
        "peel_rounds": rounds,
        "peel_layers": layers,
        "radius": radius,
    }


def valuation_histogram(speeds: tuple[int, ...], k: int) -> tuple[tuple[str, int], ...]:
    f = factor(k)
    if len(f) != 1:
        return ()
    p, a = next(iter(f.items()))
    counts: Counter[str] = Counter()
    for s in speeds:
        residue = s % k
        v = p_adic_val(residue, p, cap=a)
        label = f"v{p}>={a}" if v >= a else f"v{p}={v}"
        counts[label] += 1
    order = [f"v{p}={i}" for i in range(a)] + [f"v{p}>={a}"]
    return tuple((label, counts[label]) for label in order if counts[label])


def trienerment_from_records(records: list[dict[str, object]]) -> tuple[tuple[int, ...], tuple[tuple[int, int], ...]]:
    """States over i<j. 0=tie, 1=i beats j, -1=j beats i."""
    scores = []
    for rec in records:
        radius = rec["radius"]
        radius_num = radius.numerator if isinstance(radius, Fraction) else 0
        radius_den = radius.denominator if isinstance(radius, Fraction) else 1
        # Covered branches with many zero carriers and wide local radius are
        # more dangerous; nonempty local core would dominate if ever present.
        scores.append((
            int(rec["zero_count"] > 0),
            int(rec["local_core"]),
            int(rec["zero_count"]),
            radius_num * 10**9 // radius_den,
            -int(rec["q"]),
        ))
    states = []
    for i, j in combinations(range(len(records)), 2):
        if scores[i] == scores[j]:
            states.append(0)
        elif scores[i] > scores[j]:
            states.append(1)
        else:
            states.append(-1)
    return tuple(states), tuple((i, scores[i][2]) for i in range(len(scores)))


def states_to_oriented_adj(states: tuple[int, ...], n: int) -> list[list[int]]:
    adj = [[0] * n for _ in range(n)]
    it = iter(states)
    for i in range(n):
        for j in range(i + 1, n):
            s = next(it)
            if s == 1 or s == 0:  # tie Hamiltonian path: lower q wins ties
                adj[i][j] = 1
            else:
                adj[j][i] = 1
    return adj


def strict_directed_triangles(states: tuple[int, ...], n: int) -> int:
    mat = [[0] * n for _ in range(n)]
    it = iter(states)
    for i in range(n):
        for j in range(i + 1, n):
            s = next(it)
            mat[i][j] = s
            mat[j][i] = -s
    total = 0
    for a, b, c in combinations(range(n), 3):
        if mat[a][b] == 0 or mat[a][c] == 0 or mat[b][c] == 0:
            continue
        if (mat[a][b] == 1 and mat[b][c] == 1 and mat[c][a] == 1) or (
            mat[a][c] == 1 and mat[c][b] == 1 and mat[b][a] == 1
        ):
            total += 1
    return total


def scc_sizes(adj: list[list[int]]) -> tuple[int, ...]:
    n = len(adj)
    seen = [False] * n
    order = []

    def dfs(v: int) -> None:
        seen[v] = True
        for w in range(n):
            if adj[v][w] and not seen[w]:
                dfs(w)
        order.append(v)

    for v in range(n):
        if not seen[v]:
            dfs(v)
    rev = [[adj[j][i] for j in range(n)] for i in range(n)]
    seen = [False] * n
    sizes = []
    for v in reversed(order):
        if seen[v]:
            continue
        stack = [v]
        seen[v] = True
        size = 0
        while stack:
            x = stack.pop()
            size += 1
            for y in range(n):
                if rev[x][y] and not seen[y]:
                    seen[y] = True
                    stack.append(y)
        sizes.append(size)
    return tuple(sorted(sizes, reverse=True))


def hamiltonian_paths(adj: list[list[int]]) -> int:
    n = len(adj)
    dp = [[0] * n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for mask in range(1 << n):
        for last in range(n):
            value = dp[mask][last]
            if not value:
                continue
            for nxt in range(n):
                if mask & (1 << nxt) == 0 and adj[last][nxt]:
                    dp[mask | (1 << nxt)][nxt] += value
    return sum(dp[-1])


def trienerment_fingerprint(records: list[dict[str, object]]) -> dict[str, object]:
    states, _ = trienerment_from_records(records)
    adj = states_to_oriented_adj(states, len(records))
    tie_count = sum(1 for s in states if s == 0)
    score_hist = tuple(sorted(Counter(sum(row) for row in adj).items()))
    return {
        "vertices": len(records),
        "ties": tie_count,
        "score_hist": score_hist,
        "strict_c3": strict_directed_triangles(states, len(records)),
        "sccs": scc_sizes(adj),
        "hp": hamiltonian_paths(adj) if len(records) <= 12 else None,
    }


def summarize_case(case: Case) -> dict[str, object]:
    qnodes = prime_powers_upto(case.n)
    if case.n not in qnodes:
        # Non-prime-power n-gate is diagnostically important for n=18, but
        # marked separately by prime_power=False.
        qnodes.append(case.n)
    qnodes = sorted(set(qnodes))
    records = [branch_record(case.speeds, case.n, q) for q in qnodes]
    full_core = endpoint_core(danger_intervals(case.speeds, case.n))
    k = nstar(case.n)
    nz_count = nz_flow_count_mod(case.speeds, k)
    return {
        "records": records,
        "full_core": full_core,
        "nstar": k,
        "nstar_factor": factor(k),
        "valuation_hist": valuation_histogram(case.speeds, k),
        "nz_flow_count": nz_count,
        "trienerment": trienerment_fingerprint(records),
    }


def format_fraction(x: Fraction) -> str:
    if x == 0:
        return "0"
    return f"{x.numerator}/{x.denominator}"


def main() -> None:
    print("Prime-power zero-branch cover cores -- codex S546b")
    print("=" * 78)
    print("A q-zero branch is the set of speeds divisible by q=p^d.")
    print("If it is empty, unit times u/q are THM-369 sieve witnesses.")
    print("If it is covered, the local branch intervals are tested for a")
    print("nonpeeling endpoint-protection core.")
    print()

    all_branch_cores = []
    all_full_cores = []
    for case in CASES:
        data = summarize_case(case)
        print(f"CASE {case.name}: n={case.n}, speeds={case.speeds}")
        print(f"  n*={data['nstar']} factor={data['nstar_factor']} valuation_hist={data['valuation_hist']}")
        print(f"  NZ Z_n* full-support flow count={data['nz_flow_count']}")
        fc, fr, flayers = data["full_core"]
        all_full_cores.append(fc)
        print(f"  full danger-cover endpoint core: core={fc}, peel_rounds={fr}, layers={flayers[:8]}")
        fp = data["trienerment"]
        print(f"  branch trienerment: vertices={fp['vertices']}, ties={fp['ties']}, "
              f"score_hist={fp['score_hist']}, strict_c3={fp['strict_c3']}, "
              f"sccs={fp['sccs']}, HP={fp['hp']}")
        print("  q-branches:")
        for rec in data["records"]:
            all_branch_cores.append(int(rec["local_core"]))
            qtag = "pp" if rec["prime_power"] else "gate"
            witness = "open-witness" if rec["open_witness_if_empty"] else (
                "closed-witness" if rec["closed_witness_if_empty"] else "covered"
            )
            print(f"    q={rec['q']:>2} [{qtag}] zero={rec['zero_count']:>2} "
                  f"phi={rec['phi']:>2} status={witness:>14} "
                  f"radius={format_fraction(rec['radius'])}"
                  f" local_intervals={rec['local_intervals']:>3} "
                  f"core={rec['local_core']} peel_rounds={rec['peel_rounds']}")
        print()

    print("SYNTHESIS")
    print("-" * 78)
    print(f"Audited local branch cores: {len(all_branch_cores)}, "
          f"nonempty={sum(1 for x in all_branch_cores if x)}")
    print(f"Audited full cover cores: {len(all_full_cores)}, "
          f"nonempty={sum(1 for x in all_full_cores if x)}")
    print("Prime-power zero branches do exactly one thing locally: they kill")
    print("the unit sieve witnesses u/q.  Their local danger intervals are")
    print("nested around those unit points, so endpoint-core peeling deletes them.")
    print("Thus the branch does not store a counterexample core at q itself;")
    print("it exports endpoint debt to descendant wall/event layers.")
    print("For n=18 this matches the p-adic story: n*=9 makes zero-flow")
    print("branches abundant, but the audited cover cores remain empty.")


if __name__ == "__main__":
    main()
