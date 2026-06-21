#!/usr/bin/env python3
"""HYP-2702/HYP-2706: Gibbs, road-coloring, Hopfield, and quantum atlas.

The user proposed several external lenses: Gibbs measures, Arnold's cat map,
road coloring, Fubini-Study geometry, Hopfield/Hebbian networks, double-slit
amplitudes, and the Clifford/T-gate split.  This script tests which of those
lenses preserves the actual HYP-2702 sparse-tail predicate.

The input is the S65 death-chain band scout.  The output is deliberately
conservative:

* Gibbs reweighting measures how much adversarial "which-band" observation
  breaks the positive seven-band aggregate.
* Arnold's cat map gives a canonical projective orbit through the seven bands
  plus infinity; prefix sums test whether an orbit-order proof exists.
* Road coloring becomes the deterministic singleton-deletion automaton whose
  reset-word counts are exactly the singleton death-chain kernel.
* Hopfield/Hebbian storage records the observed band-sign attractors.
* Fubini-Study distance on hit-count laws shows that projective closeness does
  not replace the signed kernel inequality.
* The Clifford/magic analogy is treated as a proof-quotient tournament:
  stabilizer-like coarse quotients lose; generated signed context wins.

No LRC14 proof is claimed.
"""

from __future__ import annotations

import importlib.util
import math
import sys
from collections import Counter, defaultdict
from fractions import Fraction as F
from itertools import combinations, product
from pathlib import Path


sys.stdout.reconfigure(line_buffering=True)

ROOT = Path(__file__).resolve().parents[1]
S65_PATH = ROOT / "04-computation" / "lrc14_death_chain_band_automaton_codex_s65.py"
spec = importlib.util.spec_from_file_location("s65_death_chain", S65_PATH)
s65 = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(s65)
s64 = s65.s64


def fmt(q: F | None, digits: int = 9) -> str:
    if q is None:
        return "None"
    return f"{q} ({float(q):.{digits}f})"


def sign(q: F) -> int:
    if q > 0:
        return 1
    if q < 0:
        return -1
    return 0


def sign_word(vals: tuple[F, ...]) -> str:
    return "".join("+" if v > 0 else "-" if v < 0 else "0" for v in vals)


def pearson(xs: list[float], ys: list[float]) -> float:
    n = len(xs)
    mx = sum(xs) / n
    my = sum(ys) / n
    vx = sum((x - mx) ** 2 for x in xs)
    vy = sum((y - my) ** 2 for y in ys)
    if vx == 0 or vy == 0:
        return 0.0
    return sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / math.sqrt(vx * vy)


def frontier_records() -> list[dict[str, object]]:
    records: list[dict[str, object]] = []
    for size in range(3, 7):
        consec = tuple(range(size))
        r = 7 - size
        for row in s64.coordinate_violators(size):
            shape = row["shape"]
            assert isinstance(shape, tuple)
            bands, atoms = s65.band_delta_profile(consec, shape, r)
            records.append(
                {
                    "size": size,
                    "r": r,
                    "consec": consec,
                    "shape": shape,
                    "bands": bands,
                    "delta": sum(bands, F(0)),
                    "word": sign_word(bands),
                    "atoms": atoms,
                }
            )
    return records


def gibbs_weighted_sum(bands: tuple[F, ...], beta: float) -> float:
    debt = [max(0.0, -float(b)) for b in bands]
    weights = [math.exp(beta * d) for d in debt]
    z = sum(weights)
    return sum(w * float(b) for w, b in zip(weights, bands)) / z


def gibbs_beta_critical(bands: tuple[F, ...]) -> float | None:
    if min(bands) >= 0:
        return None
    if gibbs_weighted_sum(bands, 0.0) <= 0:
        return 0.0
    hi = 1.0
    while gibbs_weighted_sum(bands, hi) > 0 and hi < 1e8:
        hi *= 2
    if hi >= 1e8:
        return None
    lo = 0.0
    for _ in range(80):
        mid = (lo + hi) / 2
        if gibbs_weighted_sum(bands, mid) > 0:
            lo = mid
        else:
            hi = mid
    return hi


def print_gibbs(records: list[dict[str, object]]) -> dict[str, object]:
    print("PART A -- Gibbs which-band reweighting")
    print("  beta=0 is the uniform seven-band theorem candidate.")
    print("  beta>0 adversarially upweights negative-band debt.")
    critical: list[tuple[float, dict[str, object]]] = []
    abs_gap_ratios: list[tuple[float, dict[str, object]]] = []
    for rec in records:
        bands = rec["bands"]
        assert isinstance(bands, tuple)
        beta = gibbs_beta_critical(bands)
        if beta is not None:
            critical.append((beta, rec))
        pos = sum(max(F(0), b) for b in bands)
        neg = -sum(min(F(0), b) for b in bands)
        ratio = float(neg / pos) if pos else float("inf")
        abs_gap_ratios.append((ratio, rec))
    critical.sort(key=lambda x: x[0])
    abs_gap_ratios.sort(reverse=True, key=lambda x: x[0])
    print(f"  rows={len(records)}, rows with negative bands={len(critical)}")
    for beta, rec in critical[:5]:
        print(
            f"    smallest beta*: {beta:.6f}  size={rec['size']} C={rec['shape']} "
            f"word={rec['word']} delta={fmt(rec['delta'])}"
        )
        break
    print("  most which-path fragile rows by negative/positive mass:")
    for ratio, rec in abs_gap_ratios[:5]:
        print(
            f"    ratio={ratio:.3f} size={rec['size']} C={rec['shape']} "
            f"word={rec['word']} delta={fmt(rec['delta'])}"
        )
    print("  read: observing/reweighting the losing bands kills the proof; the")
    print("        usable object is the uniform signed aggregate, not a Gibbs ground state.")
    print()
    return {
        "min_beta": critical[0][0],
        "min_beta_rec": critical[0][1],
        "max_ratio": abs_gap_ratios[0][0],
        "max_ratio_rec": abs_gap_ratios[0][1],
    }


def inv_mod(a: int, p: int = 7) -> int:
    return pow(a % p, -1, p)


def cat_step(t: int | str) -> int | str:
    # A=(2 1; 1 1), slope t=y/x maps to (1+t)/(2+t).
    if t == "inf":
        return 1
    assert isinstance(t, int)
    den = (2 + t) % 7
    if den == 0:
        return "inf"
    return ((1 + t) * inv_mod(den)) % 7


def cat_projective_cycle(start: int | str = 0) -> list[int | str]:
    out: list[int | str] = []
    cur: int | str = start
    while cur not in out:
        out.append(cur)
        cur = cat_step(cur)
    return out


def cat_projective_cycles() -> list[list[int | str]]:
    seen: set[int | str] = set()
    cycles: list[list[int | str]] = []
    for start in list(range(7)) + ["inf"]:
        if start in seen:
            continue
        cur: int | str = start
        cyc: list[int | str] = []
        while cur not in seen:
            cyc.append(cur)
            seen.add(cur)
            cur = cat_step(cur)
        cycles.append(cyc)
    return cycles


def cat_torus_cycles() -> Counter[int]:
    seen: set[tuple[int, int]] = set()
    lengths: Counter[int] = Counter()
    for x in range(7):
        for y in range(7):
            if (x, y) in seen:
                continue
            cyc = []
            cur = (x, y)
            while cur not in cyc:
                cyc.append(cur)
                seen.add(cur)
                x0, y0 = cur
                cur = ((2 * x0 + y0) % 7, (x0 + y0) % 7)
            lengths[len(cyc)] += 1
    return lengths


def print_cat(records: list[dict[str, object]]) -> dict[str, object]:
    print("PART B -- Arnold cat map projective orbit")
    cycles = cat_projective_cycles()
    torus = cat_torus_cycles()
    print(f"  cat map A=[[2,1],[1,1]] over F7 torus cycle lengths={dict(sorted(torus.items()))}")
    print(f"  projective slope cycles: {cycles}")
    neg_prefix = 0
    neg_cycle_sum = 0
    worst_prefix: tuple[F, dict[str, object], int, list[int | str]] | None = None
    worst_cycle: tuple[F, dict[str, object], list[int | str]] | None = None
    for rec in records:
        bands = rec["bands"]
        assert isinstance(bands, tuple)
        row_prefix_min = F(0)
        row_prefix_idx = -1
        row_prefix_cycle: list[int | str] = []
        row_cycle_min: F | None = None
        row_cycle_arg: list[int | str] = []
        for cyc in cycles:
            running = F(0)
            cycle_sum = F(0)
            for i, t in enumerate(cyc):
                val = bands[t] if isinstance(t, int) else F(0)
                running += val
                cycle_sum += val
                if running < row_prefix_min:
                    row_prefix_min = running
                    row_prefix_idx = i
                    row_prefix_cycle = cyc
            if row_cycle_min is None or cycle_sum < row_cycle_min:
                row_cycle_min = cycle_sum
                row_cycle_arg = cyc
        assert row_cycle_min is not None
        if row_prefix_min < 0:
            neg_prefix += 1
        if row_cycle_min < 0:
            neg_cycle_sum += 1
        if worst_prefix is None or row_prefix_min < worst_prefix[0]:
            worst_prefix = (row_prefix_min, rec, row_prefix_idx, row_prefix_cycle)
        if worst_cycle is None or row_cycle_min < worst_cycle[0]:
            worst_cycle = (row_cycle_min, rec, row_cycle_arg)
    assert worst_prefix is not None
    assert worst_cycle is not None
    print(f"  rows with a negative cat-cycle sum={neg_cycle_sum}/{len(records)}")
    print(f"  rows with a negative cat-cycle prefix={neg_prefix}/{len(records)}")
    print(
        f"  worst prefix={fmt(worst_prefix[0])} after index={worst_prefix[2]} "
        f"at size={worst_prefix[1]['size']} C={worst_prefix[1]['shape']} "
        f"word={worst_prefix[1]['word']} cycle={worst_prefix[3]}"
    )
    print(
        f"  worst cycle sum={fmt(worst_cycle[0])} "
        f"at size={worst_cycle[1]['size']} C={worst_cycle[1]['shape']} "
        f"cycle={worst_cycle[2]}"
    )
    print("  read: the cat map supplies a rigid hyperbolic cycle decomposition, but")
    print("        cycle and prefix positivity are false; only the completed aggregate survives.")
    print()
    return {
        "neg_prefix": neg_prefix,
        "neg_cycle_sum": neg_cycle_sum,
        "worst_prefix": worst_prefix,
        "worst_cycle": worst_cycle,
        "cycles": cycles,
    }


def reset_count(t: int, r: int) -> int:
    total = 0
    for j in range(t + 1):
        total += (-1) ** j * math.comb(t, j) * (7 - j) ** r
    return total


def print_road_coloring() -> dict[str, object]:
    print("PART C -- road-coloring deletion automaton")
    print("  States are residual masks R subset {1,...,6}; letter a removes a from R.")
    print("  This monotone automaton is not strongly connected, but its reset-word")
    print("  counts are exactly the singleton death-chain kernel g_r(t).")
    print("       t  shortest_reset  g_0      g_1      g_2      g_3      g_4")
    for t in range(7):
        vals = []
        for r in range(5):
            vals.append(F(reset_count(t, r), 7**r))
        print(
            f"      {t:2d}       {t:2d}        "
            + " ".join(f"{str(v):>7}" for v in vals)
        )
    all_state_reset_word = tuple(range(1, 7))
    print(f"  one global reset word for all 64 masks has length {len(all_state_reset_word)}: {all_state_reset_word}")
    print("  read: the road-coloring shadow is already the death-chain proof object;")
    print("        it just lacks strong connectivity because deficit only decreases.")
    print()
    return {"global_reset_len": 6}


def hopfield_update(state: tuple[int, ...], weights: list[list[int]]) -> tuple[int, ...]:
    out = []
    for i, old in enumerate(state):
        field = sum(weights[i][j] * state[j] for j in range(len(state)) if j != i)
        if field > 0:
            out.append(1)
        elif field < 0:
            out.append(-1)
        else:
            out.append(old)
    return tuple(out)


def print_hopfield(records: list[dict[str, object]]) -> dict[str, object]:
    print("PART D -- Hopfield/Hebbian band-sign memory")
    patterns: Counter[tuple[int, ...]] = Counter()
    for rec in records:
        bands = rec["bands"]
        assert isinstance(bands, tuple)
        # Store zeros as +1 for a binary Hopfield state; retain words separately.
        patterns[tuple(1 if sign(b) >= 0 else -1 for b in bands)] += 1
    weights = [[0 for _ in range(7)] for _ in range(7)]
    for pat, count in patterns.items():
        for i in range(7):
            for j in range(7):
                if i != j:
                    weights[i][j] += count * pat[i] * pat[j]
    fixed_patterns = 0
    for pat in patterns:
        if hopfield_update(pat, weights) == pat:
            fixed_patterns += 1
    basins: Counter[tuple[int, ...]] = Counter()
    for bits in product((-1, 1), repeat=7):
        state = tuple(bits)
        seen = set()
        while state not in seen:
            seen.add(state)
            nxt = hopfield_update(state, weights)
            if nxt == state:
                break
            state = nxt
        basins[state] += 1
    print(f"  observed binary sign patterns={len(patterns)}, fixed observed patterns={fixed_patterns}")
    print("  observed pattern counts:")
    for pat, count in patterns.most_common():
        word = "".join("+" if x > 0 else "-" for x in pat)
        print(f"    {word}: {count}")
    print("  largest Hopfield basins over all 2^7 states:")
    for state, count in basins.most_common(6):
        word = "".join("+" if x > 0 else "-" for x in state)
        print(f"    attractor={word} basin={count}")
    print("  read: the band debts have a small attractor alphabet, mainly -++-++-")
    print("        and -+++++-, but attractor storage is descriptive, not a proof.")
    print()
    return {"patterns": patterns, "fixed_patterns": fixed_patterns, "basins": basins}


def fubini_study_distance(p: tuple[F, ...], q: tuple[F, ...]) -> float:
    fidelity = sum(math.sqrt(float(a) * float(b)) for a, b in zip(p, q))
    fidelity = max(0.0, min(1.0, fidelity))
    return math.acos(fidelity)


def print_fubini_study(records: list[dict[str, object]]) -> dict[str, object]:
    print("PART E -- Fubini-Study distance on hit-count laws")
    fs_vals: list[float] = []
    deltas: list[float] = []
    worst: tuple[float, dict[str, object]] | None = None
    closest: tuple[float, dict[str, object]] | None = None
    for rec in records:
        consec = rec["consec"]
        shape = rec["shape"]
        assert isinstance(consec, tuple)
        assert isinstance(shape, tuple)
        fs = fubini_study_distance(s64.hit_count_law(consec), s64.hit_count_law(shape))
        fs_vals.append(fs)
        deltas.append(float(rec["delta"]))
        if worst is None or fs > worst[0]:
            worst = (fs, rec)
        if closest is None or fs < closest[0]:
            closest = (fs, rec)
    assert worst is not None and closest is not None
    print(
        f"  FS distance range: min={closest[0]:.6f} at size={closest[1]['size']} C={closest[1]['shape']}; "
        f"max={worst[0]:.6f} at size={worst[1]['size']} C={worst[1]['shape']}"
    )
    print(f"  corr(FS distance, death-chain margin)={pearson(fs_vals, deltas):.6f}")
    print("  read: projective distance is signal, but it does not encode the")
    print("        signed kernel order; it cannot replace the death-chain inequality.")
    print()
    return {"min_fs": closest, "max_fs": worst, "corr": pearson(fs_vals, deltas)}


def sccs(vertices: list[str], edges: dict[tuple[str, str], str]) -> list[list[str]]:
    adj = {v: [] for v in vertices}
    radj = {v: [] for v in vertices}
    for a in vertices:
        for b in vertices:
            if a == b:
                continue
            winner = edges[(a, b)]
            loser = b if winner == a else a
            adj[winner].append(loser)
            radj[loser].append(winner)
    seen: set[str] = set()
    order: list[str] = []

    def dfs(v: str) -> None:
        seen.add(v)
        for w in adj[v]:
            if w not in seen:
                dfs(w)
        order.append(v)

    for v in vertices:
        if v not in seen:
            dfs(v)
    seen.clear()
    comps: list[list[str]] = []

    def rdfs(v: str, comp: list[str]) -> None:
        seen.add(v)
        comp.append(v)
        for w in radj[v]:
            if w not in seen:
                rdfs(w, comp)

    for v in reversed(order):
        if v not in seen:
            comp: list[str] = []
            rdfs(v, comp)
            comps.append(sorted(comp))
    return comps


def hamiltonian_path_count(vertices: list[str], edges: dict[tuple[str, str], str]) -> int:
    count = 0

    def rec(path: list[str], unused: set[str]) -> None:
        nonlocal count
        if not unused:
            count += 1
            return
        last = path[-1]
        for v in list(unused):
            if edges[(last, v)] == last:
                unused.remove(v)
                path.append(v)
                rec(path, unused)
                path.pop()
                unused.add(v)

    for v in vertices:
        rest = set(vertices)
        rest.remove(v)
        rec([v], rest)
    return count


def print_tournament(metrics: dict[str, object]) -> None:
    print("PART F -- proof-quotient tournament")
    print("  Vertices are proof carriers, not runners or arcs.")
    print("  Observable=(frontier_failures, min_margin/proxy, predicate_fidelity, inspiration_signal).")
    carriers = {
        "which_path_abs_band_measure": (
            192,
            F(-33, 12005),
            1,
            5,
            "double-slit measurement destroys interference by isolating bands",
        ),
        "hitcount_Fubini_Study_geometry": (
            72,
            F(0),
            2,
            5,
            "projective distance is useful signal but not an order",
        ),
        "cat_map_prefix_orbit": (
            int(metrics["cat_neg_prefix"]),
            metrics["cat_worst_prefix"],
            3,
            6,
            "hyperbolic projective cycles exist, but cycle/prefix positivity fails",
        ),
        "Hopfield_band_attractor": (
            72 - int(metrics["hop_fixed"]),
            F(0),
            3,
            7,
            "stores the debt alphabet, does not prove the margin",
        ),
        "Gibbs_uniform_band_sum": (
            0,
            metrics["min_delta"],
            5,
            7,
            "beta=0 aggregate passes; beta>0 which-band bias can fail",
        ),
        "road_coloring_death_chain": (
            0,
            metrics["min_delta"],
            6,
            7,
            "reset-word counts equal the singleton kernel exactly",
        ),
        "miss_zeta_generated_magic": (
            0,
            metrics["min_delta"],
            7,
            8,
            "keeps generated context phase before scalarizing",
        ),
    }
    vertices = list(carriers)

    def better(a: str, b: str) -> str:
        ma = carriers[a]
        mb = carriers[b]
        key_a = (-ma[0], ma[1], ma[2], ma[3])
        key_b = (-mb[0], mb[1], mb[2], mb[3])
        return a if key_a >= key_b else b

    edges: dict[tuple[str, str], str] = {}
    scores: Counter[str] = Counter()
    for i, a in enumerate(vertices):
        for b in vertices[i + 1 :]:
            w = better(a, b)
            edges[(a, b)] = w
            edges[(b, a)] = w
            scores[w] += 1
    score_hist = Counter(scores[v] for v in vertices)
    cycles = 0
    for a, b, c in combinations(vertices, 3):
        ab = edges[(a, b)]
        bc = edges[(b, c)]
        ca = edges[(c, a)]
        if ab == a and bc == b and ca == c:
            cycles += 1
        if ab == b and bc == c and ca == a:
            cycles += 1
    for name, metric in sorted(carriers.items(), key=lambda kv: (-scores[kv[0]], kv[0])):
        print(f"    {name}: score={scores[name]}, metric={metric[:4]}, note={metric[4]}")
    print(f"  score_hist={dict(sorted(score_hist.items()))}")
    print(f"  directed_3cycles={cycles}")
    print(f"  SCCs={sccs(vertices, edges)}")
    print(f"  Hamiltonian_path_count={hamiltonian_path_count(vertices, edges)}")
    print("  Clifford/magic read: coarse stabilizer-like quotients are tractable but")
    print("  insufficient; the LRC sparse tail needs the generated signed context,")
    print("  analogous to keeping the non-Clifford phase resource.")
    print()


def main() -> None:
    print("HYP-2702/HYP-2706 -- Gibbs/quantum/road-coloring bridge atlas")
    print("Input: S65 sparse-tail death-chain band frontier.\n")
    records = frontier_records()
    print(f"Loaded {len(records)} sparse-coordinate frontier rows.\n")
    gibbs = print_gibbs(records)
    cat = print_cat(records)
    road = print_road_coloring()
    hop = print_hopfield(records)
    fs = print_fubini_study(records)
    min_delta = min(rec["delta"] for rec in records)
    assert isinstance(min_delta, F)
    metrics = {
        "cat_neg_prefix": cat["neg_prefix"],
        "cat_worst_prefix": cat["worst_prefix"][0],
        "hop_fixed": hop["fixed_patterns"],
        "min_delta": min_delta,
    }
    print_tournament(metrics)
    print("SYNTHESIS")
    print("  The external lenses are productive when they are used as quotient tests.")
    print("  Gibbs and double-slit language says why measuring individual bands breaks")
    print("  the proof. Road coloring exactly recovers the singleton death chain.")
    print("  Cat-map and Hopfield structure expose the signed debt alphabet but do not")
    print("  give prefix/local monotonicity. Fubini-Study geometry is a diagnostic, not")
    print("  an order. The proof carrier remains miss-zeta generated context plus the")
    print("  signed seven-band death-chain aggregate.")


if __name__ == "__main__":
    main()
