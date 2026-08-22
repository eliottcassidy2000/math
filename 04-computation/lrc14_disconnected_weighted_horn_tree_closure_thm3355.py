#!/usr/bin/env python3
"""Exact review certificate for the disconnected weighted horn-tree route.

This companion verifies the analytic and finite-exact gates used by
THM-3355.  It freezes:

* the primitive centered-grid lower bound and its nonprimitive hostile;
* coefficient-positive tail polynomials for every nonhorn g=1 lane and every
  g=2,3 lane from P=421 onward;
* one internally reconstructed eight-shard exact scan for P<=420;
* independent fast/reference/literal controls and 100,000 seeded bound tests;
* the complete-multipartite graph census and its unique K_1,5 exception; and
* the exact weighted debt inequalities which discharge that exception.

All arithmetic affecting the verdict is integral or Fraction-exact.  No
``assert`` statement is used, so ordinary and ``python -O`` executions have
identical semantics and byte-identical output.
"""

from __future__ import annotations

from argparse import ArgumentParser
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction as F
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, permutations
from math import gcd, prod
from pathlib import Path
from random import Random
import ast
import os


ROOT = Path(__file__).resolve().parents[1]
THM1250 = ROOT / "01-canon/theorems/THM-1250-six-private-needles-force-fully-located-spanning-tree.md"
THM3350 = ROOT / "01-canon/theorems/THM-3350-connected-low-full-tree-atlas-dense-closure-and-uniform-tail.md"
THM3352 = ROOT / "01-canon/theorems/THM-3352-connected-low-all-head-universal-physical-forest-closure.md"
TAIL = ROOT / "04-computation/lrc14_connected_low_uniform_high_forest_tail_thm3350.py"
FAST_MASS = ROOT / "04-computation/lrc_general_reflected_pair_mass_thm3352.py"
REFERENCE_MASS = ROOT / "04-computation/lrc_general_reflected_pair_mass_reference_audit_thm3352.py"
GEOMETRY = ROOT / "04-computation/lrc14_disconnected_low_geometry_verify_20260812.py"
GEOMETRY_OUTPUT = ROOT / "05-knowledge/results/lrc14_disconnected_low_geometry_verify_20260812.out"
LARGE_RULER = ROOT / "04-computation/lrc14_disconnected_large_ruler_floor_20260812.py"
LARGE_RULER_OUTPUT = ROOT / "05-knowledge/results/lrc14_disconnected_large_ruler_floor_20260812.out"
G4_TAIL = ROOT / "04-computation/lrc14_disconnected_non35_g4_tail_20260812.py"
G4_TAIL_OUTPUT = ROOT / "05-knowledge/results/lrc14_disconnected_non35_g4_tail_20260812.out"
THREE_FIVE = ROOT / "04-computation/lrc14_disconnected_35_small_ruler_symbolic_20260812.py"
THREE_FIVE_OUTPUT = ROOT / "05-knowledge/results/lrc14_disconnected_35_small_ruler_symbolic_20260812.out"

EXPECTED_DEPENDENCIES = {
    THM1250: "4fe481ff66a5e4d5bb3b2f8b3c64a845b8ec70d887acfe93813918f06466258f",
    THM3350: "b4ddd4670c5d2f0ff3e56057d31f50e322996225f4d67b6e9cff78479d084170",
    THM3352: "b9f0a7245d3bc922f33046a4dcf0c67a922cbe492739e66e870c3c75af31752d",
    TAIL: "78daaf73966d283c0c0bafa1c0975684e6167d2ef6375a3abeece4e00cdc87f9",
    FAST_MASS: "afd417297131401254769e1ef172d89c109ad2f9a843ea55e2badc3e7891435b",
    REFERENCE_MASS: "da941a4267147d5442be81ae81880742d2f6b901bfc1d20fb667822402a2950e",
    GEOMETRY: "871d08736ab6907c7b95c0e1c9bbb62467200309039c961a82eb44bfec11fada",
    GEOMETRY_OUTPUT: "949e1269e101921bebf623232f1e0e11576d83f0b6baa2140004a2138cd0a53c",
    LARGE_RULER: "29138d28a83e83b8e5f9d5077e71cefe7402c0514cdac605ddf8f401c17c55c9",
    LARGE_RULER_OUTPUT: "54fb8d62469dc36d4362ec6850126859d4fccce4fce75dfc02fbdab3d1cd9b94",
    G4_TAIL: "d8ab54707f939e4e7d2bc3029e5b355b32c923a1ab1bc3329f5d8c5adadc2373",
    G4_TAIL_OUTPUT: "5506170d5dc8c63ace1a90ec28cc08b53ebbe31f493ead1d7a2dc0a18db45739",
    THREE_FIVE: "73cb4f0b76c4101ba60bd6b0cf4043c88c70864d90d2df3193f23d6e850aaa84",
    THREE_FIVE_OUTPUT: "89f6b255a0164e23ce856d7c9a50bf4febbeba30d08d343f40ce43dc5a2ce1bf",
}

TARGET = F(1, 294)
DEBT_MAX = F(186636088362, 11773143757375)
EXPECTED_CONTEXT = "efea6bd97522fe1c31a5a88ca9f3223f9e7a8c08e3be85c493e9f62fdfaf06e4"
EXPECTED_RANDOM = "5127b6eaf89d9a095b4c93d18a886769ee98862f3902a084c88ef4dba2208cb2"
EXPECTED_GRAPH = "f14de38dddebdafdb557df6a2afaa717c11d4367215ee921e474d74e94170564"
EXPECTED_TREE = "dff2240ec7410b545d9bc8cd22b7a7bd0b757d6501ea5cad8e281df012c7afac"
EXPECTED_DEBT_CENSUS = "aa67694894dfa788bdb19f4b00cc2d90d422aa733452c8163be98ef5f6a78b7e"
EXPECTED_POLYNOMIAL = "742222ee57bbc518eee39cfb293522736bd8f9fc997d40686eea5936f2883b5e"
EXPECTED_CONTROL = "2b77f5d05147be0b382c827b2516e9ce54c12c469f46c98618601c9c6f47a8e4"
EXPECTED_FINITE_AGGREGATE = "854238c36f2d6bdea17bfeef81baebc101dc0bd935fd2c6f1a71e8e98897dffa"

EXPECTED_SHARD_CHECKS = (530316, 936064, 897440, 898756, 1265122, 557436, 531802, 527308)
EXPECTED_SHARD_GATES = (384253, 383767, 382715, 382726, 382420, 383036, 383075, 383054)
EXPECTED_SHARD_DIGESTS = (
    "812f93d511f71ca65c1d7c94eebca238812b53818954431b3ef8fda003823534",
    "55c4c4acca42fe3f4bf5015b9c7c3aa6ad47585b504d331046302ff2f585c1e7",
    "af2f377772389ad51bb24905629f0467653415d24f0cfcfb863bdbe9e6f1cc06",
    "59331f8fee290d3a593c4e980a3e1ad00206f7d414f0477a5871a5ddddf71baf",
    "6a120e39cb94fba0e152d2c0f5565a7ca8c8346f4d1c3d0768f1e788b2618de2",
    "559a5511c9b5fd3210604f58cefaca35514d1c53bb190a4b5e5544ca494edaf0",
    "0096248ea4dc541ddb4776e95b277d21a1dc86ac0d1d1bb29c83014972581d34",
    "b8a604b36e87ad5179e9af996a52695313b1f9aa2431ebed0672070457231060",
)
EXPECTED_SHARD_BEST = (
    (F(1680, 129359), 3920, 1995, 10, 1, 2, 11, 3),
    (F(1872, 144143), 4368, 2236, 8, 1, 2, 11, 3),
    (F(1008, 77615), 2352, 1204, 8, 1, 2, 11, 3),
    (F(1584, 121967), 3696, 1892, 8, 1, 2, 11, 3),
    (F(1872, 144143), 4368, 2262, 13, 1, 2, 11, 3),
    (F(92, 7645), 168, 90, 12, 6, 4, 5, 1),
    (F(6003, 491359), 588, 315, 1, 14, 4, 5, 1),
    (F(468, 36035), 4368, 2262, 3, 4, 2, 11, 3),
)


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def normalized_hash(path: Path) -> str:
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def load(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, ("module spec", path))
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def freeze(actual: str, expected: str, name: str):
    if expected != "PENDING":
        require(actual == expected, (name, actual, expected))


def dependency_gate():
    rows = []
    for path, expected in EXPECTED_DEPENDENCIES.items():
        actual = normalized_hash(path)
        require(actual == expected, ("dependency hash", path, actual, expected))
        rows.append((path.relative_to(ROOT).as_posix(), actual))
    tree = ast.parse(Path(__file__).read_text(encoding="utf-8"), filename=__file__)
    assert_nodes = sum(isinstance(node, ast.Assert) for node in ast.walk(tree))
    require(assert_nodes == 0, ("optimization-sensitive asserts", assert_nodes))
    return tuple(rows), assert_nodes


_TAIL_MODULE = None
_FAST_MODULE = None
_REFERENCE_MODULE = None


def tail_module():
    global _TAIL_MODULE
    if _TAIL_MODULE is None:
        require(normalized_hash(TAIL) == EXPECTED_DEPENDENCIES[TAIL], "tail dependency drift")
        _TAIL_MODULE = load("weighted_horn_tail_thm3350", TAIL)
    return _TAIL_MODULE


def fast_module():
    global _FAST_MODULE
    if _FAST_MODULE is None:
        require(normalized_hash(FAST_MASS) == EXPECTED_DEPENDENCIES[FAST_MASS], "fast mass drift")
        _FAST_MODULE = load("weighted_horn_fast_mass", FAST_MASS)
    return _FAST_MODULE


def reference_module():
    global _REFERENCE_MODULE
    if _REFERENCE_MODULE is None:
        require(normalized_hash(REFERENCE_MASS) == EXPECTED_DEPENDENCIES[REFERENCE_MASS],
                "reference mass drift")
        _REFERENCE_MODULE = load("weighted_horn_reference_mass", REFERENCE_MASS)
    return _REFERENCE_MODULE


def context_atlas():
    tail = tail_module()
    body_rows = tuple(sorted((tuple(body), L) for body, L in tail.SEL.MS.body_universe()))
    require(len(body_rows) == 649, ("body count", len(body_rows)))
    rows = set()
    for body, L in body_rows:
        cell, *_ = tail.SEL.body_geometry(body, L)
        if L < 4592:
            for e in body:
                for f in body:
                    if e != f:
                        rows.add((L, cell, e, f))
    rows = tuple(sorted(rows))
    payload = "".join(f"{L} {cell} {e} {f}\n" for L, cell, e, f in rows).encode()
    require(len(rows) == 2530, ("context count", len(rows)))
    require(len({row[0] for row in rows}) == 29, "small-ruler count")
    require(sha256(payload).hexdigest() == EXPECTED_CONTEXT,
            ("context digest", sha256(payload).hexdigest()))
    lanes = tuple(sorted({(L, e, f) for L, _, e, f in rows}))
    horns = tuple(row for row in lanes if row[0] == 168 and row[1] == 12 and row[2] in (1, 2, 3, 4))
    require(len(lanes) == 1304 and horns == tuple((168, 12, f) for f in (1, 2, 3, 4)),
            ("lane/horn census", len(lanes), horns))
    return body_rows, rows, lanes, horns


def centered_bound(L: int, e: int, f: int, P: int, Q: int, g: int) -> F:
    """Lower bound at primitive pairs; elsewhere only its algebraic envelope."""
    z, w = L * g * P - e, L * g * Q - f
    require(P >= 1 and Q > P and g >= 1 and z > 0 and w > 0,
            ("centered-bound domain", L, e, f, P, Q, g))
    return (F(g * L * P, 49 * z)
            - F(2 * g * L, 7 * w)
            - F(g * L * abs(Q * e - P * f) * (P * P // 4), P * w * z))


def endpoint_bound(L: int, e: int, f: int, P: int, g: int) -> F:
    # On each side of Q=Pf/e the loss is fractional-linear, hence monotone.
    # If the breakpoint lies in [P+1,8P], its determinant loss is zero and its
    # first loss is no larger than at P+1.  Thus only these two algebraic
    # endpoints matter.  In particular Q=8P need not itself be primitive.
    return min(centered_bound(L, e, f, P, P + 1, g),
               centered_bound(L, e, f, P, 8 * P, g))


def p_add(*rows):
    size = max(map(len, rows))
    answer = [0] * size
    for row in rows:
        for index, value in enumerate(row):
            answer[index] += value
    while len(answer) > 1 and answer[-1] == 0:
        answer.pop()
    return tuple(answer)


def p_scale(row, scalar):
    return tuple(scalar * value for value in row)


def p_mul(first, second):
    answer = [0] * (len(first) + len(second) - 1)
    for i, a in enumerate(first):
        for j, b in enumerate(second):
            answer[i + j] += a * b
    return tuple(answer)


def p_derivative(row):
    return tuple(index * value for index, value in enumerate(row))[1:] or (0,)


def polynomial_record(g: int, L: int, e: int, f: int, parity: str, side: str):
    base = 421 if parity == "odd" else 422
    P = (base, 2)
    Q = p_add(P, (1,)) if side == "left" else p_scale(P, 8)
    z = p_add(p_scale(P, L * g), (-e,))
    w = p_add(p_scale(Q, L * g), (-f,))
    if side == "left":
        raw = p_add(p_scale(P, e - f), (e,))
        determinant = raw if e > f else p_scale(raw, -1)
    else:
        determinant = p_scale(P, abs(8 * e - f))
    square = p_mul(P, P)
    quarter = p_add(square, (-1,)) if parity == "odd" else square
    require(all(value % 4 == 0 for value in quarter), ("parity polynomial", parity, quarter))
    half_square = tuple(value // 4 for value in quarter)

    # 294*P*w*z*(B-1/294), a positive-denominator clearing:
    # 6gLP^2w - 84gLPz - 294gL|Qe-Pf|floor(P^2/4) - Pwz.
    numerator = p_add(
        p_scale(p_mul(p_mul(P, P), w), 6 * g * L),
        p_scale(p_mul(P, z), -84 * g * L),
        p_scale(p_mul(determinant, half_square), -294 * g * L),
        p_scale(p_mul(p_mul(P, w), z), -1),
    )
    denominator = p_scale(p_mul(p_mul(P, w), z), 294)
    require(all(value > 0 for value in denominator),
            ("positive polynomial denominator", g, L, e, f, parity, side, denominator))
    require(all(value > 0 for value in numerator),
            ("positive polynomial numerator", g, L, e, f, parity, side, numerator))
    monotonicity = p_add(
        p_mul(p_derivative(numerator), denominator),
        p_scale(p_mul(numerator, p_derivative(denominator)), -1),
    )
    require(all(value >= 0 for value in monotonicity),
            ("tail-margin monotonicity", g, L, e, f, parity, side, monotonicity))
    content = 0
    for value in numerator:
        content = gcd(content, value)
    primitive_numerator = tuple(value // content for value in numerator)
    return (g, L, e, f, parity, side, primitive_numerator, numerator, denominator)


def polynomial_tail_audit(lanes, horns):
    horn_set = set(horns)
    groups = (("g1_nonhorn", (1,), tuple(row for row in lanes if row not in horn_set)),
              ("g23_all", (2, 3), lanes))
    records = []
    group_rows = []
    for name, scales, selected in groups:
        start = len(records)
        for g in scales:
            for L, e, f in selected:
                for parity in ("odd", "even"):
                    for side in ("left", "right"):
                        records.append(polynomial_record(g, L, e, f, parity, side))
        subset = records[start:]
        group_rows.append((name, len(subset), min(min(row[6]) for row in subset)))
    require(tuple(row[1] for row in group_rows) == (5200, 10432)
            and group_rows[0][2] == 980 and group_rows[1][2] > 0,
            ("polynomial groups", group_rows))
    horn_minimum = min(min(polynomial_record(g, L, e, f, parity, side)[6])
                       for g in (2, 3) for L, e, f in horns
                       for parity in ("odd", "even") for side in ("left", "right"))
    require(horn_minimum == 30576, ("horn polynomial minimum", horn_minimum))
    semantic = sha256()
    for row in records:
        semantic.update((repr(row) + "\n").encode())
    digest = semantic.hexdigest()
    freeze(digest, EXPECTED_POLYNOMIAL, "polynomial semantic")

    candidates = []
    for name, scales, selected in groups:
        for g in scales:
            for L, e, f in selected:
                for P in (421, 422):
                    left = centered_bound(L, e, f, P, P + 1, g) - TARGET
                    right = centered_bound(L, e, f, P, 8 * P, g) - TARGET
                    candidates.extend(((left, g, L, e, f, P, "left"),
                                       (right, g, L, e, f, P, "right")))
    weakest = min(candidates)
    require(weakest[0] == F(1687815, 5745649440454)
            and weakest[1:6] == (1, 168, 1, 12, 421), ("weakest tail", weakest))

    horn_asymptotic = tuple((f, F(1, 49) - F(96 - f, 5376),
                              F(1, 49) - F(96 - f, 5376) - TARGET)
                             for f in (1, 2, 3, 4))
    require(all(row[2] < 0 for row in horn_asymptotic), horn_asymptotic)
    return tuple(group_rows), horn_minimum, len(records), digest, weakest, horn_asymptotic


def literal_mass(L: int, cell: int, e: int, p: int, f: int, q: int) -> F:
    """Definition-level clipped interval merger, independent of floor moments."""
    def intervals(endpoint, level):
        modulus = L * level - endpoint
        residue = endpoint * cell % L
        radius = F(L, 14 * modulus)
        answer = []
        for index in range(-1, level + 1):
            centre = F(residue + index * L, modulus)
            left, right = max(F(0), centre - radius), min(F(1), centre + radius)
            if left < right:
                answer.append((left, right))
        return answer

    first, second = intervals(e, p), intervals(f, q)
    i = j = 0
    answer = F(0)
    while i < len(first) and j < len(second):
        left, right = max(first[i][0], second[j][0]), min(first[i][1], second[j][1])
        if left < right:
            answer += right - left
        if first[i][1] <= second[j][1]:
            i += 1
        else:
            j += 1
    return answer


def centered_grid_controls(contexts):
    fast, reference = fast_module(), reference_module()
    route_rows = (
        ((168, 90, 12, 4, 6, 5), F(92, 7645)),
        ((336, 174, 12, 3, 3, 5), F(158, 46397)),
        ((3920, 1995, 10, 6, 1, 33), F(1680, 129359)),
        ((4368, 2236, 8, 6, 1, 33), F(1872, 144143)),
        ((2352, 1204, 8, 6, 1, 33), F(1008, 77615)),
        ((3696, 1892, 8, 6, 1, 33), F(1584, 121967)),
        ((4368, 2262, 13, 6, 1, 33), F(1872, 144143)),
        ((588, 315, 1, 4, 14, 5), F(6003, 491359)),
        ((4368, 2262, 3, 6, 4, 33), F(468, 36035)),
    )
    semantic = sha256()
    for address, expected in route_rows:
        values = (fast.mass(*address), reference.mass(*address), literal_mass(*address))
        require(values == (expected, expected, expected), ("three-route control", address, values))
        semantic.update((repr((address, expected)) + "\n").encode())

    low_hostile = (168, 90, 1, 10, 3, 2)
    low_values = (fast.mass(*low_hostile), reference.mass(*low_hostile), literal_mass(*low_hostile))
    require(low_values == (F(0), F(0), F(0)), ("low hostile", low_values))

    # Deliberately feed a nonprimitive reduction into the analytic expression.
    # Its positive answer is false for the exact overlap, proving gcd(P,Q)=1
    # is a load-bearing hypothesis rather than cosmetic normalization.
    nonprimitive = (1680, 870, 5, 1792, 10, 2688)
    nonprimitive_values = (fast.mass(*nonprimitive), reference.mass(*nonprimitive),
                           literal_mass(*nonprimitive))
    naive = centered_bound(1680, 5, 10, 896, 1344, 2)
    require(nonprimitive_values == (F(0), F(0), F(0)),
            ("nonprimitive exact hostile", nonprimitive_values))
    require(naive == F(5457530976, 271903091713) and naive > 0,
            ("nonprimitive false-positive bound", naive))
    semantic.update((repr((low_hostile, low_values, nonprimitive, nonprimitive_values, naive)) + "\n").encode())
    control_digest = semantic.hexdigest()
    freeze(control_digest, EXPECTED_CONTROL, "control semantic")

    rng = Random(3355)
    checks = 0
    weakest = None
    random_semantic = sha256()
    while checks < 100000:
        L, cell, e, f = contexts[rng.randrange(len(contexts))]
        P = rng.randrange(1, 1000)
        Q = rng.randrange(P + 1, 8 * P + 1)
        if gcd(P, Q) != 1 or P + Q < 8:
            continue
        g = rng.randrange(1, 4)
        bound = centered_bound(L, e, f, P, Q, g)
        actual = fast.mass(L, cell, e, g * P, f, g * Q)
        require(actual >= bound, ("centered-grid failure", L, cell, e, f, P, Q, g, bound, actual))
        row = (actual - bound, L, cell, e, f, P, Q, g, bound, actual)
        if weakest is None or row < weakest:
            weakest = row
        random_semantic.update((repr(row) + "\n").encode())
        checks += 1
    random_digest = random_semantic.hexdigest()
    require(random_digest == EXPECTED_RANDOM, ("random semantic", random_digest))
    expected_weakest = (F(549323712710, 14172465110565231), 3920, 2100, 1, 8, 981, 7667, 1,
                        F(16981445620880, 833674418268543), F(294835779070, 14446957299251))
    require(weakest == expected_weakest, ("random weakest", weakest))
    return len(route_rows), control_digest, low_hostile, nonprimitive, naive, checks, random_digest, weakest


def finite_worker(item):
    shard, contexts = item
    mass = fast_module()
    semantic = sha256()
    analytic_gates = exact_checks = 0
    best = None
    for context_index in range(shard, len(contexts), 8):
        L, cell, e, f = contexts[context_index]
        for P in range(1, 421):
            for g in (1, 2, 3):
                if endpoint_bound(L, e, f, P, g) >= TARGET:
                    analytic_gates += 1
                    continue
                for Q in range(P + 1, 8 * P):
                    if gcd(P, Q) != 1 or P + Q < 8 or (P, Q) == (3, 5):
                        continue
                    value = mass.mass(L, cell, e, g * P, f, g * Q)
                    exact_checks += 1
                    record = (value, L, cell, e, f, P, Q, g)
                    require(value >= TARGET, ("finite exact failure", shard, record))
                    if best is None or record < best:
                        best = record
                    semantic.update((repr(record) + "\n").encode())
    return shard, analytic_gates, exact_checks, best, semantic.hexdigest()


def finite_scan(contexts, workers: int):
    require(workers >= 1, ("workers", workers))
    jobs = tuple((shard, contexts) for shard in range(8))
    if workers == 1:
        results = tuple(finite_worker(job) for job in jobs)
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            results = tuple(pool.map(finite_worker, jobs))
    results = tuple(sorted(results))
    require(tuple(row[0] for row in results) == tuple(range(8)), "shard order")
    require(tuple(row[1] for row in results) == EXPECTED_SHARD_GATES,
            ("shard analytic gates", tuple(row[1] for row in results)))
    require(tuple(row[2] for row in results) == EXPECTED_SHARD_CHECKS,
            ("shard exact checks", tuple(row[2] for row in results)))
    require(tuple(row[3] for row in results) == EXPECTED_SHARD_BEST,
            ("shard minima", tuple(row[3] for row in results)))
    require(tuple(row[4] for row in results) == EXPECTED_SHARD_DIGESTS,
            ("shard digests", tuple(row[4] for row in results)))
    require(sum(row[1] for row in results) == 3065046, "analytic gate total")
    require(sum(row[2] for row in results) == 6144244, "exact check total")
    global_best = min(row[3] for row in results)
    require(global_best == (F(92, 7645), 168, 90, 12, 6, 4, 5, 1), global_best)
    margin = global_best[0] - TARGET
    require(margin == F(19403, 2247630), margin)
    digest = sha256((repr(results) + "\n").encode()).hexdigest()
    freeze(digest, EXPECTED_FINITE_AGGREGATE, "finite aggregate")
    return results, digest, global_best, margin


def set_partitions(size: int):
    rows = [()]
    for vertex in range(size):
        next_rows = []
        for partition in rows:
            next_rows.append(partition + ((vertex,),))
            for index in range(len(partition)):
                next_rows.append(partition[:index]
                                 + (partition[index] + (vertex,),)
                                 + partition[index + 1:])
        rows = next_rows
    return tuple(rows)


def component_count(edges):
    adjacency = {vertex: set() for vertex in range(6)}
    for first, second in edges:
        adjacency[first].add(second)
        adjacency[second].add(first)
    unseen = set(range(6))
    count = 0
    while unseen:
        count += 1
        frontier = [min(unseen)]
        unseen.remove(frontier[0])
        while frontier:
            vertex = frontier.pop()
            found = adjacency[vertex] & unseen
            unseen -= found
            frontier.extend(found)
    return count


def multipartite_tree_count(partition):
    parts = tuple(map(len, partition))
    return 6 ** (len(parts) - 2) * prod((6 - part) ** (part - 1) for part in parts)


def graph_tree_audit():
    weak = frozenset((index, 5) for index in range(4))
    exceptional = ((0, 1, 2, 3, 4), (5,))
    partitions = tuple(row for row in set_partitions(6) if len(row) > 1)
    require(len(partitions) == 202, ("disconnected partitions", len(partitions)))
    graph_semantic, tree_semantic = sha256(), sha256()
    rows = []
    for partition in partitions:
        owner = {vertex: index for index, block in enumerate(partition) for vertex in block}
        graph = tuple((i, j) for i in range(6) for j in range(i + 1, 6) if owner[i] != owner[j])
        regular = tuple(edge for edge in graph if edge not in weak)
        components = component_count(regular)
        star = partition == exceptional
        row = (partition, graph, regular, components, star)
        rows.append(row)
        graph_semantic.update((repr(row) + "\n").encode())
        tree_semantic.update((repr((partition, multipartite_tree_count(partition))) + "\n").encode())
    graph_digest, tree_digest = graph_semantic.hexdigest(), tree_semantic.hexdigest()
    freeze(graph_digest, EXPECTED_GRAPH, "graph semantic")
    freeze(tree_digest, EXPECTED_TREE, "tree semantic")
    exceptional_rows = tuple(row for row in rows if row[4])
    require(len(exceptional_rows) == 1 and exceptional_rows[0][3] == 5,
            ("exceptional K1,5", exceptional_rows))
    require(max(row[3] for row in rows if not row[4]) <= 2, "nonexceptional component bound")
    component_hist = tuple((count, sum(row[3] == count for row in rows))
                           for count in sorted({row[3] for row in rows}))
    forest_hist = tuple((6 - count, multiplicity) for count, multiplicity in component_hist)
    tree_counts = tuple(multipartite_tree_count(row) for row in partitions)
    return (len(rows), component_hist, forest_hist, min(tree_counts), max(tree_counts),
            sum(tree_counts), graph_digest, tree_digest)


def debt(L: int, labels, levels) -> F:
    return sum((F(label, 7 * (L * level - label)) for label, level in zip(labels, levels)), F())


def debt_audit(body_rows):
    semantic = sha256()
    maximum = None
    census = 0
    for body, L in body_rows:
        for levels in permutations(range(1, 7)):
            value = debt(L, body, levels)
            row = (value, L, body, levels)
            if maximum is None or row > maximum:
                maximum = row
            semantic.update((repr(row) + "\n").encode())
            census += 1
    require(census == 649 * 720 and maximum is not None, ("debt census", census))
    require(maximum == (DEBT_MAX, 168, (1, 2, 3, 4, 6, 12), (6, 5, 4, 3, 2, 1)),
            ("Dmax census", maximum))
    census_digest = semantic.hexdigest()
    freeze(census_digest, EXPECTED_DEBT_CENSUS, "debt census semantic")

    five_margin = 5 * TARGET - DEBT_MAX
    require(five_margin == F(570672686921, 494472037809750) and five_margin > 0, five_margin)
    d421 = debt(168, (6, 4, 3, 2, 1, 12), (1, 2, 3, 4, 5, 421))
    require(d421 == F(443767487288, 52278303328335), d421)
    three_margin = 3 * TARGET - d421
    require(three_margin == F(1255584224873, 731896246596690) and three_margin > 0,
            three_margin)

    expected_debts = (
        F(10171035532358753244424, 87499329204988335395285245),
        F(65302219886882882438, 90087463329358292024115),
        F(140545706290782894, 31861415435911397875),
        F(7350964239952, 883497095223435),
    )
    expected_margins = (
        F(12072720679782123134489227, 3674971826609510086601980290),
        F(3832763666951738490749, 630612243305508044168805),
        F(2583990888487810609, 446059816102759570250),
        F(32685830817806, 6184479666564045),
    )
    exceptional_rows = []
    labels = (1, 2, 3, 4)
    p = 421
    for b in range(4):
        lower_count = b + 1
        first_level = p // (6 ** lower_count) + 1
        lower_levels = tuple(range(first_level, first_level + lower_count))
        upper_levels = tuple(range(p + 1, p + 1 + 4 - b))
        best = None
        for chosen in combinations(labels, b):
            lower_labels = (6,) + chosen
            upper_labels = tuple(label for label in labels if label not in chosen)
            for lower_order in permutations(lower_labels):
                for upper_order in permutations(upper_labels):
                    value = (F(12, 7 * (168 * p - 12))
                             + debt(168, lower_order, lower_levels)
                             + debt(168, upper_order, upper_levels))
                    row = (value, lower_order, lower_levels, upper_order, upper_levels)
                    if best is None or row[0] > best[0]:
                        best = row
        require(best is not None and best[0] == expected_debts[b], ("exceptional debt", b, best))
        regular_edges = 1 + b
        margin = regular_edges * TARGET - best[0]
        require(margin == expected_margins[b] and margin > 0,
                ("exceptional margin", b, margin))
        exceptional_rows.append((b, regular_edges, *best, margin))
    exceptional_digest = sha256((repr(tuple(exceptional_rows)) + "\n").encode()).hexdigest()
    require(exceptional_digest == "fd0fc78eeabad929e0a93dcb562dfc7e631501b8d6ecb3e9a5c345285cf6087a",
            ("exceptional digest", exceptional_digest))
    return (census, census_digest, maximum, five_margin, d421, three_margin,
            tuple(exceptional_rows), exceptional_digest)


def main():
    parser = ArgumentParser()
    parser.add_argument("--workers", type=int, default=max(1, min(8, os.cpu_count() or 1)))
    args = parser.parse_args()

    dependencies, assert_nodes = dependency_gate()
    body_rows, contexts, lanes, horns = context_atlas()
    polynomial = polynomial_tail_audit(lanes, horns)
    controls = centered_grid_controls(contexts)
    graph = graph_tree_audit()
    debts = debt_audit(body_rows)
    finite = finite_scan(contexts, args.workers)

    print("LRC14 DISCONNECTED WEIGHTED HORN-TREE CLOSURE REVIEW CERTIFICATE")
    print("status=PROVED analytic gates + FINITE-EXACT verification for THM-3355")
    print("dependency_sha256=" + repr(dependencies))
    print(f"python_assert_nodes={assert_nodes};ordinary_and_O_semantics=identical")
    print(f"bodies={len(body_rows)};contexts={len(contexts)};rulers=29;lanes={len(lanes)};horn_lanes={horns};context_sha256={EXPECTED_CONTEXT}")
    print(f"centered_grid_formula=primitive P<Q: gLP/(49z)-2gL/(7w)-gL*abs(Qe-Pf)*floor(P^2/4)/(Pwz)")
    print(f"three_route_controls={controls[0]};control_sha256={controls[1]};low_zero={controls[2]}")
    print(f"nonprimitive_hostile={controls[3]};exact_mass=0;naive_false_positive={controls[4]}")
    print(f"seeded_centered_checks={controls[5]};random_sha256={controls[6]};weakest_slack={controls[7]}")
    print(f"tail_polynomial_groups={polynomial[0]};horn_min_coefficient={polynomial[1]};records={polynomial[2]};polynomial_sha256={polynomial[3]}")
    print(f"weakest_tail_margin={polynomial[4]};g1_horn_asymptotic={polynomial[5]}")
    print(f"finite_shard_gates={tuple(row[1] for row in finite[0])};total={sum(row[1] for row in finite[0])}")
    print(f"finite_shard_checks={tuple(row[2] for row in finite[0])};total={sum(row[2] for row in finite[0])};failures=0")
    print(f"finite_shard_sha256={tuple(row[4] for row in finite[0])}")
    print(f"finite_aggregate_sha256={finite[1]};global_minimum={finite[2]};margin={finite[3]}")
    print(f"multipartite_partitions={graph[0]};regular_component_hist={graph[1]};regular_forest_edge_hist={graph[2]}")
    print(f"multipartite_tree_count_min_max_sum={(graph[3], graph[4], graph[5])};graph_sha256={graph[6]};tree_sha256={graph[7]}")
    print(f"debt_census={debts[0]};debt_sha256={debts[1]};Dmax_row={debts[2]}")
    print(f"five_regular_minus_Dmax={debts[3]};D421={debts[4]};three_regular_minus_D421={debts[5]}")
    print(f"exceptional_K15_weighted_rows={debts[6]};exceptional_sha256={debts[7]}")
    print("verdict=PASS: all regular lanes close by exact head plus coefficient-positive tail; the four g=1 horns close by weighted multipartite-tree/debt compatibility in THM-3355")


if __name__ == "__main__":
    main()
