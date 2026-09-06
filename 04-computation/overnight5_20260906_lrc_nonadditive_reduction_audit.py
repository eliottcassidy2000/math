#!/usr/bin/env python3
"""Exact independent cutoff/universe audit; no physical or native census.

The frozen fourth-round polygon referee supplies inherited coefficient data.
This program checks its bytes/semantic digest and audits the new consumer.
The all-height nonadditive theorem now also follows from incoming THM-4437,
THM-4441 and the audited norm-four sharp bound; this route is independent
confirmation, not a priority claim or a necessary new theorem dependency.
"""
from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
import importlib.util
from itertools import combinations, combinations_with_replacement
from math import comb, gcd
from pathlib import Path
import re
import sys

sys.stdout.reconfigure(newline="\n")
BASE = Path(__file__).resolve().parents[1] / "05-knowledge/results"
REPO = Path(__file__).resolve().parents[1]
GATES = 0


def need(ok, label):
    global GATES
    GATES += 1
    if not ok:
        raise RuntimeError(label)


def coefficient_consumer():
    path = REPO / "04-computation/overnight4_20260906_lrc_parityfree_audit.py"
    expected = "b5f52380e28bd3883f95090ff0c06a89e667d7aafd0e4c241b8fca68c76c7c00"
    need(sha256(path.read_bytes()).hexdigest() == expected, "frozen polygon referee bytes")
    spec = importlib.util.spec_from_file_location("frozen_polygon_referee", path)
    inherited = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(inherited)
    excluded = {(0, 1, 1), (1, 1, 1), (1, 1, 2), (1, 2, 2)}
    patterns = [p for p in combinations_with_replacement(range(19), 3)
                if sum(x != 0 for x in p) >= 2 and gcd(*p) == 1
                and sum(x % 3 == 0 for x in p) <= 1 and p not in excluded]
    need(len(patterns) == 747, "inherited exact pattern universe")
    need(Counter(sum(x != 0 for x in p) for p in patterns) == {2: 49, 3: 698},
         "actual-zero patterns retained")
    digest = sha256()
    rows = []
    empty = []
    for p in patterns:
        defects, alpha, beta = inherited.polygon_pattern(p)
        digest.update(repr((p, defects, alpha, beta)).encode())
        need(alpha < Q(11, 60), ("strict consumer slope gap", p))
        if not defects:
            need(alpha == beta == 0, "empty means no carriers, not 0<0")
            empty.append(p)
        else:
            ratio = beta / (Q(11, 60) - alpha)
            need(alpha * 178 + beta < Q(11 * 178, 60), ("small-M tail", p))
            rows.append((ratio, p, alpha, beta))
    semantic = digest.hexdigest()
    need(semantic == "cf808062354debbefc1d8ead8ad0d10e9da5427cb42b8f083b6af24d0059c87c",
         "all inherited records frozen")
    need(empty == [(0, 1, 2)], "empty-defect hostile")
    top = max(row[0] for row in rows)
    maximizers = [row for row in rows if row[0] == top]
    need(top == Q(35280, 199) < 178 and len(maximizers) == 1
         and maximizers[0][1] == (7, 13, 18), "exact optimized coefficient cutoff")
    print("INHERITED_POLYGON_SHA256", expected)
    print("BOX_SEMANTIC_SHA256", semantic)
    print("SMALL_M patterns=747 nonempty=746 max_ratio", top, "maximizer", maximizers)
    print("EMPTY_DEFECT (0,1,2): N=0; strict width-count inequality is not used")


def arithmetic():
    target = Q(11, 140)
    count_target = target / Q(3, 7)
    need(count_target == Q(11, 60), "physical and every projection consumer q=3/(7c)")
    alpha = Q(6, 49) + Q(4, 7 * 19)
    gap = count_target - alpha
    A, B = Q(2, 7) / gap, Q(4, 3) / gap
    need((alpha, gap, A, B) == (Q(142, 931), Q(1721, 55860),
                                Q(15960, 1721), Q(74480, 1721)), "large-M constants")
    need(A * 53 + B == Q(920360, 1721) < 535, "short-S integer cutoff")
    g54 = Q(3 * 54**2, 16) - A * 54 - B
    derivative54 = Q(3 * 54, 8) - A
    need(g54 == Q(18547, 6884) > 0 and derivative54 > 0,
         "long-S quadratic positive and increasing on [54,infinity)")
    need(A * 53 + B > 534, "535 cutoff is not silently rounded down")
    norm5bulk = Q(3, 49) * Q(121, 128)
    norm5cutoff = Q(4, 7) / (target - norm5bulk)
    need(norm5bulk == Q(363, 6272) and norm5cutoff == Q(17920, 649) < 28,
         "inherited norm-five continuum consumer")
    need(Q(46, 665) < target and Q(6, 77) < target,
         "incoming sharper sector bounds imply strict separation")
    norm4gap = target - (Q(3, 49) + Q(4, 7 * 33))
    need(norm4gap == Q(1, 32340) > 0, "norm-four frozen tail control")
    print("LARGE_M alpha", alpha, "gap", gap, "A", A, "B", B)
    print("S_LE_53 A*S+B", A * 53 + B, "<535; g(54)", g54, "g_prime(54)", derivative54)
    print("STRICT_GENERAL_TAIL c>=535; head c<=534 suffices; H535 census has one redundant layer")
    print("NORM5_INHERITED bulk", norm5bulk, "cutoff", norm5cutoff, "head<=27")
    print("NEW_INCOMING_ROUTE generic<=6/77, norm5<=46/665, norm4<=11/140")


def mobius(H):
    # Divisor recurrence, independent of any prime sieve or native loop.
    mu = [0] * (H + 1)
    mu[1] = 1
    for d in range(1, H + 1):
        for multiple in range(2 * d, H + 1, d):
            mu[multiple] -= mu[d]
    for n in range(1, H + 1):
        need(sum(mu[d] for d in range(1, n + 1) if n % d == 0) == (n == 1),
             ("Mobius inverse divisor identity", n))
    return mu


def units(H):
    return H - H // 3


def raw_additive(H):
    # a=r+3i, b=r+3j, i<j and i+j<=floor((H-2r)/3).
    out = 0
    for r in (1, 2):
        J = (H - 2 * r) // 3
        if J >= 0:
            out += (J + 1)**2 // 4
    return out


def primitive_counts(H, mu):
    total = sum(mu[d] * comb(units(H // d), 3)
                for d in range(1, H + 1) if d % 3)
    additive = sum(mu[d] * raw_additive(H // d)
                   for d in range(1, H + 1) if d % 3)
    return total, additive, total - additive


def universe():
    H = 535
    mu = mobius(H)
    # Direct small-height count controls use no gcd-inclusion formula.
    by_height = Counter()
    additive_by_height = Counter()
    for a, b, c in combinations([x for x in range(1, 81) if x % 3], 3):
        if gcd(a, b, c) == 1:
            by_height[c] += 1
            if a + b == c:
                additive_by_height[c] += 1
    for K in range(1, 81):
        direct = sum(by_height[c] for c in range(K + 1))
        additive = sum(additive_by_height[c] for c in range(K + 1))
        need(primitive_counts(K, mu) == (direct, additive, direct - additive),
             ("literal count-only small head", K))
    # Separate one-dimensional residue count of raw additive pairs.
    for K in range(1, H + 1):
        count = 0
        for c in range(1, K + 1):
            if c % 3:
                top = (c - 1) // 2
                first = (-c) % 3
                count += 0 if top < first else (top - first) // 3 + 1
        need(raw_additive(K) == count, ("raw additive residue identity", K))
    total, additive, nonadditive = primitive_counts(H, mu)
    independent_additive = sum(1 for a in range(1, H + 1) if a % 3
                              for b in range(a + 1, H - a + 1)
                              if b % 3 and (a + b) % 3 and gcd(a, b) == 1)
    need(additive == independent_additive, "full-height additive direct coprimality count")
    need((total, additive, nonadditive) == (6514176, 10908, 6503268),
         "independent H535 universe result")
    native = BASE / "overnight5_20260906_lrc_nonadditive_native_h535.out"
    data = native.read_text(encoding="utf-8")
    match = re.search(r"UNIVERSE H=(\d+) primitive_unit=(\d+) additive=(\d+) nonadditive=(\d+)", data)
    need(match is not None and tuple(map(int, match.groups())) == (H, total, additive, nonadditive),
         "native reported universe agrees; numerical masses not independently checked here")
    print("MOBIUS_H535 primitive_unit", total, "additive", additive, "nonadditive", nonadditive)
    print("MOBIUS_H534", primitive_counts(534, mu), "redundant_H535_layer",
          tuple(x-y for x,y in zip(primitive_counts(535,mu),primitive_counts(534,mu))))
    print("COUNT_CONTROLS literal all heights1..80; raw additive identities1..535; full-height additive gcd count")
    print("NATIVE_OUTPUT_SHA256", sha256(native.read_bytes()).hexdigest())


def main():
    print("SCOPE inherited coefficient consumer, analytic tails, independent universe count; no new physical census")
    coefficient_consumer()
    arithmetic()
    universe()
    print("PASS_OPTIMIZATION_LIVE_GATES", GATES)


if __name__ == "__main__":
    main()
