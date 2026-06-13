#!/usr/bin/env python3
"""S599: two-block determinant Helly certificates for the Cprime residual.

S581/S595 turn every safe component of G(S') into a bounded language in the
multiplier w of v=nw.  The global cover attempt is the intersection of those
languages.  This script asks the next proof question:

    when the bounded language intersection is empty, how many component
    determinant rows are needed to witness emptiness?

The audit intentionally enumerates only the dominance-bounded window in w.  It
does not try to materialize the full CRT modulus.  This makes the output a
small certificate extractor rather than another big residue-set run.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import combinations
from math import ceil, floor, gcd
import random


def dist(x: F) -> F:
    x %= 1
    return min(x, 1 - x)


@dataclass(frozen=True)
class Component:
    a: F
    b: F
    length: F
    left_owners: tuple[tuple[int, int, int], ...]
    right_owners: tuple[tuple[int, int, int], ...]

    def left_safe_owners(self) -> tuple[tuple[int, int, int], ...]:
        out = tuple(o for o in self.left_owners if o[2] == 1)
        return out if out else self.left_owners

    def right_safe_owners(self) -> tuple[tuple[int, int, int], ...]:
        out = tuple(o for o in self.right_owners if o[2] == -1)
        return out if out else self.right_owners

    def small_left(self, n: int) -> tuple[tuple[int, int, int], ...]:
        return tuple(o for o in self.left_safe_owners() if o[0] < n)

    def small_right(self, n: int) -> tuple[tuple[int, int, int], ...]:
        return tuple(o for o in self.right_safe_owners() if o[0] < n)

    def owner_size_profile(self, n: int) -> str:
        left = "s" if self.small_left(n) else "L"
        right = "s" if self.small_right(n) else "L"
        return left + right


@dataclass(frozen=True)
class CoverWitness:
    w: int
    j: int
    left_owner: tuple[int, int, int]
    right_owner: tuple[int, int, int]
    r_left: int
    r_right: int
    det: int


def components(sp: tuple[int, ...], n: int) -> list[Component]:
    threshold = F(1, n)
    pts: dict[F, list[tuple[int, int, int]]] = {}
    for u in sp:
        for k in range(u + 1):
            for eps in (1, -1):
                pts.setdefault(F(k * n + eps, n * u) % 1, []).append((u, k, eps))

    order = sorted(pts)
    out: list[Component] = []
    for i, a in enumerate(order):
        b = order[(i + 1) % len(order)]
        length = b - a if b > a else b - a + 1
        mid = (a + length / 2) % 1
        if all(dist(u * mid) > threshold for u in sp):
            out.append(Component(a, b, length, tuple(pts[a]), tuple(pts[b])))
    return out


def primitive(vs: tuple[int, ...]) -> tuple[int, ...]:
    g = 0
    for v in vs:
        g = gcd(g, v)
    return tuple(sorted(v // g for v in vs))


def endpoint_pinned(owner: tuple[int, int, int], w: int, n: int) -> bool:
    u, k, eps = owner
    return (w * (k * n + eps)) % u == 0


def split_interval(a: F, length: F) -> list[tuple[F, F]]:
    b = a + length
    if b <= 1:
        return [(a, b)]
    return [(a, F(1)), (F(0), b - 1)]


def overlap(lo1: F, hi1: F, lo2: F, hi2: F) -> F:
    lo = max(lo1, lo2)
    hi = min(hi1, hi2)
    return hi - lo if hi > lo else F(0)


def n_cell_loads(comps: list[Component], n: int) -> tuple[F, ...]:
    loads = [F(0) for _ in range(n)]
    for comp in comps:
        for lo, hi in split_interval(comp.a, comp.length):
            for r in range(n):
                loads[r] += overlap(lo, hi, F(r, n), F(r + 1, n))
    return tuple(loads)


def dual_cell_cap_exit(comps: list[Component], n: int) -> bool:
    cap = F(2, n * n)
    return any(load > cap for load in n_cell_loads(comps, n))


def right_cell_of(x: F, n: int) -> int:
    y = (x % 1) * n
    return y.numerator % n if y.denominator == 1 else floor(y) % n


def left_cell_of(x: F, n: int) -> int:
    y = (x % 1) * n
    return (y.numerator - 1) % n if y.denominator == 1 else floor(y) % n


def origin_bisection_exit(comps: list[Component], n: int, w: int) -> bool:
    upper = [F(0) for _ in range(n)]
    lower = [F(0) for _ in range(n)]
    for comp in comps:
        pinned_left = [o for o in comp.small_left(n) if endpoint_pinned(o, w, n)]
        pinned_right = [o for o in comp.small_right(n) if endpoint_pinned(o, w, n)]
        if pinned_left:
            upper[right_cell_of(comp.a, n)] += comp.length
        if pinned_right:
            lower[left_cell_of(comp.b, n)] += comp.length
    cap = F(1, n * n)
    return any(load > cap for load in upper) or any(load > cap for load in lower)


def local_owner_exit(comp: Component, n: int, w: int) -> str | None:
    full_radius = F(2, n * n * w)
    half_radius = F(1, n * n * w)
    if comp.length > full_radius:
        return "Bprime_long_component"
    small_left = comp.small_left(n)
    small_right = comp.small_right(n)
    if small_left and small_right:
        return "LemmaC_two_small_owners"
    for owner in small_left + small_right:
        if not endpoint_pinned(owner, w, n):
            return "LemmaE_small_pin_off_lattice"
    if (small_left or small_right) and comp.length > half_radius:
        return "LemmaF_half_cap"
    return None


def cover_witness_for_pair(
    comp: Component,
    left_owner: tuple[int, int, int],
    right_owner: tuple[int, int, int],
    w: int,
    n: int,
) -> CoverWitness | None:
    ua, ka, ea = left_owner
    ub, kb, eb = right_owner
    aa = w * (ka * n + ea)
    bb = w * (kb * n + eb)
    lo = max(F(aa, ua) - F(1, n), F(bb, ub) - F(1, n))
    hi = min(F(aa, ua) + F(1, n), F(bb, ub) + F(1, n))

    for j in range(floor(lo), ceil(hi) + 1):
        if lo < j < hi:
            r_left = aa - j * ua
            r_right = bb - j * ub
            det = ua * r_right - ub * r_left
            return CoverWitness(w, j, left_owner, right_owner, r_left, r_right, det)
    return None


def component_language(
    comp: Component, n: int, w_bound: int
) -> tuple[frozenset[int], dict[int, CoverWitness]]:
    allowed: set[int] = set()
    witness: dict[int, CoverWitness] = {}
    for w in range(1, w_bound + 1):
        hit: CoverWitness | None = None
        for left in comp.left_safe_owners():
            for right in comp.right_safe_owners():
                hit = cover_witness_for_pair(comp, left, right, w, n)
                if hit is not None:
                    break
            if hit is not None:
                break
        if hit is not None:
            allowed.add(w)
            witness[w] = hit
    return frozenset(allowed), witness


def minimal_empty_subset(
    languages: list[frozenset[int]], max_size: int = 4
) -> tuple[int, ...] | None:
    for i, lang in enumerate(languages):
        if not lang:
            return (i,)
    for size in range(2, max_size + 1):
        for idxs in combinations(range(len(languages)), size):
            common = set(languages[idxs[0]])
            for idx in idxs[1:]:
                common &= languages[idx]
                if not common:
                    return idxs
    return None


def row_from_sample(n: int, rng: random.Random) -> tuple[int, ...] | None:
    m = n - 1
    universe = [x for x in range(1, n + 10) if x % n != 0]
    others = set(rng.sample(universe, m - 1))
    w = rng.randint(1, 4)
    v = n * w
    if v in others:
        return None
    row = primitive(tuple(sorted(others | {v})))
    if len(row) != m or not any(x % n == 0 for x in row):
        return None
    return row


def classify_row(row: tuple[int, ...], n: int, regime: str) -> dict[str, object]:
    mults = [v for v in row if v % n == 0]
    if not mults:
        return {"class": "no_multiple"}
    v = mults[0]
    w = v // n
    sp = tuple(x for x in row if x != v)
    comps = components(sp, n)
    if not comps:
        return {"class": "no_components"}

    if regime == "full_stack":
        if dual_cell_cap_exit(comps, n):
            return {"class": "dual_total_cell_cap"}
        if origin_bisection_exit(comps, n, w):
            return {"class": "origin_bisection_cap"}
        for comp in comps:
            local = local_owner_exit(comp, n, w)
            if local is not None:
                return {"class": local}
    elif regime == "BC_only":
        for comp in comps:
            local = local_owner_exit(comp, n, w)
            if local in ("Bprime_long_component", "LemmaC_two_small_owners"):
                return {"class": local}
    else:
        raise ValueError(regime)

    w_bound = ((n - 1) * max(sp)) // n
    languages: list[frozenset[int]] = []
    witnesses: list[dict[int, CoverWitness]] = []
    profiles = Counter()
    for comp in comps:
        lang, wit = component_language(comp, n, w_bound)
        languages.append(lang)
        witnesses.append(wit)
        profiles[comp.owner_size_profile(n)] += 1

    common = set(range(1, w_bound + 1))
    for lang in languages:
        common &= set(lang)
    subset = minimal_empty_subset(languages, max_size=4) if not common else None

    if common:
        cls = "bounded_live"
    elif subset is None:
        cls = "high_order_empty"
    elif len(subset) == 1:
        cls = "singleton_empty"
    elif len(subset) == 2:
        cls = "pair_empty"
    elif len(subset) == 3:
        cls = "triple_empty"
    else:
        cls = "quad_empty"

    return {
        "class": cls,
        "row": row,
        "w": w,
        "sp": sp,
        "components": comps,
        "w_bound": w_bound,
        "languages": languages,
        "witnesses": witnesses,
        "common": frozenset(common),
        "subset": subset,
        "profiles": profiles,
    }


def sample_regime(n: int, regime: str, trials: int, seed: int) -> tuple[Counter[str], dict[str, dict[str, object]]]:
    rng = random.Random(seed + 1009 * n + (17 if regime == "full_stack" else 0))
    totals: Counter[str] = Counter()
    examples: dict[str, dict[str, object]] = {}
    for _ in range(trials):
        row = row_from_sample(n, rng)
        if row is None:
            continue
        info = classify_row(row, n, regime)
        cls = str(info["class"])
        totals[cls] += 1
        examples.setdefault(cls, info)
    return totals, examples


def fmt_lang(lang: frozenset[int], limit: int = 12) -> str:
    vals = sorted(lang)
    if len(vals) <= limit:
        return "{" + ",".join(map(str, vals)) + "}"
    return "{" + ",".join(map(str, vals[:limit])) + ",...}"


def describe_example(label: str, info: dict[str, object], n: int) -> None:
    if "row" not in info:
        return
    row = info["row"]
    subset = info.get("subset")
    languages = info["languages"]
    comps = info["components"]
    print(f"  {label}: row={row} w={info['w']} w_bound={info['w_bound']} components={len(comps)}")
    print(f"    owner_profiles={dict(info['profiles'])} common={fmt_lang(info['common'])}")
    if subset is not None:
        print(f"    minimal_empty_subset={subset}")
        for idx in subset:
            comp = comps[idx]
            print(
                f"      C{idx}: len={comp.length} owners={comp.owner_size_profile(n)} "
                f"allowed={fmt_lang(languages[idx])}"
            )


def tournament_fingerprint(route_totals: Counter[str]) -> str:
    vertices = [
        ("singleton_empty", 6, route_totals["singleton_empty"]),
        ("pair_empty", 5, route_totals["pair_empty"]),
        ("triple_empty", 4, route_totals["triple_empty"]),
        ("quad_empty", 3, route_totals["quad_empty"]),
        ("high_order_empty", 2, route_totals["high_order_empty"]),
        ("bounded_live", 1, route_totals["bounded_live"]),
        ("preempted_gate", 0, route_totals["preempted_gate"]),
    ]

    def beats(a: tuple[str, int, int], b: tuple[str, int, int]) -> bool:
        return (a[1], a[2], a[0]) > (b[1], b[2], b[0])

    scores = Counter({v[0]: 0 for v in vertices})
    c3 = 0
    for a, b in combinations(vertices, 2):
        scores[(a if beats(a, b) else b)[0]] += 1
    for a, b, c in combinations(vertices, 3):
        ab = beats(a, b)
        bc = beats(b, c)
        ca = beats(c, a)
        if (ab and bc and ca) or ((not ab) and (not bc) and (not ca)):
            c3 += 1
    path = " > ".join(v[0] for v in sorted(vertices, key=lambda x: (x[1], x[2], x[0]), reverse=True))
    return (
        f"score_hist={dict(sorted(Counter(scores.values()).items()))}; "
        f"directed_3_cycles={c3}; SCCs={(1,) * len(vertices)}; "
        f"Hamiltonian_paths=1; tie_HP={path}"
    )


def main() -> None:
    print("S599 two-block determinant Helly audit for Cprime")
    print("bounded language: w in [1, floor((n-1) max(S')/n)]")
    print("component row: owner/slack 2x2 determinant window, unioned over endpoint owners")
    print()

    trials = 1200
    combined_helly: Counter[str] = Counter()
    for regime in ("BC_only", "full_stack"):
        print(f"REGIME {regime}")
        print(" n  rows  preempt  sing  pair  triple  quad  high  live")
        interesting: dict[str, tuple[int, dict[str, object]]] = {}
        for n in (6, 8, 10, 12, 14):
            totals, examples = sample_regime(n, regime, trials, 599)
            preempt = sum(v for k, v in totals.items() if k not in {
                "singleton_empty",
                "pair_empty",
                "triple_empty",
                "quad_empty",
                "high_order_empty",
                "bounded_live",
            })
            combined_helly.update(
                {
                    "singleton_empty": totals["singleton_empty"],
                    "pair_empty": totals["pair_empty"],
                    "triple_empty": totals["triple_empty"],
                    "quad_empty": totals["quad_empty"],
                    "high_order_empty": totals["high_order_empty"],
                    "bounded_live": totals["bounded_live"],
                    "preempted_gate": preempt,
                }
            )
            rows = sum(totals.values())
            print(
                f"{n:2d} {rows:5d} {preempt:8d} {totals['singleton_empty']:5d} "
                f"{totals['pair_empty']:5d} {totals['triple_empty']:7d} "
                f"{totals['quad_empty']:5d} {totals['high_order_empty']:5d} "
                f"{totals['bounded_live']:5d}"
            )
            if n == 14 and regime == "BC_only":
                for label in ("singleton_empty", "pair_empty", "triple_empty", "bounded_live"):
                    if label in examples:
                        describe_example(label, examples[label], n)
            for label in ("singleton_empty", "pair_empty", "triple_empty", "bounded_live"):
                if label in examples and label not in interesting:
                    interesting[label] = (n, examples[label])
        if interesting:
            print("  first witnesses in this regime")
            for label in ("singleton_empty", "pair_empty", "triple_empty", "bounded_live"):
                if label in interesting:
                    ex_n, ex_info = interesting[label]
                    describe_example(f"{label}@n={ex_n}", ex_info, ex_n)
        print()

    print("INTERPRETATION")
    print("- BC_only tests the S581/S595 large-owner automaton after Bprime and Lemma C.")
    print("- full_stack first removes S598 total/one-sided caps and Lemmas E/F.")
    print("- singleton_empty is often a local one-component determinant wall; pair_empty")
    print("  is the first genuinely Helly-style two-row obstruction.")
    print("- bounded_live rows are not counterexamples; they are places where this")
    print("  bounded-window certificate did not yet kill the S' family.")
    print()
    print("TOURNAMENT ANALYSIS")
    print("vertices: determinant certificate sizes, not runners or raw arcs")
    print("pair observable: (certificate_rank, sampled_route_count, name)")
    print("switch/gauge: smaller empty determinant subfamily beats larger/global residue burden")
    print(tournament_fingerprint(combined_helly))
    print()
    print("ASSUMPTION CHALLENGE")
    print("Alternative vertices considered: runners, gaps, cap centers, components,")
    print("endpoint owners, owner pairs, residue classes, prime powers, and proof obligations.")
    print("Chosen vertices are determinant component languages and minimal empty subsets.")
    print("The quotient preserves bounded cover feasibility for v=nw before dominance;")
    print("it destroys phase order beyond the bounded w window and the full CRT modulus.")
    print("Challenged assumption: automaton emptiness need not be proved globally first;")
    print("a small Helly witness may be the human proof object.")


if __name__ == "__main__":
    main()
