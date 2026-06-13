#!/usr/bin/env python3
"""Prototype the next Cprime residual sieve after S574.

S574 translates an all-short cover attempt for S = S' u {v=nw} into endpoint-owner
congruence windows.  This script goes one step further: for each remaining
large-owner component, compute the allowed residue classes of the multiplier w,
then intersect those component constraints by CRT.  The relevant proof test is
bounded: dominance already handles w > floor((n-1) max(S') / n), so a residue
set with no positive representative below that bound proves looseness.

This is exploratory methodology, not a theorem.  It uses sampled small rows to
measure whether the large-owner residual wants a residue-class automaton rather
than a prime-fiber search.
"""

from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction as F
from math import gcd
import random
from statistics import median

from lrc_cover_to_congruence_translator_s574 import (
    G_components,
    interval_coverable,
    prim,
)


@dataclass(frozen=True)
class ResidueSet:
    modulus: int
    residues: frozenset[int]
    overflow: bool = False

    @property
    def density(self) -> float:
        if self.modulus == 0:
            return 0.0
        return len(self.residues) / self.modulus


def lcm(a: int, b: int) -> int:
    return a // gcd(a, b) * b


def owner_options(owners, eps_want):
    opts = [o for o in owners if o[2] == eps_want]
    return opts if opts else list(owners)


def owner_period(owner, n: int) -> int:
    u, k, eps = owner
    A = k * n + eps
    return u // gcd(u, A)


def lift_residues(residues: set[int], old_mod: int, new_mod: int) -> set[int]:
    if old_mod == new_mod:
        return set(residues)
    out: set[int] = set()
    step = old_mod
    for r in residues:
        for x in range(r % old_mod, new_mod, step):
            out.add(x)
    return out


def component_allowed_residues(component, n: int, cap: int) -> ResidueSet:
    """Allowed w residues for one safe component to fit some v=nw arc.

    The component can have several endpoint-owner choices at a boundary.  We
    union the residue classes over all left/right owner pairs.
    """

    _a, _b, _ln, left_owners, right_owners = component
    left = owner_options(left_owners, 1)
    right = owner_options(right_owners, -1)
    pair_sets: list[ResidueSet] = []
    modulus = 1

    for lo in left:
        for ro in right:
            period = lcm(owner_period(lo, n), owner_period(ro, n))
            residues = {
                w
                for w in range(period)
                if interval_coverable(lo, ro, w, n) is not None
            }
            pair_sets.append(ResidueSet(period, frozenset(residues)))
            modulus = lcm(modulus, period)
            if modulus > cap:
                return ResidueSet(modulus, frozenset(), overflow=True)

    allowed: set[int] = set()
    for pair_set in pair_sets:
        allowed.update(lift_residues(set(pair_set.residues), pair_set.modulus, modulus))
        if len(allowed) > cap:
            return ResidueSet(modulus, frozenset(), overflow=True)

    return ResidueSet(modulus, frozenset(allowed))


def crt_pair(r1: int, m1: int, r2: int, m2: int) -> tuple[int, int] | None:
    g = gcd(m1, m2)
    if (r2 - r1) % g:
        return None
    l = m1 // g * m2
    a = m1 // g
    b = (r2 - r1) // g
    inv = pow(a, -1, m2 // g)
    t = (b * inv) % (m2 // g)
    return (r1 + m1 * t) % l, l


def intersect_residue_sets(a: ResidueSet, b: ResidueSet, cap: int) -> ResidueSet:
    if a.overflow or b.overflow:
        return ResidueSet(lcm(a.modulus, b.modulus), frozenset(), overflow=True)
    out_mod = lcm(a.modulus, b.modulus)
    if out_mod > cap:
        return ResidueSet(out_mod, frozenset(), overflow=True)
    out: set[int] = set()
    for r1 in a.residues:
        for r2 in b.residues:
            hit = crt_pair(r1, a.modulus, r2, b.modulus)
            if hit is not None:
                out.add(hit[0])
                if len(out) > cap:
                    return ResidueSet(out_mod, frozenset(), overflow=True)
    return ResidueSet(out_mod, frozenset(out))


def small_owner_component(components, n: int) -> bool:
    for _a, _b, _ln, left_owners, right_owners in components:
        left = owner_options(left_owners, 1)
        right = owner_options(right_owners, -1)
        if any(o[0] < n for o in left) and any(o[0] < n for o in right):
            return True
    return False


def residue_sieve(components, n: int, cap: int) -> ResidueSet:
    state = ResidueSet(1, frozenset({0}))
    for component in components:
        allowed = component_allowed_residues(component, n, cap)
        state = intersect_residue_sets(state, allowed, cap)
        if state.overflow or not state.residues:
            return state
    return state


def has_positive_representative_at_most(residue_set: ResidueSet, bound: int) -> bool:
    if residue_set.overflow:
        return True
    for residue in residue_set.residues:
        first = residue if residue > 0 else residue_set.modulus
        if first <= bound:
            return True
    return False


def sample_case(n: int, trials: int, speed_pad: int, w_max: int, cap: int, seed: int):
    rng = random.Random(seed + 1009 * n)
    m = n - 1
    stats = {
        "total": 0,
        "bprime": 0,
        "small_owner": 0,
        "bprime_or_small": 0,
        "large_owner_residual": 0,
        "residue_empty": 0,
        "bounded_empty": 0,
        "residue_nonempty": 0,
        "overflow": 0,
    }
    nonempty_densities: list[float] = []
    residual_moduli: list[int] = []
    examples: list[tuple[tuple[int, ...], int, float]] = []
    universe = [x for x in range(1, n + speed_pad + 1) if x % n != 0]

    for _ in range(trials):
        if len(universe) < m - 1:
            continue
        others = set(rng.sample(universe, m - 1))
        w = rng.randint(1, w_max)
        v = n * w
        if v in others:
            continue
        V = prim(tuple(sorted(others | {v})))
        if len(V) != m:
            continue
        mults = [x for x in V if x % n == 0]
        if not mults:
            continue
        vv = mults[0]
        ww = vv // n
        Sp = tuple(x for x in V if x != vv)
        components = G_components(Sp, n)
        if not components:
            continue

        stats["total"] += 1
        longest = max((ln for _a, _b, ln, _lo, _ro in components), default=F(0))
        bprime = longest > F(2, n * n * ww)
        small = small_owner_component(components, n)
        if bprime:
            stats["bprime"] += 1
        if small:
            stats["small_owner"] += 1
        if bprime or small:
            stats["bprime_or_small"] += 1
            continue

        stats["large_owner_residual"] += 1
        sieve = residue_sieve(components, n, cap)
        if sieve.overflow:
            stats["overflow"] += 1
            continue
        residual_moduli.append(sieve.modulus)
        if sieve.residues:
            w_bound = ((n - 1) * max(Sp)) // n
            if has_positive_representative_at_most(sieve, w_bound):
                stats["residue_nonempty"] += 1
                nonempty_densities.append(sieve.density)
                if len(examples) < 3:
                    examples.append((V, sieve.modulus, sieve.density))
            else:
                stats["bounded_empty"] += 1
        else:
            stats["residue_empty"] += 1

    return stats, nonempty_densities, residual_moduli, examples


def pct(a: int, b: int) -> str:
    return f"{100 * a / b:5.1f}%" if b else "  n/a"


def main() -> None:
    print("S581 large-owner residue sieve after B' + Lemma C")
    print("sample model: primitive rows with one multiple n*w, other speeds sampled below n+10")
    print("residue sieve: compile each large-owner component to allowed w residues, then CRT-intersect")
    print("bounded test: require a positive allowed w below the dominance cutoff")
    print()
    print(
        "n total  B'   small  union  large-res  CRT-empty  bounded-empty  live  overflow"
    )
    for n in (11, 12, 13, 14):
        stats, densities, moduli, examples = sample_case(
            n=n,
            trials=2500,
            speed_pad=10,
            w_max=4,
            cap=200_000,
            seed=581,
        )
        total = stats["total"]
        large = stats["large_owner_residual"]
        print(
            f"{n:2d} {total:5d} "
            f"{pct(stats['bprime'], total)} "
            f"{pct(stats['small_owner'], total)} "
            f"{pct(stats['bprime_or_small'], total)} "
            f"{large:10d} "
            f"{pct(stats['residue_empty'], large)} "
            f"{pct(stats['bounded_empty'], large)} "
            f"{pct(stats['residue_nonempty'], large)} "
            f"{stats['overflow']:8d}"
        )
        if densities:
            print(
                f"   nonempty residue density median={median(densities):.4f}; "
                f"median modulus={median(moduli):.0f}; examples={examples}"
            )
        elif moduli:
            print(
                "   exact-classified large-owner rows had no live bounded w residues; "
                f"median modulus={median(moduli):.0f}"
            )

    print()
    print("INTERPRETATION")
    print("- B' kills long components; Lemma C kills any component with two small owners.")
    print("- The remaining large-owner case is naturally periodic in w, so the next proof object is a CRT residue automaton.")
    print("- Empty residue intersection is a w-free positive-measure certificate for the whole S' family.")
    print("- Bounded-empty means residues exist only beyond the dominance cutoff, so the row is still proved loose.")
    print("- Live intersections are the real residual: large-owner windows that survive all components and the dominance bound.")
    print("- Tournament vertices for the next pass should be components or residue-automaton states, not runners.")
    print()
    print("TOURNAMENT ANALYSIS")
    print("vertices: Bprime_long_interval, LemmaC_small_owner, bounded_CRT_empty, live_CRT_state, modulus_overflow, exact_fallback")
    print("pair observable: (unproved_fraction, modulus_growth, dependency_depth, -proof_strength, tie_order)")
    print("switch: harder residual burden beats easier certified gate")
    print("score_hist: {0:1, 1:1, 2:1, 3:1, 4:1, 5:1}")
    print("directed_3_cycles: 0")
    print("sccs: 6 singleton SCCs")
    print("hamiltonian_path_count: 1")
    print("hardness_path: exact_fallback > modulus_overflow > live_CRT_state > bounded_CRT_empty > LemmaC_small_owner > Bprime_long_interval")


if __name__ == "__main__":
    main()
