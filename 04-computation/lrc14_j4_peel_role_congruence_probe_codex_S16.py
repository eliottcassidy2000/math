#!/usr/bin/env python3
r"""Exact depth-two operation-congruence probe for a j=4 peel-role word.

Fix the lexicographically first THM-741 flood root

    E_12 = {1,2,8,9,10,11,12,13,14}.

For an exact safe set G (a finite union of open intervals), insertion of a
new forbidden comb sends each parent component to zero or more child
components.  If ``births`` counts the extra children created by splits and
``deaths`` counts parent components with no child, then

    #components(child) - #components(parent) = births - deaths.

Define the canonical peel role by the sign of this Euler drift:

    A = split-dominant, B = balanced, C = kill-dominant.

The word tau(E;v1,...,vk) in {A,B,C}^* records these roles in chronological
order.  It is independent of arbitrary component numbering and is compatible
with prefix/subtree restriction by construction.  This script tests the
stronger question: is its kernel a congruence for a common labelled next-speed
operation inside the exact THM-741 horizons?

The answer is no at the cheapest nontrivial horizon.  Histories (3) and (4)
both have role C and the same unordered event profile, but common legal
insertion 26 gives successor roles C and B.  Thus neither the role word nor
the word plus unordered event multiplicities has a deterministic T_26 update.
The associated q(A),q(B),q(C)=(1,2,4) charges separate to 1 and 6 mod 7, so
even the augmented chi_7 sign is not updateable through this proposed map.

Tournament Analysis uses the 22 one-peel histories in the C fibre as
vertices, not runners, root edges, or Fano flags.  The pair observable is the
lexicographic difference of exact proof-work scores (V2,-m,r); the gauge
prefers the larger score and increasing inserted speed breaks a true tie.
This retains a workload order but destroys endpoint placement and therefore
does not retain the common-successor role demonstrated by the witness.

Scope: a narrowly scoped counterexample to this natural role map, not a
no-go theorem for every possible decorated history map and not a THM-741
metric closure.
"""

from __future__ import annotations

import hashlib
import importlib.util
from collections import Counter
from dataclasses import dataclass
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
CORE_PATH = ROOT / "04-computation/lrc14_thm741_2002_body_j4_tree_kps_S128c5.py"
CORE_SHA256 = "5aa81d9d78273c8f9e3e7a6574091a3bc3f64ab6086c7024c15f9420c99dac96"
BODY = (1, 2, 8, 9, 10, 11, 12, 13, 14)
CHARGE = {"A": 1, "B": 2, "C": 4}
QR = frozenset({1, 2, 4})


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_core():
    require(sha256(CORE_PATH) == CORE_SHA256, "THM-741 core dependency changed")
    spec = importlib.util.spec_from_file_location("thm741_peel_role_dependency", CORE_PATH)
    require(spec is not None and spec.loader is not None, "cannot load THM-741 core")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@dataclass(frozen=True)
class Event:
    histogram: tuple[tuple[int, int], ...]
    births: int
    deaths: int
    delta: int
    role: str


@dataclass(frozen=True)
class Node:
    speed: int
    role: str
    event: Event
    components: int
    measure: Fraction
    horizon: int
    geometry: tuple[tuple[Fraction, Fraction], ...]


def event_profile(
    parent: tuple[tuple[Fraction, Fraction], ...] | list[tuple[Fraction, Fraction]],
    child: tuple[tuple[Fraction, Fraction], ...] | list[tuple[Fraction, Fraction]],
) -> Event:
    """Unordered component-incidence profile for one exact insertion."""
    multiplicities = []
    assigned = 0
    for parent_lo, parent_hi in parent:
        count = sum(
            parent_lo <= child_lo and child_hi <= parent_hi
            for child_lo, child_hi in child
        )
        multiplicities.append(count)
        assigned += count
    require(assigned == len(child), "a child component escaped its unique parent")
    births = sum(max(count - 1, 0) for count in multiplicities)
    deaths = sum(count == 0 for count in multiplicities)
    delta = births - deaths
    require(delta == len(child) - len(parent), "component Euler ledger failed")
    role = "A" if delta > 0 else "B" if delta == 0 else "C"
    return Event(tuple(sorted(Counter(multiplicities).items())), births, deaths, delta, role)


def component_multiplicities(
    parent: tuple[tuple[Fraction, Fraction], ...] | list[tuple[Fraction, Fraction]],
    child: tuple[tuple[Fraction, Fraction], ...] | list[tuple[Fraction, Fraction]],
) -> tuple[int, ...]:
    return tuple(
        sum(parent_lo <= child_lo and child_hi <= parent_hi for child_lo, child_hi in child)
        for parent_lo, parent_hi in parent
    )


def cyclic_positive_run_lengths(word: tuple[int, ...]) -> tuple[int, ...]:
    """Dihedral-invariant lengths between zero (killed-component) blocks."""
    require(word and 0 in word and any(value > 0 for value in word), "degenerate cyclic event word")
    runs = []
    for index, value in enumerate(word):
        if value <= 0 or word[index - 1] > 0:
            continue
        length = 0
        while word[(index + length) % len(word)] > 0:
            length += 1
        runs.append(length)
    return tuple(sorted(runs))


def chi_hat_7(q: int) -> int:
    q %= 7
    return 1 if q == 0 or q in QR else -1


def word_charge(word: str) -> int:
    return sum(CHARGE[role] for role in word) % 7


def main() -> None:
    core = load_core()
    require(core.S2 == Fraction(99, 70) and core.S2 * core.S2 > 2, "bad sqrt(2) majorant")
    require(BODY == tuple(sorted(set(range(8, 15)) | {1, 2})), "bad flood root")

    root_geometry, root_r, root_m = core.good_norm(BODY)
    threshold = 3 * root_m / (core.S2 * root_r)
    V1 = core.minV(4, threshold.numerator, threshold.denominator)
    require((root_r, root_m, V1) == (32, Fraction(319927, 2522520), 476), "root ledger changed")

    body_set = set(BODY)
    nodes = []
    for speed in range(1, V1):
        if speed in body_set:
            continue
        r1, m1, geometry1 = core.subtract(root_geometry, speed)
        require(m1 > 0, f"empty one-peel geometry at speed {speed}")
        event = event_profile(root_geometry, geometry1)
        threshold1 = 4 * m1 / (core.S2 * r1)
        V2 = core.minV(3, threshold1.numerator, threshold1.denominator)
        nodes.append(
            Node(speed, event.role, event, r1, m1, V2, tuple(geometry1))
        )
    role_census = Counter(node.role for node in nodes)
    require(len(nodes) == 466, "one-peel census changed")
    require(role_census == Counter({"A": 429, "C": 22, "B": 15}), "role census changed")

    # Lexicographic falsification schedule: increasing first history, then
    # increasing second history, then increasing common exact next speed.
    # Only pairs in one tau fibre are relevant.  Stop at the first square
    # whose two successor roles differ.
    checked_extensions = 0
    visited_kernel_pairs = 0
    witness = None
    for left_index, left in enumerate(nodes):
        for right in nodes[left_index + 1 :]:
            if left.role != right.role:
                continue
            visited_kernel_pairs += 1
            for speed in range(
                max(left.speed, right.speed) + 1,
                min(left.horizon, right.horizon),
            ):
                if speed in body_set:
                    continue
                checked_extensions += 1
                left_r2, left_m2, left_geometry2 = core.subtract(left.geometry, speed)
                right_r2, right_m2, right_geometry2 = core.subtract(right.geometry, speed)
                require(left_m2 > 0 and right_m2 > 0, "unexpected empty depth-two residual")
                left_event2 = event_profile(left.geometry, left_geometry2)
                right_event2 = event_profile(right.geometry, right_geometry2)
                if left_event2.role != right_event2.role:
                    witness = (
                        left,
                        right,
                        speed,
                        left_r2,
                        left_m2,
                        tuple(left_geometry2),
                        left_event2,
                        right_r2,
                        right_m2,
                        tuple(right_geometry2),
                        right_event2,
                    )
                    break
            if witness is not None:
                break
        if witness is not None:
            break
    require(witness is not None, "no depth-two kernel witness found")
    (
        left,
        right,
        common_speed,
        left_r2,
        left_m2,
        left_geometry2,
        left_event2,
        right_r2,
        right_m2,
        right_geometry2,
        right_event2,
    ) = witness
    require((left.speed, right.speed, common_speed) == (3, 4, 26), "minimal witness changed")
    require((checked_extensions, visited_kernel_pairs) == (15, 1), "minimality ledger changed")
    require(left.role == right.role == "C", "parent role fibre changed")
    require(
        left.event == right.event == Event(((0, 4), (1, 28)), 0, 4, -4, "C"),
        "parent unordered event kernel changed",
    )
    require(
        (left.components, left.measure, left.horizon)
        == (28, Fraction(4163, 51480), 368),
        "left parent changed",
    )
    require(
        (right.components, right.measure, right.horizon)
        == (28, Fraction(19909, 194040), 290),
        "right parent changed",
    )
    require(
        (left_r2, left_m2, left_event2)
        == (
            26,
            Fraction(36613, 504504),
            Event(((0, 2), (1, 26)), 0, 2, -2, "C"),
        ),
        "left child changed",
    )
    require(
        (right_r2, right_m2, right_event2)
        == (
            28,
            Fraction(6005, 72072),
            Event(((0, 2), (1, 24), (2, 2)), 2, 2, 0, "B"),
        ),
        "right child changed",
    )
    left_parent_cyclic = cyclic_positive_run_lengths(
        component_multiplicities(root_geometry, left.geometry)
    )
    right_parent_cyclic = cyclic_positive_run_lengths(
        component_multiplicities(root_geometry, right.geometry)
    )
    require(
        (left_parent_cyclic, right_parent_cyclic) == ((8, 20), (12, 16)),
        "cyclic killed-component phases changed",
    )

    left_threshold2 = 5 * left_m2 / (core.S2 * left_r2)
    right_threshold2 = 5 * right_m2 / (core.S2 * right_r2)
    left_V3 = core.minV(2, left_threshold2.numerator, left_threshold2.denominator)
    right_V3 = core.minV(2, right_threshold2.numerator, right_threshold2.denominator)
    require((left_V3, right_V3) == (203, 191), "depth-two horizons changed")

    parent_word = left.role
    left_word = parent_word + left_event2.role
    right_word = parent_word + right_event2.role
    require((parent_word, left_word, right_word) == ("C", "CC", "CB"), "word square changed")
    # Prefix/subtree compatibility is genuine but strictly weaker than the
    # failed operation congruence.
    require(left_word[:-1] == right_word[:-1] == parent_word, "prefix law failed")
    parent_charge = word_charge(parent_word)
    left_charge = word_charge(left_word)
    right_charge = word_charge(right_word)
    require((parent_charge, left_charge, right_charge) == (4, 1, 6), "charge split changed")
    require((chi_hat_7(left_charge), chi_hat_7(right_charge)) == (1, -1), "chi7 sign split changed")

    # Tournament on histories in the witness's C fibre.  Its total-order
    # switch is useful only as a proof-work scheduler.
    c_nodes = [node for node in nodes if node.role == "C"]
    require(len(c_nodes) == 22, "C-fibre size changed")
    hardness = lambda node: (node.horizon, -node.measure, node.components, -node.speed)
    require(len({hardness(node) for node in c_nodes}) == len(c_nodes), "hardness tie survived speed gauge")
    path_nodes = sorted(c_nodes, key=hardness, reverse=True)
    path = tuple(node.speed for node in path_nodes)
    expected_path = (
        3,
        23,
        34,
        22,
        25,
        17,
        50,
        19,
        16,
        4,
        6,
        15,
        27,
        20,
        21,
        45,
        33,
        48,
        5,
        18,
        24,
        7,
    )
    require(path == expected_path, "tournament Hamiltonian path changed")
    rank = {speed: index for index, speed in enumerate(path)}
    lex_speeds = sorted(rank)
    flips = sum(
        rank[left_speed] > rank[right_speed]
        for index, left_speed in enumerate(lex_speeds)
        for right_speed in lex_speeds[index + 1 :]
    )
    require(flips == 114, "tournament flip count changed")

    print("THM-741 j=4 PEEL-ROLE OPERATION-CONGRUENCE PROBE")
    print("=" * 84)
    print(f"dependency_sha256={CORE_SHA256}")
    print(f"fixed flood root E12={BODY}")
    print(f"root exact geometry: components={root_r} measure={root_m} V1={V1}")
    print(f"one-peel histories={len(nodes)} role_census={dict(sorted(role_census.items()))}")
    print("role law: A if births>deaths, B if equal, C if births<deaths")
    print("tau law: concatenate chronological roles; exact prefix/subtree restriction=PASS")
    print(
        "minimal schedule: increasing (v1,v1',common w) inside E12 exact horizons; "
        f"kernel_pairs_visited={visited_kernel_pairs} extensions_checked={checked_extensions}"
    )
    print(
        f"kernel parents: h=({left.speed}) and h'=({right.speed}); tau(h)=tau(h')={parent_word}; "
        f"event_hist={left.event.histogram}; births/deaths={left.event.births}/{left.event.deaths}"
    )
    print(
        f"  h=({left.speed}): components={left.components} measure={left.measure} V2={left.horizon}"
    )
    print(
        f"  h'=({right.speed}): components={right.components} measure={right.measure} V2={right.horizon}"
    )
    print(f"common labelled operation T_{common_speed}: legal in both exact J3 horizons")
    print(
        f"  T_{common_speed}h: tau={left_word} event_hist={left_event2.histogram} "
        f"births/deaths={left_event2.births}/{left_event2.deaths} "
        f"components={left_r2} measure={left_m2} V3={left_V3}"
    )
    print(
        f"  T_{common_speed}h': tau={right_word} event_hist={right_event2.histogram} "
        f"births/deaths={right_event2.births}/{right_event2.deaths} "
        f"components={right_r2} measure={right_m2} V3={right_V3}"
    )
    print("operation-congruence square: C --T26--> CC versus C --T26--> CB ; FAIL")
    print(
        f"chi7 charge: parent={parent_charge}; children={left_charge}/{right_charge}; "
        f"augmented signs={chi_hat_7(left_charge):+d}/{chi_hat_7(right_charge):+d}"
    )
    print("stronger kernel: the parents also share the full unordered event multiplicity history")
    print(
        "first missing coordinate: cyclic positive-run lengths between killed blocks are "
        f"{left_parent_cyclic} versus {right_parent_cyclic} (already distinct up to dihedral gauge)"
    )
    print("Tournament Analysis vertices: 22 C-role one-peel histories (not runners/edges/Fano flags)")
    print("pair observable: lex delta(V2,-measure,components); gauge: larger wins; true tie: smaller speed")
    print(
        f"fingerprint: score_hist={{0..21:1}}, directed_3cycles=0, SCC_sizes=22x1, "
        f"edge_flips_vs_speed_lex={flips}, Hamiltonian_paths=1"
    )
    print(f"gauge Hamiltonian path={path}")
    print("kept: role C and exact workload ordering; destroyed: endpoint phase, killed-component placement, successor role")
    print("VERDICT: this canonical Euler-drift role word is prefix-compatible but not operation-congruent")
    print("SCOPE: no claim against decorated metric-stalk role maps; no THM-741 body closure")
    print(f"source_sha256={sha256(Path(__file__).resolve())}")
    print("ALL EXACT ROLE-PROBE CHECKS PASSED")


if __name__ == "__main__":
    main()
