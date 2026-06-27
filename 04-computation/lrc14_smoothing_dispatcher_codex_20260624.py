#!/usr/bin/env python3
"""HYP-2985 smoothing dispatcher for the LRC14 proof stack.

This script is a proof-design artifact, not a proof of LRC14.  It merges the
S162/S163 analytic-sieve, Kaczynski-boundary, Ramanujan, and Fejer lanes into a
small typed dispatcher:

    packet family -> admissible smoothing / quotient -> retained labels
                  -> missing side channels -> next proof obligation.

Tournament Analysis uses smoothing policies as vertices, not runners.  The
pairwise observable is how much LRC predicate payload a policy keeps:

    exact scale, endpoint owner, strict-open interval, exact-period q,
    prime-power q, off-resonance family control, resonance/Kaczynski approach,
    interval formalizability, state-lift handoff.

The point is to prevent a scalar analytic quotient such as mu^2/phi or
sum mu(n)/n from pretending to be a final certificate when it only supplies a
tail or family normalizer.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from itertools import combinations


LABELS = [
    "lrc_predicate",
    "exact_scale",
    "endpoint_owner",
    "strict_open_interval",
    "exact_period_q",
    "prime_power_q",
    "c27_k33_state",
    "smoothing_kernel",
    "exceptional_boundary",
    "formal_interval",
    "state_lift_handoff",
]


@dataclass(frozen=True)
class Policy:
    name: str
    labels: frozenset[str]
    coordinates: tuple[int, ...]
    note: str


@dataclass(frozen=True)
class PacketFamily:
    name: str
    requires: frozenset[str]
    preferred: frozenset[str]
    forbidden_blindness: frozenset[str]
    note: str


POLICIES = [
    Policy(
        "raw_scalar_density",
        frozenset({"lrc_predicate"}),
        (1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
        "one number; useful only as a warning light",
    ),
    Policy(
        "mertens_mu_over_n_tail",
        frozenset({"lrc_predicate", "exact_scale", "smoothing_kernel"}),
        (2, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0),
        "signed cancellation tail, not a positive lonely interval",
    ),
    Policy(
        "selberg_large_sieve_squarefree",
        frozenset(
            {
                "lrc_predicate",
                "exact_scale",
                "exact_period_q",
                "smoothing_kernel",
                "exceptional_boundary",
            }
        ),
        (3, 2, 0, 0, 2, 0, 0, 2, 1, 0, 0),
        "family upper bound over squarefree primitive units; blind at prime powers",
    ),
    Policy(
        "ramanujan_exact_period_projector",
        frozenset(
            {
                "lrc_predicate",
                "exact_scale",
                "exact_period_q",
                "prime_power_q",
                "smoothing_kernel",
            }
        ),
        (3, 3, 0, 1, 4, 3, 1, 2, 0, 0, 0),
        "primitive phase side channel, including prime-power denominators",
    ),
    Policy(
        "fejer_toeplitz_interval",
        frozenset(
            {
                "lrc_predicate",
                "exact_scale",
                "endpoint_owner",
                "strict_open_interval",
                "exact_period_q",
                "c27_k33_state",
                "smoothing_kernel",
                "formal_interval",
                "state_lift_handoff",
            }
        ),
        (4, 4, 3, 5, 3, 1, 3, 4, 1, 5, 3),
        "finite interval certificate for positive-open packets",
    ),
    Policy(
        "kaczynski_abel_boundary",
        frozenset(
            {
                "lrc_predicate",
                "exact_scale",
                "endpoint_owner",
                "exact_period_q",
                "c27_k33_state",
                "smoothing_kernel",
                "exceptional_boundary",
                "state_lift_handoff",
            }
        ),
        (4, 4, 3, 1, 3, 1, 4, 4, 5, 1, 4),
        "true-wide approach-class ledger; off resonance decorrelates, resonance reduces",
    ),
    Policy(
        "hybrid_labelled_packet_sheaf",
        frozenset(LABELS),
        (5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5),
        "allowed final object: Fejer/Ramanujan/Kaczynski/state labels glued",
    ),
]


PACKETS = [
    PacketFamily(
        "AP_GW_boundary_atom",
        frozenset(
            {
                "lrc_predicate",
                "exact_scale",
                "endpoint_owner",
                "exact_period_q",
                "c27_k33_state",
                "exceptional_boundary",
            }
        ),
        frozenset({"kaczynski_abel_boundary", "hybrid_labelled_packet_sheaf"}),
        frozenset({"strict_open_interval"}),
        "zero-open equality atom at q=14; do not demand Fejer negativity",
    ),
    PacketFamily(
        "K33_near_12_to_36",
        frozenset(
            {
                "lrc_predicate",
                "exact_scale",
                "endpoint_owner",
                "strict_open_interval",
                "exact_period_q",
                "c27_k33_state",
                "formal_interval",
                "state_lift_handoff",
            }
        ),
        frozenset({"fejer_toeplitz_interval", "hybrid_labelled_packet_sheaf"}),
        frozenset(),
        "strict-open K33 packet; interval Fejer degree 159 is the current certificate shape",
    ),
    PacketFamily(
        "P10_plus_GW_splice",
        frozenset(
            {
                "lrc_predicate",
                "exact_scale",
                "endpoint_owner",
                "strict_open_interval",
                "exact_period_q",
                "prime_power_q",
                "c27_k33_state",
                "formal_interval",
            }
        ),
        frozenset(
            {
                "fejer_toeplitz_interval",
                "ramanujan_exact_period_projector",
                "hybrid_labelled_packet_sheaf",
            }
        ),
        frozenset({"selberg_large_sieve_squarefree"}),
        "q=27/C27 splice; squarefree large-sieve weights erase the prime-power address",
    ),
    PacketFamily(
        "petal_q27_unit_strip",
        frozenset(
            {
                "lrc_predicate",
                "exact_scale",
                "strict_open_interval",
                "exact_period_q",
                "prime_power_q",
                "c27_k33_state",
            }
        ),
        frozenset(
            {
                "ramanujan_exact_period_projector",
                "fejer_toeplitz_interval",
                "hybrid_labelled_packet_sheaf",
            }
        ),
        frozenset({"selberg_large_sieve_squarefree"}),
        "unit petal packets first seen at prime-power q=27",
    ),
    PacketFamily(
        "covering_q41_or_q63_front",
        frozenset(
            {
                "lrc_predicate",
                "exact_scale",
                "strict_open_interval",
                "exact_period_q",
                "endpoint_owner",
                "formal_interval",
            }
        ),
        frozenset({"fejer_toeplitz_interval", "hybrid_labelled_packet_sheaf"}),
        frozenset(),
        "covering rows are positive-open but can look all-covered in one chart",
    ),
    PacketFamily(
        "late_prime_power_denominator_wall",
        frozenset(
            {
                "lrc_predicate",
                "exact_scale",
                "exact_period_q",
                "prime_power_q",
                "smoothing_kernel",
                "exceptional_boundary",
            }
        ),
        frozenset(
            {
                "ramanujan_exact_period_projector",
                "kaczynski_abel_boundary",
                "hybrid_labelled_packet_sheaf",
            }
        ),
        frozenset({"selberg_large_sieve_squarefree"}),
        "HYP-2901 style late denominators: useful large-sieve normalizer, unsafe scalar endpoint",
    ),
    PacketFamily(
        "truewide_off_resonance_far_packet",
        frozenset(
            {
                "lrc_predicate",
                "exact_scale",
                "endpoint_owner",
                "smoothing_kernel",
                "exceptional_boundary",
            }
        ),
        frozenset({"kaczynski_abel_boundary", "hybrid_labelled_packet_sheaf"}),
        frozenset({"raw_scalar_density"}),
        "Fatou/decorrelated boundary value plus signed Abel tail",
    ),
    PacketFamily(
        "truewide_resonant_freiman_packet",
        frozenset(
            {
                "lrc_predicate",
                "exact_scale",
                "endpoint_owner",
                "c27_k33_state",
                "exceptional_boundary",
                "state_lift_handoff",
            }
        ),
        frozenset({"kaczynski_abel_boundary", "hybrid_labelled_packet_sheaf"}),
        frozenset({"raw_scalar_density", "mertens_mu_over_n_tail"}),
        "small relation-height far tuples; scale/Freiman reduce to a bounded atlas",
    ),
]


def label_score(policy: Policy, packet: PacketFamily) -> tuple[int, list[str], list[str]]:
    missing = sorted(packet.requires - policy.labels)
    retained = sorted(packet.requires & policy.labels)
    bonus = 2 if policy.name in packet.preferred else 0
    penalty = 2 if policy.name in packet.forbidden_blindness else 0
    return len(retained) + bonus - penalty, retained, missing


def dispatcher_rows() -> list[tuple[PacketFamily, Policy, Policy, int, int, list[str], list[str]]]:
    rows = []
    for packet in PACKETS:
        scored = []
        for policy in POLICIES:
            score, retained, missing = label_score(policy, packet)
            scored.append((score, len(missing), policy.name, policy, retained, missing))
        scored.sort(key=lambda row: (-row[0], row[1], row[2]))
        score, _, _, policy, retained, missing = scored[0]
        local_score, _, _, local_policy, _, _ = next(
            row for row in scored if row[3].name != "hybrid_labelled_packet_sheaf"
        )
        rows.append((packet, policy, local_policy, score, local_score, retained, missing))
    return rows


def beats(i: int, j: int) -> bool:
    a = POLICIES[i]
    b = POLICIES[j]
    wins = sum(x > y for x, y in zip(a.coordinates, b.coordinates))
    losses = sum(x < y for x, y in zip(a.coordinates, b.coordinates))
    if wins != losses:
        return wins > losses
    return i > j


def tournament_fingerprint() -> tuple[dict[int, int], int, int, list[str]]:
    n = len(POLICIES)
    scores = [0] * n
    for i, j in combinations(range(n), 2):
        if beats(i, j):
            scores[i] += 1
        else:
            scores[j] += 1
    score_hist = dict(sorted(Counter(scores).items()))
    c3 = 0
    for a, b, c in combinations(range(n), 3):
        if beats(a, b) and beats(b, c) and beats(c, a):
            c3 += 1
        if beats(a, c) and beats(c, b) and beats(b, a):
            c3 += 1

    dp: dict[tuple[int, int], int] = {}
    nmask = 1 << n
    for i in range(n):
        dp[(1 << i, i)] = 1
    for mask in range(nmask):
        for last in range(n):
            val = dp.get((mask, last), 0)
            if not val:
                continue
            for nxt in range(n):
                if mask & (1 << nxt):
                    continue
                if beats(last, nxt):
                    key = (mask | (1 << nxt), nxt)
                    dp[key] = dp.get(key, 0) + val
    full = nmask - 1
    hpaths = sum(dp.get((full, last), 0) for last in range(n))
    path = [POLICIES[i].name for i in sorted(range(n), key=lambda idx: -scores[idx])]
    return score_hist, c3, hpaths, path


def main() -> None:
    print("HYP-2985: LRC14 smoothing dispatcher proof pass")
    print("=" * 68)
    print("")
    print("Dispatcher readout")
    print("------------------")
    for packet, policy, local_policy, score, local_score, retained, missing in dispatcher_rows():
        missing_text = ", ".join(missing) if missing else "-"
        print(f"{packet.name}")
        print(f"  final policy: {policy.name}  score={score}")
        print(f"  best local policy: {local_policy.name}  local_score={local_score}")
        print(f"  retained: {', '.join(retained)}")
        print(f"  missing: {missing_text}")
        print(f"  packet note: {packet.note}")
        print(f"  local policy note: {local_policy.note}")
    print("")

    print("Quotient guardrail")
    print("------------------")
    print(
        "A quotient may forget a label only if the LRC predicate is constant on the "
        "fiber, the label is reconstructible, a dual/orthogonality certificate "
        "annihilates it, or the missing label is routed to a named residual packet."
    )
    print(
        "In particular, mu^2/phi is a Selberg/large-sieve normalizer, not a final "
        "LRC14 certificate: it sees q=14 and prime q=41 but erases prime-power "
        "or repeated-prime packets such as q=25,27,36,63,84,98,168,280,4312."
    )
    print("")

    print("Four-regime theorem skeleton")
    print("----------------------------")
    skeleton = [
        (
            "boundary equality",
            "AP/GW are endpoint-owner/Kaczynski boundary atoms; no Fejer negative "
            "certificate should be expected.",
        ),
        (
            "positive finite packet",
            "K33, petals, splices, coverings, and few-apex rows should receive "
            "Fejer interval certificates or Ramanujan exact-period witnesses.",
        ),
        (
            "late denominator minor arc",
            "after small packets are killed, Selberg/large-sieve weights can bound "
            "families only with exact-period and boundary labels attached.",
        ),
        (
            "true-wide boundary",
            "off-resonance far packets decorrelate by signed Abel/Kaczynski control; "
            "resonant far packets reduce by Freiman/scale to bounded finite atlases.",
        ),
    ]
    for name, text in skeleton:
        print(f"- {name}: {text}")
    print("")

    score_hist, c3, hpaths, path = tournament_fingerprint()
    print("Tournament Analysis")
    print("-------------------")
    print("vertices: smoothing/quotient policies, not runners")
    print(
        "pairwise observable: retained LRC predicate payload, exact-period data, "
        "boundary approach, interval formalizability, and state-lift handoff"
    )
    print(f"score_hist={score_hist}")
    print(f"directed_3cycles={c3}")
    print(f"hamiltonian_paths={hpaths}")
    print("Hamiltonian path:")
    print("  " + " > ".join(path))
    print("")

    print("New proof target")
    print("----------------")
    print(
        "Prove an admissible-smoothing lemma on HYP-2963 packet fibers: every "
        "primitive packet is either AP/GW boundary equality, a labelled Fejer "
        "interval certificate, a Ramanujan exact-period prime-power handoff, an "
        "off-resonance Kaczynski/Abel decay packet, a Freiman-reduced resonant "
        "finite atlas row, or a HYP-2908/THM-572 state-lift obligation."
    )


if __name__ == "__main__":
    main()
