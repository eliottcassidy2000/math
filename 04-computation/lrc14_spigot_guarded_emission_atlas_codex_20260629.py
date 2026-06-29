"""HYP-3525: spigot-style guarded emission atlas for LRC14 proof packets.

This is a synthesis scout, not a numerical LRC proof.  It translates the
spigot-algorithm pattern into an exact proof-carrier discipline for recent
random031 work:

  head state + tail bound -> safe digit emission

becomes

  visible packet + hidden sidecar bound -> safe route/owner token emission.

The main target is to make "emit named debt" less slogan-like.  A quotient may
emit a theorem-facing token only when every hidden tail compatible with the
visible state lands in the same route/owner class; otherwise it must hold a
predigit sidecar such as HYP-3513 route R, HYP-3520 owner word, or HYP-3522
owner-filtration layer.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass
from itertools import combinations


@dataclass(frozen=True)
class Carrier:
    name: str
    spigot_motif: str
    lrc_translation: str
    preserved: tuple[str, ...]
    destroyed: tuple[str, ...]
    next_hook: str
    payload: int
    guard: int
    random031_fit: int
    scalar_risk: int

    @property
    def score(self) -> int:
        return 4 * self.payload + 3 * self.guard + 2 * self.random031_fit - 2 * self.scalar_risk


CARRIERS = (
    Carrier(
        name="guarded_route_emission_R",
        spigot_motif="emit a digit only after carry ambiguity is impossible",
        lrc_translation=(
            "HYP-3513 says private status can be compressed, but the full HYP-3490 "
            "route is safely emitted only with route sidecar R or a route theorem"
        ),
        preserved=("h3490_route", "private_firewall_status", "random031_named_route"),
        destroyed=("row_name", "raw_gate_order"),
        next_hook="Prove route reconstruction, or keep R in every random031 terminal packet.",
        payload=10,
        guard=10,
        random031_fit=9,
        scalar_risk=1,
    ),
    Carrier(
        name="owner_filtration_digit",
        spigot_motif="predigit holdback until the next carry resolves",
        lrc_translation=(
            "HYP-3522 splits seam owners into transport (23,93,113), "
            "branch-boundary lift (147,169), and residual (45,173)"
        ),
        preserved=("transport_word", "branch_boundary_lift", "residual_owner_pair", "route_R"),
        destroyed=("raw_seven_owner_shadow",),
        next_hook="Formalize the residual (45,173) lemma after transport, bracket, and route sidecars.",
        payload=9,
        guard=10,
        random031_fit=10,
        scalar_risk=1,
    ),
    Carrier(
        name="safe_seam_sheaf_quotient",
        spigot_motif="tail-error bound certifies the current digit",
        lrc_translation=(
            "HYP-3520's canary allows flow_class, allowed_exit, owner_union, "
            "and sheet_pgf_bucket as safe component quotients"
        ),
        preserved=("flow_class", "allowed_exit", "owner_union", "sheet_pgf_bucket"),
        destroyed=("component_size", "branch_hist", "endpoint_rank_scalar"),
        next_hook="Use these quotients as the allowed head state for random031 terminal emission.",
        payload=8,
        guard=9,
        random031_fit=10,
        scalar_risk=2,
    ),
    Carrier(
        name="terminal_split_digit",
        spigot_motif="bounded spigot preallocates enough terms before output",
        lrc_translation=(
            "HYP-3521 emits the terminal split 282=230 ordinary + 40 free-hole "
            "+ 12 pure bypass only after HYP-3486/HYP-3511/HYP-3510/HYP-3490 are attached"
        ),
        preserved=("ordinary_route", "free_hole_bracket", "pure_bypass", "firewall_guard"),
        destroyed=("raw_242_gate_routed_slogan",),
        next_hook="Turn the split into a finite theorem interface with three emission branches.",
        payload=8,
        guard=8,
        random031_fit=9,
        scalar_risk=2,
    ),
    Carrier(
        name="bbp_random_access_route",
        spigot_motif="extract one digit by modular head plus bounded tail",
        lrc_translation=(
            "For a target packet, compute only the local route token needed: "
            "owner residual, endpoint route, private firewall bit, or terminal class"
        ),
        preserved=("target_digit_address", "local_tail_bound", "proof_token"),
        destroyed=("global_enumeration_order",),
        next_hook="Build target-local certificate extractors for residual packets instead of whole-bank scans.",
        payload=7,
        guard=8,
        random031_fit=7,
        scalar_risk=2,
    ),
    Carrier(
        name="unbounded_coinductive_spigot",
        spigot_motif="produce an infinite stream by repeatedly proving the next digit safe",
        lrc_translation=(
            "A family proof can emit one packet token at a time if each token has a "
            "uniform sidecar bound and a tail theorem for all later rows"
        ),
        preserved=("family_tail_theorem", "token_order", "holdback_state"),
        destroyed=("finite_bank_completeness_slogan",),
        next_hook="State random031-style guards as reusable emission rules for future exception families.",
        payload=7,
        guard=7,
        random031_fit=6,
        scalar_risk=3,
    ),
    Carrier(
        name="continued_fraction_lft_state",
        spigot_motif="linear-fractional state transformation keeps interval enclosures exact",
        lrc_translation=(
            "Use LFT/interval heads for Farey, u=2t, and owner-current transports "
            "rather than decimalizing to floating proof scores"
        ),
        preserved=("exact_interval_head", "farey_address", "two_adic_branch"),
        destroyed=("floating_score",),
        next_hook="Attach exact interval heads to route-emission tests before Fejer or Haar smoothing.",
        payload=6,
        guard=7,
        random031_fit=6,
        scalar_risk=2,
    ),
    Carrier(
        name="raw_digit_stream_shadow",
        spigot_motif="digits without carry guard",
        lrc_translation=(
            "Raw counts such as 12, 40, 79, 242, or 282 are output shadows; "
            "they are not legal proof digits until sidecars certify the route"
        ),
        preserved=("count_shadow",),
        destroyed=("route", "owner_word", "free_hole_type", "terminal_exit"),
        next_hook="Use raw counts only as checksums after guarded emission succeeds.",
        payload=1,
        guard=1,
        random031_fit=3,
        scalar_risk=10,
    ),
)


EMISSION_TESTS = {
    "private_firewall_status": {
        "emit_if_any": ("C", "F", "N", "T", "I", "Q", "R"),
        "hold_if_missing": "colored axis or incidence/projection sidecar",
    },
    "h3490_route": {
        "emit_if_any": ("R", "route_reconstruction_theorem"),
        "hold_if_missing": "route sidecar R",
    },
    "random031_terminal_class": {
        "emit_if_all": ("flow_class", "allowed_exit", "sheet_pgf_bucket"),
        "hold_if_missing": "HYP-3520 safe seam-sheaf quotient",
    },
    "owner_residual_pair": {
        "emit_if_all": ("transport_word", "branch_boundary_lift", "residual_pair", "R"),
        "hold_if_missing": "HYP-3522 owner filtration plus HYP-3513 route sidecar",
    },
    "raw_count": {
        "emit_if_any": ("checksum_only",),
        "hold_if_missing": "never theorem-facing by itself",
    },
}


SAFE_QUOTIENTS = ("flow_class", "allowed_exit", "owner_union", "sheet_pgf_bucket")
UNSAFE_QUOTIENTS = ("owner_union_size", "endpoint_ranks", "branch_hist", "size", "mirror_closed")


def beats(left: Carrier, right: Carrier) -> bool:
    left_key = (left.score, left.payload, left.guard, left.random031_fit, -left.scalar_risk, left.name)
    right_key = (right.score, right.payload, right.guard, right.random031_fit, -right.scalar_risk, right.name)
    return left_key > right_key


def directed_three_cycles() -> int:
    total = 0
    for a, b, c in combinations(CARRIERS, 3):
        ab = beats(a, b)
        bc = beats(b, c)
        ca = beats(c, a)
        ba = beats(b, a)
        cb = beats(c, b)
        ac = beats(a, c)
        if (ab and bc and ca) or (ba and cb and ac):
            total += 1
    return total


def emit_or_hold(visible: set[str], target: str) -> str:
    rule = EMISSION_TESTS[target]
    if "emit_if_all" in rule:
        needed = set(rule["emit_if_all"])
        if needed <= visible:
            return f"emit:{target}"
        missing = tuple(sorted(needed - visible))
        return f"hold:{target}:missing={missing}"
    allowed = set(rule["emit_if_any"])
    if visible & allowed:
        return f"emit:{target}"
    return f"hold:{target}:missing={rule['hold_if_missing']}"


def main() -> None:
    print("HYP-3525 SPIGOT-STYLE GUARDED EMISSION ATLAS")
    print("status=SYNTHESIS / executable proof-carrier atlas; not an LRC14 proof")
    print("seed=https://en.wikipedia.org/wiki/Spigot_algorithm")
    print()

    print("## Spigot -> LRC Dictionary")
    dictionary = (
        ("digit", "route/owner/terminal proof token"),
        ("head state", "visible quotient fields retained by the packet"),
        ("tail bound", "sidecar theorem proving all hidden completions agree"),
        ("carry/predigit", "unresolved route or owner ambiguity"),
        ("safe emission", "quotient fiber is constant for the target predicate"),
        ("bounded spigot", "finite audit with named family theorem debt"),
        ("unbounded spigot", "coinductive packet stream with a repeated guard"),
        ("BBP-style extraction", "target-local certificate extraction without whole-bank enumeration"),
    )
    for left, right in dictionary:
        print(f"{left:<22} -> {right}")
    print()

    print("## Random031 Emission Tests")
    scenarios = {
        "colored_private_only": {"C"},
        "route_sidecar": {"R"},
        "safe_sheaf_head": {"flow_class", "allowed_exit", "sheet_pgf_bucket"},
        "owner_filtration_ready": {"transport_word", "branch_boundary_lift", "residual_pair", "R"},
        "raw_count_shadow": {"checksum_only"},
        "danger_endpoint_rank": {"endpoint_ranks"},
    }
    for name, visible in scenarios.items():
        emitted = tuple(emit_or_hold(visible, target) for target in EMISSION_TESTS)
        print(f"{name}: visible={tuple(sorted(visible))}")
        for item in emitted:
            print(f"  {item}")
    print()

    print("## Safe / Unsafe Quotient Canary")
    print(f"safe_quotients={SAFE_QUOTIENTS}")
    print(f"unsafe_quotients={UNSAFE_QUOTIENTS}")
    print(
        "reading=HYP-3520 already supplies the tail-error guard: safe quotients "
        "can emit terminal class; unsafe scalar quotients must hold sidecars."
    )
    print()

    print("## Carrier Tournament")
    ordered = sorted(CARRIERS, key=lambda item: (-item.score, item.name))
    print(f"score_hist={dict(sorted(Counter(c.score for c in CARRIERS).items()))}")
    print(f"directed_3cycles={directed_three_cycles()}")
    print("hamiltonian_path=" + " -> ".join(c.name for c in ordered))
    print()
    for carrier in ordered:
        print(f"- {carrier.name}")
        print(f"  score={carrier.score} payload={carrier.payload} guard={carrier.guard} random031_fit={carrier.random031_fit} scalar_risk={carrier.scalar_risk}")
        print(f"  spigot_motif={carrier.spigot_motif}")
        print(f"  lrc_translation={carrier.lrc_translation}")
        print(f"  preserved={carrier.preserved}")
        print(f"  destroyed={carrier.destroyed}")
        print(f"  next_hook={carrier.next_hook}")
    print()

    print("## Proof Pull")
    print("P1: State a GuardedEmission lemma: a proof token is legal only when the visible quotient is constant over all hidden tails, or a named holdback sidecar is carried.")
    print("P2: Apply it to random031: private status may emit through C/F/N/T or I/Q, but full route emission needs R or route reconstruction.")
    print("P3: Apply it to HYP-3520/HYP-3522: owner residual (45,173) may emit only after transport word, branch-boundary lift, residual pair, and route R are visible.")
    print("P4: Treat raw counts 12/40/79/242/282 as checksums, never as theorem-facing digits before guarded emission.")


if __name__ == "__main__":
    main()
