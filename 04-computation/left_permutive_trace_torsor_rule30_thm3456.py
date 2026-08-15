#!/usr/bin/env python3
"""Exact stdlib companion for THM-3456.

The load-bearing statements are elementary and are proved in the companion
reflection.  This program independently freezes their finite controls:

* the sixteen binary left-permutive elementary rules are classified;
* after fixing ``x_0,...,x_n``, every length-``n`` center trace is realized by
  exactly one left wing ``x_-1,...,x_-n`` (all sixteen rules through depth
  six, and the Rule 30/Rule 60 hostile pair through depth eight);
* two adjacent temporal columns recover every finite-cylinder state;
* a single temporal column first has a nonconstant periodic-state alias at
  width seven (exhaustive over widths three through seven);
* a width-``n`` single-seed cylinder agrees with the infinite single-seed
  center through every time ``t<n``, but not uniformly at time ``n``; and
* exact single-seed cylinder orbit and trace periods are recorded through
  width twenty-four.

All state evolution is Boolean/integer exact.  Explicit exceptions rather
than ``assert`` keep normal and optimized runs identical.
"""

from __future__ import annotations

from collections import defaultdict, deque
from hashlib import sha256
import json
from pathlib import Path


RULE30 = 30
RULE60 = 60
LEFT_PERMUTIVE_RULES = (
    15, 30, 45, 60, 75, 90, 105, 120,
    135, 150, 165, 180, 195, 210, 225, 240,
)
TRIANGULAR_DEPTH = 8
ALL_RULE_TRIANGULAR_DEPTH = 6
ALIAS_MAX_WIDTH = 7
PREFIX_MAX_WIDTH = 64
ORBIT_MAX_WIDTH = 24
EXPECTED_SEMANTIC_SHA256 = "1fda893b46126cb9653f758a88a9f0266fe93c6e12bee6d4de527fa09f5b42ca"


def require(condition: bool, detail: object) -> None:
    if not condition:
        raise RuntimeError(detail)


def local(rule: int, left: int, center: int, right: int) -> int:
    neighborhood = (left << 2) | (center << 1) | right
    return (rule >> neighborhood) & 1


def is_left_permutive(rule: int) -> bool:
    return all(local(rule, 0, center, right) != local(rule, 1, center, right)
               for center in (0, 1) for right in (0, 1))


def finite_cone_trace(
    rule: int,
    left_wing: tuple[int, ...],
    right_block: tuple[int, ...],
) -> tuple[int, ...]:
    """Return center values at times 1..n from initial cells -n..n.

    ``left_wing[j-1]`` is x_-j and ``right_block[j]`` is x_j.  At each
    update the stored interval shrinks, so no artificial boundary value is
    ever queried.
    """

    depth = len(left_wing)
    require(len(right_block) == depth + 1, "right block has wrong length")
    state = list(reversed(left_wing)) + list(right_block)
    table = tuple((rule >> neighborhood) & 1 for neighborhood in range(8))
    trace: list[int] = []
    for time in range(1, depth + 1):
        state = [
            table[(state[index] << 2) | (state[index + 1] << 1) | state[index + 2]]
            for index in range(len(state) - 2)
        ]
        trace.append(state[depth - time])
    return tuple(trace)


def bits(code: int, length: int) -> tuple[int, ...]:
    return tuple((code >> index) & 1 for index in range(length))


def invert_trace(
    rule: int,
    right_block: tuple[int, ...],
    target: tuple[int, ...],
) -> tuple[int, ...]:
    """Invert the triangular left-wing-to-trace map one coordinate at a time."""

    depth = len(target)
    recovered: list[int] = []
    for time in range(1, depth + 1):
        candidates: list[int] = []
        for new_bit in (0, 1):
            trial = tuple(recovered + [new_bit] + [0] * (depth - time))
            if finite_cone_trace(rule, trial, right_block)[time - 1] == target[time - 1]:
                candidates.append(new_bit)
        require(len(candidates) == 1, ("triangular inverse failure", rule, time, candidates))
        recovered.append(candidates[0])
    return tuple(recovered)


def audit_triangular_bijections() -> dict[str, object]:
    require(all(
        local(RULE30, left, center, right) == (left ^ (center | right))
        for left in (0, 1) for center in (0, 1) for right in (0, 1)
    ), "Rule 30 algebra mismatch")
    require(all(
        local(RULE60, left, center, right) == (left ^ center)
        for left in (0, 1) for center in (0, 1) for right in (0, 1)
    ), "Rule 60 algebra mismatch")
    actual_left_permutive = tuple(rule for rule in range(256) if is_left_permutive(rule))
    require(actual_left_permutive == LEFT_PERMUTIVE_RULES,
            ("binary left-permutive ECA list", actual_left_permutive))

    fixed_right_maps = 0
    trace_evaluations = 0
    inversion_evaluations = 0
    depth_rows: list[tuple[int, int, int]] = []
    for depth in range(1, TRIANGULAR_DEPTH + 1):
        maps_at_depth = 0
        evaluations_at_depth = 0
        rules_at_depth = LEFT_PERMUTIVE_RULES if depth <= ALL_RULE_TRIANGULAR_DEPTH else (RULE30, RULE60)
        for rule in rules_at_depth:
            full_language = set(bits(code, depth) for code in range(1 << depth))
            for right_code in range(1 << (depth + 1)):
                right_block = bits(right_code, depth + 1)
                trace_to_left: dict[tuple[int, ...], tuple[int, ...]] = {}
                for left_code in range(1 << depth):
                    left_wing = bits(left_code, depth)
                    trace = finite_cone_trace(rule, left_wing, right_block)
                    require(trace not in trace_to_left,
                            ("trace collision", rule, depth, right_block, trace))
                    trace_to_left[trace] = left_wing
                    trace_evaluations += 1
                    evaluations_at_depth += 1
                require(set(trace_to_left) == full_language,
                        ("trace language not full", rule, depth, right_block))
                # The language census above is exhaustive for every fixed
                # right block.  Independently replay the constructive inverse
                # on the two extreme right blocks for every rule and depth.
                if right_code in (0, (1 << (depth + 1)) - 1):
                    for target, expected_left in trace_to_left.items():
                        actual_left = invert_trace(rule, right_block, target)
                        require(actual_left == expected_left,
                                ("inverse mismatch", rule, depth, right_block, target))
                        inversion_evaluations += 1
                fixed_right_maps += 1
                maps_at_depth += 1
        depth_rows.append((depth, maps_at_depth, evaluations_at_depth))

    return {
        "depth_rows": depth_rows,
        "left_permutive_rules": LEFT_PERMUTIVE_RULES,
        "all_rule_depth": ALL_RULE_TRIANGULAR_DEPTH,
        "distinguished_rule_depth": TRIANGULAR_DEPTH,
        "fixed_right_maps": fixed_right_maps,
        "trace_evaluations": trace_evaluations,
        "inversion_evaluations": inversion_evaluations,
    }


def cylinder_step(state: int, width: int, rule: int = RULE30) -> int:
    result = 0
    for position in range(width):
        left = (state >> ((position - 1) % width)) & 1
        center = (state >> position) & 1
        right = (state >> ((position + 1) % width)) & 1
        result |= local(rule, left, center, right) << position
    return result


def state_word(state: int, width: int) -> str:
    """Use the mathematical order x_0,x_1,...,x_(n-1)."""

    return "".join(str((state >> position) & 1) for position in range(width))


def word_state(word: str) -> int:
    return sum(int(value) << position for position, value in enumerate(word))


def divisors(value: int) -> tuple[int, ...]:
    small: list[int] = []
    large: list[int] = []
    candidate = 1
    while candidate * candidate <= value:
        if value % candidate == 0:
            small.append(candidate)
            if candidate * candidate != value:
                large.append(value // candidate)
        candidate += 1
    return tuple(small + list(reversed(large)))


def minimal_period(word: tuple[object, ...]) -> int:
    length = len(word)
    for period in divisors(length):
        if all(word[index] == word[index % period] for index in range(length)):
            return period
    raise RuntimeError(("no period", word))


def orbit_from_state(state: int, width: int) -> tuple[int, int, tuple[int, ...]]:
    seen: dict[int, int] = {}
    path: list[int] = []
    current = state
    while current not in seen:
        seen[current] = len(path)
        path.append(current)
        current = cylinder_step(current, width)
    transient = seen[current]
    cycle = tuple(path[transient:])
    return transient, len(cycle), cycle


def periodic_states(width: int) -> tuple[int, ...]:
    size = 1 << width
    successors = tuple(cylinder_step(state, width) for state in range(size))
    indegree = [0] * size
    for successor in successors:
        indegree[successor] += 1
    queue = deque(state for state, degree in enumerate(indegree) if degree == 0)
    while queue:
        state = queue.popleft()
        successor = successors[state]
        indegree[successor] -= 1
        if indegree[successor] == 0:
            queue.append(successor)
    return tuple(state for state, degree in enumerate(indegree) if degree > 0)


def periodic_trace_key(state: int, width: int, columns: tuple[int, ...]) -> tuple[tuple[int, ...], ...]:
    transient, period, cycle = orbit_from_state(state, width)
    require(transient == 0 and cycle[0] == state, ("state is not periodic", state, width))
    symbols = tuple(tuple((configuration >> column) & 1 for column in columns)
                    for configuration in cycle)
    trace_period = minimal_period(symbols)
    return symbols[:trace_period]


def audit_sideways_and_aliases() -> dict[str, object]:
    require(all(
        left == (local(RULE30, left, center, right) ^ (center | right))
        for left in (0, 1) for center in (0, 1) for right in (0, 1)
    ), "sideways inverse mismatch")

    width_rows: list[tuple[int, int, int, int]] = []
    first_nonconstant_alias_width: int | None = None
    target_alias: tuple[int, ...] | None = None
    target_states: tuple[int, ...] | None = None
    short_trace_state: int | None = None
    for width in range(3, ALIAS_MAX_WIDTH + 1):
        states = periodic_states(width)
        one_column: dict[tuple[tuple[int, ...], ...], list[int]] = defaultdict(list)
        for state in states:
            one_column[periodic_trace_key(state, width, (0,))].append(state)

            transient, full_period, cycle = orbit_from_state(state, width)
            require(transient == 0, ("periodic-state audit", width, state))
            pair_symbols = tuple(
                ((configuration >> 0) & 1, (configuration >> 1) & 1)
                for configuration in cycle
            )
            require(minimal_period(pair_symbols) == full_period,
                    ("adjacent traces lost full period", width, state))

        aliases = [
            (trace, tuple(group))
            for trace, group in one_column.items()
            if len(group) > 1 and len(trace) > 1
        ]
        width_rows.append((width, len(states), len(one_column), len(aliases)))
        if aliases and first_nonconstant_alias_width is None:
            first_nonconstant_alias_width = width

        if width == 7:
            desired = tuple((value,) for value in (0, 1, 1, 0))
            require(desired in one_column, "missing width-seven target alias")
            target_alias = tuple(value[0] for value in desired)
            target_states = tuple(one_column[desired])
            require(tuple(state_word(state, width) for state in target_states)
                    == ("0100000", "0111000", "0000111"),
                    ("width-seven alias changed", target_states))
            for state in states:
                transient, full_period, cycle = orbit_from_state(state, width)
                center_symbols = tuple(((configuration >> 0) & 1,) for configuration in cycle)
                if full_period == 4 and minimal_period(center_symbols) == 2:
                    short_trace_state = state
                    break

    require(first_nonconstant_alias_width == 7, first_nonconstant_alias_width)
    require(target_alias == (0, 1, 1, 0) and target_states is not None, "target alias")
    require(short_trace_state is not None, "missing proper-period trace hostile")
    require(state_word(short_trace_state, 7) == "0100110", state_word(short_trace_state, 7))

    alias_cycles = []
    for state in target_states:
        _, period, cycle = orbit_from_state(state, 7)
        require(period == 4, ("alias full period", state, period))
        alias_cycles.append(tuple(state_word(configuration, 7) for configuration in cycle))

    _, short_full_period, short_cycle = orbit_from_state(short_trace_state, 7)
    short_center = tuple((configuration >> 0) & 1 for configuration in short_cycle)
    short_pair = tuple(
        ((configuration >> 0) & 1, (configuration >> 1) & 1)
        for configuration in short_cycle
    )
    require((short_full_period, minimal_period(short_center), minimal_period(short_pair)) == (4, 2, 4),
            ("proper-period hostile changed", short_full_period, short_center, short_pair))

    return {
        "width_rows": width_rows,
        "first_nonconstant_alias_width": first_nonconstant_alias_width,
        "alias_trace": target_alias,
        "alias_states": tuple(state_word(state, 7) for state in target_states),
        "alias_cycles": tuple(alias_cycles),
        "proper_period_state": state_word(short_trace_state, 7),
        "proper_period_full_center_pair": (short_full_period, minimal_period(short_center), minimal_period(short_pair)),
    }


def infinite_single_seed_trace(rule: int, horizon: int) -> tuple[int, ...]:
    state = {0}
    trace: list[int] = []
    for _time in range(horizon + 1):
        trace.append(int(0 in state))
        lower = min(state, default=0) - 1
        upper = max(state, default=0) + 1
        state = {
            position
            for position in range(lower, upper + 1)
            if local(rule, int(position - 1 in state), int(position in state), int(position + 1 in state))
        }
    return tuple(trace)


def cylinder_center_trace(width: int, horizon: int) -> tuple[int, ...]:
    state = 1
    trace: list[int] = []
    for _time in range(horizon + 1):
        trace.append(state & 1)
        state = cylinder_step(state, width)
    return tuple(trace)


def audit_prefix_boundary() -> dict[str, object]:
    infinite = infinite_single_seed_trace(RULE30, 3 * PREFIX_MAX_WIDTH)
    first_differences: list[tuple[int, int | None]] = []
    for width in range(1, PREFIX_MAX_WIDTH + 1):
        cylinder = cylinder_center_trace(width, 3 * width)
        require(cylinder[:width] == infinite[:width], ("prefix theorem", width))
        first_difference = next(
            (time for time in range(min(len(cylinder), len(infinite))) if cylinder[time] != infinite[time]),
            None,
        )
        first_differences.append((width, first_difference))
    require(dict(first_differences)[3] == 3, ("sharp width-three hostile", first_differences[2]))
    require(dict(first_differences)[8] == 10, ("delayed width-eight hostile", first_differences[7]))
    exceptions = tuple((width, first_difference) for width, first_difference in first_differences
                       if first_difference != width)
    require(exceptions == ((8, 10),), ("first-difference exceptions changed", exceptions))
    return {
        "checked_widths": PREFIX_MAX_WIDTH,
        "uniform_prefix": "times 0..n-1",
        "sharp_hostile": (3, 3),
        "first_difference_exceptions": exceptions,
    }


def audit_single_seed_orbits() -> dict[str, object]:
    rows: list[tuple[int, int, int, int, int]] = []
    for width in range(1, ORBIT_MAX_WIDTH + 1):
        transient, full_period, cycle = orbit_from_state(1, width)
        center = tuple((configuration >> 0) & 1 for configuration in cycle)
        pair = tuple(
            ((configuration >> 0) & 1, (configuration >> (1 % width)) & 1)
            for configuration in cycle
        )
        center_period = minimal_period(center)
        pair_period = minimal_period(pair)
        require(pair_period == full_period, ("seed pair/full mismatch", width, pair_period, full_period))
        rows.append((width, transient, full_period, center_period, pair_period))
    require(all(full == center for width, _mu, full, center, _pair in rows if width >= 4),
            "bounded seed center/full equality changed")
    return {"rows": rows}


def main() -> None:
    triangular = audit_triangular_bijections()
    aliases = audit_sideways_and_aliases()
    prefix = audit_prefix_boundary()
    orbits = audit_single_seed_orbits()

    rule30_seed = infinite_single_seed_trace(RULE30, 63)
    rule60_seed = infinite_single_seed_trace(RULE60, 63)
    require(rule60_seed == (1,) * 64, "Rule 60 single-seed center is not constant")
    require(rule30_seed != rule60_seed, "single-seed hostile collapsed")

    semantic = {
        "rules": {
            "left_permutive_eca": LEFT_PERMUTIVE_RULES,
            "30": "left xor (center or right)",
            "60": "left xor center",
        },
        "triangular": triangular,
        "aliases": aliases,
        "prefix": prefix,
        "orbits": orbits,
        "rule30_seed_64": "".join(map(str, rule30_seed)),
        "rule60_seed_64": "".join(map(str, rule60_seed)),
        "scope": (
            "finite exact controls and elementary reconstruction only; "
            "no Rule 30 prize problem, SOP_2/SOP_3 interpretation, randomness, "
            "frequency limit, or time lower bound"
        ),
    }
    semantic_sha256 = sha256(
        json.dumps(semantic, sort_keys=True, separators=(",", ":")).encode("ascii")
    ).hexdigest()
    script_sha256 = sha256(Path(__file__).read_bytes().replace(b"\r\n", b"\n")).hexdigest()
    if EXPECTED_SEMANTIC_SHA256 != "TO_BE_PINNED":
        require(semantic_sha256 == EXPECTED_SEMANTIC_SHA256,
                ("semantic hash mismatch", semantic_sha256, EXPECTED_SEMANTIC_SHA256))

    print("Rule 30 left-permutive trace/cylinder exact obstruction")
    print("status=FINITE-EXACT controls for an elementary proof candidate; THM-3456 RESERVED pending audit")
    print("local_rule30=left XOR (center OR right)")
    print("sideways_inverse=left=next_center XOR (center OR right)")
    print("local_rule60=left XOR center")
    print("binary_left_permutive_ECA=" + ",".join(map(str, LEFT_PERMUTIVE_RULES)))
    print("free_trace_theorem=fixed x_0..x_n makes (x_-1..x_-n)->(X_1(0)..X_n(0)) a triangular bijection")
    print("infinite_corollary=fixed right half-line gives a compatible homeomorphism from left wings to all binary temporal traces")
    print(f"triangular_all_16_rules_depth=1..{ALL_RULE_TRIANGULAR_DEPTH}")
    print(f"triangular_rule30_rule60_depth=1..{TRIANGULAR_DEPTH}")
    print(f"triangular_fixed_right_maps={triangular['fixed_right_maps']}")
    print(f"triangular_trace_evaluations={triangular['trace_evaluations']}")
    print(f"triangular_inverse_evaluations={triangular['inversion_evaluations']}")
    print("triangular_rows=depth:fixed_right_maps:trace_evaluations")
    for depth, maps, evaluations in triangular["depth_rows"]:
        print(f"  {depth}:{maps}:{evaluations}")
    print(f"rule30_single_seed_center_64={semantic['rule30_seed_64']}")
    print(f"rule60_single_seed_center_64={semantic['rule60_seed_64']}")
    print("free_language_hostile=Rule30 and Rule60 both realize every free trace, but selected single-seed centers differ")
    print("finite_cylinder_two_trace_theorem=two adjacent temporal columns reconstruct the whole cylinder")
    print("finite_cylinder_period_corollary=on a periodic orbit, adjacent-pair trace period equals full orbit period")
    print("periodic_alias_rows=width:periodic_states:one_column_traces:nonconstant_alias_classes")
    for width, states, traces, alias_classes in aliases["width_rows"]:
        print(f"  {width}:{states}:{traces}:{alias_classes}")
    print(f"first_nonconstant_one_column_alias_width={aliases['first_nonconstant_alias_width']}")
    print(f"width7_alias_trace={''.join(map(str, aliases['alias_trace']))}")
    print(f"width7_alias_states={','.join(aliases['alias_states'])}")
    for index, cycle in enumerate(aliases["alias_cycles"], start=1):
        print(f"width7_alias_cycle_{index}={'|'.join(cycle)}")
    print(f"width7_proper_period_state={aliases['proper_period_state']}")
    print("width7_proper_period_full:center:pair="
          + ":".join(map(str, aliases["proper_period_full_center_pair"])))
    print("cylinder_prefix_theorem=width n seed cylinder equals infinite single-seed center for 0<=t<n")
    print(f"cylinder_prefix_widths_checked=1..{prefix['checked_widths']}")
    print("cylinder_prefix_sharp_hostile=width 3 first differs at t=3")
    print("cylinder_prefix_first_difference_exception_to_n=width 8 first differs at t=10")
    print("seed_orbits=width:transient:full_period:center_period:adjacent_pair_period")
    for row in orbits["rows"]:
        print("  " + ":".join(map(str, row)))
    print("bounded_seed_observation=center_period=full_period for widths 4..24 (FINITE-EXACT only)")
    print("SOP_tree_boundary=trace cylinders form a full binary prefix tree as a set system; model-theoretic SOP_2 needs a specified uniform first-order coding")
    print("prize1_boundary=finite cylinders are eventually periodic and retain only an n-step guaranteed seed prefix; no infinite-tail inference")
    print("prize2_boundary=free trace balance is universal under left permutivity and says nothing about selected-seed density")
    print("prize3_boundary=existence/unique inversion is not an algorithmic lower bound for the selected nth bit")
    print(f"script_sha256_lf={script_sha256}")
    print(f"semantic_sha256={semantic_sha256}")
    print("scope=no_Rule30_prize_solution;no_eventual_nonperiodicity;no_density_limit;no_O(n)_lower_bound;no_natural_SOP_interpretation")


if __name__ == "__main__":
    main()
