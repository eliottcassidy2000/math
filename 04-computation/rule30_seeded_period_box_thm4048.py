#!/usr/bin/env python3
"""Finite exact Prize-1 controls for THM-4048.

This audits the packed/direct center convention, the first inverse-boundary
defect, every candidate period through 32768 on the first 65536 center bits,
and the inherited dyadic wrap words.  It proves no all-time prize statement.
"""

from hashlib import sha256


MAX_TIME = 65535
MAX_PERIOD = 32768
DIRECT_TIME = 512


def require(condition, label):
    if not condition:
        raise RuntimeError(label)


def packed_prefix(max_time):
    packed = 1
    center_word = 0
    power_rows = {}
    for time in range(max_time + 1):
        center_word |= ((packed >> time) & 1) << time
        if time and time & (time - 1) == 0:
            power_rows[time] = packed
        packed = packed ^ ((packed << 1) | (packed << 2))
    return center_word, power_rows


def direct_center(max_time):
    state = {0}
    answer = 0
    for time in range(max_time + 1):
        answer |= int(0 in state) << time
        lower = min(state) - 1
        upper = max(state) + 1
        state = {
            site
            for site in range(lower, upper + 1)
            if int(site - 1 in state)
            ^ (int(site in state) | int(site + 1 in state))
        }
    return answer


def evolve_center(initial, time):
    state = set(initial)
    for _ in range(time):
        lower = min(state) - 1 if state else -1
        upper = max(state) + 1 if state else 1
        state = {
            site
            for site in range(lower, upper + 1)
            if int(site - 1 in state)
            ^ (int(site in state) | int(site + 1 in state))
        }
    return int(0 in state)


def compile_left_boundary(target):
    """THM-3456 binary triangular compiler with positive boundary zero."""
    initial = {0} if target[0] else set()
    boundary = []
    for depth in range(1, len(target)):
        choices = []
        for bit in (0, 1):
            trial = set(initial)
            if bit:
                trial.add(-depth)
            if evolve_center(trial, depth) == target[depth]:
                choices.append(bit)
        require(len(choices) == 1, ("unique inverse bit", depth))
        boundary.append(choices[0])
        if choices[0]:
            initial.add(-depth)
    return boundary


def main():
    center, power_rows = packed_prefix(MAX_TIME)
    direct = direct_center(DIRECT_TIME)
    require(
        center & ((1 << (DIRECT_TIME + 1)) - 1) == direct,
        "packed/direct center indexing",
    )
    center_bits = tuple((center >> index) & 1 for index in range(MAX_TIME + 1))

    # If e is the first failed relation c_(e+p)=c_e after N, the seeded
    # periodic candidate first differs at d=e+p.  Audit that load-bearing
    # index directly through the triangular inverse compiler.
    for onset in range(9):
        for period in range(1, 9):
            relation_failure = next(
                index
                for index in range(onset, MAX_TIME - period + 1)
                if center_bits[index] != center_bits[index + period]
            )
            defect_depth = relation_failure + period
            target = list(center_bits[:onset])
            block = center_bits[onset : onset + period]
            while len(target) <= defect_depth:
                target.append(block[(len(target) - onset) % period])
            boundary = compile_left_boundary(target)
            require(
                all(bit == 0 for bit in boundary[: defect_depth - 1]),
                ("zero boundary before defect", onset, period),
            )
            require(
                boundary[defect_depth - 1] == 1,
                ("first boundary defect", onset, period),
            )

    last_mismatch = []
    first_mismatch = []
    maximum_equal_run = []
    for period in range(1, MAX_PERIOD + 1):
        comparison_count = MAX_TIME - period + 1
        comparison_mask = (1 << comparison_count) - 1
        difference = (center ^ (center >> period)) & comparison_mask
        require(difference, ("finite mismatch", period))
        first_mismatch.append((difference & -difference).bit_length() - 1)
        last_mismatch.append(difference.bit_length() - 1)
        equal = (~difference) & comparison_mask
        run = 0
        while equal:
            equal &= (equal << 1) & comparison_mask
            run += 1
        maximum_equal_run.append(run)

    box_floor = min(last_mismatch)
    require(box_floor == 32767, "period-box onset floor")
    require(all(last >= box_floor for last in last_mismatch), "period box")
    global_equal_run = max(maximum_equal_run)
    require(global_equal_run == 29, "maximum equality run")
    global_equal_periods = tuple(
        period
        for period, run in enumerate(maximum_equal_run, start=1)
        if run == global_equal_run
    )

    innovation_rows = []
    for power, packed in sorted(power_rows.items()):
        displacement = packed - 1
        valuation = (displacement & -displacement).bit_length() - 1
        wrap_length = max(valuation - power + 1, 0)
        innovation_rows.append((power.bit_length() - 1, power, valuation, wrap_length))
        if wrap_length:
            forced = (center >> power) & ((1 << wrap_length) - 1)
            require(forced == 1 << (wrap_length - 1), ("wrap word", power))

    center_bytes = center.to_bytes((MAX_TIME + 8) // 8, "little")
    period_payload = ";".join(
        f"{period}:{first_mismatch[period-1]}:{last_mismatch[period-1]}:"
        f"{maximum_equal_run[period-1]}"
        for period in range(1, MAX_PERIOD + 1)
    ).encode("ascii")
    dyadic_periods = tuple(1 << exponent for exponent in range(16))
    dyadic_last = tuple((period, last_mismatch[period - 1]) for period in dyadic_periods)
    minimum_periods = tuple(
        period
        for period, last in enumerate(last_mismatch, start=1)
        if last == box_floor
    )

    print("RULE30_SEEDED_PERIOD_BOX_THM4048")
    print("universe=time_0..65535;candidate_periods=1..32768")
    print("center_index=c_t=bit_t(R_t);R_0=1;R_(t+1)=R_t_xor((2R_t)_or_(4R_t))")
    print(f"direct_CA_crosscheck=0..{DIRECT_TIME};PASS=True")
    print("inverse_boundary_index_crosscheck=N_0..8,p_1..8;first_beta_defect=e+p;PASS=True")
    print(f"center_ones={center.bit_count()};center_prefix_sha256={sha256(center_bytes).hexdigest()}")
    print(
        "period_box=no_pair_(N,p)_with_0<=N<=32767_and_1<=p<=32768;"
        f"minimum_last_mismatch={box_floor};attained_by={minimum_periods}"
    )
    print(
        "shift_equality=max_equal_run_29_over_every_tested_shift;"
        f"attained_by_periods={global_equal_periods}"
    )
    print(f"dyadic_period_last_mismatches={dyadic_last}")
    print(f"innovation_rows_(m,q,v_m,wrap_length)={tuple(innovation_rows)}")
    print(f"period_mismatch_stream_sha256={sha256(period_payload).hexdigest()}")
    print("scope=FINITE-EXACT only;no eventual nonperiodicity or prize claim")
    print("RESULT=PASS")


if __name__ == "__main__":
    main()
