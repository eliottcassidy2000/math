#!/usr/bin/env python3
"""Exact small controls and complete bounded O/Z image of the H13 switch."""
import hashlib
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
COMPILER = HERE.parents[2] / "session251_creative_completion_20261111/ternary/compiler.py"
spec = importlib.util.spec_from_file_location("inherited_compiler", COMPILER)
cp = importlib.util.module_from_spec(spec)
spec.loader.exec_module(cp)
Q = [18, 8, 3, 5, 14]


def delete(rows, v, cut):
    vs = [u for u in range(len(rows)) if u != v]
    out = [sum(((rows[u] >> w) & 1) << j for j, w in enumerate(vs)) for u in vs]
    newcut = sum(((cut >> w) & 1) << j for j, w in enumerate(vs))
    return out, newcut


def word_state(h, B):
    state = (h, 7, 45, 13)
    for i in range(B.bit_length()):
        D, E, U, d = state
        state = (D + U - E, d, 2 * U - d, d) if B >> i & 1 else (D, 0, 2 * U, d)
    return state


def main():
    signatures = {str(cut): cp.signature(Q, cut) for cut in (6, 3)}
    table = {}
    for cut in (6, 3):
        table[str(cut)] = []
        for v in range(5):
            rows, newcut = delete(Q, v, cut)
            d, kp, km, ip, im, op, om, jp, jm = cp.signature(rows, newcut)
            table[str(cut)].append([d, op, -om, jp + jm])
    controls = []
    for B in range(16):
        record = {"B": B, "sides": []}
        for cut, seed in ((6, 31), (3, 33)):
            rows = cp.construct(Q, cut, 0, B)
            hp = cp.hp_subsets(rows)
            full = (1 << len(rows)) - 1
            deletions = [hp[full ^ (1 << v)] for v in range(len(rows))]
            predicted = 13 * cp.g(B) + 32 * B + 6 * (B & 1) + (18 if cut == 6 else 20)
            if hp[-1] != predicted or word_state(seed, B)[0] != predicted:
                raise RuntimeError((B, cut, "count"))
            if B >= 2 and min(deletions) != 13:
                raise RuntimeError((B, cut, "deletion floor"))
            record["sides"].append({"cut": cut, "H": hp[-1], "deletions": deletions,
                                     "root_state": word_state(seed, B)})
        controls.append(record)
    for A in range(4):
        for B in range(4):
            delta = cp.count(signatures["3"], A, B) - cp.count(signatures["6"], A, B)
            if delta != 3 * A - (A & 1) + 2:
                raise RuntimeError((A, B, "two-flank switch"))
    bound = 1_000_000
    cap = (bound - 31) // 32
    image = set()
    for B in range(cap + 1):
        lower = 13 * cp.g(B) + 32 * B + 6 * (B & 1) + 18
        for value in (lower, lower + 2):
            if value <= bound:
                image.add(value)
    runs = []
    for value in sorted(image):
        if runs and value == runs[-1][-1] + 2:
            runs[-1].append(value)
        else:
            runs.append([value])
    longest = max(map(len, runs))
    out = {"status": "Finite controls; all-height claims proved separately in PROOF.md.",
           "compiler_sha256": hashlib.sha256(COMPILER.read_bytes()).hexdigest(),
           "signatures": signatures, "core_deletion_coefficients_d_o_beta_c": table,
           "small_graph_controls_B_0_to_15": controls,
           "complete_OZ_image_through_target": {"bound": bound, "address_cap": cap,
             "values": len(image), "fraction_of_positive_odds": len(image) / (bound // 2),
             "longest_consecutive_odd_run": longest,
             "first_longest_run": next(run for run in runs if len(run) == longest)}}
    (HERE / "receipt.json").write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps({k: v for k, v in out.items() if k != "small_graph_controls_B_0_to_15"}, indent=2))


if __name__ == "__main__":
    main()
