"""Exact one-exit context compiler and nested source audit. Standard library only.

Run from any directory; --output PATH selects the JSON receipt. Every validation
uses require(), so python -O performs the same checks.
"""
from __future__ import annotations

import argparse
import hashlib
import itertools
import json
from fractions import Fraction
from pathlib import Path

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parent.parent


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def validate(rows):
    n = len(rows)
    for i, row in enumerate(rows):
        require(0 <= row < (1 << n) and not row >> i & 1, "bad row")
        for j in range(i):
            require(((row >> j) ^ (rows[j] >> i)) & 1, "not a tournament")


def induced(rows, vertices):
    return [sum(1 << j for j, v in enumerate(vertices) if rows[u] >> v & 1)
            for u in vertices]


def reverse(rows):
    full = (1 << len(rows)) - 1
    return [full ^ row ^ (1 << i) for i, row in enumerate(rows)]


def endpoint_dp(rows):
    n = len(rows)
    full = (1 << n) - 1
    ends = [[0] * n for _ in range(full + 1)]
    totals = [0] * (full + 1)
    totals[0] = 1
    for v in range(n):
        ends[1 << v][v] = 1
    for mask in range(1, full + 1):
        subtotal = 0
        for v in range(n):
            count = ends[mask][v]
            subtotal += count
            if count:
                choices = rows[v] & (full ^ mask)
                while choices:
                    bit = choices & -choices
                    choices -= bit
                    ends[mask | bit][bit.bit_length() - 1] += count
        totals[mask] = subtotal
    return totals, ends


def profile(rows, w):
    h, ends = endpoint_dp(rows)
    _, starts = endpoint_dp(reverse(rows))
    full = len(h) - 1
    return (h[full], starts[full][w], h[full ^ (1 << w)],
            sum(ends[mask][w] * h[full ^ mask] for mask in range(1, full + 1)))


def one_exit_coefficients(exterior, exits):
    """Return coefficients in the basis (h,s,q,F), with signs included."""
    h, ends = endpoint_dp(exterior)
    _, starts = endpoint_dp(reverse(exterior))
    full = len(h) - 1
    e = [sum(row[v] for v in exits) for row in ends]
    s = [sum(row[v] for v in exits) for row in starts]
    p = sum(h[mask] * s[full ^ mask] for mask in range(full + 1))
    q = sum(e[mask] * s[full ^ mask] for mask in range(full + 1))
    return (h[full], -e[full], -q, p)


def dot(a, b):
    return sum(x * y for x, y in zip(a, b))


def joint_end_start(rows, left, right):
    h, ends = endpoint_dp(rows)
    _, starts = endpoint_dp(reverse(rows))
    full = len(h) - 1
    return sum(ends[mask][left] * starts[full ^ mask][right]
               for mask in range(1, full + 1))


def flip_named(rows, names, left, right):
    result = list(rows)
    i, j = names.index(left), names.index(right)
    result[i] ^= 1 << j
    result[j] ^= 1 << i
    return result


def join_context(core, w, exterior, exits):
    n, m = len(core), len(exterior)
    result = list(core) + [row << n for row in exterior]
    for j in range(m):
        for i in range(n):
            if i == w and j in exits:
                result[i] |= 1 << (n + j)
            else:
                result[n + j] |= 1 << i
    return result


def all_tournaments(n):
    pairs = list(itertools.combinations(range(n), 2))
    for code in range(1 << len(pairs)):
        rows = [0] * n
        for k, (a, b) in enumerate(pairs):
            u, v = (a, b) if code >> k & 1 else (b, a)
            rows[u] |= 1 << v
        yield rows


def sccs(rows):
    n = len(rows)
    reach = [row | (1 << i) for i, row in enumerate(rows)]
    for k in range(n):
        for i in range(n):
            if reach[i] >> k & 1:
                reach[i] |= reach[k]
    unused = set(range(n))
    parts = []
    while unused:
        a = min(unused)
        part = [b for b in unused if reach[a] >> b & 1 and reach[b] >> a & 1]
        unused.difference_update(part)
        parts.append(part)
    original = list(parts)
    parts.sort(key=lambda p: -sum(bool(rows[p[0]] >> q[0] & 1)
                                 for q in original if q is not p))
    return parts


def zero_ready_decompositions(rows, names):
    h, ends = endpoint_dp(rows)
    full = len(h) - 1
    found = []
    for r in range(len(rows)):
        f = sum(ends[mask][r] * h[full ^ mask] for mask in range(1, full + 1))
        if f != h[full]:
            continue
        incoming = [i for i in range(len(rows)) if rows[i] >> r & 1]
        if len(incoming) != 1:
            continue
        w = incoming[0]
        others = [i for i in range(len(rows)) if i != r]
        parts = sccs(induced(rows, others))
        if len(parts) < 2 or others.index(w) not in parts[-1]:
            continue
        b = [others[i] for i in parts[-1]]
        a = [others[i] for part in parts[:-1] for i in part]
        br = induced(rows, b)
        ar = induced(rows, a)
        ah, _ = endpoint_dp(ar)
        p2 = sum(ah[m] * ah[(len(ah) - 1) ^ m] for m in range(len(ah)))
        bh, _, _, bf = profile(br, b.index(w))
        require(h[full] == ah[-1] * bh + p2 * bf, "zero-ready identity")
        found.append({"root": names[r], "inport": names[w],
                      "prefix": [names[i] for i in a],
                      "suffix": [names[i] for i in b],
                      "H": h[full], "suffix_H": bh, "suffix_F": bf,
                      "suffix_rows": br})
    return found


def nesting(rows, names):
    result = []
    while True:
        decomps = zero_ready_decompositions(rows, names)
        if not decomps:
            break
        require(len(decomps) == 1, "nested source has multiple zero-ready roots")
        d = decomps[0]
        result.append(d)
        rows, names = d["suffix_rows"], d["suffix"]
    return result, rows, names


def context_from_literal(rows, names, core_names):
    core_indices = [names.index(v) for v in core_names]
    exterior_indices = [i for i, v in enumerate(names) if v not in core_names]
    exterior_names = [names[i] for i in exterior_indices]
    exits = {name: [names[j] for j in exterior_indices if rows[names.index(name)] >> j & 1]
             for name in core_names}
    return induced(rows, core_indices), induced(rows, exterior_indices), exterior_names, exits


def rank(matrix):
    a = [[Fraction(x) for x in row] for row in matrix]
    row = 0
    for col in range(len(a[0]) if a else 0):
        pivot = next((i for i in range(row, len(a)) if a[i][col]), None)
        if pivot is None:
            continue
        a[row], a[pivot] = a[pivot], a[row]
        scale = a[row][col]
        a[row] = [x / scale for x in a[row]]
        for i in range(len(a)):
            if i != row and a[i][col]:
                scale = a[i][col]
                a[i] = [x - scale * y for x, y in zip(a[i], a[row])]
        row += 1
    return row


def independent_rows(candidates):
    selected = []
    for value in candidates:
        if rank(selected + [value]) > len(selected):
            selected.append(value)
        if len(selected) == 4:
            break
    return selected


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=HERE / "verification.json")
    args = parser.parse_args()
    literal_path = RESULTS / "zero_ready_strong_switch_20261103/literal/realized_depth3.json"
    data = json.loads(literal_path.read_text())
    source = data[0]
    rows, names = source["source_rows"], source["source_names"]
    validate(rows)
    ladder, bottom, bottom_names = nesting(rows, names)
    require([d["root"] for d in ladder] == [39, 32, 21, 33, 22], "wrong ladder")
    require(bottom_names == [6, 7, 8, 14, 31], "wrong bottom")
    source_core, exterior, exterior_names, exits = context_from_literal(rows, names, bottom_names)
    require(source_core == bottom, "core mismatch")
    require({v: e for v, e in exits.items() if e} == {31: [22, 39]}, "not one-exit")
    coefficients = one_exit_coefficients(exterior, {exterior_names.index(v) for v in exits[31]})
    bottom_profile = profile(bottom, bottom_names.index(31))
    require(dot(coefficients, bottom_profile) == 170807, "literal context evaluation")
    source_h = endpoint_dp(rows)[0][-1]
    require(source_h == 170807, "direct source count")

    suffix_names, suffix_rows = ladder[0]["suffix"], ladder[0]["suffix_rows"]
    _, be, ben, bx = context_from_literal(suffix_rows, suffix_names, bottom_names)
    bh_coeff = one_exit_coefficients(be, {ben.index(v) for v in bx[31]})
    differences = [a - b for a, b in zip(coefficients, bh_coeff)]
    require(all(x % 4 == 0 for x in differences), "nonintegral F transfer")
    bf_coeff = tuple(x // 4 for x in differences)
    require(dot(bh_coeff, bottom_profile) == 2195, "suffix h transfer")
    require(dot(bf_coeff, bottom_profile) == 42153, "suffix F transfer")
    require(bh_coeff == (53, 0, 0, 124), "suffix H compiler changed")
    require(bf_coeff == (0, 0, -487, 3280), "suffix F compiler changed")

    # All orders of the arbitrary bottom C and ALL outer prefix lengths are
    # covered by the proved inequalities 1 <= q <= h <= F and oddness of q,h.
    # The lower bound H_s >= (177+2793*2**s)*F makes the integer audit finite.
    arithmetic = {}
    for carry in (2, 4):
        target = 3 * source_h + carry
        smax = 0
        while 177 + 2793 * 2 ** (smax + 1) <= target:
            smax += 1
        integer_candidates = []
        odd_candidates = []
        boxes = []
        tested = 0
        for s in range(1, smax + 1):
            a, b = 124 + 3280 * 2 ** s, 487 * 2 ** s
            fmax = target // (a - b + 53)
            boxes.append({"s": s, "F_max": fmax})
            for f in range(1, fmax + 1):
                for h in range(1, f + 1):
                    tested += 1
                    numerator = 53 * h + a * f - target
                    if numerator % b:
                        continue
                    q = numerator // b
                    if 1 <= q <= h:
                        candidate = {"s": s, "h": h, "F": f, "q": q}
                        integer_candidates.append(candidate)
                        if h % 2 and q % 2:
                            odd_candidates.append(candidate)
        require(not odd_candidates, "unexpected realizable cone candidate")
        expected = [] if carry == 2 else [{"s": 1, "h": 73, "F": 79, "q": 20}]
        require(integer_candidates == expected, "integer-cone candidate audit changed")
        fixed_s_residues = [
            {"F": f, "h_residue_mod_1948": ((target - 13244 * f) * pow(53, -1, 1948)) % 1948}
            for f in range(39, 46)]
        require(all(r["h_residue_mod_1948"] > r["F"] for r in fixed_s_residues),
                "fixed s2 congruence exclusion failed")
        arithmetic[str(carry)] = {"target": target, "s_max": smax,
                                 "finite_boxes": boxes, "h_F_pairs_tested": tested,
                                 "integer_candidates_before_oddness": integer_candidates,
                                 "candidates_after_oddness": odd_candidates,
                                 "s2_h_residues": fixed_s_residues}

    # Exhaustive independent vertex-DP controls: all labeled C,E up to order3,
    # every named core root and every exterior exit set. Empty E is included.
    tests = 0
    coefficient_bank = []
    profile_bank = []
    for n in range(1, 4):
        for core in all_tournaments(n):
            for w in range(n):
                p = profile(core, w)
                profile_bank.append(p)
                require(1 <= p[2] <= p[0] <= p[3], "profile inequalities")
                require(p[0] % 2 == p[2] % 2 == 1, "Hamilton path parity")
                for m in range(4):
                    for e in all_tournaments(m):
                        for lmask in range(1 << m):
                            l = {v for v in range(m) if lmask >> v & 1}
                            c = one_exit_coefficients(e, l)
                            coefficient_bank.append(c)
                            actual = endpoint_dp(join_context(core, w, e, l))[0][-1]
                            require(actual == dot(c, p), "exhaustive one-exit mismatch")
                            tests += 1
    # Enumerate marked cores through five, retaining minimal same(h,F) collision
    # under this particular prescribed exterior. This is a hostile quotient test.
    seen = {}
    collision = None
    marked_cores = 0
    for n in range(1, 6):
        for core in all_tournaments(n):
            for w in range(n):
                p = profile(core, w)
                profile_bank.append(p)
                marked_cores += 1
                key = (p[0], p[3])
                value = dot(coefficients, p)
                old = seen.get(key)
                record = {"rows": core, "root": w, "profile": p, "context_H": value}
                if collision is None and old is not None and old["context_H"] != value:
                    collision = [old, record]
                seen.setdefault(key, record)
    require(collision is not None, "expected h,F quotient collision absent")
    profile_basis = independent_rows(profile_bank)
    coefficient_basis = independent_rows(coefficient_bank)
    require(len(profile_basis) == len(coefficient_basis) == 4, "four-coordinate rank lost")
    for witness in collision:
        observed = endpoint_dp(join_context(witness["rows"], witness["root"], exterior,
                              {exterior_names.index(v) for v in exits[31]}))[0][-1]
        require(observed == witness["context_H"], "hostile direct verification")

    branches = []
    maximal_run_paths = {
        0: [1, 6, 22, 7, 31, 4, 39, 50, 51, 55, 56, 32, 40, 33, 21, 54, 2, 8, 14],
        1: [2, 6, 54, 7, 1, 39, 50, 51, 14, 31, 22, 4, 33, 32, 40, 21, 8],
        4: [1, 6, 7, 31, 22, 14, 33, 32, 40, 39, 50, 51, 55, 21, 54, 4, 2, 8],
    }
    for index, branch in enumerate(data):
        tr, tn = branch["target_rows"], branch["target_names"]
        validate(tr)
        tc, te, ten, tx = context_from_literal(tr, tn, bottom_names)
        live_exits = {v: e for v, e in tx.items() if e}
        core_changes = [[bottom_names[i], bottom_names[j]]
                        for i in range(len(bottom_names)) for j in range(i + 1, len(bottom_names))
                        if ((bottom[i] ^ tc[i]) >> j) & 1]
        suffix_zero_roots = [d["root"] for d in zero_ready_decompositions(
            branch["suffix_rows"], branch["suffix_names"])]
        branch_record = {"index": index, "offset": branch["offset"],
                         "target_H": branch["target_H"], "core_edge_changes": core_changes,
                         "exit_relations": live_exits,
                         "max_core_runs_from_exit_width": len(live_exits) + 1,
                         "suffix_zero_ready_roots": suffix_zero_roots}
        if index in maximal_run_paths:
            path = maximal_run_paths[index]
            require(len(path) == len(tn) and set(path) == set(tn), "path misses names")
            require(all(tr[tn.index(u)] >> tn.index(v) & 1 for u, v in zip(path, path[1:])),
                    "non-directed maximal-run path")
            core_runs = []
            previous_core = False
            for v in path:
                if v in bottom_names:
                    if not previous_core:
                        core_runs.append([])
                    core_runs[-1].append(v)
                previous_core = v in bottom_names
            require(len(core_runs) == len(live_exits) + 1, "exit-width bound not sharp")
            branch_record.update({"sharpness_path": path, "core_runs": core_runs})
        branches.append(branch_record)

    # This is a physical-name hostile, not an equal-count replacement of the
    # actual source. The TWO sources and their SAME named raw operation are
    # both kept, with independently counted outputs. Native path guards are
    # not transplanted from one source to the other.
    hostile_bottom = flip_named(bottom, bottom_names, 6, 8)
    hostile_source = flip_named(rows, names, 6, 8)
    hostile_target = flip_named(data[0]["target_rows"], data[0]["target_names"], 6, 8)
    require(profile(hostile_bottom, bottom_names.index(31)) == bottom_profile, "one-port hostile drift")
    require(endpoint_dp(hostile_source)[0][-1] == 170807, "hostile source count")
    require(endpoint_dp(data[0]["target_rows"])[0][-1] == 512423, "positive target count")
    hostile_target_h = endpoint_dp(hostile_target)[0][-1]
    require(hostile_target_h == 654021, "hostile two-exit target count")
    two_exit_hostile = {
        "source_core_names": bottom_names,
        "source_core_rows": bottom, "hostile_core_rows": hostile_bottom,
        "same_root31_profile": bottom_profile,
        "source_H_pair": [170807, 170807], "target_H_pair": [512423, hostile_target_h],
        "same_raw_operation": {"named_edge_flips": data[0]["edits"],
                               "target_prefix_names": data[0]["prefix_names"],
                               "target_inport": data[0]["target_inport"]},
        "J31_6_pair": [joint_end_start(c, bottom_names.index(31), bottom_names.index(6))
                       for c in (bottom, hostile_bottom)],
        "J6_31_pair": [joint_end_start(c, bottom_names.index(6), bottom_names.index(31))
                       for c in (bottom, hostile_bottom)],
        "source6_profile_pair": [profile(c, bottom_names.index(6)) for c in (bottom, hostile_bottom)],
        "scope": "Same raw named edge/prefix operation; native Hamiltonian guards must be recompiled on hostile source.",
    }

    result = {
        "status": "PROVED identity with FINITE-EXACT controls; reusable carry closure OPEN",
        "coordinate_order": ["H(C)", "start_w(C)", "H(C-w)", "E_w(C,2)"],
        "source_H": source_h, "nested_ladder": ladder,
        "bottom_names": bottom_names, "bottom_rows": bottom,
        "bottom_profile": bottom_profile, "source_context_coefficients": coefficients,
        "suffix_h_coefficients": bh_coeff, "suffix_F_coefficients": bf_coeff,
        "exhaustive_context_controls": tests, "marked_cores_checked": marked_cores,
        "same_h_F_hostile": collision, "profile_rank_basis": profile_basis,
        "coefficient_rank_basis": coefficient_basis,
        "all_core_all_prefix_exclusion": arithmetic,
        "two_exit_named_port_hostile": two_exit_hostile,
        "literal_branch_exit_width_audit": branches,
        "input_sha256": {str(literal_path.relative_to(RESULTS)): hashlib.sha256(literal_path.read_bytes()).hexdigest()},
        "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: result[k] for k in ["status", "bottom_profile", "source_context_coefficients",
                                          "suffix_h_coefficients", "suffix_F_coefficients",
                                          "exhaustive_context_controls", "marked_cores_checked",
                                          "same_h_F_hostile", "literal_branch_exit_width_audit"]}, indent=2))


if __name__ == "__main__":
    main()
