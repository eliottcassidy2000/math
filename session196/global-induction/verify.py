#!/usr/bin/env python3
"""Independent controls for the Session 196 equal-response theorem."""
from collections import deque
from itertools import permutations
from pathlib import Path
import hashlib
import json


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def tournament(n, code):
    rows = [0] * n
    bit = 0
    for a in range(n):
        for b in range(a + 1, n):
            if code & (1 << bit):
                rows[a] |= 1 << b
            else:
                rows[b] |= 1 << a
            bit += 1
    return rows


def table(rows):
    n = len(rows)
    size = 1 << n
    end = [[0] * n for _ in range(size)]
    total = [0] * size
    total[0] = 1
    for mask in range(1, size):
        for v in range(n):
            if not (mask >> v) & 1:
                continue
            old = mask ^ (1 << v)
            if not old:
                end[mask][v] = 1
            else:
                end[mask][v] = sum(end[old][u] for u in range(n)
                                   if ((old >> u) & 1) and ((rows[u] >> v) & 1))
        total[mask] = sum(end[mask])
    return total, end


def profile(rows, root):
    n = len(rows)
    full = (1 << n) - 1
    for i in range(n):
        require(not ((rows[i] >> i) & 1), "loop")
        for j in range(i):
            require(((rows[i] >> j) & 1) + ((rows[j] >> i) & 1) == 1,
                    "not a tournament")
    total, end = table(rows)
    converse = [full ^ (1 << i) ^ rows[i] for i in range(n)]
    rev_total, start = table(converse)
    require(total == rev_total, "converse Hamilton counts")
    response = sum(start[mask][root] * total[full ^ mask]
                   for mask in range(1, full + 1) if (mask >> root) & 1)
    return (total[full], end[full][root], response,
            total[full ^ (1 << root)])


def literal_profile(rows, root):
    """Separate permutation/cut implementation, used only on small graphs."""
    n = len(rows)
    h = e = response = 0
    for p in permutations(range(n)):
        legal = [bool((rows[p[i]] >> p[i + 1]) & 1) for i in range(n - 1)]
        if all(legal):
            h += 1
            e += p[-1] == root
        cut = p.index(root)
        if all(legal[:max(0, cut - 1)]) and all(legal[cut:]):
            response += 1
    others = [x for x in range(n) if x != root]
    d = sum(all((rows[p[i]] >> p[i + 1]) & 1 for i in range(len(p) - 1))
            for p in permutations(others))
    return h, e, response, d


def near_birth(rows, root, move=False):
    n = len(rows)
    out = [row | (0 if i == root else (1 << n)) for i, row in enumerate(rows)]
    out.append(1 << root)
    return out, n if move else root


def sink_birth(rows, root):
    n = len(rows)
    return [row | (1 << n) for row in rows] + [0], root


def akd(state):
    h, e, R, d = state
    return ((h + R - e, h - e, h + R, h),
            (h + R - e, d, 2 * R - d, d),
            (h, 0, 2 * R, d))


def gadget(state, signature):
    h, e, R, d = state
    g, s, L, E, Q, B, delta, dg = signature
    return g * h - s * e + L * R, E * h - delta * e, Q * h + B * R, dg * h


def closure(seed, bound, bank):
    seen = {seed}
    q = deque([seed])
    values = {seed[0]}
    while q:
        state = q.popleft()
        for out in [*akd(state), *(gadget(state, b) for b in bank)]:
            h, e, R, d = out
            require(0 <= e <= d <= h <= R, "nonphysical transition inequalities")
            if h <= bound:
                values.add(h)
            if R <= bound and out not in seen:
                seen.add(out)
                q.append(out)
    return sorted(values), len(seen)


def expected(state):
    h, e, R, d = state
    require((e, R) == (d, h) and h >= 3 and h > d, "bad theorem input")
    require(h >= 3 * d, "physical equal-response ratio failed")
    result = {h, 2 * h - d, 3 * h}
    if h - 3 * d <= 4:
        result.add(4 * h - 3 * d)
    return sorted(result)


def main():
    base = Path(__file__).resolve().parent
    bank_bytes = (base / "bank67.txt").read_bytes()
    bank = [tuple(map(int, line.split())) for line in bank_bytes.decode().splitlines()]
    require(len(bank) == 67 and all(len(b) == 8 for b in bank), "frozen bank shape")
    require(len(set(bank)) == 67, "duplicate bank signature")
    state_graphs = {}
    graphs_checked = roots_checked = elementary_checks = 0
    C5_surpluses = {}
    for n in range(1, 6):
        for code in range(1 << (n * (n - 1) // 2)):
            rows = tournament(n, code)
            graphs_checked += 1
            for root in range(n):
                state = profile(rows, root)
                roots_checked += 1
                if n <= 4:
                    require(state == literal_profile(rows, root), "literal profile disagreement")
                if n <= 4:
                    A, K, D = akd(state)
                    for physical, want in ((near_birth(rows, root, True), A),
                                           (near_birth(rows, root), K),
                                           (sink_birth(rows, root), D)):
                        require(profile(*physical) == want, "physical elementary law")
                        elementary_checks += 1
                h, e, R, d = state
                if h == 5 and R - h in (1, 2):
                    C5_surpluses.setdefault(R - h, (rows, root))
                if h >= 3 and h > d and (e, R) == (d, h):
                    state_graphs.setdefault(state, (rows, root))

    # A root of H5 with response 6 or 7, plus a terminal singleton and
    # a new vertex beating only that root, gives equal-response H17/H19.
    positive = []
    for surplus in (1, 2):
        require(surplus in C5_surpluses, "missing physical response seed")
        rows, root = C5_surpluses[surplus]
        joined, oldroot = sink_birth(rows, root)
        out, newroot = near_birth(joined, oldroot, True)
        state = profile(out, newroot)
        require(state == (15 + 2 * surplus, 5, 15 + 2 * surplus, 5),
                "positive ternary-gap seed")
        state_graphs[state] = out, newroot
        kk, mark = near_birth(out, newroot)
        kk, mark = near_birth(kk, mark)
        kk_state = profile(kk, mark)
        require(kk_state[0] == 3 * state[0] + 2 * surplus, "KK positive escape")
        positive.append({"gap": 2 * surplus, "rows": out, "root": newroot,
                         "state": state, "KK_rows": kk, "KK_root": mark,
                         "KK_state": kk_state})

    hostile_rows = [122, 112, 67, 102, 108, 4, 32]
    hostile = profile(hostile_rows, 5)
    require(hostile == (49, 13, 49, 13), "H49 hostile profile")
    require(hostile == literal_profile(hostile_rows, 5), "H49 literal profile")
    state_graphs[hostile] = hostile_rows, 5

    infinite_family_controls = []
    for m in range(3, 9):
        n = m + 1
        chain = [sum(1 << j for j in range(i + 1, n)) for i in range(n)]
        rows, root = near_birth(chain, 0, True)
        state = profile(rows, root)
        require(state == (2**m + 1, 1, 2**m + 1, 1), "primitive hostile family")
        require(3 * state[0] + 2 not in expected(state), "lower child hostile")
        require(3 * state[0] + 4 not in expected(state), "upper child hostile")
        if m == 4:
            state_graphs[state] = rows, root
        infinite_family_controls.append(state)

    results = []
    for state, (rows, root) in sorted(state_graphs.items()):
        h, e, R, d = state
        require(h >= 3 * d, "ratio")
        actual, retained = closure(state, 3 * h + 4, bank)
        require(actual == expected(state), "complete finite-bank support disagrees")
        full_actual, full_retained = closure(state, 5 * h, bank)
        require(full_actual == sorted({h, 2*h-d, 3*h, 4*h-3*d, 5*h}),
                "five-times window support disagrees")
        kk, mark = near_birth(rows, root)
        kk, mark = near_birth(kk, mark)
        require(profile(kk, mark)[0] == 4 * h - 3 * d, "literal KK law")
        results.append({"seed": state, "ceiling": 3 * h + 4,
                        "values": actual, "retained_states": retained,
                        "five_times_ceiling": 5*h,
                        "five_times_values": full_actual,
                        "five_times_retained_states": full_retained})

    diagonal_controls = []
    for original in ([0], [2, 4, 1], [2, 12, 9, 17, 7]):
        rows, root = sink_birth(original, 0)
        root = len(rows) - 1
        state = profile(rows, root)
        h = state[0]
        require(state == (h, h, h, h), "physical diagonal state")
        values, retained = closure(state, 5*h, bank)
        require(values == [h, 3*h, 5*h], "diagonal full window")
        diagonal_controls.append({"state": state, "values": values,
                                  "retained_states": retained})

    # Exact symbolic-integer corner controls cover very large parameters.
    symbolic = 0
    for d in (1, 3, 5, 101, 10**20 + 1):
        for gap in (0, 2, 4, 6, 100):
            h = 3 * d + gap
            if h < 5:
                continue
            state = (h, d, h, d)
            A, K, D = akd(state)
            require(A[0] == K[0] == 2 * h - d, "first symbolic step")
            require(akd(A)[0][0] == akd(A)[1][0] == 3 * h, "AA/AK symbolic")
            require(akd(K)[0][0] == akd(K)[1][0] == 4 * h - 3 * d,
                    "KA/KK symbolic")
            for second in (*akd(A)[:2], *akd(K)[:2]):
                require(akd(second)[0][0] >= 5 * h, "third-step exclusion")
            symbolic += 1

    receipt = {"status": "PASS", "scope": "all-order proof with finite controls only",
               "labelled_graphs_through_five": graphs_checked,
               "rooted_graphs_through_five": roots_checked,
               "physical_elementary_checks": elementary_checks,
               "bank_sha256": hashlib.sha256(bank_bytes).hexdigest(),
               "equal_response_states": results,
               "positive_escape_witnesses": positive,
               "infinite_primitive_family_controls": infinite_family_controls,
               "physical_diagonal_controls": diagonal_controls,
               "symbolic_parameter_controls": symbolic,
               "claim": "Exact descendants through 5h from a fixed equal-response input; not global spectrum coverage"}
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
