#!/usr/bin/env python3
"""Literal controls for the W5 two-slot bridge; no class enumeration."""
import json
from pathlib import Path

INHERITED_REF = "d0b91cae1d58fc72140f14e305a054c9d22d29c5"


def require(value, message):
    if not value:
        raise AssertionError(message)


def tt(n):
    return tuple(((1 << n) - 1) ^ ((1 << (i + 1)) - 1) for i in range(n))


def table(rows):
    n = len(rows)
    full = (1 << n) - 1
    for i in range(n):
        require(rows[i] & (1 << i) == 0, "loop")
        for j in range(i):
            require(((rows[i] >> j) ^ (rows[j] >> i)) & 1, "not a tournament")
    ends = [[0] * n for _ in range(1 << n)]
    total = [0] * (1 << n)
    total[0] = 1
    for mask in range(1, full + 1):
        if mask & (mask - 1) == 0:
            ends[mask][mask.bit_length() - 1] = 1
        else:
            for j in range(n):
                if mask >> j & 1:
                    before = mask ^ (1 << j)
                    ends[mask][j] = sum(ends[before][i] for i in range(n)
                                           if rows[i] >> j & 1)
        total[mask] = sum(ends[mask])
    return total, ends


def state(rows, root):
    total, ends = table(rows)
    n, full = len(rows), (1 << len(rows)) - 1
    converse = tuple(full ^ (1 << i) ^ row for i, row in enumerate(rows))
    _, opposite_ends = table(converse)
    f = sum(ends[mask][root] * total[full ^ mask]
            for mask in range(1, full + 1) if mask >> root & 1)
    return (total[full], opposite_ends[full][root], total[full ^ (1 << root)], f)


def p2(rows):
    total, _ = table(rows)
    full = len(total) - 1
    return sum(total[m] * total[full ^ m] for m in range(full + 1))


def strong(rows):
    n, full = len(rows), (1 << len(rows)) - 1
    for v in range(n):
        reached = 1 << v
        while True:
            after = reached
            for u in range(n):
                if reached >> u & 1:
                    after |= rows[u]
            if after == reached:
                break
            reached = after
        if reached != full:
            return False
    return True


def bridge(module):
    """A(singleton) -> M -> C(TT2) -> D(TT2), with z -> A,C."""
    blocks = (tt(1), module, tt(2), tt(2))
    rows, offset = [], 0
    order = sum(map(len, blocks))
    block_names = []
    for block in blocks:
        block_names.append(tuple(range(offset, offset + len(block))))
        later = ((1 << order) - 1) ^ ((1 << (offset + len(block))) - 1)
        rows.extend((row << offset) | later for row in block)
        offset += len(block)
    apex = len(rows)
    selected = set(block_names[0] + block_names[2])
    rows.append(sum(1 << i for i in selected))
    for i in range(apex):
        if i not in selected:
            rows[i] |= 1 << apex
    return tuple(rows), apex, block_names


def birth(rows, root, retain_old):
    n = len(rows)
    result = list(rows)
    result[root] |= 1 << n
    result.append(((1 << n) - 1) ^ (1 << root))
    return tuple(result), root if retain_old else n


def run():
    modules = [("TT1", tt(1)), ("TT2", tt(2)), ("TT3", tt(3)),
               ("TT4", tt(4)), ("C3", (2, 4, 1)), ("Q4", (2, 4, 9, 3))]
    controls = []
    for name, module in modules:
        h = table(module)[0][-1]
        p = p2(module)
        rows, root, blocks = bridge(module)
        expected = (13*p + 11*h, h, h, 26*p - h)
        actual = state(rows, root)
        require(actual == expected, (name, actual, expected))
        require(strong(rows), name)
        trace = [{"command": "source", "root": root, "state": actual}]
        for command in "NNR":
            rows, root = birth(rows, root, command == "R")
            trace.append({"command": command, "root": root, "state": state(rows, root)})
        expected = (117*p + 29*h, 39*p + 9*h, 39*p + 9*h, 117*p + 29*h)
        require(trace[-1]["state"] == expected, (name, "NNR"))
        require(strong(rows), (name, "NNR strong"))
        controls.append({"module": name, "module_rows": module, "h": h, "p2": p,
                         "trace": trace, "final_rows": rows})
    core_rows, root, blocks = bridge(tt(2))
    require(state(core_rows, root) == (63, 1, 1, 103), "63 apex")
    embeddings = []
    for j in range(1, 4):
        old, new = controls[j-1], controls[j]
        old_rows, new_rows = old["final_rows"], new["final_rows"]
        mapping = [i if i <= j else i+1 for i in range(len(old_rows))]
        for i in range(len(old_rows)):
            for k in range(len(old_rows)):
                require((old_rows[i] >> k & 1) == (new_rows[mapping[i]] >> mapping[k] & 1),
                        ("induced tail extension", j, i, k))
        old_root, new_root = old["trace"][-1]["root"], new["trace"][-1]["root"]
        require(mapping[old_root] == new_root, "retained balanced root")
        b, t, _, _ = old["trace"][-1]["state"]
        require(new["trace"][-1]["state"] == (2*b-29, 2*t-9, 2*t-9, 2*b-29),
                "one-vertex balanced recurrence")
        embeddings.append({"old_module_order": j, "new_vertex": j+1,
                           "old_vertex_embedding": mapping, "retained_root": [old_root,new_root]})
    # Direct seven-adic finite controls; the proof is valuation-based.
    clocks = []
    for e in range(1, 5):
        modulus = 7**e
        cofactor_modulus = 7 * modulus
        values = [((468 * pow(8, m, cofactor_modulus) + 29) % cofactor_modulus) // 7
                  for m in range(modulus)]
        require(sorted(values) == list(range(modulus)), ("residue clock", e))
        clocks.append({"exponent": e, "modulus": modulus, "distinct": len(set(values))})
    return {"inherited_ref": INHERITED_REF,
            "scope": "six named module controls, exact rooted graphs through order13; residue controls through7^4",
            "core63": {"rows": core_rows, "apex": root, "blocks": blocks},
            "controls": controls, "tail_growth_embeddings": embeddings, "seven_adic_controls": clocks}


if __name__ == "__main__":
    result = run()
    destination = Path(__file__).with_name("receipt.json")
    destination.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"receipt": str(destination), "module_controls": len(result["controls"]),
                      "max_residue_modulus": result["seven_adic_controls"][-1]["modulus"]}))
