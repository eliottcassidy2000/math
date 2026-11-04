"""All-clock sibling guards, with a separate elementary native replay.

No asserts: normal and -O execute the same validity checks.
"""
from pathlib import Path
from itertools import permutations, combinations
import argparse
import hashlib
import json
import sys

HERE = Path(__file__).resolve().parent
RESULTS = HERE.parents[1]
OLD = RESULTS / "multiport_defect_grammar_20261103/transport"
sys.path.insert(0, str(OLD))
from check_transport import EndpointEvaluator, strong
from check_two_digit_entry import source, add_named, signature


def need(ok, why):
    if not ok:
        raise RuntimeError(why)


def legal(rows, names, word):
    return sorted(word) == sorted(names) and all(
        rows[names.index(a)] >> names.index(b) & 1 for a, b in zip(word, word[1:]))


def flip(rows, names, a, b):
    rows = list(rows)
    i, j = names.index(a), names.index(b)
    rows[i] ^= 1 << j
    rows[j] ^= 1 << i
    return rows


def native(rows, names, word):
    """Elementary pairwise interpreter, independent of the guard compiler."""
    need(legal(rows, names, word), "actual spanning directed guard")
    kept = {frozenset((a, b)) for a, b in zip(word, word[1:])}
    answer = list(rows)
    for a, b in combinations(names, 2):
        if frozenset((a, b)) not in kept:
            answer = flip(answer, names, a, b)
    need(legal(answer, names, word), "native guard survives its own move")
    return answer


def frame(word, v):
    j = word.index(v)
    return (word[j-1] if j else None, word[j+1] if j+1 < len(word) else None)


def find_guard(rows, names, a, b):
    for word in permutations(names):
        if (a, b) in zip(word, word[1:]) and legal(rows, names, word):
            return list(word)
    return None


def branch2(t):
    rows, names, prefix = source(t)
    rows, names = add_named(rows, names, 4, prefix+[0, 1, 99])
    rows, names = add_named(rows, names, 3, prefix+[0, 1, 2, 99])
    return rows, names, prefix


def branch4(t):
    rows, names, prefix = source(t)
    rows, names = add_named(rows, names, 4, prefix+[0, 1, 2])
    rows, names = add_named(rows, names, 3, prefix+[0, 1, 2, 99])
    return flip(rows, names, 1, 99), names


def padded_program(rows, names, word, a, b, helper_base):
    """Universal terminal pad makes every non-incident old frame common."""
    need(legal(rows, names, word), "source guard")
    need((a, b) in zip(word, word[1:]), "guard contains requested oriented arc")
    z, h = helper_base, helper_base+1
    need(z not in names and h not in names, "fresh helper names")
    k = word.index(a)
    return {
        "old_arc": [a, b],
        "births": [{"name": z, "incoming": list(names)},
                   {"name": h, "incoming": [v for v in names+[z] if v != a]}],
        "guards": [word+[z, h], word[:k+1]+[h]+word[k+1:]+[z]],
        "deletions": [h, z],
    }


def replay(rows, names, program):
    old_rows, old_names = list(rows), list(names)
    for item in program["births"]:
        rows, names = add_named(rows, names, item["name"], item["incoming"])
    p, q = program["guards"]
    first = native(rows, names, p)
    second = native(first, names, q)
    a, b = program["old_arc"]
    need(all(frame(p, v) == frame(q, v) for v in old_names if v not in (a, b)),
         "all non-incident old frames agree, including every prefix vertex")
    need(program["deletions"] == list(reversed(names[len(old_names):])), "only fresh helpers deleted")
    full = (1 << len(old_names))-1
    result = [r & full for r in second[:len(old_names)]]
    need(result == flip(old_rows, old_names, a, b), "exactly the requested old pair changes")
    return result


def sidecar(rows, names):
    q = EndpointEvaluator(rows)
    ports = [99, 0, 1, 2, 3, 4]
    h = q.h[-1]
    golden = {}
    for port in ports:
        r = names.index(port)
        s2 = sum(q.s[r][m]*q.h[q.full ^ m] for m in range(q.full+1))
        e = q.e[r][-1]
        need((s2-e) % 2 == 0, "integral golden half-difference")
        x, y = (s2-e)//2, (s2+e)//2
        a, b, c = (x+4*(h-y)) % 7, y % 7, (h-y) % 7
        field343 = [h % 7, x % 7, y % 7]
        golden[str(port)] = {"H_X_Y": [h, x, y], "F49_plus_charge": [a, b, c],
                             "F49_frobenius_plus_charge": [(a+4*b) % 7, -b % 7, c],
                             "F343_vector": field343,
                             "F343_frobenius": [(h+4*x+6*y) % 7, 2*x % 7, (3*x+4*y) % 7]}
    pairs = {f"{a},{b}": q.J(names.index(a), names.index(b))
             for a in ports for b in ports if a != b}
    return {"H": h, "E_root_2": q.F(names.index(99)), "golden": golden,
            "ordered_two_component_endpoint_J": pairs}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=HERE/"receipt.json")
    args = parser.parse_args()
    arcs = [(99, 4), (1, 99), (4, 2)]
    suffix_rows = [30, 28, 8, 16, 4]
    suffixes = [(suffix_rows, [0, 1, 2]), (suffix_rows, [0, 1, 2, 4]),
                (suffix_rows, [0, 2, 4]), ([30, 28, 24, 16, 0], [0, 2, 4])]
    sigma = [signature(rows, ports) for rows, ports in suffixes]
    need([[s[0]-s[4], s[2], s[4]] for s in sigma] ==
         [[-13, 21, 16], [-25, 33, 28], [-11, 27, 14], [-11, 21, 12]],
         "all four exact component signatures")
    quotient, qnames, _ = branch2(1)
    templates = []
    for a, b in arcs:
        word = find_guard(quotient, qnames, a, b)
        need(word is not None, "forward quotient adjacency exists")
        templates.append(word)
        quotient = flip(quotient, qnames, a, b)
    records = []
    for t in range(1, 33):
        rows, names, prefix = branch2(t)
        item = {"t": t, "names": names, "source_rows": rows, "steps": []}
        if t <= 6:
            item["source_sidecar"] = sidecar(rows, names)
            need(item["source_sidecar"]["H"] == -13+21*2**t, "source all-clock formula")
        for number, ((a, b), template) in enumerate(zip(arcs, templates)):
            word = [v for x in template for v in (prefix if x == 10 else [x])]
            program = padded_program(rows, names, word, a, b, 1000+2*number)
            rows = replay(rows, names, program)
            need(strong(rows), "completed sibling state is strong")
            record = {"program": program, "rows": rows}
            if t <= 6:
                record["sidecar"] = sidecar(rows, names)
                h, _, f, _, j, _, _ = sigma[number+1]
                need(record["sidecar"]["H"] == h-j+2**t*f and
                     record["sidecar"]["E_root_2"]-record["sidecar"]["H"] == j,
                     "every completed state's exact clock formula")
            item["steps"].append(record)
        target, target_names = branch4(t)
        need(rows == target and names == target_names, "actual same-name digit4 target")
        if t <= 6:
            end = item["steps"][-1]["sidecar"]
            need(end["H"] == -11+21*2**t, "target all-clock formula")
            need(end["E_root_2"]-item["source_sidecar"]["E_root_2"] == -2,
                 "endpoint response decreases by two")
        records.append(item)
    # Intrinsic forward guarding is not automatic: the TT3 long arc fails.
    tt3, ns = [6, 4, 0], [0, 1, 2]
    need(find_guard(tt3, ns, 0, 2) is None, "TT3 long-arc hostile")
    # The padded formula also covers an endpoint arc and n=2 without exceptions.
    for small_rows, small_names, word, arc in (([2, 0], [0, 1], [0, 1], (0, 1)),
                                             (tt3, ns, [0, 1, 2], (1, 2))):
        replay(small_rows, small_names, padded_program(small_rows, small_names, word, *arc, 9000))
    pins = [Path(__file__), OLD/"check_two_digit_entry.py", OLD/"check_transport.py"]
    report = {"status": "PROVED uniform native sibling transport; FINITE-EXACT literal replay",
              "quotient_guard_templates": templates,
              "suffix_signatures": sigma,
              "finite_universe": {"clock_values": list(range(1, 33)), "arc_gadgets": 96,
                                  "native_complements": 192, "fresh_sidecar_clocks": list(range(1, 7))},
              "scope": "No output closure, C4 repair, fixed-order black connectivity, or physical Frobenius transport claimed.",
              "records": records,
              "sha256": {str(p.relative_to(RESULTS)): hashlib.sha256(p.read_bytes()).hexdigest() for p in pins}}
    args.output.write_text(json.dumps(report, indent=2)+"\n")
    print(json.dumps({"templates": templates, "finite_universe": report["finite_universe"],
                      "first_source_target_H": [[r["source_sidecar"]["H"], r["steps"][-1]["sidecar"]["H"]]
                                                for r in records[:6]]}, indent=2))


if __name__ == "__main__":
    main()
