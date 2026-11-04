"""Small literal check of the two uniform letters entering different suffix types."""
from pathlib import Path
from itertools import combinations
import argparse
import hashlib
import json
from check_transport import EndpointEvaluator, literal_birth, literal_flip, strong, require

HERE = Path(__file__).resolve().parent


def signature(rows, ports):
    q = EndpointEvaluator(rows)
    full = q.full
    return [q.h[-1], sum(q.h[m]*q.h[full ^ m] for m in range(full+1)),
            sum(q.F(v) for v in ports),
            sum(q.s[v][m]*q.h[full ^ m] for v in ports for m in range(full+1)),
            sum(q.J(u, v) for u in ports for v in ports),
            sum(q.e[v][-1] for v in ports), sum(q.s[v][-1] for v in ports)]


def induced(rows, names, selected):
    indices = [names.index(v) for v in selected]
    return [sum(1 << j for j, v in enumerate(indices) if rows[u] >> v & 1) for u in indices]


def source(t):
    prefix = list(range(10, 10+t))
    names = prefix + [0, 1, 2, 99]
    rows = [0]*len(names)
    for i, a in enumerate(names):
        for j in range(i+1, len(names)):
            b = names[j]
            a_beats_b = not (b == 99 and a in prefix)
            u, v = (i, j) if a_beats_b else (j, i)
            rows[u] |= 1 << v
    return rows, names, prefix


def add_named(rows, names, new, incoming):
    indices = {names.index(v) for v in incoming}
    return list(literal_birth(rows, set(range(len(rows))), indices, False)), names + [new]


def inspect(rows, names):
    q = EndpointEvaluator(rows)
    rr = names.index(99)
    start_defects = {}
    for name in (2, 4):
        if name in names:
            r = names.index(name)
            s2 = sum(q.s[r][m]*q.h[q.full ^ m] for m in range(q.full+1))
            start_defects[str(name)] = s2-q.h[-1]
    return {"H": q.h[-1], "outer_root_defect": q.F(rr)-q.h[-1],
            "start_defects": start_defects, "strong": strong(rows)}


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=HERE/"two_digit_entry.json")
    args = parser.parse_args()
    core = []
    for mask in range(8):
        ports = [i for i in range(3) if mask >> i & 1]
        k = len(ports)
        actual = signature([2, 4, 1], ports)
        require(actual == [3, 12, 3*k, 3*k, k*(k-1), k, k], "marked triangle signature")
        core.append({"ports": ports, "signature": actual})
    suffixes = {
        "source": {"rows": [6, 4, 0], "ports": [0, 1, 2],
                   "expected": [1, 8, 7, 7, 6, 1, 1]},
        "digit2": {"rows": [30, 28, 8, 16, 4], "ports": [0, 1, 2],
                   "expected": [3, 48, 21, 39, 16, 1, 3]},
        "digit4": {"rows": [30, 28, 24, 16, 0], "ports": [0, 2, 4],
                   "expected": [1, 32, 21, 21, 12, 1, 1]},
    }
    for value in suffixes.values():
        value["signature"] = signature(value["rows"], value["ports"])
        require(value["signature"] == value["expected"], "literal suffix signature")
    families = []
    for t in range(1, 9):
        rows, names, prefix = source(t)
        observed = inspect(rows, names)
        require(observed["H"] == -5+7*2**t and observed["outer_root_defect"] == 6,
                "actual source family")
        require(observed["strong"] and observed["start_defects"]["2"] == 1,
                "source start-ready port2")
        record = {"t": t, "source_rows": rows, "source_names": names, "source": observed, "branches": []}
        for digit in (2, 4):
            aa, ns = list(rows), list(names)
            first_in = prefix + [0, 1, 99] if digit == 2 else prefix + [0, 1, 2]
            aa, ns = add_named(aa, ns, 4, first_in)
            stages = [{"operation": "birth4", "incoming_names": first_in,
                       "rows": aa, "names": ns, **inspect(aa, ns)}]
            aa, ns = add_named(aa, ns, 3, prefix + [0, 1, 2, 99])
            stages.append({"operation": "birth3", "incoming_names": prefix+[0, 1, 2, 99],
                           "rows": aa, "names": ns, **inspect(aa, ns)})
            if digit == 4:
                aa = list(literal_flip(aa, ns.index(1), ns.index(99)))
                stages.append({"operation": "reverse1_to_root99", "rows": aa, "names": ns,
                               **inspect(aa, ns)})
            require(all(s["strong"] for s in stages), "all completed states strong")
            require(stages[-1]["H"] == 3*observed["H"]+digit, "uniform digit")
            require(stages[-1]["outer_root_defect"] == (16 if digit == 2 else 12), "transported defect")
            selected = induced(aa, ns, [0, 1, 2, 3, 4])
            require(selected == suffixes[f"digit{digit}"]["rows"], "exact named suffix")
            old = induced(aa, ns, names)
            changed = [[names[i], names[j]] for i, j in combinations(range(len(names)), 2)
                       if ((old[i] ^ rows[i]) >> j) & 1]
            require(changed == ([] if digit == 2 else [[1, 99]]), "old name/edge retention")
            if digit == 2:
                require(stages[-1]["start_defects"]["4"] == 0, "start-zero output port4")
            record["branches"].append({"digit": digit, "stages": stages, "changed_old_pairs": changed})
        plus2 = record["branches"][0]["stages"][-1]
        plus4 = record["branches"][1]["stages"][-1]
        sibling_rows, sibling_names = list(plus2["rows"]), plus2["names"]
        sibling = []
        for a, b in ((99, 4), (1, 99), (4, 2)):
            require(sibling_rows[sibling_names.index(a)] >> sibling_names.index(b) & 1,
                    "sibling old arc orientation")
            sibling_rows = list(literal_flip(sibling_rows, sibling_names.index(a), sibling_names.index(b)))
            state = inspect(sibling_rows, sibling_names)
            require(state["strong"], "sibling strong intermediate")
            sibling.append({"old_arc": [a, b], "rows": sibling_rows, **state})
        require(sibling_rows == plus4["rows"] and sibling_names == plus4["names"],
                "sibling final graph is literal plus4 output")
        require(sibling[-1]["H"]-plus2["H"] == 2, "sibling H shift")
        require((sibling[-1]["H"]+sibling[-1]["outer_root_defect"])-
                (plus2["H"]+plus2["outer_root_defect"]) == -2, "sibling endpoint shift")
        record["literal_sibling_path"] = sibling
        families.append(record)
    # Small symmetry controls: a terminal marked singleton after the triangle.
    singleton = []
    for k in range(4):
        sig = signature([10, 12, 9, 0], list(range(k))+[3])
        require(sig == [3, 24, 12+3*k, 3+6*k, k*(k+3), 3, k], "terminal singleton ports")
        value = sig[0]-sig[4]+2*sig[2]
        require(value % 3 == (0 if k in (0, 3) else 2), "partial ports break triple divisibility")
        singleton.append({"triangle_port_count": k, "signature": sig, "TT1_attachment_H": value})
    result = {"status": "PROVED two uniform letters on one infinite entry family; output closure OPEN",
              "triangle_all_port_subsets": core, "suffixes": suffixes,
              "literal_t1_to8": families, "terminal_singleton_controls": singleton,
              "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    args.output.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({"suffix_signatures": {k:v["signature"] for k,v in suffixes.items()},
                      "first_four": [[r["source"]["H"], *[b["stages"][-1]["H"] for b in r["branches"]]]
                                     for r in families[:4]], "terminal_singleton_controls": singleton}, indent=2))


if __name__ == "__main__":
    main()
