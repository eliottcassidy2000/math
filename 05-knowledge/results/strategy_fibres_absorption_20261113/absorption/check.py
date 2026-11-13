"""Independent witness checks plus strict finite-census receipt checks."""
from itertools import permutations
from pathlib import Path
import hashlib
import json

HERE = Path(__file__).resolve().parent

def need(condition, message):
    if not condition:
        raise RuntimeError(message)

def rows(n, code):
    out = [0] * n
    bit = 0
    for j in range(1, n):
        for i in range(j):
            a, b = (i, j) if code >> bit & 1 else (j, i)
            out[a] |= 1 << b
            bit += 1
    return out

def literal_h(out):
    n = len(out)
    return sum(all(out[p[i]] >> p[i+1] & 1 for i in range(n-1))
               for p in permutations(range(n)))

def sccs(out):
    n = len(out)
    reach = [row | 1 << i for i, row in enumerate(out)]
    for k in range(n):
        for i in range(n):
            if reach[i] >> k & 1:
                reach[i] |= reach[k]
    todo = set(range(n))
    answer = []
    while todo:
        i = min(todo)
        cell = sorted(j for j in todo if reach[i] >> j & 1 and reach[j] >> i & 1)
        answer.append(cell)
        todo.difference_update(cell)
    return answer

def blocks(spectrum):
    ans = []
    for value in spectrum:
        if ans and value == ans[-1][1] + 2:
            ans[-1][1] = value
        else:
            ans.append([value, value])
    return ans

def main():
    names = ["corridor.json", "corridor-neutral2.json", "corridor-pair-unit.json",
             "corridor-triple-unit.json"]
    packets = [json.loads((HERE / name).read_text(encoding="utf-8-sig")) for name in names]
    expected_classes = [1, 1, 2, 4, 12, 56, 456]
    expected_components = [146, 40, 17, 10]
    expected_maximum = [25, 35, 59, 61]
    for name, packet, components, maximum in zip(names, packets, expected_components, expected_maximum):
        need([r["isomorphism_classes"] for r in packet["orders"]] == expected_classes, name + " orbit totals")
        need(sum(r["labelled"] for r in packet["orders"]) == 2131019, name + " universe size")
        need(packet["orders"][-1]["cumulative_components"] == components, name + " component total")
        need(packet["orders"][-1]["anchor23_max"] == maximum, name + " anchored interval")
        need(packet["orders"][-1]["anchor23_max_gap"] == 1, name + " no skipped count")
        need(packet["independent_permutation_checks"] == 1361, name + " count controls")
    for row in packets[-1]["orders"]:
        intervals = blocks(row["spectrum"])
        need(row["maximal_component_hulls"] == intervals, "three-arc graph covers precisely spectral blocks")
        need(row["cumulative_components"] == len(intervals), "no duplicate components inside a block")
    need(not any(a <= 35 < 37 <= b for a, b in packets[0]["orders"][-1]["maximal_component_hulls"]),
         "single-arc 35-to-37 hostile control")
    need(not any(a <= 59 < 61 <= b for a, b in packets[2]["orders"][-1]["maximal_component_hulls"]),
         "two-arc 59-to-61 hostile control")
    witness_codes = [1162, 33931, 1163, 33930, 32897, 37065]
    witnesses = [{"code": c, "rows": rows(7, c), "H_permutations": literal_h(rows(7, c)),
                  "sccs": sccs(rows(7, c))} for c in witness_codes]
    need([w["H_permutations"] for w in witnesses[:4]] == [25, 25, 15, 49], "neutral bridge square")
    need({w["H_permutations"] for w in witnesses[4:]} == {59, 61}, "three-arc boundary bridge")
    need((1162 ^ 33931).bit_count() == 2, "two-arc support")
    need((32897 ^ 37065).bit_count() == 3, "three-arc support")
    need(sorted(map(len, witnesses[0]["sccs"])) == [1, 6], "before singleton absorption")
    need(list(map(len, witnesses[1]["sccs"])) == [7], "after singleton absorption")
    # A separate brute-force spectrum control does not use the DP or DSU.
    brute_spectra = []
    for n in range(1, 6):
        spectrum = sorted({literal_h(rows(n, c)) for c in range(1 << (n * (n-1) // 2))})
        need(spectrum == packets[0]["orders"][n-1]["spectrum"], "independent full small-order spectrum")
        brute_spectra.append({"n": n, "spectrum": spectrum})
    debug_path = HERE / "corridor-triple-unit.unoptimized.json"
    need(debug_path.exists(), "unoptimized receipt missing")
    need(json.loads(debug_path.read_text(encoding="utf-8-sig")) == packets[-1], "optimized/unoptimized equality")
    receipt = {"status": "VERIFIED", "census_scope": "FINITE-EXACT orders 1 through 7",
               "optimized_unoptimized_identical": True, "witnesses": witnesses,
               "independent_literal_spectra": brute_spectra,
               "census_source_sha256": hashlib.sha256((HERE / "corridor.cpp").read_bytes()).hexdigest()}
    (HERE / "check.json").write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": "VERIFIED", "witness_counts": [w["H_permutations"] for w in witnesses],
                      "components": expected_components, "anchor_maxima": expected_maximum}))

if __name__ == "__main__":
    main()
