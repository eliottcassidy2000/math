"""Exact transitive-clone spectra on the inherited ordinary catalog.

Run coefficients use squarefree products; independent substitution counts
use a direct Held--Karp path DP on the expanded graph.
"""
from pathlib import Path
from collections import Counter
from math import factorial
import argparse
import hashlib
import json

HERE = Path(__file__).resolve().parent
PACKET = HERE.parents[1]


def need(ok, why):
    if not ok:
        raise RuntimeError(why)


def endpoint_dp(rows, allowed_starts=None):
    n = len(rows)
    size = 1 << n
    table = [[0]*n for _ in range(size)]
    allowed = range(n) if allowed_starts is None else allowed_starts
    for v in allowed:
        table[1 << v][v] = 1
    for mask in range(1, size):
        for v, count in enumerate(table[mask]):
            if count:
                available = rows[v] & ~mask
                while available:
                    bit = available & -available
                    available -= bit
                    w = bit.bit_length()-1
                    table[mask | bit][w] += count
    return table


def hamilton(rows):
    return sum(endpoint_dp(rows)[-1]) if rows else 1


def product(left, right):
    result = [0]*len(left)
    for mask in range(len(left)):
        sub = mask
        while True:
            result[mask] += left[sub]*right[mask ^ sub]
            if sub == 0:
                break
            sub = (sub-1) & mask
    return result


def run_coefficients(rows, root):
    old = [v for v in range(len(rows)) if v != root]
    sub = [sum(1 << j for j, b in enumerate(old) if rows[a] >> b & 1) for a in old]
    incoming = [i for i, a in enumerate(old) if rows[a] >> root & 1]
    outgoing = [i for i, a in enumerate(old) if rows[root] >> a & 1]
    ends = endpoint_dp(sub)
    from_out = endpoint_dp(sub, outgoing)
    size = 1 << len(old)
    left = [sum(ends[m][v] for v in incoming) for m in range(size)]
    right = [sum(from_out[m]) for m in range(size)]
    middle = [sum(from_out[m][v] for v in incoming) for m in range(size)]
    left[0] = right[0] = 1
    state = left
    coeffs = []
    for r in range(1, len(old)//2+2):
        coeffs.append(product(state, right)[-1])
        state = product(state, middle)
    need(not any(state), "two vertices consumed by every middle run")
    need(all(a >= 0 and a % factorial(r) == 0 for r, a in enumerate(coeffs)),
         "nonnegative coefficients and middle-permutation divisibility")
    return coeffs


def clones(rows, root, t):
    old = [v for v in range(len(rows)) if v != root]
    n = len(old)+t
    result = [0]*n
    for i, a in enumerate(old):
        for j, b in enumerate(old):
            if rows[a] >> b & 1:
                result[i] |= 1 << j
        for j in range(len(old), n):
            if rows[a] >> root & 1:
                result[i] |= 1 << j
            else:
                result[j] |= 1 << i
    for i in range(len(old), n):
        for j in range(i+1, n):
            result[i] |= 1 << j
    return result


def surjections(t, r):
    if r > t:
        return 0
    stirling = [1]+[0]*r
    for _ in range(t):
        stirling = [0]+[stirling[k-1]+k*stirling[k] for k in range(1, r+1)]
    return factorial(r)*stirling[r]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--catalog", type=Path, default=PACKET/"regeneration/catalog8_f126.json")
    parser.add_argument("--output", type=Path, default=HERE/"receipt.json")
    args = parser.parse_args()
    data = json.loads(args.catalog.read_text())
    need(data["counts"] == [1, 1, 2, 4, 12, 56, 456, 6880][:data["max_order"]],
         "inherited complete ordinary-class counts")
    candidates = [r for r in data["representatives"] if r["h"] in (91, 97)]
    need(len({(r["n"], r["code"]) for r in candidates}) == len(candidates), "distinct catalog keys")
    records, hits, a2hits = [], [], []
    census = Counter()
    for record in candidates:
        rows = record["rows"]
        need(hamilton(rows) == record["h"], "independent base Hamilton count")
        census[f'{record["n"]},{record["h"]},{record["strong"]}'] += 1
        for root in range(record["n"]):
            aa = run_coefficients(rows, root)
            need(aa[0] == record["h"], "one-run coefficient equals base H")
            observed2 = hamilton(clones(rows, root, 2))
            need(observed2 == aa[0]+2*aa[1], "independent literal t2 clone count")
            item = {"n": record["n"], "code": record["code"], "strong": record["strong"],
                    "rows": rows, "root": root, "A": aa, "literal_H_t2": observed2}
            if aa[1] == 63:
                for t in range(3, len(aa)+1):
                    literal = hamilton(clones(rows, root, t))
                    need(literal == sum(a*surjections(t, r) for r, a in enumerate(aa, 1)),
                         "independent full positive run spectrum")
                    item[f"literal_H_t{t}"] = literal
                a2hits.append(item)
                if not any(aa[2:]):
                    hits.append(item)
            records.append(item)
    # Inherited hostile with a genuine third run: the base-3 mode cannot be ignored.
    control = run_coefficients([20, 1, 2, 7, 14], 0)
    need(control == [9, 10, 2], "inherited genuine third-run hostile")
    for t in range(1, 6):
        need(hamilton(clones([20, 1, 2, 7, 14], 0, t)) == 2*3**t+4*2**t-5,
             "independent hostile substitution formula")
    # Accepted22 control: TT3 ports{0,2} followed by C3 ports{3,4},
    # with one prefix vertex6 and outer root7. It has the requested slope63.
    positive_rows = [190, 60, 184, 144, 160, 8, 63, 98]
    positive = run_coefficients(positive_rows, 6)
    need(positive == [89, 63, 0, 0], "attained slope63 positive control")
    for t in range(1, 6):
        need(hamilton(clones(positive_rows, 6, t)) == 63*2**t-37,
             "independent attained slope63 substitution formula")
    result = {"status": "FINITE-EXACT marked ordinary catalog scan",
              "catalog": str(args.catalog.relative_to(PACKET)), "catalog_sha256": hashlib.sha256(args.catalog.read_bytes()).hexdigest(),
              "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
              "census": dict(census), "marked_cases": len(records), "slope63_cases": a2hits,
              "uniform_hits": hits, "all_marked_coefficients": records,
              "positive_control": {"rows": positive_rows, "root": 6, "A": positive},
              "higher_run_hostile": {"rows": [20,1,2,7,14], "root": 0, "A": control},
              "scope": "Complete only through the inherited catalog order; no order9 or all-order exclusion asserted here."}
    args.output.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({"census": dict(census), "marked_cases": len(records),
                      "slope63_cases": a2hits, "uniform_hits": hits}, indent=2))


if __name__ == "__main__":
    main()
