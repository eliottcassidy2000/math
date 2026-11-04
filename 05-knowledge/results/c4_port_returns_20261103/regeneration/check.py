"""Reproduce complete searches and independently replay literal controls.

No assert statements: the same checks execute under python -O.
"""
import hashlib
import itertools
import json
from pathlib import Path
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def literal_profile(rows, ports):
    """Enumerate permutations and the unique possible deleted path edge."""
    n = len(rows)
    h = p2 = f = g = j = e = s = 0
    for word in itertools.permutations(range(n)):
        bad = [i for i in range(n-1) if not (rows[word[i]] >> word[i+1]) & 1]
        if not bad:
            h += 1
            e += word[-1] in ports
            s += word[0] in ports
            cuts = range(n+1)
        elif len(bad) == 1:
            cuts = [bad[0]+1]
        else:
            continue
        for k in cuts:
            p2 += 1
            if k:
                f += word[k-1] in ports
                g += word[0] in ports
            if 0 < k < n:
                j += word[k-1] in ports and word[k] in ports
    return [h, p2, f, g, j, e, s]


def birth(rows, cut):
    n = len(rows)
    return [row | (1 << n if cut >> i & 1 else 0) for i, row in enumerate(rows)] + [((1 << n)-1) ^ cut]


def attachment(rows, ports, clock):
    n = len(rows)
    root = n + clock
    out = list(rows) + [0] * (clock+1)
    for i in range(n):
        if i in ports:
            out[i] |= 1 << root
        else:
            out[root] |= 1 << i
    for i in range(n, root):
        out[i] = ((1 << n)-1) | sum(1 << j for j in range(i+1, root))
        out[root] |= 1 << i
    return out


def hamilton(rows):
    """Forward extension of partial paths, separate from deletion recurrence."""
    n = len(rows)
    full = (1 << n)-1
    dp = [[0] * n for _ in range(1 << n)]
    for i in range(n):
        dp[1 << i][i] = 1
    for mask in range(1, full+1):
        for end, count in enumerate(dp[mask]):
            choices = rows[end] & (full ^ mask)
            while choices:
                bit = choices & -choices
                choices -= bit
                dp[mask | bit][bit.bit_length()-1] += count
    return sum(dp[full])


def main():
    controls = [
        ([6, 4, 0], {0, 1, 2}, [1, 8, 7, 7, 6, 1, 1]),
        ([30, 28, 8, 16, 4], {0, 1, 2}, [3, 48, 21, 39, 16, 1, 3]),
        ([30, 28, 24, 16, 0], {0, 2, 4}, [1, 32, 21, 21, 12, 1, 1]),
    ]
    for rows, ports, expected in controls:
        require(literal_profile(rows, ports) == expected, "Source signature failed")
    catalogs = [(7, 63, None), (8, 126, None), (7, 63, -7)]
    run_summary = []
    with tempfile.TemporaryDirectory(prefix="regeneration_audit_") as tmp:
        small = Path(tmp) / "small.exe"
        catalog = Path(tmp) / "catalog.exe"
        for source, destination in [("small_bank.cpp", small), ("catalog.cpp", catalog)]:
            subprocess.run(["g++", "-O3", "-std=c++20", str(HERE/source), "-o", str(destination)], check=True)
        for mode in ("local", "full", "two", "complete7"):
            observed = json.loads(subprocess.check_output([str(small), mode], text=True))
            saved = json.loads((HERE / f"{mode}_small.json").read_text())
            for hit in saved["hits"]:
                hit.pop("birth_cuts", None)
            require(observed == saved, f"Search receipt changed: {mode}")
            run_summary.append({"query": mode, "counts": observed["counts"]})
        for order, slope, intercept in catalogs:
            args = [str(catalog), str(order), str(slope)]
            suffix = "" if intercept is None else f"_c{intercept}"
            if intercept is not None:
                args.append(str(intercept))
            observed = json.loads(subprocess.check_output(args, text=True))
            saved = json.loads((HERE / f"catalog{order}_f{slope}{suffix}.json").read_text())
            require(observed == {key: saved[key] for key in observed}, f"Catalog changed: {args}")
            require(not observed["hits"], f"Unexpected target hit: {args}")
            run_summary.append({"order": order, "slope": slope, "intercepts": observed["target_intercepts"], "class_counts": observed["counts"], "hits": 0})
    literal_hits = 0
    for menu in ("local", "two"):
        saved = json.loads((HERE / f"{menu}_small.json").read_text())
        for hit in saved["hits"]:
            src = hit["source"]
            old = controls[1 if src == 2 else 2][0]
            built = list(old)
            for cut in hit["birth_cuts"]:
                built = birth(built, cut)
            require(built == hit["rows"], "Literal birth reconstruction failed")
            old_mask = (1 << len(old))-1
            require([row & old_mask for row in built[:len(old)]] == old, "Old suffix arcs changed")
            ports = {i for i in range(len(built)) if hit["ports_mask"] >> i & 1}
            require(literal_profile(built, ports) == hit["signature"], "Literal split profile failed")
            for clock in (1, 2):
                require(hamilton(attachment(built, ports, clock)) == hit["intercept"] + 63 * 2**clock, "Full attachment count failed")
            literal_hits += 1
    hashes = {}
    for path in sorted(HERE.iterdir()):
        if path.suffix in (".cpp", ".py", ".md") or path.name.endswith("_small.json") or path.name.startswith("catalog") and path.suffix == ".json":
            hashes[path.name] = hashlib.sha256(path.read_bytes()).hexdigest()
    receipt = {"status": "PASS", "optimized_python": not __debug__, "queries": run_summary, "literal_source_controls": len(controls), "literal_positive_birth_controls": literal_hits, "source_sha256": hashes}
    destination = HERE / ("receipt_optimized.json" if not __debug__ else "receipt.json")
    destination.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: v for k, v in receipt.items() if k != "source_sha256"}, indent=2))


if __name__ == "__main__":
    main()
