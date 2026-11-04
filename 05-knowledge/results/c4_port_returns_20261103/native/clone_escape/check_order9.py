"""Bounded strong-deletion extension, then independent marked clone audit."""
from pathlib import Path
from collections import Counter
import argparse
import hashlib
import json
import subprocess
import tempfile
from check_clone import need, hamilton, run_coefficients, clones, surjections

HERE = Path(__file__).resolve().parent
PACKET = HERE.parents[1]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=HERE/"order9.json")
    args = parser.parse_args()
    path = PACKET/"regeneration/catalog8_f126.json"
    data = json.loads(path.read_text())
    need(data["counts"] == [1,1,2,4,12,56,456,6880], "inherited complete catalog")
    parents = [r for r in data["representatives"] if r["n"] == 8 and r["strong"] and r["h"] <= 97]
    need(len(parents) == 135, "complete eligible parent bank")
    payload = str(len(parents))+"\n"+"\n".join(" ".join(map(str,r["rows"])) for r in parents)+"\n"
    with tempfile.TemporaryDirectory(prefix="clone_escape_") as tmp:
        exe = Path(tmp)/"extend9.exe"
        subprocess.run(["g++", "-O3", "-std=c++20", str(HERE/"extend9.cpp"), "-o", str(exe)], check=True)
        result = json.loads(subprocess.check_output([str(exe)], input=payload, text=True))
    need(result["births"] == 135*256 and result["parents"] == 135, "all256 cuts per parent")
    records, hits, a2hits = [], [], []
    census = Counter()
    for rep in result["representatives"]:
        rows = rep["rows"]
        need(hamilton(rows) == rep["h"] and rep["h"] in (91,97), "independent base H")
        census[rep["h"]] += 1
        for root in range(9):
            aa = run_coefficients(rows, root)
            need(aa[0] == rep["h"], "one-run coefficient")
            h2 = hamilton(clones(rows, root, 2))
            need(h2 == aa[0]+2*aa[1], "independent t2 clone count")
            record = {"code": rep["code"], "rows": rows, "root": root, "A": aa, "literal_H_t2": h2}
            if aa[1] == 63:
                for t in range(3, len(aa)+1):
                    actual = hamilton(clones(rows, root, t))
                    need(actual == sum(a*surjections(t,r) for r,a in enumerate(aa,1)), "independent full clone spectrum")
                    record[f"literal_H_t{t}"] = actual
                a2hits.append(record)
                if not any(aa[2:]): hits.append(record)
            records.append(record)
    result.update(status="FINITE-EXACT complete strong order9 H91/97 via proved deletion reduction",
                  census=dict(census), marked_cases=len(records), slope63_cases=a2hits, uniform_hits=hits,
                  all_marked_coefficients=records,
                  source_sha256={str(p.relative_to(PACKET)):hashlib.sha256(p.read_bytes()).hexdigest()
                                 for p in [path, Path(__file__), HERE/"extend9.cpp", HERE/"check_clone.py"]})
    args.output.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps({k:result[k] for k in ("parents","births","strong_births","target_births","census","marked_cases","slope63_cases","uniform_hits")},indent=2))


if __name__ == "__main__":
    main()
