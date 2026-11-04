"""Independent bounded-ear and literal-clone interpolation audit."""
from collections import Counter
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile

HERE = Path(__file__).resolve().parent
PACKET = HERE.parent
SOURCE = PACKET / "native/clone_escape"


def need(ok, why):
    if not ok:
        raise RuntimeError(why)


def main():
    catalog_path = PACKET / "regeneration/catalog8_f126.json"
    catalog = json.loads(catalog_path.read_text())
    base = json.loads((SOURCE / "receipt.json").read_text())
    order9 = json.loads((SOURCE / "order9.json").read_text())
    parents = [r for r in catalog["representatives"] if r["n"] == 8 and r["strong"] and r["h"] <= 97]
    need(len(parents) == 135, "Parent count")
    records = base["all_marked_coefficients"] + order9["all_marked_coefficients"]
    expected_catalog_cases = {(r["n"], r["code"], root) for r in catalog["representatives"] if r["h"] in (91, 97) for root in range(r["n"])}
    need({(r["n"], r["code"], r["root"]) for r in base["all_marked_coefficients"]} == expected_catalog_cases, "Complete all-mark universe through8")
    need(len(records) == 740 and len(base["all_marked_coefficients"]) == 659, "Marked universe count")
    expected9 = {(r["code"], root) for r in order9["representatives"] for root in range(9)}
    need({(r["code"], r["root"]) for r in order9["all_marked_coefficients"]} == expected9, "Complete all-mark universe at9")
    controls = [base["positive_control"], base["higher_run_hostile"]]
    payload = [str(len(parents))] + [" ".join(map(str, r["rows"])) for r in parents]
    payload += [str(len(order9["representatives"]))] + [" ".join(map(str, r["rows"])) for r in order9["representatives"]]
    payload += [str(len(records)+len(controls))]
    for item in records + controls:
        payload += [f'{len(item["rows"])} {item["root"]} ' + " ".join(map(str, item["rows"]))]
    with tempfile.TemporaryDirectory(prefix="clone_review_") as tmp:
        exe = Path(tmp) / "review.exe"
        subprocess.run(["g++", "-std=c++20", "-O3", str(HERE / "clone_review.cpp"), "-o", str(exe)], check=True)
        result = json.loads(subprocess.check_output([str(exe)], input="\n".join(payload)+"\n", text=True))
    for key in ("parents", "births", "strong_births", "target_births"):
        need(result[key] == order9[key], f"Independent ear tally: {key}")
    need(result["class_counts"] == [2, 7], "Independent quotient classes")
    inferred = result.pop("records")
    need(len(inferred) == len(records)+2, "Returned marked cases")
    for source, observed in zip(records+controls, inferred):
        need(observed["A"] == source["A"], "Independent direct-clone interpolation")
        need(observed["h"] == source["A"][0], "Independent base count")
        if "strong" in source:
            need(observed["strong"] == source["strong"], "Independent strongness")
        if "literal_H_t2" in source:
            need(observed["h2"] == source["literal_H_t2"], "Independent literal two-clone count")
    need(all(r["A"][1] != 63 for r in inferred[:-2]), "Target slope appears")
    need(inferred[-2]["A"] == [89, 63, 0, 0], "Positive slope63 control")
    need(inferred[-1]["A"] == [9, 10, 2], "Third-run hostile")
    for name in ("receipt", "order9"):
        need((SOURCE / f"{name}.json").read_bytes() == (SOURCE / f"{name}_optimized.json").read_bytes(), "Producer normal/-O equality")
    hashes = {}
    for path in [catalog_path, SOURCE / "PROOF.md", SOURCE / "extend9.cpp", SOURCE / "check_clone.py", SOURCE / "check_order9.py", SOURCE / "receipt.json", SOURCE / "order9.json", HERE / "clone_review.cpp", Path(__file__)]:
        hashes[str(path.relative_to(PACKET))] = hashlib.sha256(path.read_bytes()).hexdigest()
    result.update(status="PASS independent clone review", marked_cases=740,
                  marked_strong_cases=sum(r["strong"] for r in inferred[:-2]),
                  full_run_vectors_recomputed=742, target_slope_cases=0,
                  producer_codes_over_32_bits=sum(r["code"] >= 2**32 for r in order9["representatives"]),
                  source_sha256=hashes)
    output = HERE / ("clone_review_optimized.json" if not __debug__ else "clone_review.json")
    output.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({key: value for key, value in result.items() if key != "source_sha256"}, indent=2))


if __name__ == "__main__":
    main()
