"""Replay the incoming all-order D2 residual proof without editing its packet.

Uses the incoming independent verifier, not a third independent implementation.
Reconstructs both affine spaces and checks every compressed producer certificate.
Builds and raw streams live in a fresh system temporary directory.
"""
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil
import subprocess
import tempfile

HERE = Path(__file__).resolve().parent
PACKET = HERE.parent / "deficit_two_centers_20261020"
SOURCE = PACKET / "session15-deficit-centers"
AUDIT = PACKET / "session15-deficit-census-audit"


def need(value, message):
    if not value:
        raise RuntimeError(message)


def digest(data):
    return hashlib.sha256(data).hexdigest()


def main():
    module_path = AUDIT / "reconstruct_spaces.py"
    spec = importlib.util.spec_from_file_location("incoming_affine_audit", module_path)
    audit = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(audit)
    compiler = shutil.which("g++")
    need(compiler, "g++ is required to replay the C++ certificates")
    scratch = Path(tempfile.mkdtemp(prefix="math-scc-d2-replay-"))
    executable = scratch / "independent_scan.exe"
    cpp = AUDIT / "independent_scan.cpp"
    subprocess.run([compiler, "-O3", "-std=c++20", str(cpp), "-o", str(executable)],
                   check=True, capture_output=True, text=True)
    specs = json.loads((SOURCE / "residual_frame_specs.json").read_text())
    by_family = {(s["n"], s["q"], s["low_out_orbit"]): s for s in specs}
    index = json.loads((SOURCE / "witness_certificate_index.json").read_text())
    by_name = {Path(s["file"]).name: s for s in index}
    rows = []
    total = 0
    for n, q, low in audit.FAMILIES:
        model = audit.affine(n, q, low)
        if model is None:
            need((n, q, low) not in by_family, "inconsistent family retained")
            rows.append(dict(n=n, q=q, low=low, consistent=False))
            continue
        pairs, equations, base, directions, free = model
        producer = by_family[n, q, low]
        need(producer["anti"] == (q in (3, 4)), "incorrect reflection assumption")
        pbase, pdirs = audit.producer_model(producer, pairs, equations, len(free))
        stem = f"n{n}_q{q}_low{low}"
        output_model = scratch / (stem + ".model")
        audit.file_model(output_model, n, q, low, pairs, pbase, pdirs)
        compressed_path = SOURCE / f"frames_{stem}.txt.witnesses.bin.gz"
        compressed = compressed_path.read_bytes()
        raw = gzip.decompress(compressed)
        retained = by_name[compressed_path.name]
        count = 1 << len(free)
        need(retained["records"] == count and len(raw) == 5 * count,
             "certificate count mismatch")
        need(digest(raw) == retained["raw_sha256"], "raw certificate hash mismatch")
        need(digest(compressed) == retained["gzip_sha256"], "gzip hash mismatch")
        input_path = scratch / (stem + ".input.bin")
        output_path = scratch / (stem + ".replayed.bin")
        input_path.write_bytes(raw)
        completed = subprocess.run(
            [str(executable), str(output_model), str(output_path), str(input_path)],
            check=True, capture_output=True, text=True)
        result = json.loads(completed.stdout)
        need(result["checked"] == count and result["survivors"] == 0,
             "incomplete certificate replay")
        need(output_path.read_bytes() == raw, "replayed certificate differs")
        result.update(consistent=True, affine_space_equality=True,
                      raw_sha256=digest(raw))
        rows.append(result)
        total += count
        print(f"{stem}: {count:,} certificates checked", flush=True)
    need(total == 2121984, "incorrect complete residual universe")
    need(len(specs) == sum(r["consistent"] for r in rows),
         "extra or missing presentation family")
    result = dict(
        status="PASS", total_presentations=total,
        provenance="Rerun of incoming independent reconstruction and verifier",
        compiler_version=subprocess.run([compiler, "--version"], check=True,
                                        capture_output=True, text=True).stdout.splitlines()[0],
        sources={str(path.relative_to(HERE.parent)): digest(path.read_bytes())
                 for path in (module_path, cpp, SOURCE / "residual_frame_specs.json",
                              SOURCE / "witness_certificate_index.json")},
        families=rows)
    HERE.joinpath("incoming_replay.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(f"PASS: all {total:,} residual presentations; zero survivors.", flush=True)


if __name__ == "__main__":
    main()
