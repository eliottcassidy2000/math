#!/usr/bin/env python3
"""Verify frozen input bytes, independently replay normal/-O, and compare receipts.

Run from any working directory. --record writes only this audit's receipts;
without it the program is read-only. No network, Git, or producer imports.
"""
from pathlib import Path
import argparse
import hashlib
import json
import subprocess
import sys

HERE = Path(__file__).resolve().parent
PINS = {
    "session19-ellipsoid-optimum": ("050fc733781c929b8d2eb36a128cec3bc3d0f7b85c95f491941dcb3dbf849dad", 71504),
    "session19-triangle-ellipsoid-repair": ("43ec88a26d022a7e10f9968a1a4e84f2634325c42550ab7defe34f9151ab17fe", 31344),
}


def require(ok, why):
    if not ok:
        raise RuntimeError(why)


def sha(data):
    return hashlib.sha256(data).hexdigest()


def verify_sources(root):
    result = []
    for name, (pin, gates) in PINS.items():
        folder = root/name
        manifest_bytes = (folder/"MANIFEST.json").read_bytes()
        require(sha(manifest_bytes) == pin, name+": manifest hash mismatch")
        manifest = json.loads(manifest_bytes)
        require(manifest["total_exact_gates"] == gates, name+": declared gate mismatch")
        require((folder/"MANIFEST.sha256").read_text().split()[0] == pin,
                name+": manifest sidecar mismatch")
        expected = {"MANIFEST.json", "MANIFEST.sha256"}
        for entry in manifest["files"]:
            relative = Path(entry["path"])
            require(not relative.is_absolute() and ".." not in relative.parts, "unsafe source path")
            data = (folder/relative).read_bytes()
            require(len(data) == entry["bytes"] and sha(data) == entry["sha256"],
                    name+": input mismatch: "+str(relative))
            expected.add(str(relative))
        actual = {str(p.relative_to(folder)) for p in folder.rglob("*") if p.is_file()}
        require(actual == expected, name+": unexpected input inventory")
        result.append({"packet": name, "manifest_sha256": pin, "pinned_files": len(expected),
                       "producer_exact_gates": gates})
    return result


def verify_audit_manifest():
    path = HERE/"MANIFEST.json"
    if not path.exists():
        return
    data = path.read_bytes()
    require((HERE/"MANIFEST.sha256").read_text().split()[0] == sha(data),
            "Audit manifest sidecar mismatch")
    expected = {"MANIFEST.json", "MANIFEST.sha256"}
    for entry in json.loads(data)["files"]:
        relative = Path(entry["path"])
        require(not relative.is_absolute() and ".." not in relative.parts, "unsafe audit path")
        content = (HERE/relative).read_bytes()
        require(len(content) == entry["bytes"] and sha(content) == entry["sha256"],
                "Audit artifact mismatch: "+str(relative))
        expected.add(str(relative))
    actual = {str(p.relative_to(HERE)) for p in HERE.rglob("*") if p.is_file()}
    require(actual == expected, "Audit inventory mismatch")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--record", action="store_true")
    parser.add_argument("--originals", type=Path, help="Also verify unchanged original packets under this parent")
    args = parser.parse_args()
    if not args.record:
        verify_audit_manifest()
    snapshots = HERE/"frozen_sources"
    inputs = verify_sources(snapshots)
    if args.originals:
        require(verify_sources(args.originals) == inputs, "Original/snapshot mismatch")
    jobs = [(name, snapshots/name/"check_exact.py", snapshots/name/"check_exact.out", gates)
            for name, (_, gates) in PINS.items()]
    jobs.append(("independent", HERE/"check_independent.py", None, None))
    outputs, details = {}, []
    for name, script, expected_path, expected_gates in jobs:
        runs = []
        for mode, flags in (("normal", []), ("optimized", ["-O"])):
            p = subprocess.run([sys.executable, "-B"]+flags+[str(script)], cwd=HERE,
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
            require(p.returncode == 0 and not p.stderr, name+" "+mode+": unsuccessful run: "+p.stderr.decode())
            report = json.loads(p.stdout)
            require(report["status"] == "PASS", name+": result not PASS")
            require(sum(report["checks_by_label"].values()) == report["total_exact_gates"],
                    name+": inconsistent exact-gate sum")
            if expected_gates is not None:
                require(report["total_exact_gates"] == expected_gates, name+": gate mismatch")
            if expected_path:
                require(p.stdout == expected_path.read_bytes(), name+": frozen output mismatch")
            key = name+"."+mode+".out"
            outputs[key] = p.stdout
            runs.append({"mode": mode, "returncode": p.returncode, "stderr_bytes": len(p.stderr),
                         "stdout_bytes": len(p.stdout), "stdout_sha256": sha(p.stdout),
                         "total_exact_gates": report["total_exact_gates"], "receipt": "receipts/"+key})
        require(outputs[name+".normal.out"] == outputs[name+".optimized.out"], name+": normal/-O differ")
        details.append({"checker": name, "script_sha256": sha(script.read_bytes()),
                        "normal_optimized_byte_identical": True, "runs": runs})
    require(verify_sources(snapshots) == inputs, "Snapshot changed during replay")
    if args.originals:
        require(verify_sources(args.originals) == inputs, "Original changed during replay")
    receipt = {"status": "PASS", "python": sys.version.split()[0], "input_packets": inputs,
               "replays": details, "source_hashes_rechecked_after_execution": True,
               "scope": "Frozen-artifact validation and finite exact replay, not a formal proof assistant"}
    outputs["replay.json"] = (json.dumps(receipt, indent=2, sort_keys=True)+"\n").encode()
    for name, content in outputs.items():
        path = HERE/"receipts"/name
        if args.record:
            path.parent.mkdir(exist_ok=True)
            path.write_bytes(content)
        else:
            if name == "replay.json":
                # Interpreter patch versions may differ on a portable replay.
                previous, current = json.loads(path.read_bytes()), json.loads(content)
                previous.pop("python", None)
                current.pop("python", None)
                require(previous == current, "Retained replay summary mismatch")
            else:
                require(path.read_bytes() == content, "Retained output mismatch: "+name)
    print(json.dumps(receipt, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
