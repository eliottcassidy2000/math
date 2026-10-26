#!/usr/bin/env python3
"""Replay this frozen packet, its historical source pins and bounded checks."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys


PACKET = Path(__file__).resolve().parent
ROOT = PACKET.parents[2]


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def digest(data):
    return hashlib.sha256(data).hexdigest()


def local(relative):
    path = PACKET / relative
    require(not path.is_symlink(), f"symlink: {relative}")
    require(path.resolve().is_relative_to(PACKET), f"outside packet: {relative}")
    return path


def git_blob(revision, path):
    result = subprocess.run(
        ["git", "show", f"{revision}:{path}"], cwd=ROOT,
        env={**os.environ, "GIT_NO_LAZY_FETCH": "1"},
        stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False,
        timeout=60,
    )
    require(result.returncode == 0,
            f"historical source unavailable {revision}:{path}: "
            + result.stderr.decode(errors="replace"))
    return result.stdout


def replay_patches(relative, live):
    record = json.loads(local(relative).read_text())
    replacements = 0
    for entry in record["files"]:
        before = git_blob(record["source_revision"], entry["path"])
        require(digest(before) == entry["before_sha256"],
                f"patch before hash: {entry['path']}")
        text = before.decode("utf-8")
        for change in entry["replacements"]:
            require(text.count(change["old"]) == 1,
                    f"patch match not unique: {entry['path']}")
            text = text.replace(change["old"], change["new"], 1)
            replacements += 1
        after = text.encode("utf-8")
        require(digest(after) == entry["after_sha256"],
                f"patch after hash: {entry['path']}")
        if live:
            require((ROOT / entry["path"]).read_bytes() == after,
                    f"live maintained file advanced: {entry['path']}")
    return {"files": len(record["files"]), "replacements": replacements}


def replay_suite(job):
    suite, optimized = job
    command = [sys.executable, "-I", "-B"]
    if optimized:
        command.append("-O")
    command.append(str(local(suite["program"])))
    result = subprocess.run(command, cwd=ROOT, stdout=subprocess.PIPE,
                            stderr=subprocess.PIPE, check=False, timeout=600)
    mode = "optimized" if optimized else "normal"
    require(result.returncode == 0,
            f"{suite['program']} {mode}: " + result.stderr.decode(errors="replace"))
    require(not result.stderr, f"unexpected stderr: {suite['program']} {mode}")
    expected = local(suite[f"{mode}_receipt"]).read_bytes()
    require(result.stdout == expected, f"receipt mismatch: {suite['program']} {mode}")
    data = json.loads(expected)
    require(data[suite["count_key"]] == suite["checks"], f"gate count: {suite['program']}")
    require(data["status"] == suite["receipt_status"], f"receipt status: {suite['program']}")
    return {"program": suite["program"], "mode": mode, "checks": suite["checks"]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--integrity-only", action="store_true")
    parser.add_argument("--live", action="store_true")
    args = parser.parse_args()
    manifest = json.loads(local("MANIFEST.json").read_text())
    actual = set()
    for path in PACKET.rglob("*"):
        require(not path.is_symlink(), f"packet symlink: {path}")
        if path.is_file():
            actual.add(path.relative_to(PACKET).as_posix())
    expected = set(manifest["files"]) | {"MANIFEST.json"}
    require(actual == expected,
            f"packet universe mismatch: extra={sorted(actual-expected)}, "
            f"missing={sorted(expected-actual)}")
    for relative, sha in manifest["files"].items():
        require(digest(local(relative).read_bytes()) == sha, f"packet hash: {relative}")
    for binding in manifest["review_bindings"]:
        require(digest(local(binding["path"]).read_bytes()) == binding["sha256"],
                f"review binding: {binding['source']} -> {binding['path']}")
    for path, sha in manifest["inherited_sources"].items():
        require(digest(git_blob(manifest["source_revision"], path)) == sha,
                f"historical source hash: {path}")
    patch_counts = {p: replay_patches(p, args.live) for p in manifest["patch_records"]}
    suites = manifest["suites"]
    for suite in suites:
        require(local(suite["normal_receipt"]).read_bytes()
                == local(suite["optimized_receipt"]).read_bytes(),
                f"normal/optimized receipts differ: {suite['program']}")
    runs = []
    if not args.integrity_only:
        jobs = [(suite, optimized) for suite in suites for optimized in (False, True)]
        with ThreadPoolExecutor(max_workers=3) as pool:
            runs = list(pool.map(replay_suite, jobs))
    print(json.dumps({"status": "PASS", "packet_files": len(expected),
                      "review_bindings": len(manifest["review_bindings"]),
                      "historical_sources": len(manifest["inherited_sources"]),
                      "live_maintained_files_checked": args.live,
                      "patches": patch_counts, "runs": runs},
                     indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
