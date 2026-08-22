#!/usr/bin/env python3
"""Inventory branch-only research strata against a chosen current main ref.

This is a provenance scanner, not a mathematical verifier.  It reports two
explicit universes:

* nonsymbolic local-head and remote-tracking refs (the branch universe);
* every Git ref, including a private stash (the all-ref comparison).

The motif census is deliberately filename-only.  It is a routing aid for the
manual truth audit recorded in the companion reflection; a keyword hit is
never evidence that a historical claim is correct or novel.
"""

from __future__ import annotations

import argparse
import collections
import pathlib
import re
import subprocess
import sys
from dataclasses import dataclass


ROOT = pathlib.Path(__file__).resolve().parents[1]

RESEARCH_PREFIXES = (
    "00-navigation/",
    "01-canon/",
    "04-computation/",
    "05-knowledge/",
    "06-writeups/",
    "07-reflections/",
    "10-conjectures/",
)

MOTIFS = {
    "continuation_transfer": re.compile(
        r"continuation|transfer[-_ ]spectrum|temperature[-_ ]ladder|automaton|recurrence",
        re.I,
    ),
    "polarized_correlation": re.compile(
        r"autocorrel|correlation|gram|phase[-_ ]transport|homometr", re.I
    ),
    "event_flow_first_passage": re.compile(
        r"nowhere[-_ ]zero[-_ ]flow|crossing[-_ ]flow|first[-_ ]passage|friendliness|survival",
        re.I,
    ),
    "valuation_unit_log": re.compile(
        r"momentum[-_ ]twistor|discrete[-_ ]log|divisor[-_ ]lattice|valuation|unit[-_ ]log",
        re.I,
    ),
    "interaction_cross_effect": re.compile(
        r"inclusion[-_ ]exclusion|mobius|m.bius|interaction|multitoggle|plaquette",
        re.I,
    ),
    "dual_certificate": re.compile(
        r"delsarte|krawtchouk|farkas|dual[-_ ]certificate|linear[-_ ]program",
        re.I,
    ),
}

POSITIVE_CONTROLS = (
    "05-knowledge/hypotheses/HYP-2355-autocorrelation-operator-unification.md",
    "05-knowledge/hypotheses/HYP-2381-transfer-spectrum-ledger-temperature-frozen-parity.md",
    "01-canon/theorems/THM-445-lrc-momentum-twistor-discrete-log.md",
    "07-reflections/lrc-nowhere-zero-flow-attack-s521.md",
    "05-knowledge/hypotheses/HYP-2325-friendliness-never-lonely-yet-first-passage.md",
)

HOSTILE_CURRENT_CONTROL = (
    "01-canon/theorems/THM-2355-component-deletion-gram-and-twist-energy-phase-transport.md"
)


@dataclass(frozen=True)
class RefRow:
    scope: str
    ref: str
    object_name: str
    symref: str


def git(*args: str) -> str:
    return subprocess.check_output(
        ["git", *args], cwd=ROOT, text=True, errors="surrogateescape"
    )


def live_ref_rows() -> tuple[list[RefRow], list[str]]:
    rows: list[RefRow] = []
    for line in git(
        "for-each-ref",
        "--format=%(refname)%09%(objectname)%09%(symref)",
    ).splitlines():
        ref, object_name, symref = line.split("\t")
        scope = "branch" if ref.startswith(("refs/heads/", "refs/remotes/")) else "other"
        rows.append(RefRow(scope, ref, object_name, symref))
    return rows, []


def manifest_ref_rows(path: pathlib.Path) -> tuple[list[RefRow], list[str]]:
    rows: list[RefRow] = []
    headers: list[str] = []
    for line_number, raw in enumerate(path.read_text().splitlines(), 1):
        if not raw:
            continue
        if raw.startswith("#"):
            headers.append(raw[1:].strip())
            continue
        fields = raw.split("\t")
        if len(fields) != 4:
            raise ValueError(f"{path}:{line_number}: expected four tab-separated fields")
        scope, ref, object_name, symref = fields
        if scope not in {"branch", "other"}:
            raise ValueError(f"{path}:{line_number}: invalid scope {scope!r}")
        rows.append(RefRow(scope, ref, object_name, "" if symref == "-" else symref))
    return rows, headers


def header_value(headers: list[str], key: str) -> str | None:
    prefix = f"{key}="
    return next((line[len(prefix) :] for line in headers if line.startswith(prefix)), None)


def revs_outside(main_object: str, objects: list[str]) -> set[str]:
    if not objects:
        return set()
    return set(git("rev-list", *objects, "--not", main_object).splitlines())


def touched_paths(main_object: str, objects: list[str]) -> set[str]:
    if not objects:
        return set()
    args = ["log", *objects, "--not", main_object, "--format=", "--name-only"]
    return {line for line in git(*args).splitlines() if line}


def print_examples(paths: list[str], limit: int = 8) -> None:
    for path in paths[:limit]:
        print(f"    {path}")
    if len(paths) > limit:
        print(f"    ... {len(paths) - limit} more")


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--main-ref", default="origin/main")
    parser.add_argument(
        "--ref-manifest",
        type=pathlib.Path,
        help="freeze ref names and object IDs; main_sha in the manifest takes precedence",
    )
    args = parser.parse_args()

    if args.ref_manifest:
        try:
            rows, headers = manifest_ref_rows(args.ref_manifest)
        except (OSError, ValueError) as exc:
            print(f"error: {exc}", file=sys.stderr)
            return 2
        main_ref = header_value(headers, "main_ref") or args.main_ref
        main_object = header_value(headers, "main_sha")
        if not main_object:
            print("error: frozen manifest has no main_sha header", file=sys.stderr)
            return 2
        mode = "frozen_manifest"
        head_at_capture = header_value(headers, "head_sha") or "UNKNOWN"
    else:
        rows, headers = live_ref_rows()
        main_ref = args.main_ref
        main_object = ""
        mode = "live_refs"
        head_at_capture = git("rev-parse", "HEAD").strip()

    try:
        main_sha = git("rev-parse", main_object or main_ref).strip()
    except subprocess.CalledProcessError:
        print(f"error: cannot resolve main object {main_object or main_ref!r}", file=sys.stderr)
        return 2

    branch_rows = [row for row in rows if row.scope == "branch"]
    refs = [(row.ref, row.object_name) for row in branch_rows if not row.symref]
    branch_objects = [object_name for _, object_name in refs]
    all_objects = [row.object_name for row in rows if not row.symref]
    outside_counts = {
        ref: int(git("rev-list", "--count", object_name, "--not", main_sha).strip())
        for ref, object_name in refs
    }
    live_branch_refs = [ref for ref, count in outside_counts.items() if count]
    branch_commits = revs_outside(main_sha, branch_objects)
    all_commits = revs_outside(main_sha, all_objects)

    main_paths = set(git("ls-tree", "-r", "--name-only", main_sha).splitlines())
    branch_touched = touched_paths(main_sha, branch_objects)
    all_touched = touched_paths(main_sha, all_objects)
    branch_absent = branch_touched - main_paths
    all_absent = all_touched - main_paths
    research_absent = sorted(
        path for path in branch_absent if path.startswith(RESEARCH_PREFIXES)
    )

    print("ALL-REF CONCEPT ARCHAEOLOGY INVENTORY")
    print("=" * 76)
    print(f"inventory_mode={mode}")
    if args.ref_manifest:
        print(f"ref_manifest={args.ref_manifest.name}")
    print(f"main_ref={main_ref}")
    print(f"main_sha={main_sha}")
    print(f"head_sha_at_capture={head_at_capture}")
    print()
    print("BRANCH / REMOTE-TRACKING UNIVERSE")
    print(f"raw_ref_entries={len(branch_rows)}")
    print(f"symbolic_ref_entries={sum(bool(row.symref) for row in branch_rows)}")
    print(f"nonsymbolic_refs={len(refs)}")
    print(f"refs_with_commits_outside_main={len(live_branch_refs)}")
    print(f"unique_commits_outside_main={len(branch_commits)}")
    print(f"unique_paths_touched_outside_main={len(branch_touched)}")
    print(f"touched_paths_absent_from_main={len(branch_absent)}")
    print(f"absent_research_paths={len(research_absent)}")
    print()
    print("ALL GIT REFS (INCLUDES PRIVATE STASH)")
    print(f"unique_commits_outside_main={len(all_commits)}")
    print(f"unique_paths_touched_outside_main={len(all_touched)}")
    print(f"touched_paths_absent_from_main={len(all_absent)}")
    print(f"extra_commits_beyond_branch_universe={len(all_commits - branch_commits)}")
    print(f"extra_absent_paths_beyond_branch_universe={len(all_absent - branch_absent)}")

    print()
    print("ABSENT RESEARCH-PATH MOTIF CENSUS (FILENAME ROUTING ONLY)")
    claimed: collections.Counter[str] = collections.Counter()
    for name, regex in MOTIFS.items():
        hits = [path for path in research_absent if regex.search(path)]
        claimed.update(hits)
        print(f"{name}: {len(hits)}")
        print_examples(hits)
    overlap_count = sum(count > 1 for count in claimed.values())
    print(f"paths_matching_multiple_motifs={overlap_count}")

    print()
    print("CONTROLS")
    for path in POSITIVE_CONTROLS:
        print(f"positive_absent={path in branch_absent}  {path}")
    print(
        "hostile_current_present="
        f"{HOSTILE_CURRENT_CONTROL in main_paths}  {HOSTILE_CURRENT_CONTROL}"
    )
    print(
        "hostile_current_absent="
        f"{HOSTILE_CURRENT_CONTROL in branch_absent}  {HOSTILE_CURRENT_CONTROL}"
    )

    print()
    print("INTERPRETATION GUARDRAIL")
    print("A path or motif hit proves only historical reachability and current absence.")
    print("Truth/status/novelty require a manual audit against current canon and MISTAKES.md.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
