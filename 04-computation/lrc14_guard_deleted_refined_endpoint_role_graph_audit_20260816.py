#!/usr/bin/env python3
"""Cheap independent graph-factor audit for the refined U_full response."""

import ast
from collections import Counter
import hashlib
import importlib.util
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
ROLE_NAME = "lrc_relation_role_chart_weighted_closure_probe_20260815.py"
ROLE_SHA256 = "e8eea838da1b4636c9796b71382e8a182e7ecfeb4ea17fef7eb265289889c502"
PRIMARY_OUTPUT = "lrc14_guard_deleted_refined_endpoint_role_probe_20260816.out"
PRIMARY_OUTPUT_SHA256 = "056093c45a05e50028f959a1a92ade136fa435abddea41266b30d92380e2552c"
PRIME = 572252886246508880869
VALUES = {
    "c1": 405336876493642499425,
    "c2": 503604956476841920373,
    "c3": 518539850465495448196,
    "H": 320618948602619577408,
    "q2": 15703541686881447885,
    "q3": 503604956476841920373,
    "q4": 503604956476841920373,
    "q5": 503604956476841920373,
}
EXPECTED_BRIDGE = 389266878372286537904
EXPECTED_CHART_SHA256 = "b7d8c2c9860e4f1aa542b1c85fdb7b65cf4985aba5a81a84ff3a324834d51c51"


def require(condition, payload):
    if not condition:
        raise RuntimeError(payload)


def lf_hash(path):
    return hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def load_role_module():
    path = ROOT / "04-computation" / ROLE_NAME
    actual = lf_hash(path)
    require(actual == ROLE_SHA256, (actual, ROLE_SHA256))
    spec = importlib.util.spec_from_file_location("refined_role_graph_audit", path)
    require(spec is not None and spec.loader is not None, "bad import")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def digest(value):
    data = json.dumps(value, separators=(",", ":"), sort_keys=True).encode("ascii")
    return hashlib.sha256(data).hexdigest()


def main():
    primary_output = ROOT / "05-knowledge/results" / PRIMARY_OUTPUT
    require(lf_hash(primary_output) == PRIMARY_OUTPUT_SHA256, "primary output drift")
    role = load_role_module()
    charts = role.role_charts()
    require(len(charts) == 72, len(charts))
    rows = []
    for chart in charts:
        weights = role.edge_weights(VALUES, chart)
        bridge = weights[role.EDGE_INDEX[(role.HUB, role.LEAF)]] % PRIME
        left = role.k4_tree_sum(
            VALUES, chart, role.WINGS[0] + (role.HUB,),
        ) % PRIME
        right = role.k4_tree_sum(
            VALUES, chart, role.WINGS[1] + (role.HUB,),
        ) % PRIME
        rows.append((chart, bridge, left, right, bridge * left % PRIME * right % PRIME))
    rows = tuple(rows)
    zeros = tuple(sum(row[index] == 0 for row in rows) for index in range(1, 5))
    require(zeros == (0, 0, 0, 0), zeros)
    require(set(row[1] for row in rows) == {EXPECTED_BRIDGE}, "bridge drift")
    require(digest(rows) == EXPECTED_CHART_SHA256, digest(rows))
    flat_rows = []
    flat = {label: 1 for label in VALUES}
    for chart in charts:
        weights = role.edge_weights(flat, chart)
        bridge = weights[role.EDGE_INDEX[(role.HUB, role.LEAF)]] % PRIME
        left = role.k4_tree_sum(flat, chart, role.WINGS[0] + (role.HUB,)) % PRIME
        right = role.k4_tree_sum(flat, chart, role.WINGS[1] + (role.HUB,)) % PRIME
        flat_rows.append((bridge, left, right, bridge * left % PRIME * right % PRIME))
    require(set(flat_rows) == {(0, 0, 0, 0)}, "flat hostile")
    source = Path(__file__).read_text(encoding="utf-8")
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(ast.parse(source))),
            "assert node")
    print("LRC REFINED ENDPOINT ROLE U_FULL INDEPENDENT GRAPH AUDIT")
    print("status=FINITE-EXACT graph post-audit only; THM-3479 RESERVED; LRC(14) OPEN")
    print("dependencies=%s:%s;%s:%s" % (
        ROLE_NAME, ROLE_SHA256, PRIMARY_OUTPUT, PRIMARY_OUTPUT_SHA256,
    ))
    print("universe=the five frozen U_full refined response residues; all 72 canonical role charts")
    print("implementation=separate canonical role-chart edge/K4 determinant engine")
    print("factor_zero_counts=(bridge,left_K4,right_K4,product)=%s" % (zeros,))
    print("bridge=%d multiplicity=72" % EXPECTED_BRIDGE)
    print("factor_histograms=(left=%s,right=%s)" % (
        tuple(sorted(Counter(row[2] for row in rows).items())),
        tuple(sorted(Counter(row[3] for row in rows).items())),
    ))
    print("chart_record_sha256=%s" % digest(rows))
    print("flat_response_hostile=all 72 factor quadruples are (0,0,0,0)")
    print("scope=no endpoint-bank rederivation, physical current, grouped coefficient, all-unit projector, U_clock statement, scalar exclusion, or LRC(14)")
    print("STATUS=PASS")


if __name__ == "__main__":
    main()
