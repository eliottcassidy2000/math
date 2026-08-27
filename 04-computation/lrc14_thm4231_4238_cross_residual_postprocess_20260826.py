"""Intersect THM-4238's proved small-ray region with THM-4231's residual."""

import contextlib
import hashlib
import io
import re
import runpy
from collections import defaultdict

BASE_POSTPROCESS = (
    "04-computation/"
    "lrc14_literal_boundary_45_residual_postprocess_thm4231.py"
)
LITERAL_LEDGER = (
    "05-knowledge/results/"
    "lrc14_small_label_hybrid_literal_thm4238.out"
)
P = {
    8, 10, 15, 16, 20, 30, 40, 42, 60, 63, 80, 84, 85, 88, 95,
    120, 126, 132, 143, 145, 168, 170, 176, 190, 193, 240, 252,
    264, 286, 290,
}
Q_SMALL = set(range(2, 50)) - P


sink = io.StringIO()
with contextlib.redirect_stdout(sink):
    base = runpy.run_path(BASE_POSTPROCESS)

base_residual = base["remainder"]
edge_fnv = base["edge_fnv"]
assert len(base_residual) == 181_194
assert edge_fnv(base_residual) == 0x3874FECAC4ECBD8A

row_pattern = re.compile(
    r"^ROW Q (\d+) R (\d+) .* MARGIN63 (\d+)$"
)
literal_edges = set()
with open(LITERAL_LEDGER, encoding="utf-8") as handle:
    for raw_line in handle:
        match = row_pattern.match(raw_line.rstrip("\n"))
        if match:
            q, r, margin = map(int, match.groups())
            assert margin > 0
            literal_edges.add((q, r))

expected_bridge = (
    {(6, r) for r in range(590, 614)}
    | {(25, r) for r in range(590, 598)}
)
assert literal_edges == expected_bridge

base_set = set(base_residual)
proved_intersection = {
    edge for edge in base_residual
    if edge[0] in Q_SMALL and edge[1] >= 590
}
assert proved_intersection == literal_edges
assert proved_intersection <= base_set

aggregate_residual = [
    edge for edge in base_residual if edge not in proved_intersection
]
layers = defaultdict(list)
for edge in aggregate_residual:
    layers[edge[1]].append(edge)
max_endpoint = max(layers)
top_layer = layers[max_endpoint]
encoded = b"".join(
    f"{a},{b}\n".encode("ascii") for a, b in aggregate_residual
)
edge_sha256 = hashlib.sha256(encoded).hexdigest()

assert len(aggregate_residual) == 181_162
assert edge_fnv(sorted(proved_intersection)) == 0x883C56217D9E6F35
assert edge_fnv(aggregate_residual) == 0x7E5F6AF58A370E3A
assert edge_sha256 == (
    "3a21737c3b7794f1f9faeae8c6683e16e7055877eb20fa05855ebfe8aa467c6c"
)
assert max_endpoint == 769
assert top_layer == [(616, 769), (721, 769)]

print("LRC14_THM4231_THM4238_CROSS_RESIDUAL")
print(
    f"BASE THM4231 COUNT {len(base_residual)} "
    f"FNV {edge_fnv(base_residual):016x}"
)
print(
    f"THM4238_INTERSECTION COUNT {len(proved_intersection)} "
    f"FNV {edge_fnv(sorted(proved_intersection)):016x}"
)
print("INTERSECTION_SHAPE Q6_R590_613 Q25_R590_597")
print("LITERAL_LEDGER ROWS 32 NONPOSITIVE 0")
print(
    f"AGGREGATE_RESIDUAL COUNT {len(aggregate_residual)} "
    f"FNV {edge_fnv(aggregate_residual):016x} SHA256 {edge_sha256}"
)
print(
    f"TOP_LAYER ENDPOINT {max_endpoint} COUNT {len(top_layer)} EDGES "
    + " ".join(f"{a},{b}" for a, b in top_layer)
)
print("CUTOFF 770 UNCHANGED LITERAL_MINIMALITY_NOT_CLAIMED")
print("VERDICT EXACT_CROSS_THEOREM_PROOF_GRAPH LRC14_OPEN")
