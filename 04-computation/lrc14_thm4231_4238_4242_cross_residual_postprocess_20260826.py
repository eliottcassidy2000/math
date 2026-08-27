"""Remove THM-4242's fixed-50 r>=590 ray from the THM-4231/4238 residual.

This is an exact proof-graph postprocess.  It does not recompute Haar masses;
THM-4242 owns the universal body statement on the removed ray.
"""

import contextlib
import hashlib
import io
import runpy
from collections import defaultdict


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


BASE_POSTPROCESS = (
    "04-computation/"
    "lrc14_thm4231_4238_cross_residual_postprocess_20260826.py"
)

sink = io.StringIO()
with contextlib.redirect_stdout(sink):
    base = runpy.run_path(BASE_POSTPROCESS)

base_residual = base["aggregate_residual"]
edge_fnv = base["edge_fnv"]
require(len(base_residual) == 181_162, "base residual count changed")
require(
    edge_fnv(base_residual) == 0x7E5F6AF58A370E3A,
    "base residual fingerprint changed",
)

# THM-4242 proves every pool body safe on the entire fixed-50 ray r>=590.
# Its intersection with the current finite proof residual is exactly this
# consecutive block; the residual itself certifies that no r>625 remains.
proved_intersection = [
    edge
    for edge in base_residual
    if 50 in edge and (edge[1] if edge[0] == 50 else edge[0]) >= 590
]
expected_intersection = [(50, r) for r in range(590, 626)]
require(proved_intersection == expected_intersection, "proved ray shape changed")
require(
    edge_fnv(proved_intersection) == 0xF50A5ABB6075F4ED,
    "proved ray fingerprint changed",
)

proved_set = set(proved_intersection)
aggregate_residual = [
    edge for edge in base_residual if edge not in proved_set
]
require(len(aggregate_residual) == 181_126, "aggregate residual count changed")
require(
    edge_fnv(aggregate_residual) == 0xBDF59726990A6C92,
    "aggregate residual fingerprint changed",
)

fixed50_slice = [edge for edge in aggregate_residual if 50 in edge]
fixed50_max_other = max(
    b if a == 50 else a for a, b in fixed50_slice
)
require(len(fixed50_slice) == 556, "fixed-50 slice count changed")
require(fixed50_max_other == 589, "fixed-50 slice endpoint changed")

layers = defaultdict(list)
for edge in aggregate_residual:
    layers[edge[1]].append(edge)
max_endpoint = max(layers)
top_layer = layers[max_endpoint]
require(max_endpoint == 769, "top residual endpoint changed")
require(top_layer == [(616, 769), (721, 769)], "top residual layer changed")

encoded = b"".join(
    f"{a},{b}\n".encode("ascii") for a, b in aggregate_residual
)
edge_sha256 = hashlib.sha256(encoded).hexdigest()
require(
    edge_sha256
    == "c0e2fe1c69cfe8cfe6e633a1eca0d8d37ca991ecdaa04b98d7c595a99b9be6bf",
    "aggregate residual SHA-256 changed",
)

removed_body_cases = len(proved_intersection) * 14_307_150
require(removed_body_cases == 515_057_400, "removed body-case count changed")

print("LRC14_THM4231_THM4238_THM4242_CROSS_RESIDUAL")
print(
    f"BASE THM4231_4238 COUNT {len(base_residual)} "
    f"FNV {edge_fnv(base_residual):016x}"
)
print(
    f"THM4242_INTERSECTION COUNT {len(proved_intersection)} "
    f"FNV {edge_fnv(proved_intersection):016x} "
    "SHAPE Q50_R590_625"
)
print(
    f"REMOVED_BODY_CASES {removed_body_cases} "
    "THM4242_LITERAL_BRIDGE_ROWS 258"
)
print(
    f"AGGREGATE_RESIDUAL COUNT {len(aggregate_residual)} "
    f"FNV {edge_fnv(aggregate_residual):016x} SHA256 {edge_sha256}"
)
print(
    f"FIXED50_REMAINDER COUNT {len(fixed50_slice)} "
    f"MAX_OTHER {fixed50_max_other}"
)
print(
    f"TOP_LAYER ENDPOINT {max_endpoint} COUNT {len(top_layer)} EDGES "
    + " ".join(f"{a},{b}" for a, b in top_layer)
)
print("CUTOFF 770 UNCHANGED LITERAL_MINIMALITY_NOT_CLAIMED")
print("VERDICT EXACT_CROSS_THEOREM_PROOF_GRAPH LRC14_OPEN")
