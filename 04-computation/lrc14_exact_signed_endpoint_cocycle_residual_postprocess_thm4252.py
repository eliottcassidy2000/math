"""Exact proof-graph update for THM-4252.

The three edge theorems are consumed only after their full 30-pool/9-body
certificates have been verified by the C++ companions.  This script performs
no Haar integration; it removes exactly those named proved edges from the
current THM-4231/4238/4242 residual and freezes the resulting graph.

Every gate is an explicit require(), so ``python -O`` cannot disable it.
"""

import contextlib
import hashlib
import io
import runpy


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


BASE_POSTPROCESS = (
    "04-computation/"
    "lrc14_thm4231_4238_4242_cross_residual_postprocess_20260826.py"
)

sink = io.StringIO()
with contextlib.redirect_stdout(sink):
    base = runpy.run_path(BASE_POSTPROCESS)

base_residual = base["aggregate_residual"]
edge_fnv = base["edge_fnv"]

require(len(base_residual) == 181_126, "base residual count changed")
require(
    edge_fnv(base_residual) == 0xBDF59726990A6C92,
    "base residual FNV changed",
)

proved_edges = [(466, 699), (616, 769), (721, 769)]
require(
    all(edge in base_residual for edge in proved_edges),
    "a THM-4252 edge is absent from the inherited residual",
)
require(len(set(proved_edges)) == 3, "proved edge list contains a duplicate")

proved_set = set(proved_edges)
residual = [edge for edge in base_residual if edge not in proved_set]
require(len(residual) == 181_123, "updated residual count changed")
require(edge_fnv(residual) == 0x6EC03ED4C4DC841B, "updated residual FNV changed")

encoded = b"".join(f"{a},{b}\n".encode("ascii") for a, b in residual)
residual_sha256 = hashlib.sha256(encoded).hexdigest()
require(
    residual_sha256
    == "9a9b6fbe14db00e9d7f8f08ecddaa1e3d263fd063c6b3c003e18c210b3334ef8",
    "updated residual SHA-256 changed",
)

maximum_endpoint = max(b for _, b in residual)
top_layer = [edge for edge in residual if edge[1] == maximum_endpoint]
expected_top_layer = [
    (616, 768),
    (721, 768),
    (744, 768),
    (750, 768),
    (765, 768),
    (766, 768),
]
require(maximum_endpoint == 768, "updated maximum endpoint changed")
require(top_layer == expected_top_layer, "updated top layer changed")

removed_body_cases = len(proved_edges) * 14_307_150
require(removed_body_cases == 42_921_450, "removed body-case count changed")

print("THM4252_EXACT_RESIDUAL_POSTPROCESS_V1")
print(
    f"BASE COUNT {len(base_residual)} FNV {edge_fnv(base_residual):016x} "
    "SHA256 c0e2fe1c69cfe8cfe6e633a1eca0d8d37ca991ecdaa04b98d7c595a99b9be6bf"
)
print(
    "PROVED_EDGES COUNT 3 EDGES "
    + " ".join(f"{a},{b}" for a, b in proved_edges)
)
print(f"REMOVED_BODY_CASES {removed_body_cases}")
print(
    f"UPDATED COUNT {len(residual)} FNV {edge_fnv(residual):016x} "
    f"SHA256 {residual_sha256}"
)
print(
    f"TOP_LAYER ENDPOINT {maximum_endpoint} COUNT {len(top_layer)} EDGES "
    + " ".join(f"{a},{b}" for a, b in top_layer)
)
print("CUTOFF 769 CERTIFICATE_EXACT LITERAL_MINIMALITY_NOT_CLAIMED")
print("VERDICT PASS THM4252_THREE_EDGE_PROOF_GRAPH_UPDATE")
