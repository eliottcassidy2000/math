"""Exact proof-graph postprocess for the THM-4261 semantic endpoint band.

The inherited THM-4256 postprocessor reconstructs the current residual.  This
script selects exactly the residual pairs with second endpoint 733 through
754, freezes their ledgers, and removes them.  Pass ``--pairs`` to emit the
same semantic universe as whitespace-separated input for the two C++ audits.
"""

import contextlib
import hashlib
import io
import os
from pathlib import Path
import runpy
import sys


def require(condition, message):
    if not condition:
        raise RuntimeError(message)


def fnv_add(value, word):
    for shift in range(0, 64, 8):
        value ^= (word >> shift) & 0xFF
        value = (value * 0x100000001B3) & ((1 << 64) - 1)
    return value


def edge_fnv(edges):
    value = 0xCBF29CE484222325
    for left, right in edges:
        value = fnv_add(value, left)
        value = fnv_add(value, right)
    return value


def edge_sha(edges):
    raw = b"".join(f"{left},{right}\n".encode("ascii") for left, right in edges)
    return hashlib.sha256(raw).hexdigest()


arguments = list(sys.argv[1:])
emit_pairs = False
if "--pairs" in arguments:
    arguments.remove("--pairs")
    emit_pairs = True
require(len(arguments) <= 1, "usage: postprocess [--pairs] [REPO]")
repo = Path(arguments[0] if arguments else ".").resolve()
script = Path(
    "04-computation/"
    "lrc14_two_three_outsider_ray_residual_postprocess_thm4256.py"
)
require((repo / script).is_file(), "missing inherited THM-4256 postprocess")

prior_cwd = Path.cwd()
prior_argv = sys.argv
try:
    os.chdir(repo)
    sys.argv = [str(script), str(repo)]
    with contextlib.redirect_stdout(io.StringIO()):
        namespace = runpy.run_path(str(script))
finally:
    sys.argv = prior_argv
    os.chdir(prior_cwd)

inherited = namespace["updated"]
require(len(inherited) == 180_991, "post-THM-4256 residual count changed")
require(edge_fnv(inherited) == 0x021BF0ED1581657F,
        "post-THM-4256 residual FNV changed")
require(edge_sha(inherited) ==
        "9192c5d73aa5f123ddd10f0115dcaf7231fa518980610042e4cd3f8e73afd44f",
        "post-THM-4256 residual SHA changed")

band = [edge for edge in inherited if 733 <= edge[1] <= 754]
expected_layer_counts = {
    733: 23, 734: 23, 735: 17, 736: 19, 737: 13, 738: 18,
    739: 17, 740: 16, 741: 15, 742: 16, 743: 12, 744: 16,
    745: 8, 746: 9, 747: 10, 748: 14, 749: 12, 750: 12,
    751: 9, 752: 6, 753: 8, 754: 4,
}
require(len(band) == 297, "THM-4261 semantic band count changed")
require(edge_fnv(band) == 0xE923D1494185B820,
        "THM-4261 semantic band FNV changed")
require(edge_sha(band) ==
        "745ef7c8809335e6d6e9623314beff917edc71cfaaaa88e7210ede9dcd97d11b",
        "THM-4261 semantic band SHA changed")
require(
    {endpoint: sum(right == endpoint for _, right in band)
     for endpoint in range(733, 755)} == expected_layer_counts,
    "THM-4261 endpoint-layer profile changed",
)

band_set = set(band)
updated = [edge for edge in inherited if edge not in band_set]
require(len(updated) == 180_694, "post-THM-4261 residual count changed")
require(edge_fnv(updated) == 0x50C911CF48E3F50A,
        "post-THM-4261 residual FNV changed")
require(edge_sha(updated) ==
        "19906a8f773517f0b29767cb16b3b64502fa38c7d03dfbbbfaeff87ba71c702c",
        "post-THM-4261 residual SHA changed")

maximum = max(right for _, right in updated)
top = [edge for edge in updated if edge[1] == maximum]
expected_top = [
    (530, 732), (540, 732), (542, 732), (550, 732), (616, 732),
    (620, 732), (626, 732), (640, 732), (650, 732), (665, 732),
    (670, 732), (672, 732), (676, 732), (690, 732), (698, 732),
    (700, 732), (703, 732), (704, 732), (711, 732), (714, 732),
    (715, 732), (718, 732), (721, 732), (726, 732),
]
require(maximum == 732 and top == expected_top,
        "post-THM-4261 top layer changed")
require(edge_fnv(top) == 0xCBF337A8EC1F6F40,
        "post-THM-4261 top-layer FNV changed")
require((542, 732) in updated, "hostile boundary pair left residual")

if emit_pairs:
    for left, right in band:
        print(left, right)
else:
    print("LRC14_PREFIX_UNION_SEMANTIC_BAND_POSTPROCESS_THM4261")
    print(f"POST_THM4256 COUNT {len(inherited)} "
          f"FNV {edge_fnv(inherited):016x} SHA256 {edge_sha(inherited)}")
    print(f"THM4261_BAND COUNT {len(band)} FNV {edge_fnv(band):016x} "
          f"SHA256 {edge_sha(band)} RANGE SECOND_ENDPOINT_733_754")
    for endpoint in range(754, 732, -1):
        layer = [edge for edge in band if edge[1] == endpoint]
        print(f"BAND_LAYER {endpoint} COUNT {len(layer)} EDGES " +
              " ".join(f"{left},{right}" for left, right in layer))
    print("CERTIFICATE_CANDIDATES PREFIX_UNION 4675 "
          "PAIR_MASK_TESTS 1388475")
    print("EXACT_BODY_CASES 4249223550")
    print("CUMULATIVE_REMOVED_FROM_POST_THM4242 432")
    print(f"UPDATED_RESIDUAL COUNT {len(updated)} "
          f"FNV {edge_fnv(updated):016x} SHA256 {edge_sha(updated)}")
    print(f"TOP_LAYER ENDPOINT {maximum} COUNT {len(top)} FNV "
          f"{edge_fnv(top):016x} EDGES " +
          " ".join(f"{left},{right}" for left, right in top))
    print("BOUNDARY_HOSTILE 542,732 REMAINS")
    print("VERDICT EXACT_PROOF_GRAPH_POSTPROCESS LRC14_OPEN")
