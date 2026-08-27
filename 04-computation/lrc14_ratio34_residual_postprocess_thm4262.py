"""Exact proof-graph postprocess for the uniform (3g,4g) outsider ray.

The inherited THM-4256 postprocessor revalidates every upstream residual hash.
This script then selects the new ray contribution semantically, rather than by
copying an edge list, and freezes the resulting proof-graph residual.
"""

import contextlib
import hashlib
import io
import os
from pathlib import Path
import runpy


ROOT = Path(__file__).resolve().parents[1]
INHERITED = ROOT / (
    "04-computation/"
    "lrc14_two_three_outsider_ray_residual_postprocess_thm4256.py"
)


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


prior_cwd = Path.cwd()
try:
    os.chdir(ROOT)
    with contextlib.redirect_stdout(io.StringIO()):
        namespace = runpy.run_path(str(INHERITED))
finally:
    os.chdir(prior_cwd)

post4256 = namespace["updated"]
require(len(post4256) == 180_991, "post-THM-4256 residual count changed")
require(edge_fnv(post4256) == 0x021BF0ED1581657F,
        "post-THM-4256 residual FNV changed")
require(edge_sha(post4256) ==
        "9192c5d73aa5f123ddd10f0115dcaf7231fa518980610042e4cd3f8e73afd44f",
        "post-THM-4256 residual SHA changed")

ray = []
scales = []
for left, right in post4256:
    if left % 3 != 0:
        continue
    scale = left // 3
    if scale >= 97 and right == 4 * scale:
        ray.append((left, right))
        scales.append(scale)

expected_scales = [*range(97, 167), 172, 180]
require(scales == expected_scales, "post-THM-4256 3:4 ray universe changed")
require(len(ray) == 72, "3:4 ray contribution count changed")
require(edge_fnv(ray) == 0x512CBBA28E2235FD,
        "3:4 ray contribution FNV changed")
require(edge_sha(ray) ==
        "b1c89073d9b82351b663e97c18c807f03f3fd2d40ddcfafe038d8cad0535cb2c",
        "3:4 ray contribution SHA changed")

ray_set = set(ray)
updated = [edge for edge in post4256 if edge not in ray_set]
require(len(updated) == 180_919, "updated residual count changed")
require(edge_fnv(updated) == 0x9FAE8A515EA17DB3,
        "updated residual FNV changed")
require(edge_sha(updated) ==
        "a44da1b8dd6daa484ae2133f6fdd4c4c8d753c25292280231b4d1aefd5af9a0e",
        "updated residual SHA changed")

maximum = max(right for _, right in updated)
top = [edge for edge in updated if edge[1] == maximum]
require(maximum == 754, "updated maximum endpoint changed")
require(top == [(616, 754), (698, 754), (704, 754), (721, 754)],
        "updated top layer changed")

print("LRC14_RATIO_3_4_RAY_RESIDUAL_POSTPROCESS")
print(
    f"POST_THM4256 COUNT {len(post4256)} "
    f"FNV {edge_fnv(post4256):016x} SHA256 {edge_sha(post4256)}"
)
print(
    f"RATIO_3_4_NEW_EDGES COUNT {len(ray)} "
    f"FNV {edge_fnv(ray):016x} SHA256 {edge_sha(ray)}"
)
print("RAY_SCALES " + ",".join(map(str, scales)))
print(
    f"UPDATED_RESIDUAL COUNT {len(updated)} "
    f"FNV {edge_fnv(updated):016x} SHA256 {edge_sha(updated)}"
)
print(
    f"TOP_LAYER ENDPOINT {maximum} COUNT {len(top)} EDGES "
    + " ".join(f"{left},{right}" for left, right in top)
)
print("VERDICT EXACT_PROOF_GRAPH_POSTPROCESS LRC14_OPEN")
