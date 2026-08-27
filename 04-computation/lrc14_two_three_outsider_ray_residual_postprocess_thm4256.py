"""Exact proof-graph postprocess for the uniform (2g,3g) ray.

Run from a checkout containing the THM-4231/4238/4242 residual script, or pass
that checkout as argv[1].  The program reconstructs the current post-THM-4254
residual, selects ray edges semantically, and freezes the exact contribution
of THM-4256 without treating already-closed ray points as new.
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


repo = Path(sys.argv[1] if len(sys.argv) == 2 else ".").resolve()
script = Path(
    "04-computation/"
    "lrc14_thm4231_4238_4242_cross_residual_postprocess_20260826.py"
)
require((repo / script).is_file(), "missing inherited residual postprocess")

prior_cwd = Path.cwd()
try:
    os.chdir(repo)
    with contextlib.redirect_stdout(io.StringIO()):
        namespace = runpy.run_path(str(script))
finally:
    os.chdir(prior_cwd)

inherited = namespace["aggregate_residual"]
require(len(inherited) == 181_126, "post-THM-4242 residual count changed")
require(edge_fnv(inherited) == 0xBDF59726990A6C92,
        "post-THM-4242 residual FNV changed")

thm4252 = {(466, 699), (616, 769), (721, 769)}
post4252 = [edge for edge in inherited if edge not in thm4252]
require(len(post4252) == 181_123, "post-THM-4252 residual count changed")

band4254 = [edge for edge in post4252 if 755 <= edge[1] <= 768]
require(len(band4254) == 59, "THM-4254 band count changed")
require(edge_fnv(band4254) == 0xB3D54B78BABBCAEC,
        "THM-4254 band FNV changed")
post4254 = [edge for edge in post4252 if edge not in set(band4254)]
require(len(post4254) == 181_064, "post-THM-4254 residual count changed")
require(edge_fnv(post4254) == 0x8F550DCC2E552962,
        "post-THM-4254 residual FNV changed")
require(edge_sha(post4254) ==
        "0167652b41139bfd00c52236338fdd50e3be604641fe03e71eb66c68ee497d35",
        "post-THM-4254 residual SHA changed")

ray = []
ray_scales = []
for left, right in post4254:
    if left % 2 != 0:
        continue
    scale = left // 2
    if scale >= 146 and right == 3 * scale:
        ray.append((left, right))
        ray_scales.append(scale)

expected_scales = [
    *range(146, 207), 208, 209, 210, 211, 212, 214, 215, 216, 217,
    220, 221, 230,
]
require(ray_scales == expected_scales, "post-THM-4254 ray universe changed")
require(len(ray) == 73, "THM-4256 new residual-edge count changed")
require(edge_fnv(ray) == 0x6BE23222A3A20764,
        "THM-4256 ray-edge FNV changed")
require(edge_sha(ray) ==
        "9d43ad3311533711e21c496141b05052c45ed68b8b03bcbce59737df3c7391ea",
        "THM-4256 ray-edge SHA changed")

ray_set = set(ray)
updated = [edge for edge in post4254 if edge not in ray_set]
require(len(updated) == 180_991, "post-THM-4256 residual count changed")
require(edge_fnv(updated) == 0x021BF0ED1581657F,
        "post-THM-4256 residual FNV changed")
require(edge_sha(updated) ==
        "9192c5d73aa5f123ddd10f0115dcaf7231fa518980610042e4cd3f8e73afd44f",
        "post-THM-4256 residual SHA changed")
maximum = max(right for _, right in updated)
top = [edge for edge in updated if edge[1] == maximum]
require(maximum == 754, "post-THM-4256 maximum endpoint changed")
require(top == [(616, 754), (698, 754), (704, 754), (721, 754)],
        "post-THM-4256 top layer changed")

print("LRC14_TWO_THREE_OUTSIDER_RAY_POSTPROCESS_THM4256")
print(f"POST_THM4254 COUNT {len(post4254)} FNV {edge_fnv(post4254):016x} "
      f"SHA256 {edge_sha(post4254)}")
print(f"THM4256_NEW_RAY_EDGES COUNT {len(ray)} FNV {edge_fnv(ray):016x} "
      f"SHA256 {edge_sha(ray)}")
print("RAY_SCALES " + ",".join(str(scale) for scale in ray_scales))
print("RAY_EDGES " + " ".join(f"{left},{right}" for left, right in ray))
print(f"UPDATED_RESIDUAL COUNT {len(updated)} FNV {edge_fnv(updated):016x} "
      f"SHA256 {edge_sha(updated)}")
print(f"TOP_LAYER ENDPOINT {maximum} COUNT {len(top)} EDGES " +
      " ".join(f"{left},{right}" for left, right in top))
print("VERDICT EXACT_PROOF_GRAPH_POSTPROCESS LRC14_OPEN")
