"""Exact proof-graph postprocess for the fixed endpoint ceiling band.

Run from a checkout containing the THM-4231/4238/4242 residual script, or pass
that checkout as argv[1].  This script recomputes the inherited residual,
applies the three THM-4252 edges, identifies the ceiling band by its semantic
predicate 755 <= second endpoint <= 768, and freezes both FNV-1a-64 and SHA-256.
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
    for a, b in edges:
        value = fnv_add(value, a)
        value = fnv_add(value, b)
    return value


def edge_sha(edges):
    encoded = b"".join(f"{a},{b}\n".encode("ascii") for a, b in edges)
    return hashlib.sha256(encoded).hexdigest()


repo = Path(sys.argv[1] if len(sys.argv) == 2 else ".").resolve()
script = Path("04-computation/lrc14_thm4231_4238_4242_cross_residual_postprocess_20260826.py")
require((repo / script).is_file(), "missing inherited THM-4231/4238/4242 postprocess")

prior_cwd = Path.cwd()
try:
    os.chdir(repo)
    with contextlib.redirect_stdout(io.StringIO()):
        inherited_namespace = runpy.run_path(str(script))
finally:
    os.chdir(prior_cwd)

inherited = inherited_namespace["aggregate_residual"]
require(len(inherited) == 181_126, "post-THM-4242 residual count changed")
require(edge_fnv(inherited) == 0xBDF59726990A6C92,
        "post-THM-4242 residual FNV changed")
require(edge_sha(inherited) ==
        "c0e2fe1c69cfe8cfe6e633a1eca0d8d37ca991ecdaa04b98d7c595a99b9be6bf",
        "post-THM-4242 residual SHA changed")

thm4252_edges = [(466, 699), (616, 769), (721, 769)]
require(all(edge in inherited for edge in thm4252_edges),
        "THM-4252 edge left inherited universe")
post4252 = [edge for edge in inherited if edge not in set(thm4252_edges)]
require(len(post4252) == 181_123, "post-THM-4252 residual count changed")
require(edge_fnv(post4252) == 0x6EC03ED4C4DC841B,
        "post-THM-4252 residual FNV changed")
require(edge_sha(post4252) ==
        "9a9b6fbe14db00e9d7f8f08ecddaa1e3d263fd063c6b3c003e18c210b3334ef8",
        "post-THM-4252 residual SHA changed")

band = [edge for edge in post4252 if 755 <= edge[1] <= 768]
expected_band = [
    (616,755),(616,756),(616,757),(616,758),(616,759),(616,760),(616,761),
    (616,762),(616,763),(616,764),(616,765),(616,766),(616,767),(616,768),
    (698,755),(698,757),
    (704,755),(704,757),(704,758),(704,759),(704,761),(704,762),(704,763),
    (704,764),(704,765),
    (721,755),(721,757),(721,758),(721,759),(721,761),(721,762),(721,763),
    (721,764),(721,765),(721,766),(721,767),(721,768),
    (726,755),(726,757),(726,758),(726,761),
    (732,755),(732,757),(732,761),(732,762),(732,763),
    (744,762),(744,763),(744,765),(744,766),(744,768),
    (750,762),(750,763),(750,765),(750,766),(750,768),
    (765,766),(765,768),(766,768),
]
expected_band.sort()
require(band == expected_band, "semantic ceiling-band universe changed")
require(len(band) == 59, "ceiling-band count changed")
require(edge_fnv(band) == 0xB3D54B78BABBCAEC, "ceiling-band FNV changed")
require(edge_sha(band) ==
        "6b54d8fa3b408325fc309bec3ed769f5e56ce370fa34fa7ad1bb6d7ed4cafc36",
        "ceiling-band SHA changed")

band_set = set(band)
updated = [edge for edge in post4252 if edge not in band_set]
require(len(updated) == 181_064, "updated residual count changed")
require(edge_fnv(updated) == 0x8F550DCC2E552962, "updated residual FNV changed")
require(edge_sha(updated) ==
        "0167652b41139bfd00c52236338fdd50e3be604641fe03e71eb66c68ee497d35",
        "updated residual SHA changed")
max_endpoint = max(b for _, b in updated)
top_layer = [edge for edge in updated if edge[1] == max_endpoint]
require(max_endpoint == 754, "updated residual maximum endpoint changed")
require(top_layer == [(616,754),(698,754),(704,754),(721,754)],
        "updated residual top layer changed")

require(len(thm4252_edges) + len(band) == 62, "cumulative removal changed")
repair_cases = len(band) * 5_852_925
body_cases = len(band) * 14_307_150
require(repair_cases == 345_322_575, "repair-case total changed")
require(body_cases == 844_121_850, "body-case total changed")

print("LRC14_FIXED_ENDPOINT_CEILING_BAND_POSTPROCESS_THM4254")
print(f"POST_THM4242 COUNT {len(inherited)} FNV {edge_fnv(inherited):016x} "
      f"SHA256 {edge_sha(inherited)}")
print(f"THM4252_REMOVED COUNT {len(thm4252_edges)} EDGES " +
      " ".join(f"{a},{b}" for a, b in thm4252_edges))
print(f"POST_THM4252 COUNT {len(post4252)} FNV {edge_fnv(post4252):016x} "
      f"SHA256 {edge_sha(post4252)}")
print(f"THM4254_BAND COUNT {len(band)} FNV {edge_fnv(band):016x} "
      f"SHA256 {edge_sha(band)} RANGE SECOND_ENDPOINT_755_768")
for endpoint in range(768, 754, -1):
    layer = [edge for edge in band if edge[1] == endpoint]
    print(f"BAND_LAYER {endpoint} COUNT {len(layer)} EDGES " +
          " ".join(f"{a},{b}" for a, b in layer))
print(f"EXACT_CASES REPAIRS {repair_cases} BODIES {body_cases}")
print(f"CUMULATIVE_REMOVED_FROM_POST_THM4242 62")
print(f"UPDATED_RESIDUAL COUNT {len(updated)} FNV {edge_fnv(updated):016x} "
      f"SHA256 {edge_sha(updated)}")
print(f"TOP_LAYER ENDPOINT {max_endpoint} COUNT {len(top_layer)} EDGES " +
      " ".join(f"{a},{b}" for a, b in top_layer))
print("CUTOFF 755 LITERAL_MINIMALITY_NOT_CLAIMED")
print("VERDICT EXACT_PROOF_GRAPH_POSTPROCESS LRC14_OPEN")
