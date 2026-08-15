#!/usr/bin/env python3
"""Exact post-THM3361 projected-k3 z216 wall closure audit.

Theorem-grade companion.  It reconstructs the 110 historical
post-THM3361 rows and all 12 complete families, then removes THM3378's row 94
to recover the 109-row live wall.  It replays the inherited exact
ray/common-status screen,
and closes every residual packet by a lawful common-cell terminal.  Positive
duplicate-two-high gap rows reduce to one high and use located torsion or
translation-uniform residue-support cardinality.  The unique nonpositive-gap
row is enumerated exactly and uses common-modulus support, a
multiplicity-preserving top-h fibre majorant, or the denominator-two measure
terminal.  No raw drift ratio is used; every high is quantified over every
unit, height, and local-coordinate translate on its denominator ray.
"""

from __future__ import annotations

import argparse
import ast
import hashlib
import importlib.util
import multiprocessing as mp
from array import array
from collections import Counter
from fractions import Fraction
from itertools import combinations_with_replacement
from math import gcd, lcm
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PREFIX = ROOT / "04-computation/lrc14_j7_k3_z216_sixteen_family_prefix_translated_two_high_closure_scout_20260812.py"
PREFIX_OUT = ROOT / "05-knowledge/results/lrc14_j7_k3_z216_sixteen_family_prefix_translated_two_high_closure_scout_20260812.out"
THM3361 = ROOT / "04-computation/lrc14_j7_k3_z216_L720720_one_high_translated_residue_thm3361.py"
THM3361_OUT = ROOT / "05-knowledge/results/lrc14_j7_k3_z216_L720720_one_high_translated_residue_thm3361.out"
THM3378 = ROOT / "04-computation/lrc14_j7_k3_z216_gcd24_L129360_row94_one_high_torsion_thm3378.py"
THM3378_OUT = ROOT / "05-knowledge/results/lrc14_j7_k3_z216_gcd24_L129360_row94_one_high_torsion_thm3378.out"
AUDIT = ROOT / "04-computation/lrc14_j7_k3_z216_cost_prefix_independent_audit_20260803.py"
NATURAL = ROOT / "04-computation/lrc14_j7_k3_z216_natural_wall_family_screen_descent_thm3281.py"
SOURCE3139 = ROOT / "04-computation/lrc14_j7_k3_z225_terminal_z224_screen_double_layer_descent_thm3139.py"
SOURCE3114 = ROOT / "04-computation/lrc14_j7_k3_z227_screen_z226_terminal_double_layer_descent_thm3114.py"
SOURCE3113 = ROOT / "04-computation/lrc14_j7_k3_z229_terminal_z228_screen_double_layer_descent_thm3113.py"
SOURCE3111 = ROOT / "04-computation/lrc14_j7_k3_z230_exact_screen_compressed_complete_cell_descent_thm3111.py"
SOURCE3109 = ROOT / "04-computation/lrc14_j7_k3_z231_exact_screen_complete_cell_cardinality_descent_thm3109.py"
SOURCE3106 = ROOT / "04-computation/lrc14_j7_k3_z232_exact_screen_complete_cell_cardinality_descent_thm3106.py"
TORSION = ROOT / "04-computation/lrc14_j7_k3_z297_ray_status_torsion_closure_thm2941.py"
FOURTH = ROOT / "04-computation/lrc14_j7_k3_z216_fourth_intrinsic_ruler_cost_prefix_affine_multicover_closure_scout_20260803.py"
FOURTH_OUT = ROOT / "05-knowledge/results/lrc14_j7_k3_z216_fourth_intrinsic_ruler_cost_prefix_affine_multicover_closure_scout_20260803.out"
Z378 = ROOT / "04-computation/lrc14_j7_k3_z378_ray_status_closure_thm2941.py"
ATLAS = ROOT / "04-computation/lrc14_j7_k3_projected_scalar_body_atlas_thm2941.py"
ATLAS_OUT = ROOT / "05-knowledge/results/lrc14_j7_k3_projected_scalar_body_atlas_thm2941.out"
THM3391_SOURCE = ROOT / "04-computation/lrc14_weighted_common_source_cyclic_support_capacity_thm3391.py"
THM3391_OUT = ROOT / "05-knowledge/results/lrc14_weighted_common_source_cyclic_support_capacity_thm3391.out"
ROW195_AUDIT = ROOT / "04-computation/lrc14_j7_k3_z216_row195_two_high_top_weight_majorant_independent_audit_20260814.py"
ROW195_AUDIT_OUT = ROOT / "05-knowledge/results/lrc14_j7_k3_z216_row195_two_high_top_weight_majorant_independent_audit_20260814.out"
ROW195_ACCELERATOR = ROOT / "04-computation/lrc14_j7_k3_z216_row195_three_block_top_h_torsion_accelerator_audit_20260815.py"
ROW195_ACCELERATOR_OUT = ROOT / "05-knowledge/results/lrc14_j7_k3_z216_row195_three_block_top_h_torsion_accelerator_audit_20260815.out"

EXPECTED_HASHES = {
    PREFIX: "cfb020bfc6636a52f1eaf55f82a925e70c11c90da7f87f36b0bd77ece1ec6a62",
    PREFIX_OUT: "a88646fbd28d807a0cc9671c509c4424056a539b49d04a2076ba17de57ef5ee4",
    THM3361: "52d70825de6199182c6a8dc3c2674cb0f6f92cfdec8a2b49e22c478307fc3285",
    THM3361_OUT: "22732fcf2ed2d094834dbbafa9893129856de52922892881f485e2e2b34cdede",
    THM3378: "62769c8024ccd3f85a71b858e65f909a216e4ea85a7bca2bfb5fd0aa61a43a73",
    THM3378_OUT: "26db8dc95acac0afe179aedce04543150755860215cd506b55c09cf0c872c37a",
    AUDIT: "d3581bc937597270558a77c0a02398a4cfd368306f5a302be01240cdb4cef7c3",
    NATURAL: "430dee7ba03e0d5c9ae0df72ac512500de4f7056cb4663d1c8468bfb93a49bfe",
    SOURCE3139: "92e1c22088998d37db89020ac1ebbb7d1ee17e3539152ab205f8b0cd92532e36",
    SOURCE3114: "8de1e3d03b5070a84b040ac13a173a598646107f85e7ba0defc2ca070808f162",
    SOURCE3113: "1e23ec19fa147c55fb6d38a965eedae0132f5e069b9f820bfd5c300dce4d8f89",
    SOURCE3111: "42323171481deba2371eed9947b2079976cb367dac340cf58b8f1f0c0afb5082",
    SOURCE3109: "1f74f2b2368c04f514f2c388b54c70a9ee66c9387fbc437093884b807b3eb23c",
    SOURCE3106: "f6f64ab8d8ea9b04a1a03e26fc6026efc864e44518e9cb40df4fe8471a4a7991",
    TORSION: "d062c7ac8ebf6a433c8fb1543293e941c85625e2eb40b82fcf05fc2404539b0a",
    FOURTH: "b515c70174d58ad859a08c29949cdee36a4e04122451a66237732570ab5ee213",
    FOURTH_OUT: "d1845611f27d427a1d38afe349ed07bc964590fd39a3e88cd45e6ea34a86bc38",
    Z378: "2ef5e0639354c38b13e17e41f91acb4143c7f60973295b0e2dd0f57eb8f38db2",
    ATLAS: "2af6d96882f336a409a8657070ed76a75c09a53b3789101b83103b051e864ded",
    ATLAS_OUT: "cee82237ce1f51729813b9c916edd3353204c18172abe1d71278dee2c5562eda",
    THM3391_SOURCE: "22c2ea187e3d43ca55dd61611a0f6d8a70cf7b1111b1f01cb7338bc1aef7e195",
    THM3391_OUT: "9cc8b652ae3552f970fae1b8f46f3b6c1d4316a5170d2f9a37eb5e59495e3062",
    ROW195_AUDIT: "fccc10392624bbdfc2083993ad51a423e8974c135b9bc635351304a71cb0de74",
    ROW195_AUDIT_OUT: "55ad1da385d35f5b38fdc7de2d9916f54ee874b806496d4a7b1270cb526ad30c",
    ROW195_ACCELERATOR: "00e1855dd428466c113d3d559e629baefcfdd9b0ea8ac732069fee3bb4ac74fe",
    ROW195_ACCELERATOR_OUT: "b188805a09b484466f4c36395b2b426beba8ad470704d9c28414f92ec7a43f13",
}

LEVEL = 216
THM3361_CLOSED = (191, 228, 332)
EXPECTED_ATLAS = (480, 447, 33)
EXPECTED_QUEUE = (110, 12)
EXPECTED_INDICES_SHA256 = "75f7e326e7694af41d6eedc3c520d4f7dd3c38e7fa83b34d7eb5bccc5cda460d"
EXPECTED_LIVE_QUEUE = (109, 12)
EXPECTED_LIVE_INDICES_SHA256 = "e0d24b7f8a46cfdaa4b3c2fbe3aca070b4004b7b6a81f1f1083799b8e9138d80"

# Filled after the first complete scratch census, then frozen for normal/-O.
EXPECTED_FAMILIES_SHA256 = "3f9bd7b1ff330b6a804fae4a280265b2b5791d7bd33ddeec00f03f1b377d3921"
EXPECTED_HEAVY_INDICES = (459, 463, 477)
EXPECTED_HEAVY_SCREEN_SHA256 = "3cb8e34caa7b1d9f3d3e058bd9bc5ea58fa42cf841124e03877f5b74ba2b6d58"
EXPECTED_HEAVY_TERMINAL_SHA256 = "9bbac6d9e329647d830097de5a02be273c85b785039c77bd96ee7c4a54a45ed9"
EXPECTED_NONPOSITIVE_ROWS = (195,)
EXPECTED_ROW195_TWO_CASES = (
    137,
    "5811cb585e35831a9c58b2d0af35ab6db89e4503229e6420259d6ce6d4586c3c",
)
EXPECTED_ROW195_TWO_MECHANISMS = (
    ("common-modulus-support", 76),
    ("denominator-two-measure", 1),
    ("top-h-weighted-common-source", 60),
)
EXPECTED_ROW195_WEAKEST_TOP_H = (10456, 130, (32760, 458640), 24800, 65520)
EXPECTED_ROW195_CLOSURES_SHA256 = "4cdb5acfda8b7db6e061b54ea2d71102c29a8694a7fce85f6a2c25e317a93ef9"
EXPECTED_ROW195_CELL_PACKET = (
    (351, 100776, 37987, 420652, "89236b7d22a6be06afecaaaba8f2e0be1346b1cab3c0b720314baa4014d87c36"),
)
EXPECTED_ROW195_PERIOD_OBSTRUCTION = (
    (
        50960,
        ((1, 760), (2, 7944), (3, 9848), (4, 11626), (5, 1616)),
        5,
        37987,
        88947,
        (True, False, False, True, True, False, False, False, False),
    ),
)
EXPECTED_ROW195_CONTROL_HASHES = (
    ("selfcontained_independent_source", "fccc10392624bbdfc2083993ad51a423e8974c135b9bc635351304a71cb0de74"),
    ("selfcontained_independent_output", "55ad1da385d35f5b38fdc7de2d9916f54ee874b806496d4a7b1270cb526ad30c"),
    ("selfcontained_independent_semantic", "7a70daaa1f25fc6fb9bfd7469ce3e5a18ba650296119789c81bc70fd32016b9c"),
    ("selfcontained_accelerator_source", "00e1855dd428466c113d3d559e629baefcfdd9b0ea8ac732069fee3bb4ac74fe"),
    ("selfcontained_accelerator_output", "b188805a09b484466f4c36395b2b426beba8ad470704d9c28414f92ec7a43f13"),
    ("selfcontained_accelerator_semantic", "5ad577c0a388c0021cda9506a8147a25732f36f801e7bcf746b856a5d9ff1070"),
)
EXPECTED_ROW195_ACCELERATOR_SEMANTIC = "5ad577c0a388c0021cda9506a8147a25732f36f801e7bcf746b856a5d9ff1070"
EXPECTED_ROW195_ACCELERATOR_STATES_SHA256 = "2a8584886131c177149055d208f5ccd3b356cc42dcdc84b6533773f6b800b935"
EXPECTED_ROW195_ACCELERATOR_SUFFIX_SHA256 = "7d8644411a2d78ed40a6c593a7f977979803aac2e0f2427f788ee58da7c9489c"
EXPECTED_ROW195_ACCELERATOR_FIRST_MULTIPLICITIES = ((1, 65014), (2, 711), (3, 3))
EXPECTED_ROW195_ACCELERATOR_NEEDED = (
    78,
    "6b12fcbf3b4d7415cc4006d72ed57099bd635c09a084dda231e9643dad2c7803",
)
EXPECTED_ROW195_ACCELERATOR_FIRST_SOURCE = (
    119368,
    37800,
    420839,
    "745925d13ce82fd6f85d5184de86a12a3cb2c370890d6f60ec3d95d354d13cae",
)
EXPECTED_ROW195_ACCELERATOR_TOP_H_SHA256 = "ce68e2ad70faf1bd0d5581ac490d84e345eea97479a42969927431c15bcc5985"
EXPECTED_ROW195_ACCELERATOR_TRIPLES_SHA256 = "9ab6d15791f021aaef022dfb27e0e231cc1826b2bb9ec1ac65747ea6fb6eab82"
EXPECTED_ROW195_ACCELERATOR_FAILED_SHA256 = "716c0d4435ad8e7ff3b7912acf748f17591eb2c109bdea3343870aea8cb5c320"
EXPECTED_ROW195_ACCELERATOR_ATTAINED_SHA256 = "82edd53458fef7fafbec9e10683fa8d6799827cf6e4b49138a3b142609fae6a5"
EXPECTED_ROW195_ACCELERATOR_CRUDE_SHA256 = "66bd5ec297b30d419efdb585e0df215d0de00334db89b5e539af54b7d1d89b83"
EXPECTED_ROW195_ACCELERATOR_STATUS_SHA256 = "3b96905d9f42302b530dc8978de1390697de086f709cf10086c613b345c5183f"
EXPECTED_ROW195_ACCELERATOR_TRIPLE_PARTITION = (82160, 80687, 1473)
EXPECTED_ROW195_ACCELERATOR_FAILED_PARTITION = (1473, 695, 295, 243, 157)
EXPECTED_ROW195_ACCELERATOR_RESIDUAL_PARTITION = (1565, 1408, 157)
EXPECTED_ROW195_ACCELERATOR_STRICT_RESIDUAL_SHA256 = "3e3fe58c04ccbbbc1751d48bd199b76782c37589c58625f23767f5f870b8818d"
EXPECTED_ROW195_ACCELERATOR_STRICT_CERTIFICATE_SHA256 = "3ddbcf69c2c1521e1440973cdc46be0e85ca769a6ab68a247c8d00115340e372"
EXPECTED_ROW195_ACCELERATOR_HOSTILE_SHA256 = "e997f387c89563b9d87fca78e049f4484f0712a21879a7f4b6088461f704e947"
EXPECTED_ROW195_ACCELERATOR_MIN_SUCCESS = (
    4,
    (1872, 6370, 57330, 458640),
    (1872, 57330, 458640),
    (18182, 35662, 65520),
)
EXPECTED_ROW195_ACCELERATOR_CLOSEST_HOSTILE = (
    (10192, 65520, 458640),
    0,
    (22696, 31152, 65520),
)
EXPECTED_ROW195_ACCELERATOR_ONE_HIGH = (
    599,
    156,
    "935127e640015f44655de13600eb68fe282d7c77fb7c398341d2ef5da4dcbe7e",
)
EXPECTED_ROW195_ACCELERATOR_CASE_STATE_SHA256 = "f5d7e673e8394bf35147c3ce8f08e3198fefdc8da9f1295c78af15828a0f7cbe"
EXPECTED_ROW195_ACCELERATOR_ZERO_HIGH = (
    109,
    "a1111e7538da5d8bebbc938d28db591396fb30f4504c20fe2b4a3023ccbad092",
)
EXPECTED_ROW195_ACCELERATOR_ONE_HIGH_AST_SHA256 = "001859942059725805de93d67634758a40734d8c7e55f3507ab516ddef8b94fa"
EXPECTED_ROW195_ACCELERATOR_TORSION_SUMMARY = (
    599,
    156,
    179,
    599,
    ((2, 467), (3, 89), (4, 7), (6, 36)),
    ((2, 469), (3, 124), (4, 6)),
    (
        1,
        (2, 5733, 6370, 152880),
        2,
        (237, 320),
        (2, 2, 2, 1, 37800, 37801, 0, 1, 1, 1, 86268),
    ),
    352552,
    57330,
)
EXPECTED_ROW195_ACCELERATOR_TORSION_SHA256 = "a67160cf481123639297f2e121776f6c17e552dbf04dee711bc9eb402a0201ac"
EXPECTED_ROW195_ACCELERATOR_HOSTILE_TWO_HIGH = (
    1,
    "0153f587817922eab7f9012644fe168b000233b3bded211229d3abb934eed194",
    0,
    (
        "denominator-two-measure",
        Fraction(3, 7),
        Fraction(3, 91),
        "7a70daaa1f25fc6fb9bfd7469ce3e5a18ba650296119789c81bc70fd32016b9c",
    ),
)
EXPECTED_ROW195_ACCELERATOR_TWO_HIGH_TAIL = (
    (
        (2, 2, 3920, 6370),
        (2, 2),
        351,
        Fraction(6457, 31302180),
    ),
)
EXPECTED_ROW195_ACCELERATOR_PARTITION_SHA256 = "a353502ccd9bcee61c6981948bbbba9af8a353d17dfc60e25f41c95efdbb6f38"
EXPECTED_ROW195_ACCELERATOR_PACKET_SHA256 = "8dcede7f9ee5b92eb9cc45c5b101e5bbd6f0862ca85942d8cc9ff12e73f94b2d"
EXPECTED_SCREEN_TOTALS = (596799, 279934, 298387, 18478)
EXPECTED_SCREEN_SHA256 = "748f1a4f9590d0eb72da8a03c3a14d5278a63f26eaba4c99e091f9057b61a0e0"
EXPECTED_TERMINAL_SUMMARY = (67, 18478, 16236, 113728, 17069, 4690, (('common-modulus-support', 76), ('denominator-two-measure', 1), ('located-torsion', 113727), ('top-h-weighted-common-source', 60), ('translated-cardinality', 1)), 113728, 1, 112514, 4232424, 630630, '9bbac6d9e329647d830097de5a02be273c85b785039c77bd96ee7c4a54a45ed9', 137, (('common-modulus-support', 76), ('denominator-two-measure', 1), ('top-h-weighted-common-source', 60)), (10456, 130, (32760, 458640), 24800, 65520), 403104, 1834560, ((195, '5811cb585e35831a9c58b2d0af35ab6db89e4503229e6420259d6ce6d4586c3c', '4cdb5acfda8b7db6e061b54ea2d71102c29a8694a7fce85f6a2c25e317a93ef9', ((351, 100776, 37987, 420652, '89236b7d22a6be06afecaaaba8f2e0be1346b1cab3c0b720314baa4014d87c36'),), ((50960, ((1, 760), (2, 7944), (3, 9848), (4, 11626), (5, 1616)), 5, 37987, 88947, (True, False, False, True, True, False, False, False, False)),)),), ((195, ((119368, 37800, 420839, '745925d13ce82fd6f85d5184de86a12a3cb2c370890d6f60ec3d95d354d13cae'), (78, '6b12fcbf3b4d7415cc4006d72ed57099bd635c09a084dda231e9643dad2c7803'), 'ce68e2ad70faf1bd0d5581ac490d84e345eea97479a42969927431c15bcc5985', (82160, 80687, 1473), (1473, 695, 295, 243, 157), (1565, 1408, 157), '3e3fe58c04ccbbbc1751d48bd199b76782c37589c58625f23767f5f870b8818d', '3ddbcf69c2c1521e1440973cdc46be0e85ca769a6ab68a247c8d00115340e372', 'e997f387c89563b9d87fca78e049f4484f0712a21879a7f4b6088461f704e947', (4, (1872, 6370, 57330, 458640), (1872, 57330, 458640), (18182, 35662, 65520)), ((10192, 65520, 458640), 0, (22696, 31152, 65520)), (599, 156, '935127e640015f44655de13600eb68fe282d7c77fb7c398341d2ef5da4dcbe7e', 'f5d7e673e8394bf35147c3ce8f08e3198fefdc8da9f1295c78af15828a0f7cbe', 109, 'a1111e7538da5d8bebbc938d28db591396fb30f4504c20fe2b4a3023ccbad092', '001859942059725805de93d67634758a40734d8c7e55f3507ab516ddef8b94fa'), (599, 156, 179, 599, ((2, 467), (3, 89), (4, 7), (6, 36)), ((2, 469), (3, 124), (4, 6)), (1, (2, 5733, 6370, 152880), 2, (237, 320), (2, 2, 2, 1, 37800, 37801, 0, 1, 1, 1, 86268)), 352552, 57330), 'a67160cf481123639297f2e121776f6c17e552dbf04dee711bc9eb402a0201ac', (1, '0153f587817922eab7f9012644fe168b000233b3bded211229d3abb934eed194', 0, ('denominator-two-measure', Fraction(3, 7), Fraction(3, 91), '7a70daaa1f25fc6fb9bfd7469ce3e5a18ba650296119789c81bc70fd32016b9c'), (((2, 2, 3920, 6370), (2, 2), 351, Fraction(6457, 31302180)),), (('denominator-two-measure', Fraction(3, 7), Fraction(3, 91)),), (('selfcontained_independent_source', 'fccc10392624bbdfc2083993ad51a423e8974c135b9bc635351304a71cb0de74'), ('selfcontained_independent_output', '55ad1da385d35f5b38fdc7de2d9916f54ee874b806496d4a7b1270cb526ad30c'), ('selfcontained_independent_semantic', '7a70daaa1f25fc6fb9bfd7469ce3e5a18ba650296119789c81bc70fd32016b9c'))), ('2a8584886131c177149055d208f5ccd3b356cc42dcdc84b6533773f6b800b935', '7d8644411a2d78ed40a6c593a7f977979803aac2e0f2427f788ee58da7c9489c', ((1, 65014), (2, 711), (3, 3))), (('selfcontained_accelerator_source', '00e1855dd428466c113d3d559e629baefcfdd9b0ea8ac732069fee3bb4ac74fe'), ('selfcontained_accelerator_output', 'b188805a09b484466f4c36395b2b426beba8ad470704d9c28414f92ec7a43f13'), ('selfcontained_accelerator_semantic', '5ad577c0a388c0021cda9506a8147a25732f36f801e7bcf746b856a5d9ff1070')), ('historical_nonloadbearing', 'a5cdcedc61c07714cacb2164e8793c655fe57ee80fd17ca82149b6710e254036', 314615), 'a353502ccd9bcee61c6981948bbbba9af8a353d17dfc60e25f41c95efdbb6f38')),))
EXPECTED_TERMINAL_SHA256 = "c39709a32e08e8c24b7d25259a85a85c2816d09df250a7fdab33feed1943dd89"
EXPECTED_SEMANTIC_SHA256 = "125200a407d2d055ff907097cabb9894817cae37f468b55e6a79b1c3d1a77920"


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def sha(path):
    payload = path.read_bytes().replace(b"\r\n", b"\n").replace(b"\r", b"\n")
    return hashlib.sha256(payload).hexdigest()


def digest(value):
    return hashlib.sha256(repr(value).encode("ascii")).hexdigest()


def digest_u32(values):
    state = hashlib.sha256()
    for value in values:
        state.update(int(value).to_bytes(4, "big"))
    return state.hexdigest()


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    require(spec is not None and spec.loader is not None, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def function_ast_hash(path, function_name):
    tree = ast.parse(path.read_text(encoding="utf-8"))
    matches = [
        node
        for node in tree.body
        if isinstance(node, ast.FunctionDef) and node.name == function_name
    ]
    require(len(matches) == 1, (path, function_name, len(matches)))
    node = matches[0]
    return hashlib.sha256(
        ast.dump(node, annotate_fields=True, include_attributes=False).encode("ascii")
    ).hexdigest(), node


def screen_worker(item):
    index, row = item
    prefix = load(f"wall110_screen_prefix_{index}", PREFIX)
    index2, screen, direct, legacy = prefix.screen_worker((index, row))
    require(index2 == index, (index, index2))
    screen = tuple(screen)
    require(screen[16] == screen[11], (index, "status replay"))
    require(direct + legacy == screen[11], (index, "Farkas partition"))
    require(screen[12] == len(screen[13]), (index, "residual length"))
    return index, screen, direct, legacy


def direct_cell_clean(cell, label, ruler):
    residue = (label * cell) % ruler
    return 14 * residue >= ruler and 14 * (residue + label) <= 13 * ruler


def direct_fixed_safe_cells(stream, labels):
    return tuple(
        cell
        for left, right in stream.ranges
        for cell in range(left, right)
        if all(direct_cell_clean(cell, label, stream.L) for label in (stream.first, *labels))
    )


def compact_fixed_safe_cells(stream, labels):
    """Rebuild the exact common-cell set with four bytes per cell."""
    require(stream.L < 2**32, (stream.L, "compact cell width"))
    return array(
        "I",
        (
            cell
            for left, right in stream.ranges
            for cell in range(left, right)
            if all(
                direct_cell_clean(cell, label, stream.L)
                for label in (stream.first, *labels)
            )
        ),
    )


def bit_is_set(bits, residue):
    return bool(bits[residue >> 3] & (1 << (residue & 7)))


def compact_torsion_pigeonhole(clean_cells, denominator):
    """Exact torsion replay using a denominator-bit support representation.

    This returns the identical least-order, least-residue witness selected by
    the inherited tuple/dictionary implementation, while its peak auxiliary
    support representation is exactly ceil(denominator/8) bytes.
    """
    bits = bytearray((denominator + 7) // 8)
    support = 0
    for cell in clean_cells:
        residue = cell % denominator
        byte = residue >> 3
        mask = 1 << (residue & 7)
        if not bits[byte] & mask:
            bits[byte] |= mask
            support += 1

    for qualifying_order in range(2, 8):
        if denominator % qualifying_order:
            continue
        quotient = denominator // qualifying_order
        if support <= quotient:
            continue

        # The inherited routine inserts cosets while scanning sorted residues
        # and selects the first crowded coset, then its first two residues.
        # Minimize the first present residue over quotient classes to replay
        # that ordering exactly without materializing a residue dictionary.
        selected = None
        for base in range(quotient):
            present = []
            for multiple in range(qualifying_order):
                residue = base + multiple * quotient
                if bit_is_set(bits, residue):
                    present.append(residue)
                    if len(present) == 2:
                        candidate = (present[0], present[1])
                        if selected is None or candidate[0] < selected[0]:
                            selected = candidate
                        break
        require(selected is not None, (denominator, qualifying_order, support))
        first_residue, second_residue = selected
        first_cell = None
        second_cell = None
        for cell in clean_cells:
            residue = cell % denominator
            if residue == first_residue and first_cell is None:
                first_cell = cell
            if residue == second_residue and second_cell is None:
                second_cell = cell
            if first_cell is not None and second_cell is not None:
                break
        require(first_cell is not None and second_cell is not None, selected)
        shift = (second_residue - first_residue) % denominator
        effective = denominator // gcd(denominator, shift)
        require(2 <= effective <= qualifying_order <= 7, (denominator, selected))
        primitive = tuple(
            min(unit, effective - unit)
            for unit in range(1, effective)
            if gcd(unit, effective) == 1
        )
        require(primitive and 7 * min(primitive) >= effective, (denominator, selected))
        return (
            qualifying_order,
            effective,
            support,
            quotient,
            first_cell,
            second_cell,
            first_residue,
            second_residue,
            shift,
            min(primitive),
            len(clean_cells),
        )
    return (
        None,
        None,
        support,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        len(clean_cells),
    )


def verify_torsion_witness(cells, denominator, witness, ruler, labels):
    qorder, effective, support, quotient, first, second, rfirst, rsecond, shift, phase, ncell = witness
    require(qorder is not None, witness)
    require(ncell == len(cells) and first in cells and second in cells, witness)
    require((first % denominator, second % denominator) == (rfirst, rsecond), witness)
    require((second - first) % denominator == shift, witness)
    require(effective == denominator // gcd(denominator, shift), witness)
    require(2 <= effective <= qorder <= 7, witness)
    primitive = tuple(
        min(unit, effective - unit)
        for unit in range(1, effective)
        if gcd(unit, effective) == 1
    )
    require(primitive and phase == min(primitive), witness)
    require(7 * phase >= effective, witness)
    require(support > quotient, witness)
    require(
        all(direct_cell_clean(cell, label, ruler) for cell in (first, second) for label in labels),
        (denominator, "witness cell not directly clean"),
    )
    return support - quotient


def compact_common_capacity(cells, denominators):
    """Inherited common-modulus support bound with one bit per residue."""
    modulus = lcm(*denominators)
    bits = bytearray((modulus + 7) // 8)
    support = 0
    for cell in cells:
        residue = cell % modulus
        byte = residue >> 3
        mask = 1 << (residue & 7)
        if not bits[byte] & mask:
            bits[byte] |= mask
            support += 1
    capacities = tuple(
        ((denominator + 6) // 7) * (modulus // denominator)
        for denominator in denominators
    )
    return modulus, support, capacities, support - sum(capacities), len(bits)


def compact_top_h_majorant(cells, denominator):
    """Bound every unit/window pullback by the h heaviest residue fibres."""
    weights = array("I", (0 for _ in range(denominator)))
    require(weights.itemsize == 4, weights.itemsize)
    for cell in cells:
        weights[cell % denominator] += 1
    require(sum(weights) == len(cells), (denominator, sum(weights)))
    histogram = Counter(weights)
    width = (denominator + 6) // 7
    remaining = width
    upper = 0
    for weight in sorted(histogram, reverse=True):
        take = min(remaining, histogram[weight])
        upper += take * weight
        remaining -= take
        if remaining == 0:
            break
    require(remaining == 0, (denominator, width, remaining))
    histogram_packet = tuple(sorted(histogram.items()))
    return (
        denominator,
        len(cells),
        sum(weight > 0 for weight in weights),
        width,
        upper,
        max(weights),
        digest(histogram_packet),
        weights.itemsize * len(weights),
    )


def row195_accelerator_state_partition(stream, engine, residual, needed):
    """Split row 195 into 1,408 uniform states and 157 hostile states."""
    require(
        (stream.body, stream.L, stream.first, stream.first_d)
        == ((1, 5, 8, 9, 13, 14), 458640, LEVEL, 6370),
        (stream.body, stream.L, stream.first, stream.first_d),
    )
    needed = tuple(sorted(needed))
    require(
        (len(needed), digest(needed)) == EXPECTED_ROW195_ACCELERATOR_NEEDED,
        (len(needed), digest(needed)),
    )

    _trials, states, _ray_checks, _signs = engine.ray.ray_quotient_states(stream)
    state_keys = tuple(sorted(states))
    states_hash = digest(state_keys)
    require(states_hash == EXPECTED_ROW195_ACCELERATOR_STATES_SHA256, states_hash)
    suffix_bijection = []
    for divisors in state_keys:
        suffix = tuple(sorted(engine.suffix_slots(divisors, stream.first_d)))
        require(tuple(sorted((stream.first_d, *suffix))) == divisors, (divisors, suffix))
        suffix_bijection.append((divisors, suffix))
    suffix_bijection = tuple(suffix_bijection)
    require(
        len({suffix for _divisors, suffix in suffix_bijection}) == len(states),
        "row195 suffix collision",
    )
    suffix_hash = digest(suffix_bijection)
    require(suffix_hash == EXPECTED_ROW195_ACCELERATOR_SUFFIX_SHA256, suffix_hash)
    first_multiplicities = tuple(
        sorted(Counter(ds.count(stream.first_d) for ds in state_keys).items())
    )
    require(
        first_multiplicities == EXPECTED_ROW195_ACCELERATOR_FIRST_MULTIPLICITIES,
        first_multiplicities,
    )

    first_cells = compact_fixed_safe_cells(stream, ())
    require(first_cells and first_cells.itemsize == 4, (len(first_cells), first_cells.itemsize))
    first_source = (
        len(first_cells),
        first_cells[0],
        first_cells[-1],
        digest_u32(first_cells),
    )
    require(first_source == EXPECTED_ROW195_ACCELERATOR_FIRST_SOURCE, first_source)
    top_h = tuple(compact_top_h_majorant(first_cells, denominator) for denominator in needed)
    top_h_hash = digest(top_h)
    require(top_h_hash == EXPECTED_ROW195_ACCELERATOR_TOP_H_SHA256, top_h_hash)
    upper = {record[0]: record[4] for record in top_h}
    triple_records = []
    failed_triples = []
    for suffix in combinations_with_replacement(needed, 3):
        margin = len(first_cells) - sum(upper[denominator] for denominator in suffix)
        record = (suffix, margin, tuple(upper[denominator] for denominator in suffix))
        triple_records.append(record)
        if margin <= 0:
            failed_triples.append(record)
    triple_records = tuple(triple_records)
    failed_triples = tuple(failed_triples)
    triple_partition = (
        len(triple_records),
        len(triple_records) - len(failed_triples),
        len(failed_triples),
    )
    require(
        triple_partition == EXPECTED_ROW195_ACCELERATOR_TRIPLE_PARTITION,
        triple_partition,
    )
    require(
        digest(triple_records) == EXPECTED_ROW195_ACCELERATOR_TRIPLES_SHA256,
        digest(triple_records),
    )
    require(
        digest(failed_triples) == EXPECTED_ROW195_ACCELERATOR_FAILED_SHA256,
        digest(failed_triples),
    )

    record_by_suffix = {record[0]: record for record in triple_records}
    failed_suffixes = {record[0] for record in failed_triples}
    failed_state_by_suffix = {
        tuple(sorted((stream.first_d, *record[0]))): record
        for record in failed_triples
    }
    require(len(failed_state_by_suffix) == len(failed_triples), "failed state collision")
    attained_failed = {
        divisors: states[divisors]
        for divisors in sorted(failed_state_by_suffix)
        if divisors in states
    }
    attained_packet = tuple((ds, attained_failed[ds]) for ds in sorted(attained_failed))
    require(
        digest(attained_packet) == EXPECTED_ROW195_ACCELERATOR_ATTAINED_SHA256,
        digest(attained_packet),
    )
    crude, status, screened_failed = engine.exact_common_status_screen(stream, attained_failed)
    crude_keys = set(crude)
    status_keys = set(status)
    screened_failed_keys = set(screened_failed)
    require(not (crude_keys & status_keys), "row195 crude/status overlap")
    require(not (crude_keys & screened_failed_keys), "row195 crude/residual overlap")
    require(not (status_keys & screened_failed_keys), "row195 status/residual overlap")
    require(
        set(attained_failed) == crude_keys | status_keys | screened_failed_keys,
        "row195 attained-failed partition",
    )
    require(
        digest(tuple(sorted(crude.items()))) == EXPECTED_ROW195_ACCELERATOR_CRUDE_SHA256,
        digest(tuple(sorted(crude.items()))),
    )
    require(
        digest(tuple(sorted(status.items()))) == EXPECTED_ROW195_ACCELERATOR_STATUS_SHA256,
        digest(tuple(sorted(status.items()))),
    )
    residual_failed = tuple(
        divisors
        for divisors in residual
        if tuple(sorted(engine.suffix_slots(divisors, stream.first_d))) in failed_suffixes
    )
    require(tuple(screened_failed) == residual_failed, "row195 subset/full screen mismatch")
    failed_partition = (
        len(failed_triples),
        len(attained_failed),
        len(crude),
        len(status),
        len(residual_failed),
    )
    require(
        failed_partition == EXPECTED_ROW195_ACCELERATOR_FAILED_PARTITION,
        failed_partition,
    )
    residual_failed_set = set(residual_failed)
    strict_residual = tuple(
        divisors for divisors in residual if divisors not in residual_failed_set
    )
    residual_partition = (len(residual), len(strict_residual), len(residual_failed))
    require(
        residual_partition == EXPECTED_ROW195_ACCELERATOR_RESIDUAL_PARTITION,
        residual_partition,
    )
    strict_residual_hash = digest(strict_residual)
    if EXPECTED_ROW195_ACCELERATOR_STRICT_RESIDUAL_SHA256 is not None:
        require(
            strict_residual_hash == EXPECTED_ROW195_ACCELERATOR_STRICT_RESIDUAL_SHA256,
            strict_residual_hash,
        )
    strict_certificates = tuple(
        (
            divisors,
            tuple(sorted(engine.suffix_slots(divisors, stream.first_d))),
            record_by_suffix[
                tuple(sorted(engine.suffix_slots(divisors, stream.first_d)))
            ][1],
            record_by_suffix[
                tuple(sorted(engine.suffix_slots(divisors, stream.first_d)))
            ][2],
        )
        for divisors in strict_residual
    )
    require(all(record[2] > 0 for record in strict_certificates), "nonstrict success")
    strict_certificate_hash = digest(strict_certificates)
    require(
        strict_certificate_hash == EXPECTED_ROW195_ACCELERATOR_STRICT_CERTIFICATE_SHA256,
        strict_certificate_hash,
    )
    min_success = min(
        (record[2], record[0], record[1], record[3])
        for record in strict_certificates
    )
    require(min_success == EXPECTED_ROW195_ACCELERATOR_MIN_SUCCESS, min_success)
    closest_hostile = max(failed_triples, key=lambda record: (record[1], record[0]))
    require(
        closest_hostile == EXPECTED_ROW195_ACCELERATOR_CLOSEST_HOSTILE,
        closest_hostile,
    )
    hostile_hash = digest(residual_failed)
    require(hostile_hash == EXPECTED_ROW195_ACCELERATOR_HOSTILE_SHA256, hostile_hash)

    one_high_ast_hash, one_high_ast = function_ast_hash(TORSION, "one_high_cases")
    require(
        one_high_ast_hash == EXPECTED_ROW195_ACCELERATOR_ONE_HIGH_AST_SHA256,
        one_high_ast_hash,
    )
    state_loops = [
        node
        for node in ast.walk(one_high_ast)
        if isinstance(node, ast.For)
        and isinstance(node.target, ast.Name)
        and node.target.id == "ds"
        and isinstance(node.iter, ast.Name)
        and node.iter.id == "residuals"
    ]
    case_adds = [
        node
        for node in ast.walk(one_high_ast)
        if isinstance(node, ast.Call)
        and isinstance(node.func, ast.Attribute)
        and isinstance(node.func.value, ast.Name)
        and node.func.value.id == "cases"
        and node.func.attr == "add"
    ]
    require(len(state_loops) == 1 and case_adds, "one-high locality shape")
    require(
        all(
            len(call.args) == 1
            and isinstance(call.args[0], ast.Tuple)
            and call.args[0].elts
            and isinstance(call.args[0].elts[0], ast.Name)
            and call.args[0].elts[0].id == "ds"
            for call in case_adds
        ),
        "one-high case not tagged by current state",
    )
    prefix_packet = (
        first_source,
        (len(needed), digest(needed)),
        top_h_hash,
        triple_partition,
        failed_partition,
        residual_partition,
        strict_residual_hash,
        strict_certificate_hash,
        hostile_hash,
        min_success,
        closest_hostile,
    )
    state_map_controls = (states_hash, suffix_hash, first_multiplicities)
    return residual_failed, prefix_packet, state_map_controls, one_high_ast_hash


def clipped_periodic_orbit_obstruction(cells, stream, low_label):
    """Freeze the carrier loss in the tempting row-195 L/9 orbit."""
    require(stream.L == 458640 and low_label == 351, (stream.L, low_label))
    period = stream.L // 9
    require((stream.first * period) % stream.L == 0, "first period")
    require((low_label * period) % stream.L == 0, "low period")
    counts = array("B", (0 for _ in range(period)))
    for cell in cells:
        counts[cell % period] += 1
    histogram = tuple(
        sorted(
            (multiplicity, number)
            for multiplicity, number in Counter(counts).items()
            if multiplicity
        )
    )
    require(max(counts) == 5 < 9, max(counts))
    source = 37987
    orbit = tuple(source + offset * period for offset in range(9))
    present = bytearray(stream.L)
    for cell in cells:
        present[cell] = 1
    mask = tuple(bool(present[cell]) for cell in orbit)
    require(
        mask == (True, False, False, True, True, False, False, False, False),
        mask,
    )
    target = source + period
    require(
        direct_cell_clean(target, stream.first, stream.L)
        and direct_cell_clean(target, low_label, stream.L),
        "intrinsic period target",
    )
    require(
        not any(left <= target < right for left, right in stream.ranges),
        "period target unexpectedly in carrier",
    )
    return period, histogram, max(counts), source, target, mask


def compact_two_high_closure(index, stream, engine, cases):
    """Stream scalar two-high cases through common-source certificates."""
    if not cases:
        return (), (), None, 0, 0, (), ()
    cases_by_low = {}
    for case_index, case in enumerate(cases):
        divisors, high_ds, low_label, excess = case
        require(low_label is not None and len(high_ds) == 2, (index, case_index, case))
        cases_by_low.setdefault(low_label, []).append((case_index, case))

    closure_by_index = {}
    mechanism_counts = Counter()
    weakest_top_h = None
    peak_cell_bytes = 0
    peak_aux_bytes = 0
    cell_packets = []
    orbit_obstructions = []
    for low_label, indexed_cases in sorted(cases_by_low.items()):
        cells = compact_fixed_safe_cells(stream, (low_label,))
        require(cells and cells.itemsize == 4, (index, low_label, len(cells)))
        peak_cell_bytes = max(peak_cell_bytes, cells.itemsize * len(cells))
        cell_packet = (
            low_label,
            len(cells),
            cells[0],
            cells[-1],
            digest_u32(cells),
        )
        cell_packets.append(cell_packet)
        if index == 195:
            orbit_obstructions.append(
                clipped_periodic_orbit_obstruction(cells, stream, low_label)
            )
        top_h_cache = {}
        for case_index, case in indexed_cases:
            divisors, high_ds, _low_label, excess = case
            common = compact_common_capacity(cells, high_ds)
            peak_aux_bytes = max(peak_aux_bytes, common[-1])
            if common[3] > 0:
                certificate = ("common-modulus-support", common[3])
            elif high_ds == (2, 2):
                certificate = (
                    "denominator-two-measure",
                    Fraction(3, 7),
                    Fraction(3, 91),
                )
            else:
                for denominator in high_ds:
                    if denominator not in top_h_cache:
                        top_h_cache[denominator] = compact_top_h_majorant(
                            cells, denominator
                        )
                        peak_aux_bytes = max(
                            peak_aux_bytes,
                            top_h_cache[denominator][-1],
                        )
                majorants = tuple(top_h_cache[denominator] for denominator in high_ds)
                lower = len(cells) - sum(record[4] for record in majorants)
                require(lower > 0, (index, case_index, high_ds, lower, majorants))
                certificate = (
                    "top-h-weighted-common-source",
                    lower,
                    majorants,
                )
                candidate = (
                    lower,
                    case_index,
                    high_ds,
                    majorants[0][4],
                    majorants[1][4],
                )
                if weakest_top_h is None or candidate < weakest_top_h:
                    weakest_top_h = candidate
            mechanism_counts[certificate[0]] += 1
            closure_by_index[case_index] = (
                case_index,
                divisors,
                high_ds,
                low_label,
                excess,
                common,
                certificate,
            )
    closures = tuple(closure_by_index[position] for position in range(len(cases)))
    require(len(closures) == len(cases), (index, len(closures), len(cases)))
    return (
        closures,
        tuple(sorted(mechanism_counts.items())),
        weakest_top_h,
        peak_cell_bytes,
        peak_aux_bytes,
        tuple(cell_packets),
        tuple(orbit_obstructions),
    )


def terminal_worker(item):
    index, body, residual = item
    prefix = load(f"wall110_terminal_prefix_{index}", PREFIX)
    audit = load(f"wall110_terminal_audit_{index}", prefix.AUDIT_SOURCE)
    natural = audit.load(f"wall110_terminal_natural_{index}", audit.NATURAL_SOURCE)
    engine = audit.status_engine(natural, f"wall110_terminal_engine_{index}")
    engine.FIRST = LEVEL
    engine.ray.FIRST = LEVEL
    stream = engine.ray.Stream(body)
    require(stream.first == LEVEL < stream.high_floor, (index, "high gate inactive"))
    needed = {
        denominator
        for divisors in residual
        for denominator in engine.suffix_slots(divisors, stream.first_d)
    }
    one_high_residual = residual
    accelerator_prefix = None
    accelerator_state_controls = None
    accelerator_ast_hash = None
    if index == 195:
        (
            one_high_residual,
            accelerator_prefix,
            accelerator_state_controls,
            accelerator_ast_hash,
        ) = row195_accelerator_state_partition(stream, engine, residual, needed)
    low, high, sign_census, recurrence_checks = engine.build_literal_tables(stream, needed)
    two_gap, two_witness = engine.duplicate_two_high_gap(stream, residual, low, high)
    two_cases = ()
    three_high_masks = 0
    if two_gap <= 0:
        require(index in EXPECTED_NONPOSITIVE_ROWS, (index, "unexpected nonpositive gap", two_gap))
        two_cases, three_high_masks = prefix.two_high_scalar_cases(
            engine, stream, residual, low, high
        )
        require(not three_high_masks, (index, "three-high scalar survivor"))
    zero_high = engine.zero_high_scalar_passes(stream, one_high_residual, low)
    one_high = engine.one_high_cases(stream, one_high_residual, low, high)
    hostile_two_cases = ()
    if index == 195:
        hostile_states = set(one_high_residual)
        hostile_two_cases = tuple(case for case in two_cases if case[0] in hostile_states)
        require(hostile_two_cases == EXPECTED_ROW195_ACCELERATOR_TWO_HIGH_TAIL, hostile_two_cases)
        require(not three_high_masks, three_high_masks)
        require(
            (len(zero_high), digest(zero_high))
            == EXPECTED_ROW195_ACCELERATOR_ZERO_HIGH,
            (len(zero_high), digest(zero_high)),
        )
        require(
            (len(one_high), len({case[0] for case in one_high}), digest(one_high))
            == EXPECTED_ROW195_ACCELERATOR_ONE_HIGH,
            (len(one_high), len({case[0] for case in one_high}), digest(one_high)),
        )

    witness_cache = {}
    certificate_cache = {}
    cell_counts = {}
    translated_labels = set()
    legacy_parity_checks = 0
    peak_compact_cell_bytes = 0
    peak_support_bytes = 0

    keys_by_labels = {}
    for _divisors, high_d, low_rows, _excess in one_high:
        labels = tuple(sorted(label for _d, label in low_rows))
        keys_by_labels.setdefault(labels, set()).add(high_d)
    for labels, denominators in sorted(keys_by_labels.items()):
        cells = compact_fixed_safe_cells(stream, labels)
        require(cells, (index, labels, "no fixed-safe cell"))
        require(cells.itemsize == 4, cells.itemsize)
        cell_counts[labels] = len(cells)
        peak_compact_cell_bytes = max(
            peak_compact_cell_bytes, cells.itemsize * len(cells)
        )

        # Every non-heavy row except the bounded row-195 accelerator is checked
        # byte-for-byte against the inherited tuple reconstruction and witness
        # selector.  Row 195 is controlled by the canonical self-contained
        # accelerator; the three heavy rows retain direct literal checks and
        # are pinned to their frozen family screen/terminal packets below.
        legacy_cells = None
        if stream.L < 5_000_000 and index != 195:
            legacy_cells = engine.fixed_safe_cells(stream, labels)
            require(tuple(cells) == legacy_cells, (index, labels, "cell parity"))

        for high_d in sorted(denominators):
            key = (labels, high_d)
            witness = compact_torsion_pigeonhole(cells, high_d)
            witness_cache[key] = witness
            peak_support_bytes = max(peak_support_bytes, (high_d + 7) // 8)
            if legacy_cells is not None:
                inherited_witness = engine.torsion_pigeonhole(legacy_cells, high_d)
                require(witness == inherited_witness, (index, key, "witness parity"))
                legacy_parity_checks += 1
            if witness[0] is not None:
                surplus = verify_torsion_witness(
                    cells, high_d, witness, stream.L, (stream.first, *labels)
                )
                certificate_cache[key] = ("located-torsion", surplus, witness)
                continue

            # This branch is rare (one known row-181 key).  Reconstruct the
            # tuple independently here, not in the compact producer, and pin
            # the complete translated residue support.
            direct_cells = direct_fixed_safe_cells(stream, labels)
            require(tuple(cells) == direct_cells, (index, labels, "direct cell reconstruction"))
            support = tuple(sorted({cell % high_d for cell in direct_cells}))
            capacity = (high_d + 6) // 7
            margin = len(support) - capacity
            require(margin > 0, (index, high_d, labels, margin))
            certificate_cache[key] = (
                "translated-cardinality",
                margin,
                len(support),
                capacity,
                digest(direct_cells),
                digest(support),
                direct_cells[0],
                direct_cells[-1],
            )
            translated_labels.add(labels)

    cases = []
    mechanisms = Counter()
    qualifying = Counter()
    effective = Counter()
    minimum_torsion = None
    minimum_translated = None
    for divisors, high_d, low_rows, excess in one_high:
        labels = tuple(sorted(label for _d, label in low_rows))
        key = (labels, high_d)
        witness = witness_cache[key]
        certificate = certificate_cache[key]
        if witness[0] is not None:
            surplus = certificate[1]
            qualifying[witness[0]] += 1
            effective[witness[1]] += 1
            candidate = (surplus, index, divisors, high_d, labels, witness)
            if minimum_torsion is None or candidate < minimum_torsion:
                minimum_torsion = candidate
        else:
            margin, support_size, capacity = certificate[1:4]
            candidate = (margin, index, divisors, high_d, labels, support_size, capacity)
            if minimum_translated is None or candidate < minimum_translated:
                minimum_translated = candidate
        mechanisms[certificate[0]] += 1
        cases.append(
            (divisors, high_d, low_rows, excess, labels, cell_counts[labels], certificate)
        )

    (
        two_closures,
        two_mechanisms,
        weakest_top_h,
        peak_two_cell_bytes,
        peak_two_aux_bytes,
        two_cell_packets,
        orbit_obstructions,
    ) = compact_two_high_closure(index, stream, engine, two_cases)
    for mechanism, count in two_mechanisms:
        mechanisms[mechanism] += count
    if index == 195:
        require(
            (len(two_cases), digest(two_cases)) == EXPECTED_ROW195_TWO_CASES,
            (len(two_cases), digest(two_cases)),
        )
        require(two_mechanisms == EXPECTED_ROW195_TWO_MECHANISMS, two_mechanisms)
        require(weakest_top_h == EXPECTED_ROW195_WEAKEST_TOP_H, weakest_top_h)
        require(digest(two_closures) == EXPECTED_ROW195_CLOSURES_SHA256, digest(two_closures))
        require(two_cell_packets == EXPECTED_ROW195_CELL_PACKET, two_cell_packets)
        require(
            orbit_obstructions == EXPECTED_ROW195_PERIOD_OBSTRUCTION,
            orbit_obstructions,
        )
        hostile_two_closures = tuple(
            record for record in two_closures if record[1] in set(one_high_residual)
        )
        require(len(hostile_two_closures) == 1, hostile_two_closures)
        require(
            hostile_two_closures[0][-1]
            == ("denominator-two-measure", Fraction(3, 7), Fraction(3, 91)),
            hostile_two_closures,
        )
    else:
        hostile_two_closures = ()

    case_packet = tuple(cases)
    passports = tuple(sorted({case[0] for case in one_high}))
    low_pairs = tuple(sorted({case[4] for case in case_packet}))
    require(peak_compact_cell_bytes <= 4 * stream.L, (index, peak_compact_cell_bytes))
    require(
        peak_support_bytes <= (stream.L + 7) // 8,
        (index, peak_support_bytes),
    )
    accelerator_packet = None
    if index == 195:
        case_state_counts = tuple(sorted(Counter(case[0] for case in one_high).items()))
        require(
            digest(case_state_counts) == EXPECTED_ROW195_ACCELERATOR_CASE_STATE_SHA256,
            digest(case_state_counts),
        )
        accelerator_closures = tuple(
            (
                divisors,
                high_d,
                low_rows,
                excess,
                labels,
                cell_count,
                certificate[1],
                certificate[2],
            )
            for (
                divisors,
                high_d,
                low_rows,
                excess,
                labels,
                cell_count,
                certificate,
            ) in case_packet
        )
        require(
            all(record[6] > 0 and record[7][0] is not None for record in accelerator_closures),
            "row195 non-torsion accelerated case",
        )
        accelerator_minimum = min(
            (record[6], record[0], record[1], record[4], record[7])
            for record in accelerator_closures
        )
        accelerator_torsion_summary = (
            len(one_high),
            len(passports),
            len(keys_by_labels),
            len(witness_cache),
            tuple(sorted(qualifying.items())),
            tuple(sorted(effective.items())),
            accelerator_minimum,
            peak_compact_cell_bytes,
            peak_support_bytes,
        )
        require(
            accelerator_torsion_summary == EXPECTED_ROW195_ACCELERATOR_TORSION_SUMMARY,
            accelerator_torsion_summary,
        )
        accelerator_closure_hash = digest(accelerator_closures)
        require(
            accelerator_closure_hash == EXPECTED_ROW195_ACCELERATOR_TORSION_SHA256,
            accelerator_closure_hash,
        )
        reduced_one_high_packet = (
            len(one_high),
            len(passports),
            digest(one_high),
            digest(case_state_counts),
            len(zero_high),
            digest(zero_high),
            accelerator_ast_hash,
        )
        hostile_two_summary = (
            len(hostile_two_cases),
            digest(hostile_two_cases),
            three_high_masks,
            (
                "denominator-two-measure",
                Fraction(3, 7),
                Fraction(3, 91),
                EXPECTED_ROW195_CONTROL_HASHES[2][1],
            ),
        )
        require(
            hostile_two_summary == EXPECTED_ROW195_ACCELERATOR_HOSTILE_TWO_HIGH,
            hostile_two_summary,
        )
        hostile_two_packet = (
            *hostile_two_summary,
            hostile_two_cases,
            tuple(record[-1] for record in hostile_two_closures),
            EXPECTED_ROW195_CONTROL_HASHES[:3],
        )
        accelerator_packet = (
            *accelerator_prefix,
            reduced_one_high_packet,
            accelerator_torsion_summary,
            accelerator_closure_hash,
            hostile_two_packet,
            accelerator_state_controls,
            EXPECTED_ROW195_CONTROL_HASHES[3:],
            ("historical_nonloadbearing", "a5cdcedc61c07714cacb2164e8793c655fe57ee80fd17ca82149b6710e254036", 314615),
            EXPECTED_ROW195_ACCELERATOR_PARTITION_SHA256,
        )
        if EXPECTED_ROW195_ACCELERATOR_PACKET_SHA256 is not None:
            require(
                digest(accelerator_packet) == EXPECTED_ROW195_ACCELERATOR_PACKET_SHA256,
                digest(accelerator_packet),
            )
    terminal_repair_record = (
        index,
        body,
        stream.L,
        stream.high_floor,
        stream.first_d,
        len(residual),
        tuple(sorted(needed)),
        two_gap,
        two_witness,
        len(zero_high),
        len(one_high),
        len(low_pairs),
        sign_census,
        recurrence_checks,
        len(witness_cache),
        tuple(sorted(qualifying.items())),
        tuple(sorted(effective.items())),
        digest(one_high),
        digest(tuple(sorted(witness_cache.items()))),
        len(two_cases),
        sum(count for mechanism, count in two_mechanisms if mechanism == "common-modulus-support"),
        sum(count for mechanism, count in two_mechanisms if mechanism != "common-modulus-support"),
        digest(two_cases),
        digest(tuple(record for record in two_closures if record[-1][0] == "common-modulus-support")),
        digest(tuple(record for record in two_closures if record[-1][0] != "common-modulus-support")),
        weakest_top_h,
        tuple(
            record[-1][1]
            for record in two_closures
            if record[-1][0] == "denominator-two-measure"
        ),
    )
    return (
        index,
        body,
        stream.L,
        stream.high_floor,
        stream.first_d,
        len(residual),
        digest(residual),
        two_gap,
        digest(two_witness),
        len(zero_high),
        digest(zero_high),
        len(one_high),
        len(passports),
        len(low_pairs),
        tuple(sorted(mechanisms.items())),
        len(witness_cache),
        len(translated_labels),
        minimum_torsion,
        minimum_translated,
        digest(case_packet),
        sign_census,
        recurrence_checks,
        terminal_repair_record,
        legacy_parity_checks,
        peak_compact_cell_bytes,
        peak_support_bytes,
        len(two_cases),
        two_mechanisms,
        weakest_top_h,
        digest(two_cases),
        digest(two_closures),
        peak_two_cell_bytes,
        peak_two_aux_bytes,
        two_cell_packets,
        orbit_obstructions,
        accelerator_packet,
    )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--processes", type=int, default=4)
    args = parser.parse_args()
    require(args.processes >= 1, args.processes)
    for path, expected in EXPECTED_HASHES.items():
        require(sha(path) == expected, (path, "dependency changed", sha(path)))
    require(
        "semantic_sha256=436c5160c5e9d3dcaa0f4f3dc4104450670a67582ab5e375c8751dc0ea82c93f"
        in THM3378_OUT.read_text(encoding="utf-8"),
        "THM3378 semantic pin",
    )
    require(
        "semantic_sha256=7a70daaa1f25fc6fb9bfd7469ce3e5a18ba650296119789c81bc70fd32016b9c"
        in ROW195_AUDIT_OUT.read_text(encoding="utf-8"),
        "row195 independent semantic pin",
    )
    require(
        f"semantic_sha256={EXPECTED_ROW195_ACCELERATOR_SEMANTIC}"
        in ROW195_ACCELERATOR_OUT.read_text(encoding="utf-8"),
        "row195 accelerator semantic pin",
    )
    syntax = ast.parse(Path(__file__).read_text(encoding="utf-8"))
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(syntax)), "assert node")
    require(
        not any(
            isinstance(node, ast.Constant) and isinstance(node.value, float)
            for node in ast.walk(syntax)
        ),
        "float literal",
    )

    chain_natural = load("wall110_chain_natural", NATURAL)
    wrapper_chain = (
        (chain_natural, "SOURCE_3139", SOURCE3139),
        (load("wall110_chain_3139", SOURCE3139), "SOURCE_3114", SOURCE3114),
        (load("wall110_chain_3114", SOURCE3114), "SOURCE_3113", SOURCE3113),
        (load("wall110_chain_3113", SOURCE3113), "SOURCE_3111", SOURCE3111),
        (load("wall110_chain_3111", SOURCE3111), "SOURCE_3109", SOURCE3109),
        (load("wall110_chain_3109", SOURCE3109), "SOURCE_3106", SOURCE3106),
    )
    for module, attribute, expected_path in wrapper_chain:
        require(
            Path(getattr(module, attribute)).resolve() == expected_path,
            (module.__file__, attribute, getattr(module, attribute), expected_path),
        )
    chain_base = load("wall110_chain_3106", SOURCE3106)
    require(Path(chain_base.__file__).resolve() == SOURCE3106, chain_base.__file__)

    prefix = load("wall110_prefix", PREFIX)
    rows, components, inherited, ranked, old_prefix, old_indices, old_hash = prefix.reconstruct_queue()
    census = (len(rows), sum(row[3] for row in rows), sum(not row[3] for row in rows))
    require(census == EXPECTED_ATLAS, census)
    thm3361_family = ranked[prefix.PREFIX_FAMILIES]
    require(thm3361_family[4] == THM3361_CLOSED, thm3361_family)
    families = tuple(ranked[prefix.PREFIX_FAMILIES + 1 :])
    indices = tuple(index for family in families for index in family[4])
    require((len(indices), len(families)) == EXPECTED_QUEUE, (len(indices), len(families)))
    require(digest(indices) == EXPECTED_INDICES_SHA256, digest(indices))
    live_indices = tuple(index for index in indices if index != 94)
    require((len(live_indices), len(families)) == EXPECTED_LIVE_QUEUE, (len(live_indices), len(families)))
    require(digest(live_indices) == EXPECTED_LIVE_INDICES_SHA256, digest(live_indices))
    require(tuple(index for index in indices if index not in live_indices) == (94,), "live subtraction")
    families_hash = digest(families)
    if EXPECTED_FAMILIES_SHA256 is not None:
        require(families_hash == EXPECTED_FAMILIES_SHA256, families_hash)

    tasks = tuple((index, rows[index]) for index in indices)
    if args.processes == 1:
        screens = tuple(screen_worker(task) for task in tasks)
    else:
        parallel_tasks = tuple(task for task in tasks if task[1][1] < 5_000_000)
        heavy_tasks = tuple(task for task in tasks if task[1][1] >= 5_000_000)
        with mp.get_context("spawn").Pool(min(args.processes, len(parallel_tasks))) as pool:
            parallel_screens = tuple(pool.imap(screen_worker, parallel_tasks, chunksize=1))
        heavy_screens = tuple(screen_worker(task) for task in heavy_tasks)
        by_result = {record[0]: record for record in (*parallel_screens, *heavy_screens)}
        screens = tuple(by_result[index] for index in indices)
    require(tuple(record[0] for record in screens) == indices, "screen order")
    screen_totals = tuple(sum(record[1][pos] for record in screens) for pos in (9, 10, 11, 12))
    require(screen_totals[0] == sum(screen_totals[1:]), screen_totals)
    screen_hash = digest(screens)
    heavy_family_screens = tuple(
        record for record in screens if record[0] in EXPECTED_HEAVY_INDICES
    )
    require(
        tuple(record[0] for record in heavy_family_screens) == EXPECTED_HEAVY_INDICES,
        tuple(record[0] for record in heavy_family_screens),
    )
    require(
        digest(heavy_family_screens) == EXPECTED_HEAVY_SCREEN_SHA256,
        ("heavy screen packet changed", digest(heavy_family_screens)),
    )
    if EXPECTED_SCREEN_TOTALS is not None:
        require(screen_totals == EXPECTED_SCREEN_TOTALS, screen_totals)
    if EXPECTED_SCREEN_SHA256 is not None:
        require(screen_hash == EXPECTED_SCREEN_SHA256, screen_hash)

    by_screen = {record[0]: record for record in screens}
    terminal_tasks = tuple(
        (index, rows[index][0], tuple(record[1][13]))
        for index, record in sorted(by_screen.items())
        if record[1][12]
    )
    if args.processes == 1:
        terminals = tuple(terminal_worker(task) for task in terminal_tasks)
    else:
        parallel_terminal_tasks = tuple(
            task for task in terminal_tasks if rows[task[0]][1] < 5_000_000
        )
        heavy_terminal_tasks = tuple(
            task for task in terminal_tasks if rows[task[0]][1] >= 5_000_000
        )
        with mp.get_context("spawn").Pool(
            min(args.processes, len(parallel_terminal_tasks))
        ) as pool:
            parallel_terminals = tuple(
                pool.imap(terminal_worker, parallel_terminal_tasks, chunksize=1)
            )
        heavy_terminals = tuple(terminal_worker(task) for task in heavy_terminal_tasks)
        by_result = {
            record[0]: record for record in (*parallel_terminals, *heavy_terminals)
        }
        terminals = tuple(by_result[task[0]] for task in terminal_tasks)
    require(
        tuple(record[0] for record in terminals) == tuple(task[0] for task in terminal_tasks),
        "terminal order",
    )
    terminal_hash = digest(terminals)
    heavy_inherited_packet = tuple(
        (record[0], True, record[22])
        for record in terminals
        if record[0] in EXPECTED_HEAVY_INDICES
    )
    require(
        tuple(record[0] for record in heavy_inherited_packet) == EXPECTED_HEAVY_INDICES,
        tuple(record[0] for record in heavy_inherited_packet),
    )
    heavy_inherited_hash = digest(heavy_inherited_packet)
    require(
        heavy_inherited_hash == EXPECTED_HEAVY_TERMINAL_SHA256,
        ("heavy inherited terminal packet changed", heavy_inherited_hash),
    )
    by_terminal = {record[0]: record for record in terminals}
    mechanism_totals = Counter()
    two_mechanism_totals = Counter()
    for record in terminals:
        for mechanism, count in record[14]:
            mechanism_totals[mechanism] += count
        for mechanism, count in record[27]:
            two_mechanism_totals[mechanism] += count
    nonpositive_rows = tuple(record[0] for record in terminals if record[26])
    require(nonpositive_rows == EXPECTED_NONPOSITIVE_ROWS, nonpositive_rows)
    accelerator_packets = tuple(
        (record[0], record[35]) for record in terminals if record[35] is not None
    )
    require(
        tuple(index for index, _packet in accelerator_packets) == (195,),
        accelerator_packets,
    )
    terminal_summary = (
        len(terminals),
        sum(record[5] for record in terminals),
        sum(record[9] for record in terminals),
        sum(record[11] for record in terminals),
        sum(record[12] for record in terminals),
        sum(record[13] for record in terminals),
        tuple(sorted(mechanism_totals.items())),
        sum(record[15] for record in terminals),
        sum(record[16] for record in terminals),
        sum(record[23] for record in terminals),
        max(record[24] for record in terminals),
        max(record[25] for record in terminals),
        heavy_inherited_hash,
        sum(record[26] for record in terminals),
        tuple(sorted(two_mechanism_totals.items())),
        min((record[28] for record in terminals if record[28] is not None), default=None),
        max(record[31] for record in terminals),
        max(record[32] for record in terminals),
        tuple(
            (record[0], record[29], record[30], record[33], record[34])
            for record in terminals
            if record[26]
        ),
        accelerator_packets,
    )
    if EXPECTED_TERMINAL_SUMMARY is not None:
        require(terminal_summary == EXPECTED_TERMINAL_SUMMARY, terminal_summary)
    if EXPECTED_TERMINAL_SHA256 is not None:
        require(terminal_hash == EXPECTED_TERMINAL_SHA256, terminal_hash)

    family_records = []
    for rank, family in enumerate(families, 1):
        cost, count, divisor_gcd, ruler, family_indices, packet = family
        family_screens = tuple(by_screen[index] for index in family_indices)
        family_terminals = tuple(
            by_terminal[index] for index in family_indices if index in by_terminal
        )
        totals = tuple(sum(record[1][pos] for record in family_screens) for pos in (9, 10, 11, 12))
        mechanisms = Counter()
        two_mechanisms = Counter()
        for record in family_terminals:
            for mechanism, number in record[14]:
                mechanisms[mechanism] += number
            for mechanism, number in record[27]:
                two_mechanisms[mechanism] += number
        family_records.append(
            (
                rank,
                cost,
                count,
                divisor_gcd,
                ruler,
                family_indices,
                tuple(index for index in family_indices if index != 94),
                totals,
                sum(record[1][12] > 0 for record in family_screens),
                sum(record[1][12] == 0 for record in family_screens),
                sum(record[11] for record in family_terminals),
                tuple(sorted(mechanisms.items())),
                min((record[7], record[0]) for record in family_terminals),
                min(
                    (record[17] for record in family_terminals if record[17] is not None),
                    default=None,
                ),
                min(
                    (record[18] for record in family_terminals if record[18] is not None),
                    default=None,
                ),
                digest(family_screens),
                digest(family_terminals),
                sum(record[26] for record in family_terminals),
                tuple(sorted(two_mechanisms.items())),
                min(
                    (record[28] for record in family_terminals if record[28] is not None),
                    default=None,
                ),
            )
        )
    family_records = tuple(family_records)
    require(sum(record[2] for record in family_records) == 110, family_records)
    require(sum(len(record[6]) for record in family_records) == 109, family_records)

    weakest_gap = min((record[7], record[0], record[1]) for record in terminals)
    weakest_torsion = min(record[17] for record in terminals if record[17] is not None)
    weakest_translated = min(record[18] for record in terminals if record[18] is not None)
    hostile = (14, (2, 3), (14 + 6) // 7, (2 > Fraction(3, 2), 2 < Fraction(7, 2), 3 > Fraction(3, 2), 3 < Fraction(7, 2)))
    require(hostile == (14, (2, 3), 2, (True, True, True, True)), hostile)

    semantic_packet = (
        tuple((path.name, expected) for path, expected in EXPECTED_HASHES.items()),
        EXPECTED_ROW195_CONTROL_HASHES,
        census,
        thm3361_family,
        families,
        families_hash,
        indices,
        live_indices,
        screens,
        screen_totals,
        screen_hash,
        terminals,
        terminal_summary,
        terminal_hash,
        accelerator_packets,
        family_records,
        weakest_gap,
        weakest_torsion,
        weakest_translated,
        hostile,
        (372913, 372804),
        (109, 0),
        (12, 0),
        (216, 215),
    )
    semantic_hash = digest(semantic_packet)
    if EXPECTED_SEMANTIC_SHA256 is not None:
        require(semantic_hash == EXPECTED_SEMANTIC_SHA256, semantic_hash)

    lines = [
        "LRC14 projected-k3 z216 post-THM3378 complete-wall hybrid closure",
        "status=FINITE-EXACT;scope=necessary_projected_k3_z1_216_atlas_only",
        "dependency=" + ";".join(f"{path.name}:{expected}" for path, expected in EXPECTED_HASHES.items()),
        "row195_controls=" + ";".join(f"{name}:{value}" for name, value in EXPECTED_ROW195_CONTROL_HASHES),
        f"universe=atlas:{census};historical_post_THM3361_rows:{len(indices)};current_post_THM3378_rows:{len(live_indices)};families:{len(families)};historical_indices_sha256:{digest(indices)};live_indices_sha256:{digest(live_indices)};families_sha256:{families_hash}",
        f"SCREEN_TOTAL;states={screen_totals[0]};crude={screen_totals[1]};status={screen_totals[2]};residual={screen_totals[3]};direct={sum(r[2] for r in screens)};legacy={sum(r[3] for r in screens)};screen_sha256={screen_hash}",
    ]
    for record in family_records:
        rank, cost, count, divisor_gcd, ruler, family_indices, live_family_indices, totals, residual_rows, zero_rows, one_cases, mechanisms, family_gap, family_torsion, family_translated, family_screen_hash, family_terminal_hash, two_cases, two_mechanisms, family_top_h = record
        lines.append(
            f"FAMILY;rank={rank};gcd={divisor_gcd};L={ruler};historical_rows={count};live_rows={len(live_family_indices)};cost={cost};"
            f"indices={','.join(map(str, family_indices))};states={totals[0]};crude={totals[1]};"
            f"status={totals[2]};residual={totals[3]};residual_rows={residual_rows};"
            f"screen_empty_rows={zero_rows};one_high_cases={one_cases};two_high_cases={two_cases};"
            f"mechanisms={mechanisms};two_high_mechanisms={two_mechanisms};"
            f"weakest_two_high_gap={family_gap};weakest_torsion={family_torsion};"
            f"weakest_translated={family_translated};weakest_top_h={family_top_h};"
            f"screen_sha256={family_screen_hash};terminal_sha256={family_terminal_hash}"
        )
    for record in terminals:
        if record[2] < 5_000_000:
            continue
        lines.append(
            f"HEAVY_TERMINAL;index={record[0]};E={','.join(map(str, record[1]))};"
            f"L={record[2]};high={record[3]};first_d={record[4]};residual={record[5]};"
            f"two_gap={record[7]};zero_high_hostiles={record[9]};one_high_cases={record[11]};"
            f"one_high_passports={record[12]};low_pairs={record[13]};mechanisms={record[14]};"
            f"high_keys={record[15]};translated_low_pair_keys={record[16]};"
            f"weakest_torsion={record[17]};weakest_translated={record[18]};"
            f"cases_sha256={record[19]};inherited_terminal_sha256={digest(record[22])};"
            f"nonheavy_legacy_parity_checks={record[23]};"
            f"peak_compact_cell_bytes={record[24]};peak_support_bytes={record[25]}"
        )
    for record in terminals:
        if not record[26]:
            continue
        lines.append(
            f"TWO_HIGH_TERMINAL;index={record[0]};E={','.join(map(str, record[1]))};"
            f"L={record[2]};two_gap={record[7]};cases={record[26]};"
            f"mechanisms={record[27]};weakest_top_h={record[28]};"
            f"cases_sha256={record[29]};closures_sha256={record[30]};"
            f"peak_cell_bytes={record[31]};peak_aux_bytes={record[32]};"
            f"cell_packets={record[33]};period_obstruction={record[34]}"
        )
    for record in terminals:
        packet = record[35]
        if packet is None:
            continue
        lines.append(
            f"ONE_HIGH_ACCELERATOR;index={record[0]};"
            f"first_source={packet[0]};needed={packet[1]};top_h_sha256={packet[2]};"
            f"triple_partition={packet[3]};failed_state_partition={packet[4]};"
            f"residual_state_partition={packet[5]};strict_states_sha256={packet[6]};"
            f"strict_certificates_sha256={packet[7]};hostile_states_sha256={packet[8]};"
            f"weakest_strict={packet[9]};closest_hostile={packet[10]};"
            f"reduced_one_high={packet[11]};torsion_summary={packet[12]};"
            f"torsion_closures_sha256={packet[13]};hostile_two_high={packet[14]};"
            f"state_map_controls={packet[15]};canonical_controls={packet[16]};"
            f"historical_only={packet[17]};partition_sha256={packet[18]};"
            f"accelerator_sha256={digest(packet)}"
        )
    lines.extend(
        [
            f"TERMINAL_TOTAL;rows={terminal_summary[0]};passports={terminal_summary[1]};"
            f"zero_high_hostiles={terminal_summary[2]};one_high_cases={terminal_summary[3]};"
            f"one_high_passports={terminal_summary[4]};low_pairs={terminal_summary[5]};"
            f"mechanisms={terminal_summary[6]};torsion_keys={terminal_summary[7]};"
            f"translated_low_pair_keys={terminal_summary[8]};"
            f"nonheavy_legacy_parity_checks={terminal_summary[9]};"
            f"peak_compact_cell_bytes={terminal_summary[10]};"
            f"peak_support_bytes={terminal_summary[11]};"
            f"heavy_inherited_terminal_sha256={terminal_summary[12]};"
            f"two_high_cases={terminal_summary[13]};"
            f"two_high_mechanisms={terminal_summary[14]};"
            f"weakest_top_h={terminal_summary[15]};"
            f"peak_two_high_cell_bytes={terminal_summary[16]};"
            f"peak_two_high_aux_bytes={terminal_summary[17]};"
            f"nonpositive_packets={terminal_summary[18]};"
            f"accelerator_packets={len(terminal_summary[19])};"
            f"accelerator_sha256={digest(terminal_summary[19])};"
            f"terminal_sha256={terminal_hash}",
            f"weakest_two_high_gap={weakest_gap}",
            f"weakest_torsion_support_surplus={weakest_torsion}",
            f"weakest_translated_cardinality_margin={weakest_translated}",
            "hybrid_lemma=the strict high gate gives at least one high; outside row195 every positive-gap residual is exactly-one-high and closes by located torsion or translated support; on row195 the three-block top-h atom majorant uniformly closes 1408 residual states, the remaining 599 one-high cases on 156 hostile states close by located torsion, and the unique hostile state without a one-high case routes through the full 137-case two-high audit to its denominator-two measure terminal",
            f"translated_equality_hostile={hostile}",
            "independent_control=every non-heavy non-row195 compact one-high common-cell array and torsion witness exactly equals the inherited tuple implementation; the three heavy inherited terminal records reproduce their frozen packet hash; a self-contained independent row195 two-high audit reconstructs all 137 two-high cases and exact common-source closures with source/output/semantic hashes fccc10392624bbdfc2083993ad51a423e8974c135b9bc635351304a71cb0de74/55ad1da385d35f5b38fdc7de2d9916f54ee874b806496d4a7b1270cb526ad30c/7a70daaa1f25fc6fb9bfd7469ce3e5a18ba650296119789c81bc70fd32016b9c; a separate self-contained bounded row195 accelerator freezes the 1408-state/599-case partition with source/output/semantic hashes 00e1855dd428466c113d3d559e629baefcfdd9b0ea8ac732069fee3bb4ac74fe/b188805a09b484466f4c36395b2b426beba8ad470704d9c28414f92ec7a43f13/5ad577c0a388c0021cda9506a8147a25732f36f801e7bcf746b856a5d9ff1070; every selected witness cell is rechecked by a second literal cell-clean formula; ordinary and optimized runs byte-match",
            "composition=THM3378 already removes row94;the historical 110-row scan separately and redundantly re-closes row94 and closes every one of the 109 current live rows",
            "consequence=all_109_live_rows_excluded;ledger:372913->372804;wall:109->0;families:12->0;projected_k3_cap:216->215",
            "first_surviving_row=none_at_z1_216",
            "scope=no_physical_entry_or_arbitrary_k_or_rung_or_LRC14_conclusion",
            f"semantic_sha256={semantic_hash}",
            "all_exact_controls=PASS",
        ]
    )
    print("\n".join(lines))


if __name__ == "__main__":
    main()
