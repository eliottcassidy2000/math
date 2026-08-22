#!/usr/bin/env python3
"""Lightweight hash/semantic gate for the dual exact boundary through d=10000."""

import ast
import hashlib
import json
from pathlib import Path


ARTIFACTS = (
    ("04-computation/factorial_adaptive_rho_block_6000.py",
     "b65edcf2870714ca57456b8297afdd05284a09b302ec4b84d2e57829520c94d1"),
    ("04-computation/factorial_adaptive_rho_block_6000_independent_audit.py",
     "d416cb2955fd745394cf1043ac8c2eba28a6a97beb264dd9cbe9919ed8c96724"),
    ("04-computation/factorial_adaptive_rho_block_10000.py",
     "105e62698d0c3cf0066a100e9d205a5c1f1c31e64cfdba0dc2fe23decd8f0eba"),
    ("05-knowledge/results/factorial_adaptive_rho_block_6000.out",
     "ab629edc04e31d1889741688897bfe60f5249df60b64116391217393962b1ddf"),
    ("05-knowledge/results/factorial_adaptive_rho_block_6000_independent_audit.out",
     "8d6adbcaa14c85d022f726db97a80a7bafc1288505ccde720b3f5fbf6ee2a922"),
    ("05-knowledge/results/factorial_adaptive_rho_block_10000.out",
     "18b131aed2f380b1c1bace8beeb8488ced0e24599f4b7484d66a14e5869c0d22"),
)
BLOCK_SEMANTICS = (
    (4001, 6000, "7f8ab74ae9fae027f32fd7eabaf0338c217319e274594bd603859a1bbcca28bd"),
    (6001, 10000, "d90179fdebd48dd82cd368b957c9602fbd287774287de0fecb73b4a84dca69f3"),
)
EXPECTED_COMBINED_SHA256 = "4a364bbafdfef0dc6d905063c9570b987e00eca2506506b8160ce24e453878bf"


def require(condition, payload):
    if not condition:
        raise RuntimeError(payload)


def digest(value):
    data = json.dumps(value, separators=(",", ":"), sort_keys=True).encode("ascii")
    return hashlib.sha256(data).hexdigest()


def main():
    root = Path(__file__).resolve().parents[1]
    verified = []
    for relative, expected in ARTIFACTS:
        path = root / relative
        actual = hashlib.sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()
        require(actual == expected, (relative, actual, expected))
        verified.append((relative, actual))

    output_6000 = (root / ARTIFACTS[3][0]).read_text(encoding="utf-8")
    independent_6000 = (root / ARTIFACTS[4][0]).read_text(encoding="utf-8")
    output_10000 = (root / ARTIFACTS[5][0]).read_text(encoding="utf-8")
    require("semantic_sha256=" + BLOCK_SEMANTICS[0][2] in output_6000,
            "primary d6000 semantic")
    require("semantic_sha256=" + BLOCK_SEMANTICS[0][2] in independent_6000,
            "independent d6000 semantic")
    require("semantic_sha256=" + BLOCK_SEMANTICS[1][2] in output_10000,
            "d10000 semantic")
    require("survivors=()" in independent_6000 and "survivors=()" in output_10000,
            "survivor gate")
    require("d6001=(6001, 'rho', 11" in output_10000, "d6001 boundary record")
    require("(6518, 29)" in output_10000, "d6518 hostile killer")
    require("(6518, (3087, 3430, 4802, 5145, 5488, 5831))" in output_10000,
            "d6518 hostile packet")
    combined = digest(BLOCK_SEMANTICS)
    require(combined == EXPECTED_COMBINED_SHA256, combined)
    require(not any(isinstance(node, ast.Assert)
                    for node in ast.walk(ast.parse(Path(__file__).read_text()))),
            "assert node")

    print("FACTORIAL ADAPTIVE RHO BOUNDARY THROUGH d=10000 HASH GATE")
    print("artifacts=%s" % (tuple(verified),))
    print("block_semantic_sha256=%s" % (BLOCK_SEMANTICS,))
    print("combined_basis=json.dumps(block_semantics,separators=(',',':'),sort_keys=True).encode('ascii')")
    print("combined_semantic_sha256=%s" % EXPECTED_COMBINED_SHA256)
    print("d6001=(divisor_packet_size=59,rho_p7_post=(5152,),rho_killer=11,p2_inadmissible_and_skipped)")
    print("d6518_hostile=(divisor_packet=(3087,3430,4802,5145,5488,5831);p11_post=(5831,);p13,p17,p19_retain;p23_inadmissible;p29_gap=(5669,5887)_kills_5831)")
    print("consequence=FINITE-EXACT closure through d=10000, equivalently r<=9998; no arbitrary-support FC/SFC consequence")
    print("STATUS=PASS")


if __name__ == "__main__":
    main()
