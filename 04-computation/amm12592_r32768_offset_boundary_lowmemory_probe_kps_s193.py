#!/usr/bin/env python3
"""Exact Rule-A state replay with the diagnostic junk-L1 observer elided."""

from __future__ import annotations

import ast
from functools import lru_cache
from hashlib import sha256
import json
from math import comb
import os
from pathlib import Path
import sys
import time


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "04-computation"))
import amm12592_R16384_offset_transition_thm3644 as T

ENGINE = ROOT / "04-computation/amm12592_transient_fast_junkflow_boxeph.py"
EXPECTED_ENGINE_SHA256 = (
    "8887080fc6e30760efa4a0ba76218ec97676cc717c6e76ccefbaeec6c73684ad"
)


def require(condition, payload):
    if condition is not True:
        raise RuntimeError(payload)


class ElideJunkL1(ast.NodeTransformer):
    def __init__(self):
        self.count = 0

    def visit_AugAssign(self, node):
        self.generic_visit(node)
        if (isinstance(node.target, ast.Name) and node.target.id == "junkL1"
                and isinstance(node.op, ast.Add)):
            self.count += 1
            return ast.copy_location(ast.Pass(), node)
        return node


def load_lowmemory_engine():
    require(sha256(ENGINE.read_bytes()).hexdigest() == EXPECTED_ENGINE_SHA256,
            "engine hash drift")
    tree = ast.parse(ENGINE.read_text(encoding="utf-8"), str(ENGINE))
    selected = {"two_G_coeffs", "initial_junk", "run_fast"}
    body = [node for node in tree.body
            if isinstance(node, ast.FunctionDef) and node.name in selected]
    require({node.name for node in body} == selected, "engine AST selection")
    module = ast.Module(body=body, type_ignores=[])
    transformer = ElideJunkL1()
    module = transformer.visit(module)
    ast.fix_missing_locations(module)
    require(transformer.count == 2, ("junk observer edits", transformer.count))
    namespace = {
        "__builtins__": __builtins__, "comb": comb, "json": json,
        "os": os, "sys": sys, "time": time,
        "floor_gamma_star": T.floor_gamma_star,
    }
    exec(compile(module, str(ENGINE), "exec"), namespace)
    namespace["two_G_coeffs"] = lru_cache(maxsize=None)(namespace["two_G_coeffs"])
    return namespace["run_fast"]


def main():
    offset = int(sys.argv[1]) if len(sys.argv) > 1 else 854
    low = load_lowmemory_engine()
    original = T.load_engine(ROOT)
    controls = ((64, 0), (64, 1), (512, 4), (512, 5))
    for r, d0 in controls:
        require(T.stable_result(low, r, d0) == T.stable_result(original, r, d0),
                ("observer-elision control", r, d0))
    print(f"engine={EXPECTED_ENGINE_SHA256};observer_elision=junkL1_only;controls={controls}",
          flush=True)
    started = time.perf_counter()
    row = T.stable_result(low, 32768, offset)
    print(f"R=32768;offset={offset};stable_result={row};seconds={time.perf_counter()-started:.3f}",
          flush=True)


if __name__ == "__main__":
    main()

