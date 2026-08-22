#!/usr/bin/env python3
"""Independent exact audit of the disconnected head-263 C ledger.

Run the C scanner first and pass its generated ledger.  This audit checks all
201,377 reported minimizers with the slower literal-audited THM-3352 reference
engine, then recomputes all 2,530 contexts on 220 hostile/spread channels.
"""

from argparse import ArgumentParser
from fractions import Fraction as F
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from math import gcd
from pathlib import Path
from random import Random
import ast

ROOT = Path(__file__).resolve().parents[1]
ENGINE = ROOT / "04-computation/lrc_general_reflected_pair_mass_reference_audit_thm3352.py"
CONTEXTS = ROOT / "04-computation/lrc14_disconnected_head263_contexts_20260812.txt"
CSOURCE = ROOT / "04-computation/lrc14_disconnected_head263_exact_scan_20260812.c"
GENERATOR = ROOT / "04-computation/lrc14_disconnected_head263_contexts_20260812.py"
EXPECTED_ENGINE = "da941a4267147d5442be81ae81880742d2f6b901bfc1d20fb667822402a2950e"
EXPECTED_CONTEXT = "efea6bd97522fe1c31a5a88ca9f3223f9e7a8c08e3be85c493e9f62fdfaf06e4"
EXPECTED_LEDGER = "1caa410d5e57d6085e730270091da3d8433f14d0bd74d98f6283f2ec8a4ca7a0"
EXPECTED_CSOURCE = "498ef114fbf7b3d54e62de4556cd0f669b0300f882d61d16c39f1566f5efb23f"
EXPECTED_GENERATOR = "468a88781de94cb6d0d49d500371bb234b1ddb2615e8e8a8c6203c827d1ca298"
EXPECTED_TASKS = "2771f2f901f2f052952343fd77412114ae1d1d99543f42539bc28d0f0f1948af"
EXPECTED_CONTROLS = "d24a7ea39659a271fb87f145aeaa16a98c85851dbfdb2762226f3f59655e5841"
TARGET = F(186636088362, 58865718786875)


def require(condition, description):
    if not condition:
        raise RuntimeError(description)


def file_hash(path):
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def load(name, path):
    spec = spec_from_file_location(name, path)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def main():
    parser = ArgumentParser()
    parser.add_argument("ledger", type=Path)
    args = parser.parse_args()
    ledger = args.ledger.resolve()

    require(file_hash(ENGINE) == EXPECTED_ENGINE, ("engine hash", file_hash(ENGINE)))
    require(file_hash(CONTEXTS) == EXPECTED_CONTEXT, ("context hash", file_hash(CONTEXTS)))
    require(file_hash(ledger) == EXPECTED_LEDGER, ("ledger hash", file_hash(ledger)))
    require(file_hash(CSOURCE) == EXPECTED_CSOURCE, ("C source hash", file_hash(CSOURCE)))
    require(file_hash(GENERATOR) == EXPECTED_GENERATOR, ("generator hash", file_hash(GENERATOR)))
    for path in (GENERATOR, Path(__file__)):
        tree = ast.parse(path.read_text(), filename=str(path))
        require(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)),
                ("optimization-sensitive assert", path))

    mass_engine = load("head263_reference_mass", ENGINE)
    contexts = tuple(tuple(map(int, line.split())) for line in CONTEXTS.read_text().splitlines())
    require(len(contexts) == 2530 and contexts == tuple(sorted(contexts)),
            ("contexts", len(contexts)))
    rows = []
    for line in ledger.read_text().splitlines():
        g, P, Q, numerator, denominator, L, cell, e, f = map(int, line.split())
        rows.append((g, P, Q, F(numerator, denominator), L, cell, e, f))
    require(len(rows) == 201377, ("rows", len(rows)))

    expected = tuple(
        (g, P, Q)
        for g in (1, 2, 3)
        for P in range(1, 264)
        if g * P < 264
        for Q in range(P + 1, 8 * P)
        if gcd(P, Q) == 1 and P + Q >= 8
    )
    task_payload = "".join(f"{g} {P} {Q}\n" for g, P, Q in expected).encode()
    require(sha256(task_payload).hexdigest() == EXPECTED_TASKS, "task semantic digest")
    require(tuple(row[:3] for row in rows) == expected, "task universe/order")
    require(all(row[3] > TARGET for row in rows), "ledger threshold")

    for g, P, Q, value, L, cell, e, f in rows:
        require(mass_engine.mass(L, cell, e, g * P, f, g * Q) == value,
                ("argmin mismatch", g, P, Q, value, L, cell, e, f))

    ranked = sorted(rows, key=lambda row: row[3])
    chosen = {row[:3] for row in ranked[:20]}
    rng = Random(3352)
    while len(chosen) < 220:
        chosen.add(rows[rng.randrange(len(rows))][:3])
    controls = []
    for g, P, Q in sorted(chosen):
        best = min(
            (mass_engine.mass(L, cell, e, g * P, f, g * Q), L, cell, e, f)
            for L, cell, e, f in contexts
        )
        claimed = next(row for row in rows if row[:3] == (g, P, Q))
        require(best == (claimed[3], *claimed[4:]),
                ("global min mismatch", g, P, Q, best, claimed))
        controls.append(((g, P, Q), best))
    control_semantic = sha256(repr(tuple(controls)).encode()).hexdigest()
    require(control_semantic == EXPECTED_CONTROLS, ("control semantic", control_semantic))
    by_g = tuple((g, min((row for row in rows if row[0] == g), key=lambda row: row[3]))
                 for g in (1, 2, 3))

    print("DISCONNECTED HEAD263 INDEPENDENT PYTHON AUDIT")
    print("rows", len(rows), "contexts", len(contexts), "argmin_exact_checks", len(rows),
          "full_context_controls", len(controls),
          "control_mass_checks", len(controls) * len(contexts))
    print("failures", 0, "by_g", by_g)
    print("task_semantic_sha256", sha256(task_payload).hexdigest())
    print("engine_sha256", file_hash(ENGINE), "context_sha256", file_hash(CONTEXTS),
          "ledger_sha256", file_hash(ledger))
    print("c_source_sha256", file_hash(CSOURCE), "generator_sha256", file_hash(GENERATOR),
          "python_assert_nodes", 0)
    print("control_semantic_sha256", control_semantic)


if __name__ == "__main__":
    main()
