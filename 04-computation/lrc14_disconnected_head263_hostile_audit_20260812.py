#!/usr/bin/env python3
"""Exact hostile comparison of the head-263 C kernel and THM-3352 reference."""

from argparse import ArgumentParser
from fractions import Fraction as F
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
from subprocess import run
import ast

ROOT = Path(__file__).resolve().parents[1]
ENGINE = ROOT / "04-computation/lrc_general_reflected_pair_mass_reference_audit_thm3352.py"
CONTEXTS = ROOT / "04-computation/lrc14_disconnected_head263_contexts_20260812.txt"
CSOURCE = ROOT / "04-computation/lrc14_disconnected_head263_exact_scan_20260812.c"
EXPECTED_ENGINE = "da941a4267147d5442be81ae81880742d2f6b901bfc1d20fb667822402a2950e"
EXPECTED_CONTEXT = "efea6bd97522fe1c31a5a88ca9f3223f9e7a8c08e3be85c493e9f62fdfaf06e4"
EXPECTED_CSOURCE = "498ef114fbf7b3d54e62de4556cd0f669b0300f882d61d16c39f1566f5efb23f"
EXPECTED_QUERY = "a767e1ae984c120ae615dd9cd57deba0bc6a2f15bb320c7e8c0ca128ee51e494"
EXPECTED_C_OUTPUT = "e6da2cc4eb9446ce05c9e438fed9e291ff12a6f88d9763d05947d0df58a531c2"
EXPECTED_SEMANTIC = "43a551e5efd27e31940fa452f5d044b714cd6d82304834f706c14b0536a538c5"


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
    parser.add_argument("binary", type=Path)
    parser.add_argument("--work", type=Path, default=Path("/tmp/lrc14-head263-hostiles"))
    args = parser.parse_args()
    binary = args.binary.resolve()
    work = args.work.resolve()
    work.mkdir(parents=True, exist_ok=True)
    query = work / "hostiles.query"
    c_output = work / "hostiles.c.out"

    require(file_hash(ENGINE) == EXPECTED_ENGINE, ("engine", file_hash(ENGINE)))
    require(file_hash(CONTEXTS) == EXPECTED_CONTEXT, ("contexts", file_hash(CONTEXTS)))
    require(file_hash(CSOURCE) == EXPECTED_CSOURCE, ("C source", file_hash(CSOURCE)))
    tree = ast.parse(Path(__file__).read_text(), filename=__file__)
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)),
            "optimization-sensitive assert")
    mass_engine = load("head263_hostile_reference", ENGINE)
    contexts = tuple(tuple(map(int, line.split())) for line in CONTEXTS.read_text().splitlines())
    pairs = ((1, 1), (3, 5), (6, 10), (9, 15), (1, 7),
             (263, 2103), (2103, 263), (261, 2087))
    rows = [(L, cell, e, p, f, q)
            for L, cell, e, f in contexts for p, q in pairs]
    equal_hostile = (168, 33, 3, 1, 8, 1)
    if equal_hostile not in rows:
        rows.append(equal_hostile)
    query.write_text("".join(" ".join(map(str, row)) + "\n" for row in rows),
                     encoding="ascii", newline="\n")
    require(file_hash(query) == EXPECTED_QUERY, ("query hash", file_hash(query)))
    process = run([str(binary), "--query", str(query), str(c_output)],
                  capture_output=True, text=True)
    require(process.returncode == 0, ("C query", process.returncode, process.stderr))
    require(process.stderr.strip() == f"query_rows {len(rows)}", process.stderr)
    require(file_hash(c_output) == EXPECTED_C_OUTPUT, ("C output hash", file_hash(c_output)))
    got = tuple(F(*map(int, line.split())) for line in c_output.read_text().splitlines())
    expected = tuple(mass_engine.mass(L, cell, e, p, f, q)
                     for L, cell, e, p, f, q in rows)
    require(got == expected, "C/reference mismatch")
    index = rows.index(equal_hostile)
    require(got[index] == F(8, 55), ("equal-level regression", got[index]))
    semantic = sha256("".join(f"{row} {value}\n" for row, value in zip(rows, got)).encode()).hexdigest()
    require(semantic == EXPECTED_SEMANTIC, ("semantic", semantic))

    print("DISCONNECTED HEAD263 C HOSTILE AUDIT")
    print("query_rows", len(rows), "reference_equalities", len(rows), "failures", 0)
    print("pairs", pairs, "equal_level_regression", got[index])
    print("engine_sha256", file_hash(ENGINE), "c_source_sha256", file_hash(CSOURCE))
    print("query_sha256", file_hash(query), "c_output_sha256", file_hash(c_output),
          "semantic_sha256", semantic)
    print("python_assert_nodes", 0)


if __name__ == "__main__":
    main()
