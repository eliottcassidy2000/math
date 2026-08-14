#!/usr/bin/env python3
"""Replay every reported 679--698 bridge argmin with the reference engine."""
from concurrent.futures import ProcessPoolExecutor
from fractions import Fraction as F
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from math import gcd
from pathlib import Path
import ast


ENGINE = Path("/tmp/canonical_reference_engine_thm3352.py")
LEDGER = Path("/tmp/disconnected_bridge679_698_exact_scan.ledger")
EXPECTED_ENGINE = "da941a4267147d5442be81ae81880742d2f6b901bfc1d20fb667822402a2950e"
EXPECTED_LEDGER = "0915c0c3a6eaf9fa4577e4449900865992d8fa646184e4ba6600801e89a55421"
EXPECTED_TASKS = "9b2087d4a05e1ed453b05594d3c4d6ec748ff9c4f5cab95ea615d8238e002fa3"
TARGET = F(186_636_088_362, 58_865_718_786_875)
MASS = None


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def filehash(path):
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def init_worker():
    global MASS
    spec = spec_from_file_location("bridge_679_698_reference", ENGINE)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    MASS = module


def replay(row):
    p, q, n, d, L, j, e, f = row
    value = MASS.mass(L, j, e, p, f, q)
    if value != F(n, d):
        return row, value
    return None


def main():
    require(filehash(ENGINE) == EXPECTED_ENGINE, ("engine hash", filehash(ENGINE)))
    require(filehash(LEDGER) == EXPECTED_LEDGER, ("ledger hash", filehash(LEDGER)))
    tree = ast.parse(Path(__file__).read_text(), filename=__file__)
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)), "assert node")
    rows = tuple(tuple(map(int, line.split())) for line in LEDGER.read_text().splitlines())
    expected = tuple(
        (p, q)
        for p in range(679, 699)
        for q in range(p + 1, 8 * p)
        if gcd(p, q) <= 3
    )
    payload = "".join(f"{p} {q}\n" for p, q in expected).encode()
    require(tuple(row[:2] for row in rows) == expected, ("universe", len(rows), len(expected)))
    require(sha256(payload).hexdigest() == EXPECTED_TASKS, "task digest")
    require(all(F(row[2], row[3]) > TARGET for row in rows), "threshold")
    failures = []
    semantic = sha256()
    with ProcessPoolExecutor(max_workers=8, initializer=init_worker) as pool:
        for row, result in zip(rows, pool.map(replay, rows, chunksize=128)):
            if result is not None:
                failures.append(result)
            semantic.update((repr((row, result)) + "\n").encode())
    require(not failures, failures[:10])
    print("DISCONNECTED RAW BRIDGE 679--698 ALL-ARGMIN REFERENCE AUDIT")
    print("rows", len(rows), "independent_argmin_mass_checks", len(rows), "failures", len(failures))
    print("engine_sha256", filehash(ENGINE))
    print("ledger_sha256", filehash(LEDGER))
    print("task_semantic_sha256", sha256(payload).hexdigest())
    print("audit_semantic_sha256", semantic.hexdigest())
    print("python_assert_nodes", 0)
    print("status=PASS")


if __name__ == "__main__":
    main()
