#!/usr/bin/env python3
"""Exact certificate and full substitution replay for THM-4144.

Default mode verifies the frozen seventeen-shard certificate and independently
replays every shard extremum with the audited THM-4131 Python evaluator.  Pass
``--full`` to regenerate all 4,055,870 substitution rows with nauty and the
pinned standalone C++ kernel.  ``--emit-regime`` is the streaming generator
used by the full replay.
"""

from __future__ import annotations

from argparse import ArgumentParser
from concurrent.futures import ThreadPoolExecutor
from fractions import Fraction
from hashlib import sha256
import ast
import importlib.util
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
from tempfile import TemporaryDirectory


sys.dont_write_bytecode = True
sys.stdout.reconfigure(newline="\n")

ROOT = Path(__file__).resolve().parents[1]
BASE_PATH = ROOT / "04-computation/tournament_strong_centrality_through_order_eight_thm4131.py"
FAMILY_PATH = ROOT / "04-computation/tournament_strong_cyclic_substitution_centrality_counterexample_thm4133.py"
ENGINE_PATH = ROOT / "04-computation/tournament_order_eleven_large_module_centrality_thm4144_engine.cpp"

EXPECTED_BASE_SHA256 = "6b195b0379d1ae3e5d215aa1c495f7180daeecae189df86269d07ef855867881"
EXPECTED_FAMILY_SHA256 = "7bd4c518464d4c48baf9cb9c1c8c2012a9f79f029c8e07141c0e51c338ffeeb2"
EXPECTED_ENGINE_SHA256 = "03dfb092830313d71d41c35dcf7fb7d300db9850bbd52b1291ace7e674bf63bc"
EXPECTED_GENTOURNG_SHA256 = "89df605922cc574b28688248b7c256d24342cc615f887e89b2d096038970c110"
EXPECTED_SEMANTIC = "4c7ca82098f748137cc99e52a10b47c5eb59b405746c9daf20b664701bed10c9"

ALL_CLASS_COUNTS = {3: 2, 4: 4, 5: 12, 6: 56, 7: 456, 8: 6880, 9: 191536}
STRONG_CLASS_COUNTS = {3: 1, 4: 1, 5: 6, 6: 35, 7: 353, 8: 6008, 9: 178133}
REGIMES = ((3, 9, 4), (4, 8, 1), (5, 7, 1), (6, 6, 1),
           (7, 5, 1), (8, 4, 1), (9, 3, 8))

# Each certificate row is
# (q,r,shard,shards,(rows,strong,rational_fail,coset_fail),
#  (max_label,(H,W,D4,C_hd),rho,margin),
#  (min_label,(H,W,D4,C_hd),rho,margin)).
CERTIFICATE_ROWS = (
    (3, 9, 0, 4, (143652, 143652, 0, 0),
     ("1000000000111111111111111111111111110100110101101110111", (8253, 117117, 4167524088, 1129039884), (1493439, 2756299), 1812),
     ("1111111110111111111111111111111111111111111111111101111", (771, 24819, 183002157, 24649812), (5477736, 20333573), 598)),
    (3, 9, 1, 4, (143652, 143652, 0, 0),
     ("1000000000111111111110100111101011110111110111111111111", (8253, 117117, 4167524088, -1129039884), (1493439, 2756299), 1812),
     ("1000000000111111111111111111111111111111111111111111101", (771, 24819, 183002157, 24649812), (5477736, 20333573), 598)),
    (3, 9, 2, 4, (143652, 143652, 0, 0),
     ("1111111110111111110110100101101010110110110101110110101", (8253, 117117, 4167524088, 1129039884), (1493439, 2756299), 1812),
     ("1011111110111111110111111101111110111110111101110110101", (771, 24819, 183002157, -24649812), (5477736, 20333573), 598)),
    (3, 9, 3, 4, (143652, 143652, 0, 0),
     ("1111111110111111111111111111101001110101110111101111111", (8253, 117117, 4167524088, 1129039884), (1493439, 2756299), 1812),
     ("1111111110111111110111111101111110111110111101010110101", (771, 24819, 183002157, 24649812), (5477736, 20333573), 598)),
    (4, 8, 0, 1, (27520, 27520, 0, 0),
     ("1000000000111111111111111111111111110100110101101110111", (8253, 117117, 4167524088, 1129039884), (1493439, 2756299), 1812),
     ("1000000000111111111111111111111111111111111111111111101", (771, 24819, 183002157, 24649812), (5477736, 20333573), 598)),
    (5, 7, 0, 1, (13680, 13680, 0, 0),
     ("1111111101110100111110101111101111110111111111111111101", (6615, 102123, 3185303373, -907491438), (3201028, 5617819), 1494),
     ("1101111111111111111100000001111111111111111111111111101", (585, 20421, 125387919, 19657296), (1456096, 4643997), 468)),
    (6, 6, 0, 1, (11760, 11760, 0, 0),
     ("1111111101110011111111011111101111111111111111111111101", (2727, 55323, 934218519, -263172504), (175448336, 311406173), 818),
     ("1110111111111111111111111111000000111111111111111111101", (495, 17703, 96064920, 13692870), (152143, 533694), 430)),
    (7, 5, 0, 1, (29652, 29652, 0, 0),
     ("1111111101110010111110101111110111110111101111111111101", (4905, 78525, 1883312325, -519709950), (4619644, 8370277), 1194),
     ("1111011111111111111111111111111111100000111111111111101", (459, 16299, 82439328, 10556748), (879729, 3434972), 418)),
    (8, 4, 0, 1, (192256, 192256, 0, 0),
     ("1111111101100101111111011111101111101111111111111111101", (2223, 46371, 655081417, -176794806), (353589612, 655081417), 722),
     ("1011101111111111111111111111111111111111100001111111101", (405, 15633, 76033377, 10601010), (112180, 402293), 388)),
    (9, 3, 0, 8, (400800, 400800, 0, 0),
     ("1111111101101100111111001111100111111111111111111111101", (5985, 95013, 2757707790, -774391140), (51626076, 91923593), 1412),
     ("1011011111111111111111111111111111111111101111111000111", (255, 12147, 44439174, 924918), (308306, 7406529), 390)),
    (9, 3, 1, 8, (400800, 400800, 0, 0),
     ("1011111110111111111111111111111001111001110011111111111", (5319, 86835, 2300518719, 627952454), (1255904908, 2300518719), 1336),
     ("1011011111111111111111111111111111111111101111111000101", (405, 15633, 76033377, 10601010), (112180, 402293), 388)),
    (9, 3, 2, 8, (400799, 400799, 0, 0),
     ("1011111110111111111111111111011001111001110011111111111", (5985, 95013, 2757707790, 774391140), (51626076, 91923593), 1412),
     ("1101111111111111111111111111110111111111111111000111111", (225, 11625, 40926020, 551482), (275741, 10231505), 382)),
    (9, 3, 3, 8, (400799, 400799, 0, 0),
     ("1111111101111100111101101111110111110111111111111111101", (5301, 86769, 2298444159, -629903862), (419935908, 766148053), 1330),
     ("1110111111111111111111111111111111111101111111111111101", (243, 11871, 42303227, -878818), (1757636, 42303227), 380)),
    (9, 3, 4, 8, (400799, 400799, 0, 0),
     ("1111111101111100111111101111011111111111111111111111101", (5301, 86769, 2298444159, -629903862), (419935908, 766148053), 1330),
     ("1011110111111111111111111111111111111111111111000111111", (243, 11871, 42303227, 878818), (1757636, 42303227), 380)),
    (9, 3, 5, 8, (400799, 400799, 0, 0),
     ("1111111101110000111111101111111111101111111111111111101", (5301, 86769, 2298444159, -629903862), (419935908, 766148053), 1330),
     ("1110111111101111111111111111111111111101111111111111101", (405, 15489, 74011266, -9917910), (1090, 4067), 388)),
    (9, 3, 6, 8, (400799, 400799, 0, 0),
     ("1111111101110000111110001111111111111111101111111111101", (5301, 86769, 2298444159, -629903862), (419935908, 766148053), 1330),
     ("1000111101111111111111111111111111111111111111111111101", (243, 11871, 42303227, -878818), (1757636, 42303227), 380)),
    (9, 3, 7, 8, (400799, 400799, 0, 0),
     ("1111111101110001111111101111111111101111111111111111101", (4833, 80253, 1966784871, -538167240), (358778160, 655594957), 1234),
     ("1011101111111111111111111111111111111111111111101111110", (243, 11871, 42303227, 878818), (1757636, 42303227), 380)),
)

SHARP_QUOTIENT_LABEL = "1101111101"
SHARP_BLOCK_LABEL = "110100110101101110111"
SHARP_DISTINGUISHED_VERTEX = 1
DELETION_LABEL = "1111111101111100011111100111111011111111111111111111111"
BEST_REPAIR_LABEL = "1111111101111100011111100111111011111111111111111111110"

CHECKS = 0


def require(condition: bool, label: object) -> None:
    global CHECKS
    CHECKS += 1
    if not condition:
        raise RuntimeError(f"requirement failed: {label}")


def file_digest(path: Path) -> str:
    hasher = sha256()
    with path.open("rb") as source:
        for block in iter(lambda: source.read(1 << 20), b""):
            hasher.update(block)
    return hasher.hexdigest()


def digest(value: object) -> str:
    raw = json.dumps(value, sort_keys=True, separators=(",", ":")).encode("ascii")
    return sha256(raw).hexdigest()


def load_module(name: str, path: Path, expected_hash: str):
    require(file_digest(path) == expected_hash, ("pinned dependency", path.name))
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def adjacency_from_label(label: str, order: int) -> tuple[int, ...]:
    require(len(label) == order * (order - 1) // 2, ("label length", order))
    adjacency = [0] * order
    cursor = 0
    for left in range(order):
        for right in range(left + 1, order):
            require(label[cursor] in "01", ("binary label", cursor))
            if label[cursor] == "1":
                adjacency[left] |= 1 << right
            else:
                adjacency[right] |= 1 << left
            cursor += 1
    return tuple(adjacency)


def label_from_adjacency(adjacency: tuple[int, ...]) -> str:
    return "".join(
        "1" if adjacency[left] & (1 << right) else "0"
        for left in range(len(adjacency))
        for right in range(left + 1, len(adjacency))
    )


def substitute(quotient: tuple[int, ...], vertex: int,
               block: tuple[int, ...]) -> tuple[int, ...]:
    q = len(quotient)
    r = len(block)
    require(0 <= vertex < q and q + r - 1 == 11, "substitution sizes")
    sizes = [1] * q
    sizes[vertex] = r
    starts = []
    cursor = 0
    for size in sizes:
        starts.append(cursor)
        cursor += size
    adjacency = [0] * 11
    for source, row in enumerate(block):
        for target in range(r):
            if row & (1 << target):
                adjacency[starts[vertex] + source] |= 1 << (starts[vertex] + target)
    for source_block in range(q):
        for target_block in range(q):
            if not (quotient[source_block] & (1 << target_block)):
                continue
            for source in range(starts[source_block], starts[source_block] + sizes[source_block]):
                for target in range(starts[target_block], starts[target_block] + sizes[target_block]):
                    adjacency[source] |= 1 << target
    return tuple(adjacency)


def coset_margin(record: dict[str, object]) -> int:
    central = {-1, 1}
    best_central = max(
        layer["coset_floor"] for layer in record["layers"] if layer["t"] in central
    )
    best_outer = max(
        layer["coset_floor"] for layer in record["layers"] if layer["t"] not in central
    )
    return best_central - best_outer


def packet(record: dict[str, object]) -> tuple[int, int, int, int]:
    return (
        record["H"], int(record["W"]), int(record["D4"]), int(record["C_hd"])
    )


def certificate_summary(base) -> dict[str, object]:
    require(tuple((row[0], row[1], row[2], row[3]) for row in CERTIFICATE_ROWS)
            == tuple((q, r, shard, shards) for q, r, shards in REGIMES
                     for shard in range(shards)), "complete ordered shard bank")
    regime_rows = []
    for q, r, shards in REGIMES:
        selected = tuple(row for row in CERTIFICATE_ROWS if row[:2] == (q, r))
        total = sum(row[4][0] for row in selected)
        expected = q * STRONG_CLASS_COUNTS[q] * ALL_CLASS_COUNTS[r]
        require(len(selected) == shards and total == expected, ("regime count", q, r))
        require(all(row[4][0] == row[4][1] and row[4][2:] == (0, 0)
                    for row in selected), ("regime verdict", q, r))
        regime_rows.append((q, r, STRONG_CLASS_COUNTS[q], ALL_CLASS_COUNTS[r], total))
    require(sum(row[-1] for row in regime_rows) == 4_055_870, "atlas row count")

    maximum = max(
        (Fraction(*row[5][2]), row[5]) for row in CERTIFICATE_ROWS
    )[1]
    minimum_margin = min(row[6][3] for row in CERTIFICATE_ROWS)
    require(maximum[2] == (3_201_028, 5_617_819)
            and maximum[0] == "1111111101110100111110101111101111110111111111111111101",
            "sharp atlas tilt")
    require(minimum_margin == 380, "sharp strict coset margin")

    replayed = []
    for row in CERTIFICATE_ROWS:
        for role, frozen in (("max", row[5]), ("min", row[6])):
            label, expected_packet, expected_rho, expected_margin = frozen
            adjacency = adjacency_from_label(label, 11)
            require(base.is_strong(adjacency), ("extremum strong", row[:4], role))
            record = base.analyze(adjacency)
            require(packet(record) == expected_packet, ("extremum packet", row[:4], role))
            require(base.pair(record["normalized_tilt"]) == expected_rho,
                    ("extremum rho", row[:4], role))
            require(coset_margin(record) == expected_margin,
                    ("extremum margin", row[:4], role))
            replayed.append((row[:4], role, label, expected_packet, expected_rho,
                             record["rational_t"], record["coset_t"],
                             record["actual_t"], expected_margin))
    return {
        "regimes": tuple(regime_rows),
        "rows": 4_055_870,
        "strong_rows": 4_055_870,
        "failures": (0, 0),
        "maximum": maximum,
        "minimum_coset_margin": minimum_margin,
        "certificate_sha256": digest(CERTIFICATE_ROWS),
        "extremal_replays_sha256": digest(replayed),
    }


def structural_controls(base, family) -> dict[str, object]:
    quotient = adjacency_from_label(SHARP_QUOTIENT_LABEL, 5)
    block = adjacency_from_label(SHARP_BLOCK_LABEL, 7)
    require(base.is_strong(quotient), "sharp quotient strong")
    require(tuple((row.bit_count()) for row in block) == (3,) * 7,
            "sharp block regular")
    sharp = substitute(quotient, SHARP_DISTINGUISHED_VERTEX, block)
    sharp_label = label_from_adjacency(sharp)
    require(sharp_label == "1111111101110100111110101111101111110111111111111111101",
            "sharp substitution provenance")

    witness = family.compact_family_row(9)["adjacency"]
    deleted = []
    deleted_vertex = 9
    for source in range(len(witness)):
        if source == deleted_vertex:
            continue
        row = 0
        for target in range(len(witness)):
            if target == deleted_vertex:
                continue
            if witness[source] & (1 << target):
                row |= 1 << (target - (target > deleted_vertex))
        deleted.append(row)
    deleted = tuple(deleted)
    require(label_from_adjacency(deleted) == DELETION_LABEL, "deletion hostile label")
    require(not base.is_strong(deleted), "deletion hostile nonstrong")
    hostile = base.analyze(deleted)
    require(packet(hostile) == (9253, 147483, 5815378668, -4396377498),
            "deletion hostile packet")
    require(base.pair(hostile["normalized_tilt"]) == (732729583, 484614889)
            and hostile["rational_t"] == hostile["coset_t"] == (-3,),
            "deletion hostile noncentrality")

    repairs = []
    for left in range(11):
        for right in range(left + 1, 11):
            adjacency = list(deleted)
            adjacency[left] ^= 1 << right
            adjacency[right] ^= 1 << left
            adjacency = tuple(adjacency)
            if not base.is_strong(adjacency):
                continue
            record = base.analyze(adjacency)
            repairs.append((base.pair(record["normalized_tilt"]),
                            label_from_adjacency(adjacency), record["rational_t"],
                            record["coset_t"], coset_margin(record)))
    require(len(repairs) == 10, "one-edge strong repair count")
    best_repair = max(repairs, key=lambda row: Fraction(*row[0]))
    require(best_repair[:2] == ((5121379056, 12507930943), BEST_REPAIR_LABEL),
            "best strong repair")
    require(all(Fraction(*row[0]) < 1 and row[2] in ((-1,), (1,))
                and row[3] in ((-1,), (1,)) and row[4] > 0 for row in repairs),
            "all strong repairs central")
    return {
        "sharp_provenance": (SHARP_QUOTIENT_LABEL, SHARP_DISTINGUISHED_VERTEX,
                              SHARP_BLOCK_LABEL, sharp_label),
        "deletion_hostile": (DELETION_LABEL, packet(hostile),
                             base.pair(hostile["normalized_tilt"]),
                             hostile["rational_t"], hostile["coset_t"],
                             coset_margin(hostile)),
        "strong_repairs": (len(repairs), best_repair),
    }


ENGINE_PATTERN = re.compile(
    r"rows=(\d+);strong=(\d+);rational_failures=(\d+);coset_failures=(\d+)\n"
    r"max_rho_label=([01]+);packet=\((-?\d+),(-?\d+),(-?\d+),(-?\d+)\);"
    r"rho=(\d+)/(\d+);coset_margin=(-?\d+)\n"
    r"min_coset_label=([01]+);packet=\((-?\d+),(-?\d+),(-?\d+),(-?\d+)\);"
    r"rho=(\d+)/(\d+);coset_margin=(-?\d+)\n"
)


def parse_engine_output(text: str, task: tuple[int, int, int, int]):
    match = ENGINE_PATTERN.fullmatch(text)
    require(match is not None, ("engine transcript", task, text[-200:]))
    values = match.groups()
    counts = tuple(map(int, values[:4]))
    maximum = (values[4], tuple(map(int, values[5:9])),
               tuple(map(int, values[9:11])), int(values[11]))
    minimum = (values[12], tuple(map(int, values[13:17])),
               tuple(map(int, values[17:19])), int(values[19]))
    return (*task, counts, maximum, minimum)


def find_compiler() -> str:
    configured = os.environ.get("THM4144_CXX")
    candidates = ((configured,) if configured else ()) + ("clang++", "g++")
    for candidate in candidates:
        if candidate and shutil.which(candidate):
            return shutil.which(candidate)
    raise RuntimeError("no C++17 compiler found")


def emit_regime(base, q: int, r: int, shard: int, shards: int) -> None:
    require((q, r, shards) in REGIMES and 0 <= shard < shards,
            "declared emission regime")
    gentourng = Path(base.find_gentourng())
    require(file_digest(gentourng) == EXPECTED_GENTOURNG_SHA256,
            "pinned gentourng binary")
    quotients = tuple(
        base.decode_bits(line, q)
        for line in subprocess.check_output(
            [str(gentourng), "-q", "-c", str(q)], text=True
        ).splitlines()
    )
    blocks = tuple(
        base.decode_bits(line, r)
        for line in subprocess.check_output(
            [str(gentourng), "-q", str(r)], text=True
        ).splitlines()
    )
    require(len(quotients) == STRONG_CLASS_COUNTS[q]
            and len(blocks) == ALL_CLASS_COUNTS[r], "nauty class counts")
    buffer = []
    index = 0
    for quotient in quotients:
        for block in blocks:
            for vertex in range(q):
                if index % shards == shard:
                    buffer.append(label_from_adjacency(substitute(quotient, vertex, block)))
                    if len(buffer) == 4096:
                        sys.stdout.write("\n".join(buffer) + "\n")
                        buffer.clear()
                index += 1
    if buffer:
        sys.stdout.write("\n".join(buffer) + "\n")


def full_replay() -> None:
    require(file_digest(ENGINE_PATH) == EXPECTED_ENGINE_SHA256, "pinned primary engine")
    tasks = tuple((q, r, shard, shards) for q, r, shards in REGIMES
                  for shard in range(shards))
    with TemporaryDirectory(prefix="thm4144-primary-") as directory:
        engine = Path(directory) / "engine"
        subprocess.run(
            [find_compiler(), "-O3", "-std=c++17", str(ENGINE_PATH), "-o", str(engine)],
            check=True,
        )

        def run_task(task):
            generator = subprocess.Popen(
                [sys.executable, "-B", str(Path(__file__).resolve()),
                 "--emit-regime", *map(str, task)], stdout=subprocess.PIPE,
            )
            completed = subprocess.run(
                [str(engine)], stdin=generator.stdout, stdout=subprocess.PIPE,
                stderr=subprocess.PIPE, text=True,
            )
            generator.stdout.close()
            generator_status = generator.wait()
            require(generator_status == 0, ("generator status", task, generator_status))
            require(completed.returncode == 0,
                    ("engine status", task, completed.stderr[-400:]))
            return parse_engine_output(completed.stdout, task)

        workers = min(len(tasks), max(1, os.cpu_count() or 1))
        with ThreadPoolExecutor(max_workers=workers) as pool:
            rows = tuple(pool.map(run_task, tasks))
    require(rows == CERTIFICATE_ROWS, "full certificate byte semantics")


def main() -> None:
    parser = ArgumentParser()
    parser.add_argument("--full", action="store_true")
    parser.add_argument("--emit-regime", nargs=4, type=int,
                        metavar=("Q", "R", "SHARD", "SHARDS"))
    arguments = parser.parse_args()
    tree = ast.parse(Path(__file__).read_text(encoding="utf-8"))
    require(not any(isinstance(node, ast.Assert) for node in ast.walk(tree)),
            "no assertions")
    base = load_module("thm4131_base_for_4144", BASE_PATH, EXPECTED_BASE_SHA256)
    if arguments.emit_regime:
        emit_regime(base, *arguments.emit_regime)
        return
    family = load_module("thm4133_family_for_4144", FAMILY_PATH, EXPECTED_FAMILY_SHA256)
    require(file_digest(ENGINE_PATH) == EXPECTED_ENGINE_SHA256, "primary engine hash")
    summary = certificate_summary(base)
    controls = structural_controls(base, family)
    if arguments.full:
        full_replay()
    ledger = {
        "theorem": "THM-4144",
        "statement": (
            "every strong order-eleven tournament with a homogeneous module of size "
            "three through nine has only central rational and exact-coset support-floor optimizers"
        ),
        "class_counts": (tuple(sorted(STRONG_CLASS_COUNTS.items())),
                         tuple(sorted(ALL_CLASS_COUNTS.items()))),
        "certificate": summary,
        "controls": controls,
        "scope": "large homogeneous modules only; pair-module and prime order-eleven tournaments remain open",
    }
    semantic = digest(ledger)
    if EXPECTED_SEMANTIC is not None:
        require(semantic == EXPECTED_SEMANTIC, "frozen semantic digest")
    print("THM4144_ORDER11_LARGE_MODULE_CENTRALITY_20260825")
    print("status=PASS;scope=homogeneous module size 3..9")
    print(f"regimes={summary['regimes']}")
    print(f"construction_rows={summary['rows']};strong_rows={summary['strong_rows']};failures={summary['failures']}")
    print(f"sharp_rho={summary['maximum'][2]};sharp_label={summary['maximum'][0]};sharp_packet={summary['maximum'][1]}")
    print(f"minimum_strict_coset_margin={summary['minimum_coset_margin']}")
    print(f"certificate_sha256={summary['certificate_sha256']};extremal_replays_sha256={summary['extremal_replays_sha256']}")
    print(f"deletion_hostile={controls['deletion_hostile']}")
    print(f"strong_repairs={controls['strong_repairs']}")
    print("residual=homogeneous pairs plus prime tournaments; no actual-maximizer classification")
    print(f"checks={CHECKS}")
    print(f"semantic_sha256={semantic}")


if __name__ == "__main__":
    main()
