from __future__ import annotations

from bisect import bisect_left
from collections import Counter
from importlib.util import module_from_spec, spec_from_file_location
from math import comb
from pathlib import Path
import re
import subprocess
import tempfile


ROOT = Path(__file__).resolve().parents[1]
SCANNER = ROOT / "04-computation" / "sun_2468_exact_residue_class_block_scan.py"
MODULUS = 1_062_347
RESIDUE = 459_490
TARGET = 896_315_812_331_399
UINT64_MAX = 2**64 - 1
UINT128_MAX = 2**128 - 1


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def load_cpp_source() -> str:
    spec = spec_from_file_location("sun_block_scan", SCANNER)
    require(spec is not None and spec.loader is not None, "cannot import scanner")
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.CPP_SOURCE


def higher_values(k: int, high: int) -> list[int]:
    values = []
    shifted = k - 1
    while True:
        value = comb(shifted, k)
        if value >= high:
            return values
        values.append(value)
        shifted += 1


def direct_counts(start: int, high: int, modulus: int, residue: int):
    targets = [n for n in range(start, high + 1) if n % modulus == residue]
    counts = Counter({target: 0 for target in targets})
    c4 = higher_values(4, high)
    c6 = higher_values(6, high)
    c8 = higher_values(8, high)
    feasible = 0
    for z_value in c8:
        for y_value in c6:
            if z_value + y_value >= high:
                break
            for x_value in c4:
                higher = z_value + y_value + x_value
                if higher >= high:
                    break
                feasible += 1
                triangular_index = 1
                while True:
                    target = higher + triangular_index * (triangular_index + 1) // 2
                    if target > high:
                        break
                    if target >= start and target % modulus == residue:
                        counts[target] += 1
                    triangular_index += 1
    return counts, (len(c4), len(c6), len(c8)), feasible


LOW_RE = re.compile(r"^LOW target=(\d+) count=(\d+)$")


def parse_engine(stdout: str):
    fields = {}
    low = {}
    for line in stdout.splitlines():
        match = LOW_RE.match(line)
        if match:
            low[int(match.group(1))] = int(match.group(2))
        elif "=" in line and not line.startswith("interval="):
            key, value = line.split("=", 1)
            fields[key] = value
    return fields, low


def run_engine(executable: Path, start: int, high: int, modulus: int, residue: int):
    completed = subprocess.run(
        [str(executable), str(start), str(high - start + 1), str(modulus), str(residue), "65535"],
        check=False,
        text=True,
        capture_output=True,
    )
    require(completed.returncode == 0, f"engine failed: {completed.stderr}")
    return parse_engine(completed.stdout)


def check_small_exhaustive(executable: Path) -> int:
    cases = [
        (1, 1, 5, 1),
        (1, 1, 5, 2),  # no selected target
        (1, 200, 2, 0),
        (1, 200, 2, 1),
        (7, 500, 3, 0),
        (7, 500, 3, 1),
        (35, 3034, 4, 3),
        (90, 5089, 7, 0),
        (1, 10_000, 33, 20),
    ]
    checks = 0
    for start, high, modulus, residue in cases:
        direct, sizes, feasible = direct_counts(start, high, modulus, residue)
        fields, low = run_engine(executable, start, high, modulus, residue)
        require(fields.get("RESULT") == "PASS", f"PASS marker for case {(start, high, modulus, residue)}")
        require(int(fields["targets"]) == len(direct), "target count")
        require(int(fields["zeros"]) == sum(value == 0 for value in direct.values()), "zero count")
        if not direct:
            require(low == {}, "empty class LOW output")
            checks += 4
            continue
        require(low == dict(direct), f"all exact counts for case {(start, high, modulus, residue)}")
        first = min(direct)
        minimum_target, minimum = min(direct.items(), key=lambda item: (item[1], item[0]))
        require(int(fields["base"]) == first, "base target")
        require(int(fields["minimum_count"]) == minimum, "minimum count")
        require(int(fields["minimum_target"]) == minimum_target, "minimum address")
        require(
            tuple(int(fields[f"c{k}_values"]) for k in (4, 6, 8)) == sizes,
            "higher value box",
        )
        require(int(fields["feasible_higher_triples"]) == feasible, "higher triple universe")
        require(int(fields["representation_marks"]) == sum(direct.values()), "mark total")
        checks += 11
    return checks


BLOCKS = [
    (1, 1_999_999_999_999),
    (2_000_000_000_000, 11_999_999_999_999),
    (12_000_000_000_000, 101_999_999_999_999),
    (102_000_000_000_000, 201_999_999_999_999),
    (202_000_000_000_000, 301_999_999_999_999),
    (302_000_000_000_000, 401_999_999_999_999),
    (402_000_000_000_000, 601_999_999_999_999),
    (602_000_000_000_000, 801_999_999_999_999),
    (802_000_000_000_000, 1_001_999_999_999_999),
]

LAST_THREE = [
    {
        "base": 402_000_000_622_394,
        "sizes": (10_963, 868, 262),
        "triples": 2_220_090_999,
        "marks": 6_842_995_836,
        "targets": 188_262_404,
        "minimum": 4,
        "minimum_target": 406_999_014_660_698,
        "zeros": 0,
    },
    {
        "base": 602_000_000_724_582,
        "sizes": (11_778, 911, 272),
        "triples": 2_594_651_419,
        "marks": 6_944_628_451,
        "targets": 188_262_404,
        "minimum": 4,
        "minimum_target": 677_862_454_957_862,
        "zeros": 0,
    },
    {
        "base": 802_000_000_826_770,
        "sizes": (12_452, 945, 279),
        "triples": 2_927_893_995,
        "marks": 6_984_950_518,
        "targets": 188_262_404,
        "minimum": 0,
        "minimum_target": TARGET,
        "zeros": 1,
    },
]


def class_base(start: int) -> int:
    return start + (RESIDUE - start) % MODULUS


def class_count(start: int, high: int) -> int:
    base = class_base(start)
    return 0 if base > high else (high - base) // MODULUS + 1


def feasible_triple_count(high: int, c4: list[int], c6: list[int], c8: list[int]) -> int:
    total = 0
    for z_value in c8:
        for y_value in c6:
            remaining = high - z_value - y_value
            if remaining <= 0:
                break
            total += bisect_left(c4, remaining)
    return total


def check_manifest_arithmetic() -> int:
    checks = 0
    require(BLOCKS[0][0] == 1, "manifest starts at one")
    for previous, current in zip(BLOCKS, BLOCKS[1:]):
        require(previous[1] + 1 == current[0], "contiguous manifest blocks")
        checks += 1
    bound = BLOCKS[-1][1]
    global_base = class_base(1)
    total = sum(class_count(start, high) for start, high in BLOCKS)
    direct_total = (bound - global_base) // MODULUS + 1
    require(total == direct_total == 943_194_644, "manifest class-address partition")
    require(TARGET % MODULUS == RESIDUE, "counterexample class")
    require(BLOCKS[-1][0] <= TARGET <= BLOCKS[-1][1], "counterexample final block")
    require((TARGET - global_base) // MODULUS == 843_712_847, "counterexample global address")
    checks += 5

    for (start, high), expected in zip(BLOCKS[-3:], LAST_THREE):
        c4 = higher_values(4, high)
        c6 = higher_values(6, high)
        c8 = higher_values(8, high)
        require(class_base(start) == expected["base"], "late block base")
        require(class_count(start, high) == expected["targets"], "late block targets")
        require((len(c4), len(c6), len(c8)) == expected["sizes"], "late atom box")
        require(
            feasible_triple_count(high, c4, c6, c8) == expected["triples"],
            "late feasible triple universe",
        )
        require(expected["marks"] < 2**64, "mark counter uint64 safety")
        require(start <= expected["minimum_target"] <= high, "minimum interval address")
        require(expected["minimum_target"] % MODULUS == RESIDUE, "minimum class address")
        checks += 7

    # Numeric safety in the actual census, independent of the engine's guards.
    high = bound
    for k in (4, 6, 8):
        values = higher_values(k, high)
        top = len(values) + k - 1
        first_excluded = comb(top, k)
        require(first_excluded < UINT64_MAX, f"C(n,{k}) uint64")
        require(k * first_excluded < UINT128_MAX, f"C(n,{k}) multiply-before-divide")
        checks += 2
    require(8 * high + 1 < UINT128_MAX, "triangular discriminant uint128")
    require(2 * MODULUS < UINT64_MAX, "triangular root period uint64")
    require(2 * MODULUS + 50_000_000 < UINT64_MAX, "root offset/index addition")
    checks += 3
    return checks


def main() -> None:
    source = load_cpp_source()
    with tempfile.TemporaryDirectory(prefix="sun-block-hostile-audit-") as temporary:
        directory = Path(temporary)
        source_path = directory / "engine.cpp"
        executable = directory / ("engine.exe" if __import__("os").name == "nt" else "engine")
        source_path.write_text(source, encoding="utf-8", newline="\n")
        subprocess.run(
            ["g++", "-std=c++20", "-O2", "-fno-tree-vectorize", "-Wall", "-Wextra", str(source_path), "-o", str(executable)],
            check=True,
        )
        small_checks = check_small_exhaustive(executable)
    manifest_checks = check_manifest_arithmetic()
    print(f"small_exact_checks={small_checks}")
    print(f"manifest_arithmetic_checks={manifest_checks}")
    print(f"manifest_blocks={len(BLOCKS)}")
    print("manifest_union=[1,1001999999999999]")
    print("manifest_class_targets=943194644")
    print("claimed_zero_vector=0,0,0,0,0,0,0,0,1")
    print(f"claimed_unique_zero={TARGET}")
    print("all_checks=PASS")


if __name__ == "__main__":
    main()
