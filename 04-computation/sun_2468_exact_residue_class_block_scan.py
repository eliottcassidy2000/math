#!/usr/bin/env python3
"""Exact block scanner for one target residue class in Sun's 2-4-6-8 problem.

The scanner counts the canonical OEIS A306477 representations

    n = C(w+2,2) + C(x+3,4) + C(y+5,6) + C(z+7,8),  w,x,y,z >= 0,

for every n in [start, start+width) with n == residue (mod modulus).  It is
an independent, target-block-oriented complement to the one-input certificate
used by THM-4026.  Counts saturate at 65535, but the zero/nonzero predicate and
all counts below that threshold are exact.

The implementation is Python only as an orchestration/reference path.  For
large blocks it compiles the embedded C++20 engine into a requested temporary
directory and runs it.  No floating-point decision is load-bearing: the
initial long-double square-root estimate is repaired by exact __int128 tests.

Example (the five-prime class containing THM-4026's counterexample):

    python -B 04-computation/sun_2468_exact_residue_class_block_scan.py \
      --start 802000000000000 --width 100000000000000 \
      --modulus 1062347 --residue 459490 --low-cutoff 1
"""

from __future__ import annotations

import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
from math import comb
import os
from pathlib import Path
import subprocess
import sys
import tempfile


TARGET_CLASS_MODULUS = 1_062_347
TARGET_CLASS_RESIDUE = 459_490
TARGET_HOLE = 896_315_812_331_399

# These half-open interval blocks partition [1, 1_002_000_000_000_000).
# The remaining entries are exact expected controls, not search heuristics:
# start, width, base, c4, c6, c8, triples, marks, targets, min, argmin, zeros.
TARGET_CLASS_CENSUS = (
    (1, 1_999_999_999_999, 459_490, 2631, 334, 127,
     99_211_900, 51_852_255, 1_882_624, 3, 636_132_780_743, 0),
    (2_000_000_000_000, 10_000_000_000_000, 2_000_000_418_018, 4119, 451, 159,
     263_480_161, 283_941_711, 9_413_120, 2, 9_480_152_433_497, 0),
    (12_000_000_000_000, 90_000_000_000_000, 12_000_000_210_658, 7033, 645, 209,
     845_239_561, 2_789_161_302, 84_718_082, 2, 99_934_198_737_404, 0),
    (102_000_000_000_000, 100_000_000_000_000, 102_000_000_469_112, 8343, 723, 228,
     1_225_867_610, 3_272_771_677, 94_131_202, 3, 122_533_388_547_887, 0),
    (202_000_000_000_000, 100_000_000_000_000, 202_000_000_520_206, 9226, 773, 240,
     1_525_588_515, 3_321_577_932, 94_131_202, 4, 202_687_664_107_388, 0),
    (302_000_000_000_000, 100_000_000_000_000, 302_000_000_571_300, 9910, 811, 249,
     1_782_451_689, 3_390_881_393, 94_131_202, 3, 350_011_776_754_814, 0),
    (402_000_000_000_000, 200_000_000_000_000, 402_000_000_622_394, 10963, 868, 262,
     2_220_090_999, 6_842_995_836, 188_262_404, 4, 406_999_014_660_698, 0),
    (602_000_000_000_000, 200_000_000_000_000, 602_000_000_724_582, 11778, 911, 272,
     2_594_651_419, 6_944_628_451, 188_262_404, 4, 677_862_454_957_862, 0),
    (802_000_000_000_000, 200_000_000_000_000, 802_000_000_826_770, 12452, 945, 279,
     2_927_893_995, 6_984_950_518, 188_262_404, 0, TARGET_HOLE, 1),
)


CPP_SOURCE = r'''#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

using u64 = std::uint64_t;
using u128 = unsigned __int128;

static u64 choose_exact(u64 n, int k) {
    if (n < static_cast<u64>(k)) {
        return 0;
    }
    u128 value = 1;
    for (int i = 1; i <= k; ++i) {
        value = value * static_cast<u64>(n - k + i) / static_cast<u64>(i);
    }
    if (value > std::numeric_limits<u64>::max()) {
        throw std::overflow_error("binomial coefficient exceeds uint64");
    }
    return static_cast<u64>(value);
}

static u128 triangular(u64 i) {
    return static_cast<u128>(i) * static_cast<u128>(i + 1) / 2;
}

static u64 ceil_triangular_index(u64 q) {
    long double raw = std::sqrt(static_cast<long double>(8) * q + 1);
    u64 d = static_cast<u64>(raw);
    const u128 discriminant = static_cast<u128>(8) * q + 1;
    while (static_cast<u128>(d) * d < discriminant) {
        ++d;
    }
    while (d > 0 && static_cast<u128>(d - 1) * (d - 1) >= discriminant) {
        --d;
    }
    u64 i = (d - 1) / 2;
    while (triangular(i) < q) {
        ++i;
    }
    while (i > 0 && triangular(i - 1) >= q) {
        --i;
    }
    return i;
}

static u64 floor_triangular_index(u64 q) {
    const u64 i = ceil_triangular_index(q);
    return triangular(i) > q ? i - 1 : i;
}

int main(int argc, char** argv) {
    if (argc != 6) {
        std::cerr << "usage: engine START WIDTH MODULUS RESIDUE LOW_CUTOFF\n";
        return 2;
    }
    const u64 start = std::stoull(argv[1]);
    const u64 width = std::stoull(argv[2]);
    const u64 modulus = std::stoull(argv[3]);
    const u64 residue = std::stoull(argv[4]);
    const unsigned low_cutoff = static_cast<unsigned>(std::stoul(argv[5]));
    if (start == 0 || width == 0 || modulus == 0 || residue >= modulus) {
        throw std::invalid_argument("require START,WIDTH,MODULUS>0 and RESIDUE<MODULUS");
    }
    if (width - 1 > std::numeric_limits<u64>::max() - start) {
        throw std::overflow_error("target interval overflows uint64");
    }
    if (modulus > std::numeric_limits<u64>::max() / 2) {
        throw std::overflow_error("2*MODULUS overflows uint64");
    }

    const u64 high = start + width - 1;
    const u64 delta = (residue + modulus - start % modulus) % modulus;
    if (delta > std::numeric_limits<u64>::max() - start) {
        throw std::overflow_error("first target overflows uint64");
    }
    const u64 base = start + delta;
    if (base > high) {
        std::cout << "RESULT=PASS\n"
                  << "targets=0\nzeros=0\n";
        return 0;
    }

    std::vector<u64> c4;
    std::vector<u64> c6;
    std::vector<u64> c8;
    for (int k : {4, 6, 8}) {
        auto& values = (k == 4 ? c4 : (k == 6 ? c6 : c8));
        for (u64 x = 0;; ++x) {
            const u64 value = choose_exact(x + static_cast<u64>(k - 1), k);
            // Equality cannot occur in a representation because C(w+2,2)>=1.
            if (value >= high) {
                break;
            }
            values.push_back(value);
        }
    }

    const std::size_t target_count = static_cast<std::size_t>((high - base) / modulus + 1);
    std::vector<std::uint16_t> counts(target_count, 0);
    const u64 triangular_period = 2 * modulus;
    std::vector<std::vector<u64>> triangular_roots(static_cast<std::size_t>(modulus));
    for (u64 i = 0; i < triangular_period; ++i) {
        const u64 value = static_cast<u64>(triangular(i) % modulus);
        triangular_roots[static_cast<std::size_t>(value)].push_back(i);
    }

    u64 feasible_higher_triples = 0;
    u64 representation_marks = 0;
    for (u64 z_value : c8) {
        for (u64 y_value : c6) {
            if (static_cast<u128>(z_value) + y_value >= high) {
                break;
            }
            for (u64 x_value : c4) {
                const u128 sum128 = static_cast<u128>(x_value) + y_value + z_value;
                if (sum128 >= high) {
                    break;
                }
                const u64 sum = static_cast<u64>(sum128);
                ++feasible_higher_triples;
                const u64 low_remainder = start > sum ? start - sum : 1;
                const u64 high_remainder = high - sum;
                const u64 first_index = ceil_triangular_index(low_remainder);
                const u64 last_index = floor_triangular_index(high_remainder);
                if (first_index > last_index) {
                    continue;
                }
                const u64 needed = (residue + modulus - sum % modulus) % modulus;
                for (u64 root : triangular_roots[static_cast<std::size_t>(needed)]) {
                    const u64 offset = (root + triangular_period - first_index % triangular_period)
                                     % triangular_period;
                    for (u64 i = first_index + offset; i <= last_index;) {
                        const u128 target128 = static_cast<u128>(sum) + triangular(i);
                        if (target128 > high) {
                            throw std::logic_error("triangular bound escaped target interval");
                        }
                        const u64 target = static_cast<u64>(target128);
                        if (target < base || target % modulus != residue) {
                            throw std::logic_error("mark escaped selected residue class");
                        }
                        const std::size_t address = static_cast<std::size_t>((target - base) / modulus);
                        if (address >= counts.size()) {
                            throw std::logic_error("mark address escaped count array");
                        }
                        if (counts[address] != std::numeric_limits<std::uint16_t>::max()) {
                            ++counts[address];
                        }
                        ++representation_marks;
                        if (last_index - i < triangular_period) {
                            break;
                        }
                        i += triangular_period;
                    }
                }
            }
        }
    }

    const auto minimum_it = std::min_element(counts.begin(), counts.end());
    const unsigned minimum = *minimum_it;
    const u64 minimum_target = base + modulus * static_cast<u64>(minimum_it - counts.begin());
    std::size_t zeros = 0;
    std::size_t low_count = 0;
    for (std::size_t address = 0; address < counts.size(); ++address) {
        const unsigned count = counts[address];
        if (count == 0) {
            ++zeros;
        }
        if (count <= low_cutoff) {
            ++low_count;
            std::cout << "LOW target="
                      << (base + modulus * static_cast<u64>(address))
                      << " count=" << count << "\n";
        }
    }

    std::cout << "ENGINE=sun_2468_exact_residue_class_block_scan_v1\n"
              << "interval=[" << start << "," << high << "]\n"
              << "modulus=" << modulus << "\n"
              << "residue=" << residue << "\n"
              << "base=" << base << "\n"
              << "c4_values=" << c4.size() << "\n"
              << "c6_values=" << c6.size() << "\n"
              << "c8_values=" << c8.size() << "\n"
              << "feasible_higher_triples=" << feasible_higher_triples << "\n"
              << "representation_marks=" << representation_marks << "\n"
              << "targets=" << counts.size() << "\n"
              << "minimum_count=" << minimum << "\n"
              << "minimum_target=" << minimum_target << "\n"
              << "low_count_targets=" << low_count << "\n"
              << "zeros=" << zeros << "\n"
              << "RESULT=PASS\n";
}
'''


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--start", type=int)
    parser.add_argument("--width", type=int)
    parser.add_argument("--modulus", type=int)
    parser.add_argument("--residue", type=int)
    parser.add_argument("--low-cutoff", type=int, default=0)
    parser.add_argument(
        "--target-class-census",
        action="store_true",
        help="replay the exact nine-block THM-4040 target-class census",
    )
    parser.add_argument("--jobs", type=int, default=3, help="parallel census engines (default: 3)")
    parser.add_argument(
        "--build-dir",
        type=Path,
        help="directory for the generated C++ source and executable; defaults to a temporary directory",
    )
    parser.add_argument("--compiler", default=os.environ.get("CXX", "g++"))
    return parser.parse_args()


def require_inputs(args: argparse.Namespace) -> None:
    if args.target_class_census:
        if any(value is not None for value in (args.start, args.width, args.modulus, args.residue)):
            raise SystemExit("do not combine --target-class-census with single-block inputs")
        if args.jobs <= 0:
            raise SystemExit("--jobs must be positive")
        return
    if None in (args.start, args.width, args.modulus, args.residue):
        raise SystemExit("single-block mode requires START, WIDTH, MODULUS, and RESIDUE")
    if args.start <= 0 or args.width <= 0 or args.modulus <= 0:
        raise SystemExit("START, WIDTH, and MODULUS must be positive")
    if not 0 <= args.residue < args.modulus:
        raise SystemExit("require 0 <= RESIDUE < MODULUS")
    if not 0 <= args.low_cutoff <= 65535:
        raise SystemExit("require 0 <= LOW_CUTOFF <= 65535")


def require(condition: bool, message: str) -> None:
    if not condition:
        raise AssertionError(message)


def compile_engine(args: argparse.Namespace, build_dir: Path) -> Path:
    build_dir.mkdir(parents=True, exist_ok=True)
    source_path = build_dir / "sun_2468_exact_residue_class_block_scan.cpp"
    executable_name = "sun_2468_exact_residue_class_block_scan.exe" if os.name == "nt" else "sun_2468_exact_residue_class_block_scan"
    executable_path = build_dir / executable_name
    source_path.write_text(CPP_SOURCE, encoding="utf-8", newline="\n")
    digest = hashlib.sha256(CPP_SOURCE.encode("utf-8")).hexdigest()
    print(f"embedded_cpp_sha256={digest}")
    compile_command = [
        args.compiler,
        "-std=c++20",
        "-O3",
        "-Wall",
        "-Wextra",
        "-Wconversion",
        "-Wshadow",
        str(source_path),
        "-o",
        str(executable_path),
    ]
    subprocess.run(compile_command, check=True)
    return executable_path


def run_engine(
    executable_path: Path,
    start: int,
    width: int,
    modulus: int,
    residue: int,
    low_cutoff: int,
) -> str:
    run_command = [
        str(executable_path),
        str(start),
        str(width),
        str(modulus),
        str(residue),
        str(low_cutoff),
    ]
    completed = subprocess.run(run_command, check=True, text=True, capture_output=True)
    require(not completed.stderr, f"engine wrote stderr: {completed.stderr}")
    return completed.stdout


def parse_engine_output(output: str) -> tuple[dict[str, int | str], dict[int, int]]:
    fields: dict[str, int | str] = {}
    low: dict[int, int] = {}
    for line in output.splitlines():
        if line.startswith("LOW target="):
            pieces = line.replace("LOW target=", "").replace(" count=", " ").split()
            require(len(pieces) == 2, f"malformed LOW line: {line}")
            low[int(pieces[0])] = int(pieces[1])
        elif "=" in line and not line.startswith("interval="):
            key, value = line.split("=", 1)
            fields[key] = int(value) if value.isdigit() else value
    require(fields.get("RESULT") == "PASS", "engine did not return RESULT=PASS")
    return fields, low


def small_independent_control(executable_path: Path) -> None:
    start, width, modulus, residue = 1, 10_000, 33, 20
    high = start + width - 1
    targets = tuple(n for n in range(start, high + 1) if n % modulus == residue)
    direct = {n: 0 for n in targets}
    atom_lists = []
    for k, lower in ((2, 2), (4, 3), (6, 5), (8, 7)):
        values = []
        top = lower
        while (value := comb(top, k)) <= high:
            values.append(value)
            top += 1
        atom_lists.append(values)
    for a in atom_lists[0]:
        for b in atom_lists[1]:
            if a + b > high:
                break
            for c in atom_lists[2]:
                if a + b + c > high:
                    break
                for d in atom_lists[3]:
                    total = a + b + c + d
                    if total > high:
                        break
                    if total in direct:
                        direct[total] += 1
    output = run_engine(executable_path, start, width, modulus, residue, 65_535)
    fields, low = parse_engine_output(output)
    require(low == direct, "C++ block counts disagree with direct Python Cartesian control")
    require(fields["minimum_count"] == 2 and fields["minimum_target"] == 9590,
            "small-control minimum changed")
    print("small_direct_control=PASS;targets=303;minimum_count=2;minimum_target=9590;zeros=0")


def census_block(executable_path: Path, index: int, row: tuple[int, ...]):
    start, width = row[:2]
    output = run_engine(
        executable_path, start, width, TARGET_CLASS_MODULUS, TARGET_CLASS_RESIDUE, 0
    )
    fields, low = parse_engine_output(output)
    expected_keys = (
        "base", "c4_values", "c6_values", "c8_values", "feasible_higher_triples",
        "representation_marks", "targets", "minimum_count", "minimum_target", "zeros",
    )
    expected_values = row[2:]
    actual_values = tuple(fields[key] for key in expected_keys)
    require(actual_values == expected_values,
            f"census block {index} mismatch: {actual_values} != {expected_values}")
    expected_low = {TARGET_HOLE: 0} if row[-1] else {}
    require(low == expected_low, f"census block {index} LOW rows changed: {low}")
    high = start + width - 1
    return (
        f"block={index};interval=[{start},{high}];base={fields['base']};"
        f"triples={fields['feasible_higher_triples']};marks={fields['representation_marks']};"
        f"targets={fields['targets']};minimum_count={fields['minimum_count']};"
        f"minimum_target={fields['minimum_target']};zeros={fields['zeros']}"
    )


def target_class_census(executable_path: Path, jobs: int) -> None:
    small_independent_control(executable_path)
    with ThreadPoolExecutor(max_workers=jobs) as executor:
        futures = [
            executor.submit(census_block, executable_path, index, row)
            for index, row in enumerate(TARGET_CLASS_CENSUS, 1)
        ]
        summaries = [future.result() for future in futures]
    for summary in summaries:
        print(summary)
    total_targets = sum(row[-4] for row in TARGET_CLASS_CENSUS)
    require(total_targets == 943_194_644, "wrong target census cardinality")
    print(f"census_modulus={TARGET_CLASS_MODULUS};census_residue={TARGET_CLASS_RESIDUE}")
    print("census_interval=[1,1001999999999999]")
    print(f"census_targets={total_targets}")
    print(f"census_holes=[{TARGET_HOLE}]")
    print("RESULT=PASS")


def build_and_run(args: argparse.Namespace, build_dir: Path) -> int:
    executable_path = compile_engine(args, build_dir)
    if args.target_class_census:
        target_class_census(executable_path, args.jobs)
        return 0
    output = run_engine(
        executable_path, args.start, args.width, args.modulus, args.residue, args.low_cutoff
    )
    sys.stdout.write(output)
    return 0


def main() -> int:
    args = parse_args()
    require_inputs(args)
    if args.build_dir is not None:
        return build_and_run(args, args.build_dir.resolve())
    with tempfile.TemporaryDirectory(prefix="sun2468-class-scan-") as temporary:
        return build_and_run(args, Path(temporary))


if __name__ == "__main__":
    sys.exit(main())
