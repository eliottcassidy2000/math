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
import hashlib
import os
from pathlib import Path
import subprocess
import sys
import tempfile


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
    parser.add_argument("--start", type=int, required=True)
    parser.add_argument("--width", type=int, required=True)
    parser.add_argument("--modulus", type=int, required=True)
    parser.add_argument("--residue", type=int, required=True)
    parser.add_argument("--low-cutoff", type=int, default=0)
    parser.add_argument(
        "--build-dir",
        type=Path,
        help="directory for the generated C++ source and executable; defaults to a temporary directory",
    )
    parser.add_argument("--compiler", default=os.environ.get("CXX", "g++"))
    return parser.parse_args()


def require_inputs(args: argparse.Namespace) -> None:
    if args.start <= 0 or args.width <= 0 or args.modulus <= 0:
        raise SystemExit("START, WIDTH, and MODULUS must be positive")
    if not 0 <= args.residue < args.modulus:
        raise SystemExit("require 0 <= RESIDUE < MODULUS")
    if not 0 <= args.low_cutoff <= 65535:
        raise SystemExit("require 0 <= LOW_CUTOFF <= 65535")


def build_and_run(args: argparse.Namespace, build_dir: Path) -> int:
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
    run_command = [
        str(executable_path),
        str(args.start),
        str(args.width),
        str(args.modulus),
        str(args.residue),
        str(args.low_cutoff),
    ]
    completed = subprocess.run(run_command, check=False)
    return completed.returncode


def main() -> int:
    args = parse_args()
    require_inputs(args)
    if args.build_dir is not None:
        return build_and_run(args, args.build_dir.resolve())
    with tempfile.TemporaryDirectory(prefix="sun2468-class-scan-") as temporary:
        return build_and_run(args, Path(temporary))


if __name__ == "__main__":
    sys.exit(main())
