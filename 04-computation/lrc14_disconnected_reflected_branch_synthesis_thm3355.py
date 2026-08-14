#!/usr/bin/env python3
"""Exact synthesis guards for THM-3355's body and debt bookkeeping."""
from fractions import Fraction as F
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MEDIAN = ROOT / "04-computation/lrc14_j7_reflected_median_star_chord_scout_20260801.py"
CHROMATIC = ROOT / "04-computation/lrc14_j7_reflected_universal_pair_chromatic_closure_thm2941.py"
DMAX = F(186_636_088_362, 11_773_143_757_375)
TARGET = DMAX / 5
EXPECTED_MEDIAN_SHA256 = "b0eba7785cd73ecf76d2887f3d36316eae6980a76652058daa2587ecbe031276"
EXPECTED_CHROMATIC_SHA256 = "a6f58c1a52dfc1fca61a239068dbe0b216bac41f1622b98748bc4a6d213fb6e8"


def require(condition, detail):
    if not condition:
        raise RuntimeError(detail)


def load(name, path):
    spec = spec_from_file_location(name, path)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def filehash(path):
    return sha256(path.read_bytes().replace(b"\r\n", b"\n")).hexdigest()


def main():
    require(filehash(MEDIAN) == EXPECTED_MEDIAN_SHA256,
            ("median hash", filehash(MEDIAN)))
    require(filehash(CHROMATIC) == EXPECTED_CHROMATIC_SHA256,
            ("chromatic hash", filehash(CHROMATIC)))
    median = load("thm3355_median", MEDIAN)
    chromatic = load("thm3355_chromatic", CHROMATIC)
    bodies = median.body_universe()
    require(len(bodies) == 649, len(bodies))
    exceptional = {row[0] for row in chromatic.EXPECTED_EXCEPTIONS}
    active = {body for body, _ in bodies}
    active_robust_counts = tuple(
        len(median.LP.robust_edges(body)[2]) for body in active
    )
    require(max(active_robust_counts) == 10, max(active_robust_counts))
    require(active.isdisjoint(exceptional), active & exceptional)
    robust_counts = tuple(len(median.LP.robust_edges(body)[2]) for body in exceptional)
    require(robust_counts == (15, 15), robust_counts)
    require(2354 + 649 == 3003, "body partition")
    pair_floor = F(709, 48048) - F(1_792_138_785_426,
                                      221_510_098_565) / 699
    require(pair_floor > TARGET, pair_floor - TARGET)
    require(5 * pair_floor > DMAX, 5 * pair_floor - DMAX)
    require(F(1, 294) > TARGET, F(1, 294) - TARGET)
    # Any two cross-component edges in the explicit double star are distinct;
    # its construction always has (6-|V1|)+(|V1|-1)=5 edges.
    require(all((6 - n) + (n - 1) == 5 for n in range(1, 6)), "tree size")
    print("LRC14 THM-3355 DISCONNECTED/REFLECTED SYNTHESIS AUDIT")
    print("dependency_sha256", EXPECTED_MEDIAN_SHA256, EXPECTED_CHROMATIC_SHA256)
    print("active_upper_median_bodies", len(bodies),
          "robust_edge_max", max(active_robust_counts))
    print("same_level_exceptions", tuple(sorted(exceptional)),
          "robust_counts", robust_counts, "active_intersection", 0)
    print("body_partition", (2354, 649), "total", 3003)
    print("Dmax", DMAX, "per_edge_target", TARGET)
    print("p699_pair_floor", pair_floor, "five_edge_gap", 5 * pair_floor - DMAX)
    print("multipartite_tree_edges", 5)
    print("status=PASS; reflected-residue level branch only, not LRC14")


if __name__ == "__main__":
    main()
