from __future__ import annotations

import csv
import hashlib
from fractions import Fraction
from pathlib import Path


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
UNIVERSE = (
    ROOT
    / "05-knowledge/results/lrc14_rank2_wall_graph_complete_pair_closure_thm4326"
    / "thm4231_remainder181194.csv"
)
SCREEN = HERE / "pair_rank3_2over27_screen_all.csv"
EXACT = HERE / "pair_rank3_2over27_exact_v2.csv"

EXPECTED_SHA = {
    UNIVERSE: "9dfbf0a8948bf23016ae40f40f9118020c9429ac60421681a9286fc4d34041a1",
    SCREEN: "a9fed47f8aa1fbd6e3822eb6eecda15805646aebb6d541e53afa5e641a970070",
    EXACT: "c25e57f205c49828fbc93b8b92a7ff5e04eabb8731a53e0b54f224e697186738",
}


def demand(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


for path, expected in EXPECTED_SHA.items():
    demand(path.is_file(), f"missing artifact: {path}")
    demand(sha256(path) == expected, f"hash mismatch: {path.name}")

universe: list[tuple[int, int]] = []
for line in UNIVERSE.read_text(encoding="ascii").splitlines():
    q_text, r_text = line.split(",")
    q, r = int(q_text), int(r_text)
    demand(0 < q < r, "bad universe ordering")
    universe.append((q, r))
demand(len(universe) == 181_194, "universe count changed")
demand(len(set(universe)) == len(universe), "duplicate universe pair")

with SCREEN.open(newline="", encoding="ascii") as handle:
    screen = list(csv.DictReader(handle))
demand(len(screen) == len(universe), "screen count changed")
demand(
    [(int(row["q"]), int(row["r"])) for row in screen] == universe,
    "screen does not exactly preserve the universe order",
)

fallback_pairs: list[tuple[int, int]] = []
coarse_positive: list[tuple[Fraction, dict[str, str]]] = []
for row in screen:
    q, r = int(row["q"]), int(row["r"])
    grid = int(row["grid"])
    mass = int(row["coarse_mass"])
    ticks = int(row["coarse_ticks"])
    positive = row["positive"]
    body = int(row["top9"], 16)
    demand(grid > 0, "nonpositive grid")
    demand(ticks == 27 * mass - 2 * grid, "coarse tick identity failed")
    demand(body.bit_count() == 9, "coarse mask is not a nine-body")
    demand(positive in {"0", "1"}, "invalid positive flag")
    demand((ticks > 0) == (positive == "1"), "positive flag mismatch")
    if positive == "1":
        coarse_positive.append((Fraction(ticks, grid), row))
    else:
        fallback_pairs.append((q, r))

with EXACT.open(newline="", encoding="ascii") as handle:
    exact = list(csv.DictReader(handle))
demand(len(coarse_positive) == 113_996, "coarse-positive count changed")
demand(len(fallback_pairs) == 67_198, "fallback count changed")
demand(len(exact) == len(fallback_pairs), "exact count changed")
demand(
    [(int(row["q"]), int(row["r"])) for row in exact] == fallback_pairs,
    "exact ledger is not precisely the fallback ledger",
)

exact_by_pair: dict[tuple[int, int], dict[str, str]] = {}
exact_positive: list[tuple[Fraction, dict[str, str]]] = []
nodes = 0
prunes = 0
screen_by_pair = {
    (int(row["q"]), int(row["r"])): row for row in screen
}
for row in exact:
    q, r = int(row["q"]), int(row["r"])
    pair = (q, r)
    parent = screen_by_pair[pair]
    grid = int(row["grid"])
    total = int(row["rank3_total"])
    mass = int(row["min_mass"])
    ticks = int(row["min_ticks"])
    body = int(row["least_body"], 16)
    demand(grid == int(parent["grid"]), "exact grid mismatch")
    demand(total == int(parent["rank3_total"]), "exact total mismatch")
    demand(ticks == 27 * mass - 2 * grid, "exact tick identity failed")
    demand(ticks > 0, "nonpositive exact fallback")
    demand(body.bit_count() == 9, "exact mask is not a nine-body")
    demand(int(row["nodes"]) > 0, "empty exact search")
    demand(0 <= int(row["prunes"]) <= int(row["nodes"]), "bad prune count")
    nodes += int(row["nodes"])
    prunes += int(row["prunes"])
    exact_positive.append((Fraction(ticks, grid), row))
    exact_by_pair[pair] = row

demand(nodes == 75_201_674, "node total changed")
demand(prunes == 37_153_312, "prune total changed")

least_coarse_ratio, least_coarse = min(coarse_positive, key=lambda item: item[0])
least_exact_ratio, least_exact = min(exact_positive, key=lambda item: item[0])
demand((least_coarse["q"], least_coarse["r"]) == ("509", "640"),
       "least coarse certificate moved")
demand(least_coarse_ratio == Fraction(1_895_053_421, 515_819_452_388_240),
       "least coarse ratio changed")
demand((least_exact["q"], least_exact["r"]) == ("50", "70"),
       "least exact fallback moved")
demand(least_exact["least_body"] == "011cd402", "hostile body changed")
demand(least_exact_ratio == Fraction(1_565_212_028_297, 5_066_988_726_800),
       "least exact ratio changed")

# For any nine-body H over the fixed pool and pair (q,r), positive component
# count is at most sum(H)+q+r.  The largest nine pool labels sum to 2061.
# If delta = certified_mass-2/27, an added speed s loses at most
# 6*c/(49*s) from the (6/7)-scaled mass. Thus s >= c/(7*delta) suffices;
# equality reaches the closed 4/63 target.
cutoffs: list[tuple[Fraction, dict[str, str], int, int]] = []
for row in screen:
    pair = (int(row["q"]), int(row["r"]))
    if row["positive"] == "1":
        ticks = int(row["coarse_ticks"])
    else:
        ticks = int(exact_by_pair[pair]["min_ticks"])
    grid = int(row["grid"])
    component_bound = 2061 + pair[0] + pair[1]
    cutoff = Fraction(27 * component_bound * grid, 7 * ticks)
    cutoffs.append((cutoff, row, ticks, component_bound))
worst_cutoff, cutoff_row, cutoff_ticks, cutoff_components = max(
    cutoffs, key=lambda item: item[0]
)
demand((cutoff_row["q"], cutoff_row["r"]) == ("509", "640"),
       "cutoff extremizer moved")
demand(cutoff_ticks == 272_887_692_624, "cutoff ticks changed")
demand(cutoff_components == 3210, "cutoff component bound changed")
demand(worst_cutoff == Fraction(6_386_581_705_498_394_400, 1_895_053_421),
       "cutoff ratio changed")
integer_cutoff = -(-worst_cutoff.numerator // worst_cutoff.denominator)
demand(integer_cutoff == 3_370_132_808, "integer cutoff changed")

print("PAIR_RANK3_PACKET_VERIFY_V1")
print(
    "COUNTS universe=181194 coarse_positive=113996 fallback=67198 "
    "exact_positive=67198"
)
print(f"SEARCH nodes={nodes} prunes={prunes}")
print(
    "COARSE_MIN pair=509,640 ticks/grid="
    f"{least_coarse_ratio.numerator}/{least_coarse_ratio.denominator}"
)
print(
    "EXACT_FALLBACK_MIN pair=50,70 body=011cd402 ticks/grid="
    f"{least_exact_ratio.numerator}/{least_exact_ratio.denominator}"
)
print(
    "COFINAL_THIRD_OUTSIDER_CUTOFF pair=509,640 c<=3210 "
    f"s>={integer_cutoff}"
)
print("VERDICT PASS")
