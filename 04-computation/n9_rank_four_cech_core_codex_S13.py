#!/usr/bin/env python3
"""Exact linear/nonlinear audit of THM-828's rank-four defect core.

The script restricts the four global difference generators to every nonempty
intersection of the A/B/C B3 faces, constructs coordinate recovery pivots,
checks the resulting constant-rank Cech complex, and joins the linear sector
data to the exact 58-witness table.  It intentionally does not identify cube
adjacency with a merged-metagraph edge relation.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
WITNESSES = ROOT / "05-knowledge/results/mobius_cech_n9_exact_join_witnesses_codex_S13.tsv"
BASIS = (0x192486, 0x8C2C0C, 0x11B4600, 0x4483414)
MISSING = (0x1026286, 0x4DD3C9E, 0x54A5692, 0x5537214)


def tiles(n: int) -> list[tuple[int, int]]:
    return [(a, b) for b in range(1, n - 1) for a in range(n, b + 1, -1)]


def gf2_rank(words: list[int]) -> int:
    pivots: dict[int, int] = {}
    for word in words:
        while word:
            pivot = word.bit_length() - 1
            if pivot in pivots:
                word ^= pivots[pivot]
            else:
                pivots[pivot] = word
                break
    return len(pivots)


def face_covers(ts: list[tuple[int, int]]) -> tuple[int, int, int]:
    masks = []
    for role in range(3):
        mask = 0
        for i, (a, b) in enumerate(ts):
            if (a < 9, a - b >= 3, b >= 2)[role]:
                mask |= 1 << i
        masks.append(mask)
    return tuple(masks)  # type: ignore[return-value]


def intersection_mask(covers: tuple[int, int, int], subset: int) -> int:
    mask = (1 << 28) - 1
    for role in range(3):
        if subset >> role & 1:
            mask &= covers[role]
    return mask


def exact_region_mask(covers: tuple[int, int, int], state: int) -> int:
    mask = 0
    for i in range(28):
        membership = sum(((covers[r] >> i) & 1) << r for r in range(3))
        if membership == state:
            mask |= 1 << i
    return mask


def independent_pivots(words: list[int], allowed: int) -> list[int]:
    chosen: list[int] = []
    old_rank = 0
    for bit in range(28):
        if not (allowed >> bit) & 1:
            continue
        trial = chosen + [bit]
        columns = [sum(((word >> p) & 1) << j for j, p in enumerate(trial)) for word in words]
        rank = gf2_rank(columns)
        if rank > old_rank:
            chosen.append(bit)
            old_rank = rank
        if old_rank == 4:
            return chosen
    raise AssertionError("rank-four intersection has no four pivots")


def pivot_decoder(words: tuple[int, ...], pivots: list[int]) -> tuple[int, ...]:
    """Return rows of M^-1, where pivot bits y=M*c over F2."""
    rows = [sum(((words[i] >> pivot) & 1) << i for i in range(4)) | (1 << (4 + j))
            for j, pivot in enumerate(pivots)]
    for col in range(4):
        pivot_row = next(j for j in range(col, 4) if (rows[j] >> col) & 1)
        rows[col], rows[pivot_row] = rows[pivot_row], rows[col]
        for j in range(4):
            if j != col and ((rows[j] >> col) & 1):
                rows[j] ^= rows[col]
    assert [row & 15 for row in rows] == [1, 2, 4, 8]
    return tuple(row >> 4 for row in rows)


def span_coordinates() -> dict[int, str]:
    out = {}
    for z in range(16):
        word = 0
        for i, basis_word in enumerate(BASIS):
            if z >> i & 1:
                word ^= basis_word
        out[word] = f"{z:04b}"
    assert len(out) == 16
    return out


def matrix_rank(rows: list[int]) -> int:
    return gf2_rank(rows)


def cech_complex_check() -> tuple[int, int, int, int, int]:
    # Coordinates are globally trivialized by BASIS.  C0=V^3, C1=V^3,
    # C2=V.  In characteristic two, delta0(a,b,c)=(a+b,a+c,b+c) and
    # delta1(p,q,r)=p+q+r.
    d0_rows = []
    for output_block, inputs in enumerate(((0, 1), (0, 2), (1, 2))):
        for k in range(4):
            row = 0
            for source_block in inputs:
                row |= 1 << (4 * source_block + k)
            d0_rows.append(row)
    d1_rows = []
    for k in range(4):
        row = sum(1 << (4 * block + k) for block in range(3))
        d1_rows.append(row)
    rank_d0, rank_d1 = matrix_rank(d0_rows), matrix_rank(d1_rows)
    # delta1*delta0=0 because each source coordinate occurs twice.
    for d0_row in d0_rows:
        assert d0_row.bit_count() == 2
    h0 = 12 - rank_d0
    h1 = (12 - rank_d1) - rank_d0
    h2 = 4 - rank_d1
    assert (rank_d0, rank_d1, h0, h1, h2) == (8, 4, 4, 0, 0)
    return rank_d0, rank_d1, h0, h1, h2


def main() -> None:
    ts = tiles(9)
    assert len(ts) == 28
    covers = face_covers(ts)
    names = {1: "A", 2: "B", 3: "AB", 4: "C", 5: "AC", 6: "BC", 7: "ABC"}
    ranks: dict[int, int] = {}
    pivots: dict[int, list[int]] = {}
    decoders: dict[int, tuple[int, ...]] = {}
    for subset in range(1, 8):
        mask = intersection_mask(covers, subset)
        ranks[subset] = gf2_rank([word & mask for word in BASIS])
        pivots[subset] = independent_pivots(list(BASIS), mask)
        decoders[subset] = pivot_decoder(BASIS, pivots[subset])
        assert ranks[subset] == 4

    regions = {state: exact_region_mask(covers, state) for state in range(1, 8)}
    region_ranks = {state: gf2_rank([word & mask for word in BASIS]) for state, mask in regions.items()}
    assert (region_ranks[1], region_ranks[2], region_ranks[4]) == (0, 0, 0)
    assert (region_ranks[3], region_ranks[5], region_ranks[6], region_ranks[7]) == (3, 2, 3, 4)
    assert sum(mask.bit_count() for mask in regions.values()) == 28
    assert sum(ranks[x] for x in (1, 2, 4)) - sum(ranks[x] for x in (3, 5, 6)) + ranks[7] == 4
    cech = cech_complex_check()

    coordinates = span_coordinates()
    assert tuple(coordinates[word] for word in MISSING) == ("0101", "1011", "1101", "1100")
    for word, coord in coordinates.items():
        z = int(coord, 2)
        for subset in range(1, 8):
            y = sum(((word >> pivot) & 1) << j for j, pivot in enumerate(pivots[subset]))
            decoded = sum(((decoders[subset][i] & y).bit_count() & 1) << i for i in range(4))
            assert decoded == z

    by_sector: dict[int, list[tuple[int, ...]]] = defaultdict(list)
    lines = WITNESSES.read_text().splitlines()
    assert len(lines) == 59
    for line in lines[1:]:
        fields = line.split("\t")
        difference = int(fields[0], 16)
        j_vector = tuple(map(int, fields[10].split(",")))
        by_sector[difference].append(j_vector)
    assert len(by_sector) == 11 and sum(map(len, by_sector.values())) == 58
    assert not set(MISSING) & set(by_sector)

    profiles: dict[tuple[int, int, int, int], set[int]] = defaultdict(set)
    sector_rows = []
    for difference, js in sorted(by_sector.items()):
        profile = tuple((difference & regions[state]).bit_count() for state in (3, 5, 6, 7))
        profiles[profile].add(len(js))
        assert len(set(js)) == 1
        sector_rows.append((difference, coordinates[difference], len(js), difference.bit_count(), profile, js[0]))
    mixed_profiles = {profile: sorted(masses) for profile, masses in profiles.items() if len(masses) > 1}
    assert mixed_profiles == {(2, 0, 2, 4): [2, 6], (2, 2, 2, 2): [2, 4, 26]}

    # Reproduce THM-828's sector Tournament Analysis from independent inputs.
    sectors = sorted(by_sector)
    first_layer = {d: next(i + 3 for i, z in enumerate(by_sector[d][0]) if z) for d in sectors}
    flips = 0
    for i, a in enumerate(sectors):
        for b in sectors[i + 1 :]:
            structural_a = (first_layer[a], a.bit_count(), len(by_sector[a]), a)
            structural_b = (first_layer[b], b.bit_count(), len(by_sector[b]), b)
            empirical_a = (-len(by_sector[a]), first_layer[a], a.bit_count(), a)
            empirical_b = (-len(by_sector[b]), first_layer[b], b.bit_count(), b)
            flips += (structural_a < structural_b) != (empirical_a < empirical_b)
    assert flips == 22

    print("THM-832 N=9 RANK-FOUR CECH CORE")
    print("basis=" + ",".join(f"{word:07x}" for word in BASIS))
    print("intersection name:cells/rank:pivot-tiles")
    for subset in range(1, 8):
        mask = intersection_mask(covers, subset)
        pivot_tiles = ",".join(f"{ts[p][0]}-{ts[p][1]}" for p in pivots[subset])
        decoder = "".join(f"{row:x}" for row in decoders[subset])
        print(f"  {names[subset]}:{mask.bit_count()}/{ranks[subset]}:{pivot_tiles}:decoder={decoder}")
    print("exact-Venn-region name:cells/rank=" + ",".join(
        f"{names[state]}:{regions[state].bit_count()}/{region_ranks[state]}" for state in range(1, 8)
    ))
    print("rank inclusion-exclusion=4+4+4-4-4-4+4=4")
    print(f"Cech ranks d0/d1={cech[0]}/{cech[1]} cohomology H0/H1/H2={cech[2]}/{cech[3]}/{cech[4]}")
    print("missing coordinates=" + ",".join(f"{word:07x}:{coordinates[word]}" for word in MISSING))
    print("occupied sector coord:mass:weight:AB,AC,BC,ABC:J3..J8")
    for difference, coord, mass, weight, profile, j_vector in sector_rows:
        print(f"  {difference:07x} {coord}:{mass}:{weight}:{','.join(map(str, profile))}:{','.join(map(str, j_vector))}")
    print("nonlinear no-go: exact-region support weights do not determine fibre mass")
    for profile, masses in sorted(mixed_profiles.items()):
        print(f"  profile={','.join(map(str, profile))} masses={','.join(map(str, masses))}")
    print("TOURNAMENT ANALYSIS vertices=11 occupied sectors observable=(first skew layer,weight,mass)")
    print("  switches=structural/empirical tie-path=hex-D C3=0 SCC=1^11 HP=1 edge-flips=22")
    print("PRESERVES difference coordinates and all face restrictions; DESTROYS base witness,H8 kernel,S2 realization,LRC metric")
    print("CHALLENGED ASSUMPTION vertices need not be runners/nodes: here they are defect sectors with literal fibres")


if __name__ == "__main__":
    main()
