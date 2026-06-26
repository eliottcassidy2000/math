"""HYP-2937 scout: lift marked C27 shell transfers into q=3 unital blocks.

Assumption: HYP-2937 is the proposed "marked C27 transfer" frame from the
current prompt.  The local repo does not yet contain a HYP-2937 file.

Construction:
  * Build the classical Hermitian unital of order q=3 in PG(2,9).
  * Use its affine chart of 27 points plus infinity.
  * Label the 27 affine points by C27 residues via base-3 digits.
  * Attach AP/Goddyn-Wong C27 shell labels.
  * For each shell doubling transfer a -> fold(2a mod 27), list the unique
    q=3 unital 4-point blocks carrying the source-target point pairs.

The point is to test whether the q=3 unital can be a real incidence carrier
for the C27 transfer story, and where it loses the cyclic 3-adic carry data.
"""

from __future__ import annotations

from collections import Counter, defaultdict
from itertools import combinations, product
from math import gcd


MOD = 27
INF = 27


# GF(9) = F3[w]/(w^2 + 1), so w^2 = 2 in F3.
# Encode a + b*w as a + 3*b, with a,b in {0,1,2}.
def gf_a(x: int) -> int:
    return x % 3


def gf_b(x: int) -> int:
    return x // 3


def gf(a: int, b: int = 0) -> int:
    return (a % 3) + 3 * (b % 3)


def gf_add(x: int, y: int) -> int:
    return gf(gf_a(x) + gf_a(y), gf_b(x) + gf_b(y))


def gf_neg(x: int) -> int:
    return gf(-gf_a(x), -gf_b(x))


def gf_sub(x: int, y: int) -> int:
    return gf_add(x, gf_neg(y))


def gf_mul(x: int, y: int) -> int:
    a, b = gf_a(x), gf_b(x)
    c, d = gf_a(y), gf_b(y)
    return gf(a * c + 2 * b * d, a * d + b * c)


def gf_pow(x: int, n: int) -> int:
    out = 1
    base = x
    while n:
        if n & 1:
            out = gf_mul(out, base)
        base = gf_mul(base, base)
        n >>= 1
    return out


def gf_inv(x: int) -> int:
    if x == 0:
        raise ZeroDivisionError("GF(9) zero inverse")
    # GF(9)^* has order 8.
    return gf_pow(x, 7)


def gf_str(x: int) -> str:
    a, b = gf_a(x), gf_b(x)
    if b == 0:
        return str(a)
    if a == 0:
        return "w" if b == 1 else "2w"
    return f"{a}+{'' if b == 1 else '2'}w"


def normalize(v: tuple[int, int, int]) -> tuple[int, int, int]:
    for x in v:
        if x != 0:
            inv = gf_inv(x)
            return tuple(gf_mul(y, inv) for y in v)  # type: ignore[return-value]
    raise ValueError("zero projective vector")


def hermitian_value(p: tuple[int, int, int]) -> int:
    x, y, z = p
    left = gf_add(gf_mul(gf_pow(y, 3), z), gf_mul(y, gf_pow(z, 3)))
    right = gf_pow(x, 4)
    return gf_sub(left, right)


def all_projective_points() -> list[tuple[int, int, int]]:
    pts = {
        normalize((x, y, z))
        for x, y, z in product(range(9), repeat=3)
        if (x, y, z) != (0, 0, 0)
    }
    return sorted(pts)


def all_projective_lines() -> list[tuple[int, int, int]]:
    return all_projective_points()


def line_contains(line: tuple[int, int, int], point: tuple[int, int, int]) -> bool:
    a, b, c = line
    x, y, z = point
    val = gf_add(gf_add(gf_mul(a, x), gf_mul(b, y)), gf_mul(c, z))
    return val == 0


def residue_to_point(r: int) -> tuple[int, int, int]:
    d0 = r % 3
    d1 = (r // 3) % 3
    d2 = (r // 9) % 3
    x = gf(d0, d1)
    norm = gf_pow(x, 4)
    assert gf_b(norm) == 0
    # Trace(y)=y^3+y=2*real(y).  So real(y)=2*norm.
    y = gf(2 * gf_a(norm), d2)
    return normalize((x, y, 1))


def point_label(label: int) -> str:
    return "inf" if label == INF else str(label)


def build_unital() -> tuple[dict[int, tuple[int, int, int]], list[tuple[int, ...]]]:
    unital_points = [p for p in all_projective_points() if hermitian_value(p) == 0]
    assert len(unital_points) == 28

    label_to_point = {r: residue_to_point(r) for r in range(MOD)}
    label_to_point[INF] = normalize((0, 1, 0))
    assert set(label_to_point.values()) == set(unital_points)

    point_to_label = {p: label for label, p in label_to_point.items()}
    blocks: set[tuple[int, ...]] = set()
    tangent_count = 0
    for line in all_projective_lines():
        hit = [point_to_label[p] for p in unital_points if line_contains(line, p)]
        if len(hit) == 4:
            blocks.add(tuple(sorted(hit)))
        elif len(hit) == 1:
            tangent_count += 1
        else:
            raise AssertionError(f"unexpected line intersection size {len(hit)}")
    assert len(blocks) == 63
    assert tangent_count == 28
    return label_to_point, sorted(blocks)


def fold_shell(r: int) -> int:
    r %= MOD
    if r == 0:
        return 0
    return min(r, MOD - r)


def shell_points(shell: int) -> tuple[int, int]:
    assert 1 <= shell <= 13
    return (shell, MOD - shell)


def shell_stratum(shell: int) -> str:
    if shell == 0:
        return "zero"
    g = gcd(shell, MOD)
    if g == 1:
        return "unit"
    if g == 3:
        return "gcd3"
    if g == 9:
        return "gcd9"
    raise AssertionError(shell)


def point_state(label: int) -> str:
    if label == INF:
        return "inf"
    if label == 0:
        return "zero"
    shell = fold_shell(label)
    if shell == 3:
        return "GW_collision_shell3"
    if shell == 12:
        return "GW_hole_shell12"
    if shell == 9:
        return "C27_gcd9_core"
    return shell_stratum(shell)


def block_state(block: tuple[int, ...]) -> tuple[str, ...]:
    return tuple(sorted(point_state(x) for x in block))


def unique_pair_block(blocks: list[tuple[int, ...]], a: int, b: int) -> tuple[int, ...]:
    hits = [block for block in blocks if a in block and b in block]
    assert len(hits) == 1
    return hits[0]


def shell_transfer_rows(blocks: list[tuple[int, ...]]) -> list[dict[str, object]]:
    rows = []
    for shell in range(1, 14):
        target = fold_shell(2 * shell)
        pairs = []
        carrier_blocks = set()
        for a in shell_points(shell):
            for b in shell_points(target):
                if a == b:
                    continue
                block = unique_pair_block(blocks, a, b)
                pairs.append((a, b, block))
                carrier_blocks.add(block)
        rows.append(
            {
                "source": shell,
                "target": target,
                "source_stratum": shell_stratum(shell),
                "target_stratum": shell_stratum(target),
                "pair_count": len(pairs),
                "carrier_count": len(carrier_blocks),
                "carriers": tuple(sorted(carrier_blocks)),
            }
        )
    return rows


def verify_design(blocks: list[tuple[int, ...]]) -> dict[str, object]:
    point_rep = Counter()
    pair_rep = Counter()
    for block in blocks:
        assert len(block) == 4
        for p in block:
            point_rep[p] += 1
        for a, b in combinations(block, 2):
            pair_rep[(a, b)] += 1
    return {
        "blocks": len(blocks),
        "block_size_hist": dict(Counter(len(b) for b in blocks)),
        "point_rep_hist": dict(Counter(point_rep.values())),
        "pair_rep_hist": dict(Counter(pair_rep.values())),
    }


def print_unital_summary(label_to_point: dict[int, tuple[int, int, int]], blocks: list[tuple[int, ...]]) -> None:
    print("[q=3 Hermitian unital verification]")
    print(f"  points={len(label_to_point)}; blocks={len(blocks)}")
    print(f"  design={verify_design(blocks)}")
    print("  affine labels use C27 residues via base-3 digits r=d0+3*d1+9*d2:")
    for r in range(9):
        p = label_to_point[r]
        print(
            "    r={:<2d} -> (X,Y,Z)=({}, {}, {})".format(
                r, gf_str(p[0]), gf_str(p[1]), gf_str(p[2])
            )
        )
    print("  warning: this is a C27 -> F3^3 digit lift, not a cyclic action lift.")
    print()


def print_block_statistics(blocks: list[tuple[int, ...]]) -> None:
    print("[AP/Goddyn-Wong labelled block statistics]")
    stratum_hist = Counter()
    state_hist = Counter()
    for block in blocks:
        strata = []
        for x in block:
            if x == INF:
                strata.append("inf")
            elif x == 0:
                strata.append("zero")
            else:
                strata.append(shell_stratum(fold_shell(x)))
        stratum_hist[tuple(sorted(strata))] += 1
        state_hist[block_state(block)] += 1
    print("  block stratum histogram:")
    for key, val in sorted(stratum_hist.items(), key=lambda kv: (kv[0], kv[1])):
        print(f"    {key}: {val}")
    print("  block states involving GW hole/collision/core:")
    for key, val in sorted(state_hist.items(), key=lambda kv: (kv[0], kv[1])):
        if any("GW_" in item or item == "C27_gcd9_core" for item in key):
            print(f"    {key}: {val}")
    print()


def print_transfer_lift(blocks: list[tuple[int, ...]]) -> None:
    print("[Marked C27 doubling transfers lifted to unital blocks]")
    print("  transfer is shell a -> fold(2a mod 27); GW petal is 12 -> 3.")
    rows = shell_transfer_rows(blocks)
    print("   a -> b   strata        pair carriers  AP/GW mark")
    for row in rows:
        mark = ""
        if row["source"] == 12 and row["target"] == 3:
            mark = "GW 12->24 petal: hole shell12 -> collision shell3"
        elif row["source"] == 3 and row["target"] == 6:
            mark = "inside nonunit C27 orbit"
        elif row["source"] == 9 and row["target"] == 9:
            mark = "gcd9 fixed core"
        print(
            "  {source:2d} -> {target:<2d}  {source_stratum:>5s}->{target_stratum:<5s}"
            "    {pair_count:>2d}       {carrier_count:>2d}      {mark}".format(
                **row, mark=mark
            )
        )
    print()

    gw = next(row for row in rows if row["source"] == 12 and row["target"] == 3)
    print("[GW transfer carrier blocks]")
    for block in gw["carriers"]:  # type: ignore[index]
        labels = ", ".join(point_label(x) for x in block)
        states = ", ".join(point_state(x) for x in block)
        print(f"  {{{labels}}}  states=({states})")
    print()


def print_readout() -> None:
    print("[Readout]")
    print("  Positive: the q=3 unital gives a real lambda=1 block carrier: every")
    print("  marked C27 point-pair transfer lands in a unique 4-point block.")
    print("  Negative/guardrail: the natural 27-point affine chart is F3^3.  The")
    print("  C27 cyclic carry is not an automorphism of this lift, so block patterns")
    print("  are label-dependent unless a stronger C27-compatible marking is proved.")
    print("  Useful proof role: use unital blocks as a pair-incidence ledger after")
    print("  AP/GW shell labels are fixed; do not treat the unital as a standalone")
    print("  LRC invariant.")


def main() -> None:
    print("HYP-2937 MARKED C27 TRANSFERS INTO Q=3 UNITAL BLOCKS")
    print("=" * 78)
    label_to_point, blocks = build_unital()
    print_unital_summary(label_to_point, blocks)
    print_block_statistics(blocks)
    print_transfer_lift(blocks)
    print_readout()


if __name__ == "__main__":
    main()
