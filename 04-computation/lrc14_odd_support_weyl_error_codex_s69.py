#!/usr/bin/env python3
"""HYP-2718/S69: odd-support dominance in the LRC14 Weyl error.

The current carrier-product gap is the origin atom

    Q_0 = ProductCover - ActualCover
        = sum_j (-1)^j W_j,

where W_j=sum_{|R|=j}(z_product(R)-z_actual(R)) is the discrepancy of the
factorial moment E[binom(T,j)] for T=#missed inner sectors.  This script asks
whether the Weyl/origin-atom error is dominated by the odd-support side

    odd_support  = sum_{j odd} W_j,
    even_support = sum_{j even} W_j,
    Q_0 = even_support - odd_support.

No proof is claimed.  The aim is to turn "odd support dominated" into a
checkable carrier statement and to locate the exceptions that must be routed
through the HYP-2717 finite low-relation ledger.
"""

from __future__ import annotations

import importlib.util
import itertools
from collections import Counter, defaultdict
from fractions import Fraction as F
from math import gcd
from pathlib import Path


HERE = Path(__file__).resolve().parent
MISS_ZETA_PATH = HERE / "lrc14_multiblock_miss_zeta_layers_codex_20260621.py"
spec = importlib.util.spec_from_file_location("miss_zeta_layers", MISS_ZETA_PATH)
miss_zeta = importlib.util.module_from_spec(spec)
assert spec.loader is not None
spec.loader.exec_module(miss_zeta)


CAP = {
    8: F(2243, 5880),
    9: F(1979, 4004),
    10: F(55, 91),
    11: F(66, 91),
    12: F(6, 7),
}
BLOCK4 = (0, 1, 2, 3)


def fmt(q: F) -> str:
    return f"{q} ({float(q):.9f})"


def relation_height(offsets: tuple[int, ...]) -> int:
    """Primitive l1 height of the shortest obvious relation for two carriers."""
    if len(offsets) != 2:
        return 0
    a, b = offsets
    g = gcd(a, b)
    return abs(a // g) + abs(b // g)


def row_report(name: str, offset_blocks: tuple[tuple[int, tuple[int, ...]], ...]) -> dict[str, object]:
    blocks = tuple(block for _offset, block in offset_blocks)
    offsets = tuple(offset for offset, _block in offset_blocks)
    row = miss_zeta.row_from_blocks(offset_blocks)
    actual_z = miss_zeta.zeta_from_hit_profile(miss_zeta.actual_hit_profile(row))
    product_z = miss_zeta.product_miss_zeta_shared_x(blocks)
    actual = miss_zeta.cover_from_zeta(actual_z)
    product = miss_zeta.cover_from_zeta(product_z)
    delta = tuple(product_z[mask] - actual_z[mask] for mask in range(64))
    W, shadow = miss_zeta.krawtchouk_shadow(delta)
    atom = miss_zeta.atom_profile_from_factorial(W)
    q0 = product - actual
    assert shadow[6] == q0
    assert atom[0] == q0

    even_support = sum(W[j] for j in range(0, 7, 2))
    odd_support = sum(W[j] for j in range(1, 7, 2))
    assert even_support - odd_support == q0
    even_l1 = sum(abs(W[j]) for j in range(0, 7, 2))
    odd_l1 = sum(abs(W[j]) for j in range(1, 7, 2))
    signed_odd_dominates = abs(odd_support) > abs(even_support)
    sign_from_odd = q0 != 0 and (q0 > 0) != (odd_support > 0)
    cap = CAP.get(len(row))
    cap_product = cap - product if cap is not None else None
    ratio = abs(q0) / cap_product if cap_product and cap_product > 0 else None
    support_word = "".join("O" if W[j] and j % 2 else "E" if W[j] else "." for j in range(7))
    layer_signs = "".join("+" if W[j] > 0 else "-" if W[j] < 0 else "0" for j in range(7))
    return {
        "name": name,
        "row": row,
        "blocks": blocks,
        "offsets": offsets,
        "relation_height": relation_height(offsets),
        "actual": actual,
        "product": product,
        "q0": q0,
        "cap_product": cap_product,
        "ratio": ratio,
        "W": W,
        "atom": atom,
        "even_support": even_support,
        "odd_support": odd_support,
        "even_l1": even_l1,
        "odd_l1": odd_l1,
        "odd_l1_share": odd_l1 / (odd_l1 + even_l1) if odd_l1 + even_l1 else F(0),
        "signed_odd_dominates": signed_odd_dominates,
        "sign_from_odd": sign_from_odd,
        "support_word": support_word,
        "layer_signs": layer_signs,
    }


def print_report(r: dict[str, object]) -> None:
    print("=" * 92)
    print(r["name"])
    print(f"  row={r['row']}")
    print(f"  blocks={r['blocks']} offsets={r['offsets']} Hrel={r['relation_height']}")
    print(f"  actual={fmt(r['actual'])}")
    print(f"  product={fmt(r['product'])}")
    print(f"  Q0=Product-Actual={fmt(r['q0'])}")
    if r["ratio"] is not None:
        print(f"  |Q0|/(cap-product)={fmt(r['ratio'])}")
    print(f"  even_support=sum W_j even={fmt(r['even_support'])}")
    print(f"  odd_support=sum W_j odd={fmt(r['odd_support'])}")
    print(f"  check Q0=even-odd={fmt(r['even_support'] - r['odd_support'])}")
    print(f"  even_L1={fmt(r['even_l1'])}")
    print(f"  odd_L1={fmt(r['odd_l1'])}")
    print(f"  odd_L1_share={fmt(r['odd_l1_share'])}")
    print(
        "  odd-support flags: "
        f"signed_odd_dominates={r['signed_odd_dominates']} "
        f"sign_from_odd={r['sign_from_odd']}"
    )
    print(f"  layer signs W_0..W_6={r['layer_signs']} support parity={r['support_word']}")
    print("  W_j:")
    for j, v in enumerate(r["W"]):
        if v:
            print(f"    W_{j}={fmt(v)}")
    print("  Q_t atoms:")
    for t, v in enumerate(r["atom"]):
        if v:
            mark = "  <-- origin/Weyl atom" if t == 0 else ""
            print(f"    Q_{t}={fmt(v)}{mark}")


def tournament(rows: list[dict[str, object]]) -> None:
    print("\n" + "=" * 92)
    print("TOURNAMENT ANALYSIS")
    print("  vertices: tested split rows")
    print("  pairwise observable: odd_L1_share, then signed odd dominance, then cap risk")
    print("  switch/gauge: factorial-support parity before origin-atom scalarization")
    scores = Counter()
    edges = set()
    for i, j in itertools.combinations(range(len(rows)), 2):
        ai = (
            rows[i]["odd_l1_share"],
            rows[i]["signed_odd_dominates"],
            rows[i]["ratio"] or F(0),
        )
        aj = (
            rows[j]["odd_l1_share"],
            rows[j]["signed_odd_dominates"],
            rows[j]["ratio"] or F(0),
        )
        if ai >= aj:
            edges.add((i, j))
            scores[i] += 1
        else:
            edges.add((j, i))
            scores[j] += 1
    cycles = 0
    for a, b, c in itertools.combinations(range(len(rows)), 3):
        if (a, b) in edges and (b, c) in edges and (c, a) in edges:
            cycles += 1
        if (a, c) in edges and (c, b) in edges and (b, a) in edges:
            cycles += 1
    hist = Counter(scores[i] for i in range(len(rows)))
    print(f"  score_hist={dict(sorted(hist.items()))}")
    print(f"  directed_3cycles={cycles}")
    print("  Hamiltonian pressure path:")
    order = sorted(
        range(len(rows)),
        key=lambda i: (
            rows[i]["odd_l1_share"],
            rows[i]["signed_odd_dominates"],
            rows[i]["ratio"] or F(0),
        ),
        reverse=True,
    )
    for i in order:
        r = rows[i]
        print(
            f"    {r['name']}: odd_share={fmt(r['odd_l1_share'])} "
            f"odd_dom={r['signed_odd_dominates']} Q0={fmt(r['q0'])}"
        )


def aggregate(rows: list[dict[str, object]]) -> None:
    print("\n" + "=" * 92)
    print("AGGREGATE SUPPORT PARITY")
    even_l1 = sum(r["even_l1"] for r in rows)
    odd_l1 = sum(r["odd_l1"] for r in rows)
    even_signed = sum(r["even_support"] for r in rows)
    odd_signed = sum(r["odd_support"] for r in rows)
    q0_l1 = sum(abs(r["q0"]) for r in rows)
    print(f"  aggregate even_L1={fmt(even_l1)}")
    print(f"  aggregate odd_L1={fmt(odd_l1)}")
    print(f"  aggregate odd_L1_share={fmt(odd_l1 / (odd_l1 + even_l1))}")
    print(f"  aggregate even_support={fmt(even_signed)}")
    print(f"  aggregate odd_support={fmt(odd_signed)}")
    print(f"  aggregate signed even-odd={fmt(even_signed - odd_signed)}")
    print(f"  aggregate sum |Q0|={fmt(q0_l1)}")
    counts = Counter((r["signed_odd_dominates"], r["sign_from_odd"]) for r in rows)
    print(f"  flag_counts={dict(counts)}")


def main() -> None:
    print("HYP-2718/S69 odd-support Weyl-error scout")
    print("Exact miss-zeta/factorial arithmetic; no proof claimed.\n")
    cases = [
        ("two 4-blocks, ratio 2:1", ((14, BLOCK4), (28, BLOCK4))),
        ("two 4-blocks, moderate gap", ((14, BLOCK4), (30, BLOCK4))),
        ("two 4-blocks, wider gap", ((30, BLOCK4), (80, BLOCK4))),
        ("two 4-blocks, high relation phase", ((15, BLOCK4), (31, BLOCK4))),
        ("two 4-blocks, positive Q0", ((40, BLOCK4), (61, BLOCK4))),
        ("two 4-blocks, positive Q0 high", ((80, BLOCK4), (121, BLOCK4))),
        ("5+3 split", ((20, (0, 1, 2, 3, 4)), (55, (0, 1, 2)))),
        (
            "3+3+2 split",
            ((18, (0, 1, 2)), (45, (0, 1, 2)), (90, (0, 1))),
        ),
        (
            "five 2-blocks",
            (
                (15, (0, 1)),
                (30, (0, 1)),
                (46, (0, 1)),
                (63, (0, 1)),
                (81, (0, 1)),
            ),
        ),
        (
            "seven singleton carriers",
            (
                (19, (0,)),
                (31, (0,)),
                (44, (0,)),
                (58, (0,)),
                (73, (0,)),
                (89, (0,)),
                (106, (0,)),
            ),
        ),
    ]
    rows = [row_report(name, blocks) for name, blocks in cases]
    for r in rows:
        print_report(r)
    aggregate(rows)
    tournament(rows)
    print("\nSYNTHESIS")
    print("  In the factorial-origin basis, Q0 is literally even-support minus")
    print("  odd-support.  The exact evidence supports an odd-support L1 envelope,")
    print("  not a naive signed odd cone: aggregate odd_L1 is larger, but most")
    print("  negative Product-Actual rows are signed-even-led after cancellation.")
    print("  The positive-Q0 exceptions are useful, not noise: they are the rows where")
    print("  odd support wins in the signed aggregate and flips the origin atom.")
    print("  A plausible proof target is therefore an odd-support envelope plus a")
    print("  finite ledger for signed-even-led low-height packets, followed by the")
    print("  HYP-2717 high-height tail estimate for the remaining origin atom.")


if __name__ == "__main__":
    main()
