#!/usr/bin/env python3
"""Cross-check the primary and independent forty-ray ledgers."""

from pathlib import Path


Q_SMALL = {
    2, 3, 4, 5, 6, 7, 9, 11, 12, 13, 14, 17, 18, 19, 21, 22,
    23, 24, 25, 26, 27, 28, 29, 31, 32, 33, 34, 35, 36, 37,
    38, 39, 41, 43, 44, 45, 46, 47, 48, 49,
}


def rows(path: str) -> dict[int, list[str]]:
    answer = {}
    for line in Path(path).read_text(encoding="ascii").splitlines():
        fields = line.split()
        answer[int(fields[0])] = fields
    return answer


def main() -> None:
    primary = rows("/tmp/lrc14-small-label-primary-ledger.out")
    independent = rows("/tmp/lrc14-small-label-independent-ledger.out")
    assert set(primary) == Q_SMALL == set(independent)
    total_bodies = 0
    for q in sorted(Q_SMALL):
        p = primary[q]
        i = independent[q]
        assert p[1] == i[1]
        assert p[5] == i[5] == i[6] == "14307150"
        assert p[6] == "0"
        assert p[9] == i[9]
        assert p[10] == i[10]
        assert p[11] == i[11]
        assert p[12] == i[14]
        assert p[13] == i[15]
        assert p[14] == i[16]
        assert p[16] == "UNIFORM_COFINAL_TAIL"
        total_bodies += int(p[5])
    assert total_bodies == 572_286_000
    assert int(primary[6][9]) == 614
    assert int(independent[6][12]) == 587
    assert int(primary[25][9]) == 598
    assert int(independent[25][12]) == 590
    assert max(int(primary[q][9]) for q in Q_SMALL - {6, 25}) == 589
    print(
        "PASS Q_ROWS 40 BODY_ROWS 572286000 ESSENTIAL_MISMATCHES 0 "
        "ANALYTIC_MAX 614 Q6_SECOND 587 Q25_MAX 598 Q25_SECOND 590 "
        "OTHER_MAX 589 HYBRID_UNIFORM_R 590"
    )


if __name__ == "__main__":
    main()
