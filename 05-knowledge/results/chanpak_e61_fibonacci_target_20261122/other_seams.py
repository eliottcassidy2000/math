"""Replay two-sided Euclid-center flips using exact small-core path profiles."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "chanpak_gap_star_ports_20261122"))
from replay import count_from_profile, flipped, profile, star_rows, strong


def certificate(exponent, core, arc, tags, x, y):
    center = (1 << (exponent - 1)) * ((1 << exponent) - 1)
    profiles = []
    hosts = []
    for side, rows in enumerate((core, flipped(core, *arc))):
        data = profile(rows, tags)
        host = star_rows(rows, tags, x, y)
        height = count_from_profile(data, x, y)
        assert strong(host) and height == center - 1 + 2 * side
        profiles.append(data)
        hosts.append(host)
        print(exponent, len(host), height, data)
    assert flipped(hosts[0], *(x.bit_length() + i for i in arc)) == hosts[1]
    return profiles, hosts


if __name__ == "__main__":
    profiles, _ = certificate(18, (80, 21, 73, 131, 44, 15, 58, 119), (4, 6), 89, 7, 749445)
    assert profiles == [(413, 6200, (1042, 1253), (1266, 1527), 1575),
                        (413, 6200, (1043, 1252), (1266, 1527), 1579)]
    profiles, first = certificate(21, (4, 17, 34, 7, 13, 27), (2, 3), 22, 7, 923014801)
    assert profiles == [(29, 320, (59, 70), (73, 89), 113),
                        (29, 320, (58, 71), (72, 90), 111)]
    _, second = certificate(21, (240, 141, 153, 177, 194, 150, 174, 0), (0, 4), 209, 0, 1331008449)
    for a, b in zip(first, second):
        assert sum(row.bit_count() == 1 for row in a) == 1
        assert sum(row.bit_count() == 1 for row in b) == 2
