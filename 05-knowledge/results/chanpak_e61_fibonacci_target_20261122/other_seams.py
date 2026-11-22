"""Replay two-sided Euclid-center flips using exact small-core path profiles."""

import sys
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

PACKETS = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PACKETS / "chanpak_gap_star_ports_20261122"))
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
    profiles, hosts18 = certificate(18, (80, 21, 73, 131, 44, 15, 58, 119), (4, 6), 89, 7, 749445)
    assert profiles == [(413, 6200, (1042, 1253), (1266, 1527), 1575),
                        (413, 6200, (1043, 1252), (1266, 1527), 1579)]
    profiles, first = certificate(21, (4, 17, 34, 7, 13, 27), (2, 3), 22, 7, 923014801)
    assert profiles == [(29, 320, (59, 70), (73, 89), 113),
                        (29, 320, (58, 71), (72, 90), 111)]
    _, second = certificate(21, (240, 141, 153, 177, 194, 150, 174, 0), (0, 4), 209, 0, 1331008449)
    for a, b in zip(first, second):
        assert sum(row.bit_count() == 1 for row in a) == 1
        assert sum(row.bit_count() == 1 for row in b) == 2

    def scores(host):
        degrees = tuple(sorted(row.bit_count() for row in host))
        return {degrees, tuple(len(host) - 1 - d for d in reversed(degrees))}

    for host in hosts18:
        assert len(scores(host)) == 2
    for a, b in zip(first, second):
        assert len(scores(a) | scores(b)) == 4

    sys.path.insert(0, str(PACKETS / "chanpak_decimal_kernel_20261122"))
    from e15_pair import A, D, V, MARKS, EDGE, ordinal

    tags = MARKS[0] | MARKS[1] << len(A) | MARKS[2] << (len(A) + len(D))
    for middle in (D, flipped(D, *EDGE)):
        assert len(scores(star_rows(ordinal(ordinal(A, middle), V), tags, 0, 13))) == 2

    spec = spec_from_file_location("platonic_orbit_bound", PACKETS / "chanpak_fermat_solids_seams_20261122/check.py")
    rotations = module_from_spec(spec)
    spec.loader.exec_module(rotations)
    assert rotations.orbit_floor(37) == 32973313212691 > rotations.perfect_shape(21) + 1
    allowed = {}
    for exponent in range(3, 22):
        upper = rotations.perfect_shape(exponent) + 1
        primes = tuple(q for q in rotations.small_prime_factors(upper, 31)
                       if rotations.orbit_floor(q) <= upper)
        allowed.setdefault(primes, []).append(exponent)
    assert allowed == {(): [3, 4, 6, 9, 10, 12, 13, 15, 16, 18, 19, 21],
                       (7,): [5, 8, 11, 20], (11,): [7], (7, 11): [14, 17]}
    print("Upper exponent groups, by allowed automorphism primes:", allowed)
