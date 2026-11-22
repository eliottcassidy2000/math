"""Exact finite controls for the perfect seam endpoints and periodic drift formula."""

from pathlib import Path
import sys

results = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(results / "chanpak_gap_star_ports_20261122"))
from replay import flipped, p as digit_cost, profile, star_rows, subpaths
from check_perfect_band import CERTIFICATES


def one_site(u, w):
    common = 2 * u * w + digit_cost(u) + digit_cost(w) + 1
    return (common + u + u % 2, common + w + w % 2)


def check_one_site():
    for exponent in range(3, 12):
        a = 2 ** (exponent - 1)
        center = a * (2 * a - 1)
        candidates = []
        for u in range(0, 2 * a):
            for w in range(u + 1, u + 4):
                h, hp = one_site(u, w)
                if hp - h == 2:
                    candidates.append(h + 1)
        assert max(c for c in candidates if c < center) == center - a // 2
        assert min(c for c in candidates if c > center) == center + 5 * a // 2
        assert center not in candidates
        if exponent <= 5:
            for u, w in ((a - 2, a), (a - 1, a + 1)):
                left = format(u, 'b')
                right = format(w, 'b')[::-1].translate(str.maketrans('01', '10'))
                masks = tuple(star_rows((), 0, int(left + digit + right, 2), 0)
                              for digit in '01')
                assert tuple(subpaths(graph)[-1] for graph in masks) == one_site(u, w)
    print('literal one-apex witnesses through exponent 5 and address checks through exponent 11')


def check_periods():
    for width in range(1, 10):
        k = 2 ** width
        for d in range(k):
            word = format(d, f'0{width}b')
            z = k - 1 - int(word[::-1], 2)
            n = [sum(word[i] == '1' and word[(i+j) % width] == '0'
                     for i in range(width)) for j in range(width)]
            trans = sum(word[i] == '1' and word[(i+1) % width] == '0'
                        for i in range(width))
            numerator = (k-1) * (digit_cost(d) + ((word[0]=='0') * (d % 2))) - d*z
            assert numerator == trans * (k-1) - sum(n[j] * 2**(j-1) for j in range(1,width))
            assert numerator >= width * trans
            for first in (0, 5):
                w = k * first + d
                nxt = k * w + d
                assert (k-1)*(digit_cost(nxt)-digit_cost(w))-z*(nxt-w) == numerator
            assert (numerator == 0) == (d == 0 or d == k-1)
    print('drifts and cyclic-distance bounds for 1022 nonempty words of width through 9')


def check_exponent31():
    rows = (964,965,984,979,995,975,128,512,192,320)
    twin = flipped(rows, 1, 3)
    a, b = profile(rows, 513), profile(twin, 513)
    assert {a, b} == {(125, 9728, (523,548), (3769,3884), 1021),
                      (125, 9728, (523,548), (3784,3899), 1023)}
    y = 1973896132951836
    lower = 1021 + 1071*y + 125*digit_cost(y)
    assert lower + 1 == 2**30 * (2**31-1)
    from importlib.util import module_from_spec, spec_from_file_location
    spec = spec_from_file_location("exponent31_physical", results / "chanpak_e31_physical_20261122/check.py")
    source = module_from_spec(spec)
    spec.loader.exec_module(source)
    def power_of_two(n):
        return n > 0 and n & (n-1) == 0

    for (d, tags, edge, v, vt, _, order), slope in zip(source.WITNESSES, (30, 30, 66)):
        rows = source.ordinal(d, v)
        tag = tags | vt << len(d)
        lo = profile(rows, tag)
        hi = profile(flipped(rows, *edge), tag)
        assert lo[:3] == hi[:3] and hi[4] - lo[4] == 2
        assert tuple(b-a for a,b in zip(lo[3],hi[3])) == (slope//2, slope//2)
        c = sum(lo[2])
        assert not any(map(power_of_two, (c, c+lo[0])))
    for _, _, rows, edge, tags, x, y in CERTIFICATES:
        h, q, st, _, _ = profile(rows, tags)
        c = q*x+sum(st)
        assert not any(map(power_of_two, (c,c+h)))
    print('three exponent-31 pairs and all twenty fixed-left band templates fail the dyadic endpoint test')


if __name__ == '__main__':
    check_one_site()
    check_periods()
    check_exponent31()
