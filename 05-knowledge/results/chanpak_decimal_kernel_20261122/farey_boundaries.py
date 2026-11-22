"""Replay the two order-118 absolute witnesses and construct words with four prescribed boundaries."""
from math import gcd

from fermat_receiver import A, D, EDGE, TARGET, address_height, flipped, join, receiver, tagged_data
from local_word_inverse import advance, digits, steer


def boundary_words(chart, target, modulus, lp, ls, rp, rs):
    """Return canonical binary words with the given boundary strings and lower reading target mod 2M."""
    assert lp.startswith("1") and rp.startswith("1") and ls and rs
    assert modulus > 0 and target % 2
    h, q, _, t, _, _, _ = chart
    assert h % 2 and gcd(t, h) == 1
    # Put the right address in residue zero modulo h without changing either end.
    width = h.bit_length()
    right_scale = 1 << len(rs)
    right_mid = (-int(rs, 2) * pow(right_scale, -1, h) - (int(rp, 2) << width)) % h
    right = rp + digits(right_mid, width) + rs
    y = int(right, 2)
    assert y % h == 0

    scale = 1 << len(ls)
    i = int(ls, 2)
    slope = scale * (t + q * y) + h * (scale - 1 - int(ls[::-1], 2))
    baseline = address_height(chart, i, y)
    assert (target - baseline) % 2 == 0
    half = (target - baseline) // 2
    g = gcd(h, modulus)
    a = half * pow(slope, -1, g) % g
    m = modulus // g
    z = (half - slope * a) // g * pow(h // g, -1, m) % m
    word = steer(modulus, advance((0, 0), lp, modulus), (2 * a, z))
    left = lp + word + ls
    assert (address_height(chart, int(left, 2), y) - target) % (2 * modulus) == 0
    return left, right


def main():
    low = join(tagged_data(A, 12), tagged_data(flipped(D, *EDGE), 35))
    high = join(tagged_data(A, 12), tagged_data(D, 35))
    lower, upper = join(low, receiver(4)), join(high, receiver(4))
    for x, y, offset in ((999325040886328, 1159777636980293, -5378156638),
                         (1000669119751650, 1158219846718188, 464575890)):
        assert address_height(lower, x, y) == TARGET + offset
        assert address_height(upper, x, y) == TARGET + offset + 2
        order = 17 + x.bit_length() + y.bit_length()
        print("spine 4, order", order, "midpoint offset", offset, "flanks", x, y)
    for n in range(4, 33, 4):
        chart = join(low, receiver(n))
        for modulus in (2, 4, 8, 16, 25, 81, 289, 841, 2**6 * 17**2 * 29**2):
            targets = range(1, 2 * modulus, 2) if modulus <= 16 else [TARGET]
            for target in targets:
                for ls, rs in (("1"*16, "1"*20), ("0"*12, "10"), ("01", "001"), ("101", "0")):
                    left, right = boundary_words(chart, target, modulus, "11", ls, "101", rs)
                    assert left.startswith("11") and left.endswith(ls)
                    assert right.startswith("101") and right.endswith(rs)
        print("boundary controls passed at spine", n)


if __name__ == "__main__":
    main()
