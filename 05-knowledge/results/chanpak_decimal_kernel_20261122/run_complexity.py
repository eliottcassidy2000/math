"""Finite inputs and small integer controls for RUN-COMPLEXITY.md."""
from fermat_receiver import A, D, EDGE, flipped, join, p, receiver, tagged_data


def ones_runs(value):
    return (value & ~(value >> 1)).bit_count()


def ceiling_log2(value):
    return (value - 1).bit_length()


def lowbit(value):
    return value & -value


def check_charts():
    left = join(tagged_data(A, 12), tagged_data(flipped(D, *EDGE), 35))
    q_lowbits = []
    for size in range(4, 33, 4):
        h, q, c, t, _, _, _ = join(left, receiver(size))
        slopes = (c, c + h, t, t + h)
        valuations = tuple(map(lowbit, slopes))
        assert valuations == (2, 1, 16, 1)
        assert all(lowbit(slope) < min(slope, lowbit(q)) for slope in slopes)
        assert q % 5 == 0 and h % 2 == 1
        q_lowbits.append(lowbit(q))
        print("spine", size, "largest power of two dividing q", lowbit(q),
              "and the four slopes", valuations)
    assert q_lowbits == [65536, 131072, 65536, 262144, 65536, 131072, 65536, 524288]


def check_product_bound():
    cases = 0
    for x in range(1, 65):
        for y in range(1, 65):
            width = x.bit_length() + y.bit_length()
            terms = 4 * ones_runs(x) * ones_runs(y)
            for q in (5, 12, 35, 320):
                gap = ceiling_log2(2 * q * terms)
                product = q * x * y
                # Consecutive powers around the product minimize the distance.
                for j in (product.bit_length() - 1 - width, product.bit_length() - width):
                    if j >= 0:
                        numerator = abs(product - (1 << (j + width)))
                        assert (numerator << (terms * gap + 1)) >= (1 << width)
                        cases += 1
    print("sparse-product integer inequalities", cases)


def check_reversal():
    for y in range(1, 4096):
        width = y.bit_length()
        scale = 1 << width
        reversed_digits = int(format(y, "b")[::-1], 2)
        # f(y) = y*(1-radical_inverse_2(y)) + O(width), with explicit slack.
        assert abs(scale * p(y) - y * (scale - reversed_digits)) <= (width + 1) * scale
    print("reversed-phase estimates", 4095)


if __name__ == "__main__":
    check_charts()
    check_product_bound()
    check_reversal()
