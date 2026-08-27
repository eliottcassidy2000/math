"""Independent 24-coordinate F2 rank audit for THM-4264."""


N = 12
NCOLS = 2 * N
ZERO, ONE, W, W2 = 0, 1, 2, 3


def mul(left, right):
    a0, a1 = left & 1, (left >> 1) & 1
    b0, b1 = right & 1, (right >> 1) & 1
    return ((a0 * b0) ^ (a1 * b1)) | (
        ((a0 * b1) ^ (a1 * b0) ^ (a1 * b1)) << 1
    )


def operator_rows(coefficients):
    rows = []
    for index in range(N):
        pair = [0, 0]
        for offset, coefficient in enumerate(coefficients):
            site = (index + offset) % N
            for input_bit in range(2):
                image = mul(coefficient, 1 << input_bit)
                for output_bit in range(2):
                    if (image >> output_bit) & 1:
                        pair[output_bit] ^= 1 << (2 * site + input_bit)
        rows.extend(pair)
    return rows


def telescope_rows():
    return [sum(1 << (2 * index + bit) for index in range(N))
            for bit in range(2)]


def zero_site_rows(*sites):
    return [1 << (2 * site + bit) for site in sites for bit in range(2)]


def rank(rows):
    pivots = {}
    for original in rows:
        row = original
        while row:
            pivot = row.bit_length() - 1
            if pivot in pivots:
                row ^= pivots[pivot]
            else:
                pivots[pivot] = row
                break
    return len(pivots)


def nullity(rows):
    return NCOLS - rank(rows)


def main():
    p0_rows = operator_rows((ONE, W, W2, ONE))
    pv_rows = operator_rows((W2, W, W, ONE))
    ambient = p0_rows + telescope_rows()
    projected = ambient + pv_rows
    ambient_two = ambient + zero_site_rows(0, 1)
    projected_one = projected + zero_site_rows(0)
    projected_two = projected + zero_site_rows(0, 1)

    assert (rank(ambient), nullity(ambient)) == (18, 6)
    assert (rank(projected), nullity(projected)) == (20, 4)
    assert nullity(ambient_two) == 2
    assert nullity(projected_one) == 2
    assert nullity(projected_two) == 0

    print("VISIBLE-GATED TWO-EDGE INDEPENDENT BINARY-RANK AUDIT")
    print(f"ambient_rank={rank(ambient)} nullity={nullity(ambient)} "
          f"words={1 << nullity(ambient)}")
    print(f"projected_rank={rank(projected)} nullity={nullity(projected)} "
          f"words={1 << nullity(projected)}")
    print(f"ambient_two_prefix_zero nullity={nullity(ambient_two)} "
          f"words={1 << nullity(ambient_two)}")
    print(f"projected_one_prefix_zero nullity={nullity(projected_one)} "
          f"words={1 << nullity(projected_one)}")
    print(f"projected_two_prefix_zero nullity={nullity(projected_two)} "
          f"words={1 << nullity(projected_two)}")
    print("VISIBLE_TWO_EDGE_RANK_AUDIT: PASS scope=F4_module_only")


if __name__ == "__main__":
    main()
