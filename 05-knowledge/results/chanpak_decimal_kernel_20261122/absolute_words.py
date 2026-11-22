"""Two ordinary isomorphism classes realizing one lower height at order 41."""
from e15_pair import A, D, EDGE, flipped, literal_count, ordinal, star_rows, strong
from fermat_receiver import literal_receiver, nontags


def main():
    receiver, receiver_tags = literal_receiver(nontags(4))
    tags = 12 | (35 << len(A)) | (receiver_tags << (len(A) + len(D)))
    lower_count = 32_103_068_571_795
    degrees = [set(), set()]
    for x, y in ((60, 231007), (1881, 7438)):
        i, j = x.bit_length(), y.bit_length()
        sizes = (1,) * i + (len(A), len(D), len(receiver)) + (1,) * j
        hosts = []
        for orientation, middle in enumerate((flipped(D, *EDGE), D)):
            core = ordinal(ordinal(A, middle), receiver)
            host = star_rows(core, tags, x, y)
            height = literal_count(host, sizes)
            assert height == lower_count + 2 * orientation
            assert len(host) == 41 and strong(host)
            degrees[orientation].add(tuple(sorted(row.bit_count() for row in host)))
            hosts.append(host)
            print(x, y, len(host), height)
        edge = tuple(i + len(A) + v for v in EDGE)
        assert flipped(hosts[0], *edge) == hosts[1]
    assert all(len(spectra) == 2 for spectra in degrees)


if __name__ == "__main__":
    main()
