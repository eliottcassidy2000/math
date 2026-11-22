"""One physical arc across E15 and independent subset counts for its two hosts."""

from pathlib import Path
import sys

PACKETS = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(PACKETS / "chanpak_gap_star_ports_20261122"))
sys.path.insert(0, str(PACKETS / "chanpak_e31_physical_20261122"))
from e31_certificates import ordinal
from replay import flipped, p, profile, star_rows, strong, subpaths


A = (2, 8, 3, 5)
D = (16, 5, 9, 19, 38, 15)
V = (310, 492, 504, 497, 226, 448, 385, 257, 16)
MARKS = (1, 35, 451)
EDGE = (1, 4)


def tagged_data(rows, marks):
    h, q, (se, so), (te, to), k = profile(rows, marks)
    # c, t are endpoint totals; e and r are the signed parity corrections.
    return h, q, se+so, te+to, so-se, to-te, k


def join(a, b):
    """Six endpoint/slot readings and the ordinary apex height of an ordinal join."""
    h, q, c, t, e, r, k = a
    H, Q, C, T, E, R, K = b
    return (h*H, q*Q, Q*c+h*C, q*T+H*t, h*E, H*r,
            H*k+h*(K-H)+c*T+e*R)


def local_tables(rows, marks):
    """Literal paths on each subset with optional first/last apex contact."""
    all_paths = subpaths(rows)
    all_paths[0] = 1
    from_apex = subpaths(rows, marks)
    full = (1 << len(rows)) - 1
    backwards = tuple(r ^ full ^ (1 << i) for i, r in enumerate(rows))
    into_apex = subpaths(backwards, full ^ marks)
    return all_paths, from_apex, into_apex


def literal_count(rows, sizes):
    """Split a physical Hamilton path at its apex and concatenate occupied blocks."""
    apex = len(rows) - 1
    assert sum(sizes) == apex
    tables = []
    first = 0
    for size in sizes:
        mask = (1 << size) - 1
        block = tuple((r >> first) & mask for r in rows[first:first+size])
        later = ((1 << apex) - 1) ^ ((1 << (first+size)) - 1)
        assert all(rows[i] & later == later for i in range(first, first+size))
        tables.append(local_tables(block, (rows[apex] >> first) & mask))
        first += size
    total = 0
    # Empty pre-apex path has last block -1; empty post-apex path has first block b.
    b = len(tables)
    for last_pre in range(-1, b):
        for first_post in range(b+1):
            product = 1
            for j, (hp, from_apex, into_apex) in enumerate(tables):
                full = len(hp)-1
                local = 0
                for pre in range(full+1):
                    post = full ^ pre
                    if (j > last_pre and pre) or (j < first_post and post):
                        continue
                    pre_count = into_apex[pre] if j == last_pre else hp[pre]
                    post_count = from_apex[post] if j == first_post else hp[post]
                    local += pre_count*post_count
                product *= local
            total += product
    return total


def main():
    e15 = 2**14 * (2**15-1)
    a, v = tagged_data(A, MARKS[0]), tagged_data(V, MARKS[2])
    tags = MARKS[0] | MARKS[1] << len(A) | MARKS[2] << (len(A)+len(D))
    sample = star_rows(ordinal(A, D), MARKS[0] | MARKS[1] << len(A), 0, 0)
    assert literal_count(sample, (len(A), len(D))) == subpaths(sample)[-1]
    hosts = []
    for middle, expected in ((D, e15+1), (flipped(D, *EDGE), e15-1)):
        d = tagged_data(middle, MARKS[1])
        ad = join(a, d)
        h, q, c, t, e, r, k = join(ad, v)
        core = ordinal(ordinal(A, middle), V)
        hosts.append(star_rows(core, tags, 0, 13))
        literal = literal_count(hosts[-1], (len(A), len(D), len(V), 4))
        print("D", d, "AD", ad, "ADV", (h, q, c, t, e, r, k))
        print("H", literal, "strong", strong(hosts[-1]))
        assert k + e + c*13 + h*p(13) == literal == expected
        assert len(hosts[-1]) == 24 and strong(hosts[-1])
    assert hosts[1] == flipped(hosts[0], len(A)+EDGE[0], len(A)+EDGE[1])
    print("A", a, "V", v, "physical edge", (len(A)+EDGE[0], len(A)+EDGE[1]))


if __name__ == "__main__":
    main()
