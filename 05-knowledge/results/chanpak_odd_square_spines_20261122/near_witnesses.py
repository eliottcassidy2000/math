"""Construct E61 brackets from selected weighted sums of the two flank prefixes."""
from math import isqrt
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "chanpak_decimal_kernel_20261122"))
from fermat_receiver import (A, D, EDGE, TARGET, address_height, flipped, join,
                             literal_receiver, nontags, receiver, tagged_data)
from e15_pair import literal_count, ordinal, star_rows, strong


# (weights, two low bits, displacement from the leading square root, physical record)
WINNING_SLICES = (
    (15, 7, 0, 0, 6803, (735435556549430, 1575929833923500, 118, -844710)),
    (1, 20, 0, 0, 692570, (4814468855604222, 240731609083994, 118, 184814588)),
)


def last_root(bound, quadratic, linear):
    """Largest nonnegative j with quadratic*j**2 + linear*j <= bound."""
    if bound < 0:
        return -1
    return (isqrt(linear**2 + 4*quadratic*bound) - linear) // (2*quadratic)


def slice_candidates(lower, a, b, left_parity, right_parity, offset, tolerance):
    """Every address pair in this slice whose height could be within tolerance."""
    h, q, c, t, e, r, k = lower
    i, j = left_parity, right_parity
    # Extract one suffix bit. f(2X+i)=(1-i)X+f(X)+(1-i)[X odd].
    Q = 4*q
    T = 2*(t+q*j) + h*(1-i)
    C = 2*(c+q*i) + h*(1-j)
    K = k+q*i*j+t*i+c*j+r*i+e*j
    S = isqrt(4*a*b*(TARGET-K)//Q) + offset
    residue = (pow(a, -1, b)*S)%b if b > 1 else 0
    curvature = Q*a*b
    center = residue+b*((Q*S+b*T-a*C-2*a*Q*residue+curvature)//(2*curvature))
    other = (S-a*center)//b
    beta = b*T-a*C+Q*(b*other-a*center)
    excess = K + Q*center*other + T*center + C*other - TARGET
    allowance = h*(S//min(a,b)+2)
    for sign in (-1, 1):
        linear = -sign*beta
        start = max(int(sign==1), last_root(excess-tolerance-1, curvature, linear)+1)
        stop = last_root(excess+allowance+tolerance, curvature, linear)
        for step in range(start, stop+1):
            X, Y = center + sign*b*step, other - sign*a*step
            if min(X,Y)>=0:
                yield 2*X+i, 2*Y+j


def verify_vertex_rows(x, y, delta):
    """Recount Hamilton paths on each full tournament, independently of the chart."""
    recv_rows, recv_tags = literal_receiver(nontags(4))
    tags = 12 | (35 << len(A)) | (recv_tags << (len(A)+len(D)))
    hosts = []
    for middle in (flipped(D, *EDGE), D):
        core = ordinal(ordinal(A, middle), recv_rows)
        hosts.append(star_rows(core, tags, x, y))
    start = x.bit_length()+len(A)
    edge = (start+EDGE[0], start+EDGE[1])
    assert flipped(hosts[0], *edge)==hosts[1]
    sizes = (1,)*x.bit_length() + (len(A),len(D),len(recv_rows)) + (1,)*y.bit_length()
    for orientation, host in enumerate(hosts):
        assert len(host)==118 and strong(host)
        assert literal_count(host,sizes)==TARGET+delta+2*orientation
    return edge


def main():
    prefix = tagged_data(A, 12)
    lower = join(join(prefix, tagged_data(flipped(D, *EDGE), 35)), receiver(4))
    upper = join(join(prefix, tagged_data(D, 35)), receiver(4))
    for a, b, i, j, offset, expected in WINNING_SLICES:
        tolerance = 10**9
        records = []
        for x, y in slice_candidates(lower, a, b, i, j, offset, tolerance):
            low_height = address_height(lower, x, y)
            delta = low_height - TARGET
            if abs(delta)>tolerance:
                continue
            high_height = address_height(upper, x, y)
            assert high_height - low_height == 2
            order = 17 + x.bit_length() + y.bit_length()
            records.append((x,y,order,delta))
        assert records == [expected], (expected, records)
        edge = verify_vertex_rows(expected[0],expected[1],expected[-1])
        print("weighted sum offset", offset, "reconstructed", records, "physical arc", edge)


if __name__ == "__main__":
    main()
