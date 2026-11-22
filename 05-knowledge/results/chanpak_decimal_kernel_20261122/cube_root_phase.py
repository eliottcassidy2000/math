"""Find absolute, phase-matched brackets on a fixed Farey antidiagonal."""
from math import gcd, isqrt
from fermat_receiver import TARGET, address_height
from hyperbola_controls import chart


def bracket(data, target, modulus, reach=2500):
    h, q, c, t, e, r, k = data
    z = isqrt(target // q)
    B = 1 << ((z.bit_length() - 1) // 3)
    # The odd part chooses a compatible antidiagonal, as in the controller proof.
    odd = modulus // (modulus & -modulus)
    g = gcd(h, odd)
    v0 = (2*z + 2*B + 1 + B-1) // B
    v0 += ((target - k - r + t - t*B*v0) * pow(t*B, -1, g)) % g
    below = above = None
    for v in range(v0, v0 + 64*g, g):
        total = B*v - 1
        # At d=x-total/2, the smoothed height is
        # -q*x*x+(q*total+t-c)*x+c*total+k+h*(v-1)*(B-1)/2.
        linear = q*total + t-c
        disc = linear**2 + 4*q*(c*total+k-target) + 2*q*h*(v-1)*(B-1)
        if disc < 0:
            continue
        root = (linear + isqrt(disc)) // (2*q)
        for x in range(root-reach, root+reach+1):
            y = total-x
            if min(x,y) <= 0 or not x&1:
                continue
            diff = address_height(data, x, y) - target
            if diff % modulus:
                continue
            if diff < 0 and (below is None or diff > below[0]):
                below = (diff, x, y)
            if diff > 0 and (above is None or diff < above[0]):
                above = (diff, x, y)
            if not diff:
                print("exact hit",x,y)
    return below, above


if __name__ == "__main__":
    for m in (60, 510, 3600):
        pair = bracket(chart(4), TARGET, m)
        assert all(item is not None for item in pair)
        print("modulus", m, pair)
